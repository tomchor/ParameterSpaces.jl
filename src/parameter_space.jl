# parameter_space.jl
#
# The generic machinery of ParameterSpaces: a set of variables linked by equations, and a solver
# that computes the other variables from values for any sufficient subset of them.
#
# A solve has three steps:
# 1. Structure. The rank of the Jacobian of the equations with respect to the unknowns, at a generic
#    point that satisfies the equations, tells which unknowns the given values determine, how many
#    values are missing, which equations are independent, and which are redundant (they over-specify
#    the input). Zeros are handled first: a given zero can force other variables to be zero.
# 2. Solution. Equations with a single unknown are solved one at a time, and any coupled remainder
#    as a smallest block of equations, with `symbolic_solve` (equations of the form a uⁿ + b = 0 are
#    solved directly). Each variable's domain picks the root.
# 3. Consistency. The redundant equations are evaluated at the solution, with a tolerance.
# Formulas repeat step 2 with the given variables kept symbolic.

#+++ Units
const BASE_UNITS  = (:kg, :m, :s, :A, :K, :mol, :cd)        # the SI base units
const NO_UNITS    = ntuple(_ -> 0//1, length(BASE_UNITS))
const UNIT_SCALES = (2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0)   # one factor per base unit, for check_dimensions

"""
    parse_units(str)

Exponents of the SI base units (kg, m, s, A, K, mol, cd) in a unit string such as `"m/s^2"`,
`"kg/(m*s^2)"` or `"1/s"`. `"1"` is dimensionless.
"""
function parse_units(str::AbstractString)
    exponents = unit_exponents(Meta.parse(str))
    return Tuple(get(exponents, u, 0//1) for u in BASE_UNITS)
end

unit_exponents(x::Number) = isone(x) ? Dict{Symbol, Rational{Int}}() :
                                      throw(ArgumentError("units cannot contain the factor $x"))

function unit_exponents(x::Symbol)
    x in BASE_UNITS || throw(ArgumentError("unknown unit $x: use $(join(BASE_UNITS, ", "))"))
    return Dict(x => 1//1)
end

function unit_exponents(ex::Expr)
    Meta.isexpr(ex, :call) || throw(ArgumentError("cannot read the units $ex"))
    op, args = ex.args[1], ex.args[2:end]
    op === :* && return mergewith(+, map(unit_exponents, args)...)
    op === :/ && return mergewith(+, unit_exponents(args[1]), scale_exponents(unit_exponents(args[2]), -1))
    op === :^ && return scale_exponents(unit_exponents(args[1]), rational_exponent(args[2]))
    throw(ArgumentError("cannot read the units $ex"))
end

scale_exponents(exponents, p) = Dict(k => p * v for (k, v) in exponents)

rational_exponent(p::Integer) = p // 1
function rational_exponent(ex::Expr)
    Meta.isexpr(ex, :call) && ex.args[1] === :// && return ex.args[2] // ex.args[3]
    Meta.isexpr(ex, :call) && ex.args[1] === :- && return -rational_exponent(ex.args[2])
    throw(ArgumentError("the unit exponent $ex is not a rational number"))
end

#++++ Printing units
const SUPERSCRIPTS = Dict(zip("0123456789-/", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻ᐟ"))
superscript(p::Rational) = join(SUPERSCRIPTS[c] for c in (isone(denominator(p)) ? string(numerator(p)) :
                                                          string(numerator(p), "/", denominator(p))))

"Units as text, such as `m s⁻²`; empty for a dimensionless variable."
format_units(u) = join((isone(p) ? string(b) : string(b, superscript(p)) for (b, p) in zip(BASE_UNITS, u) if !iszero(p)), " ")
#----
#---

#+++ The parameter space
"""
    ParameterSpace(vars, equations; units = (;), domains = (;), optional = (), constants = Dict())

A set of variables `vars` (Symbolics variables, in priority order: the order in which missing inputs
are suggested) linked by `equations` (`lhs ~ rhs`).

- `equations` must be polynomial in the variables, without division by a variable or radicals,
  because `symbolic_solve` solves polynomial systems: write `T^2 * g ~ 4π^2 * L` for T = 2π √(L/g).
- `units` maps variable names to unit strings over the SI base units (kg, m, s, A, K, mol, cd), such
  as `"kg/(m*s^2)"`. Variables not listed are dimensionless. Every equation is checked for
  dimensional consistency.
- `domains` maps names to `:positive` (the default), `:nonnegative` or `:real`. Given values outside
  a variable's domain are rejected, and roots outside it are discarded.
- `optional` lists variables that may stay undetermined without counting as missing.
- `constants` maps symbols, such as a symbolic π, to their values. They stay symbolic in formulas
  and must be dimensionless; give a dimensional constant, such as g, as a variable instead.

# Example
```julia
@variables U D ν Re
space = ParameterSpace([U, D, ν, Re], [Re * ν ~ U * D];   # Re = U D / ν
                       units = (U = "m/s", D = "m", ν = "m^2/s"))
solve_parameters(space; Re = 1e5, D = 0.05, ν = 1e-6)[:U]   # 2.0
```
"""
struct ParameterSpace
    vars            :: Vector{Num}
    names           :: Vector{Symbol}
    index           :: Dict{Symbol, Int}
    equations       :: Vector{Equation}
    residuals       :: Vector{Num}                         # lhs - rhs of each equation
    units           :: Vector{NTuple{length(BASE_UNITS), Rational{Int}}}   # exponents of the base units
    domains         :: Vector{Symbol}
    optional        :: BitVector
    constants       :: Dict{Num, Float64}
    constant_values :: Vector{Float64}                     # in the order of the compiled functions' arguments
    lhs             :: Function                            # compiled left sides, at [variables; constants]
    rhs             :: Function                            # compiled right sides
    jacobian        :: Function                            # compiled ∂(lhs - rhs)/∂variables
    points          :: Dict{Vector{Symbol}, Vector{Float64}}   # points on the equations, by variables at zero
end

function ParameterSpace(vars::AbstractVector, equations::AbstractVector{Equation};
                        units = (;), domains = (;), optional = (), constants = Dict{Num, Float64}())
    #++++ Names and keyword arguments
    vars  = collect(Num, vars)
    names = Symbolics.tosymbol.(vars)
    allunique(names) || throw(ArgumentError("the variable names must be unique"))
    index = Dict(n => i for (i, n) in enumerate(names))
    for n in (keys(units)..., keys(domains)..., optional...)
        haskey(index, n) || throw(ArgumentError("$n is not one of the variables"))
    end
    constants = Dict{Num, Float64}(constants)
    #----

    #++++ Every symbol in the equations must be a variable or a constant
    residuals = [Num(eq.lhs) - Num(eq.rhs) for eq in equations]
    allowed = Set(Symbolics.unwrap.(vcat(vars, collect(keys(constants)))))
    for (eq, r) in zip(equations, residuals), s in Symbolics.get_variables(r)
        s in allowed || throw(ArgumentError("$eq contains $s, which is neither a variable nor a constant"))
    end
    #----

    #++++ Units, domains, optional variables
    unit_list   = [haskey(units, n) ? parse_units(units[n]) : NO_UNITS for n in names]
    domain_list = [Symbol(get(domains, n, :positive)) for n in names]
    for (n, d) in zip(names, domain_list)
        d in (:positive, :nonnegative, :real) ||
            throw(ArgumentError("the domain of $n is $d: use :positive, :nonnegative or :real"))
    end
    #----

    #++++ Compiled sides and Jacobian, as functions of [variables; constants]
    constant_symbols = collect(keys(constants))
    compile(ex) = first(build_function(ex, vcat(vars, constant_symbols); expression = Val(false)))
    space = ParameterSpace(vars, names, index, collect(equations), residuals, unit_list, domain_list,
                           BitVector(n in optional for n in names), constants,
                           [constants[c] for c in constant_symbols],
                           compile([Num(eq.lhs) for eq in equations]), compile([Num(eq.rhs) for eq in equations]),
                           compile(Symbolics.jacobian(residuals, vars)), Dict{Vector{Symbol}, Vector{Float64}}())
    #----
    check_dimensions(space)
    return space
end

"""
    degrees_of_freedom(space)

How many values determine every variable of `space`: the number of variables minus the rank of the
equations' Jacobian at a generic point that satisfies them.
"""
function degrees_of_freedom(space::ParameterSpace)
    x = point_on_equations(space, Symbol[])
    return length(space.vars) - rank(space.jacobian(arguments_at(space, x)) .* transpose(x); rtol = 1e-9)
end

function Base.show(io::IO, space::ParameterSpace)
    counted(n, word) = string(n, " ", word, n == 1 ? "" : "s")
    print(io, "ParameterSpace with ", counted(length(space.vars), "variable"), ", ",
              counted(length(space.equations), "equation"), " and ",
              counted(degrees_of_freedom(space), "degree"), " of freedom")
end

function Base.show(io::IO, ::MIME"text/plain", space::ParameterSpace)
    show(io, space)
    print(io, ":")
    for eq in space.equations
        print(io, "\n  ", eq)
    end
end
#---

#+++ Helpers
"The arguments of the compiled functions for the variable values `x`."
arguments_at(space::ParameterSpace, x) = vcat(x, space.constant_values)

"Values of order one for every variable, drawn from a seeded generator."
generic_values(space::ParameterSpace; seed = 1) = (rng = Xoshiro(seed); [exp(randn(rng) / 2) for _ in space.vars])

"Whether the expression `x` is zero once expanded."
is_zero(x::Number) = iszero(x)
is_zero(x) = (v = Symbolics.value(expand(x)); v isa Number && iszero(v))

"The names of the variables of `space` that appear in the expression `ex`."
names_in(space::ParameterSpace, ex) = [Symbolics.tosymbol(v) for v in Symbolics.get_variables(ex)
                                       if haskey(space.index, Symbolics.tosymbol(v))]

"Exact (rational) substitutions for the values in `known` (a Dict from names to numbers) and for the constants."
function substitutions(space::ParameterSpace, known)
    subs = Dict{Num, Any}(space.vars[space.index[n]] => exact(v) for (n, v) in known)
    for (c, v) in space.constants
        subs[c] = exact(v)
    end
    return subs
end

exact(v::Integer) = big(v) // 1
exact(v::Rational) = big(numerator(v)) // big(denominator(v))
exact(v::AbstractFloat) = rationalize(BigInt, v)

"A root returned by the solver as a complex number."
function to_number(x)
    v = Symbolics.value(Symbolics.unwrap(x))
    v isa Number && return ComplexF64(v)
    v = try
        Symbolics.symbolic_to_float(x)
    catch
        nothing
    end
    v isa Number || throw(ArgumentError("cannot evaluate the root $x numerically"))
    return ComplexF64(v)
end

"Whether `z` is a real number inside the domain of variable `j`."
function admissible(space::ParameterSpace, j, z)
    abs(imag(z)) <= 1e-9 * abs(z) || return false
    x, domain = real(z), space.domains[j]
    return domain === :real || (domain === :nonnegative && x >= 0) || (domain === :positive && x > 0)
end
#---

#+++ Dimensional consistency
"""
    check_dimensions(space)

Throw if an equation is dimensionally inconsistent: at a random point, scaling each base unit by
its own factor (masses by 2, lengths by 3, times by 5, …) must change both sides of each equation
by the same factor.
"""
function check_dimensions(space::ParameterSpace)
    x  = generic_values(space)
    xs = x .* [prod(s^p for (s, p) in zip(UNIT_SCALES, u)) for u in space.units]
    l,  r  = space.lhs(arguments_at(space, x)),  space.rhs(arguments_at(space, x))
    ls, rs = space.lhs(arguments_at(space, xs)), space.rhs(arguments_at(space, xs))
    for (i, eq) in enumerate(space.equations)
        isapprox(ls[i] * r[i], rs[i] * l[i]; rtol = 1e-10) ||
            throw(ArgumentError("$eq is dimensionally inconsistent with the declared units"))
    end
end
#---

#+++ Zeros
"""
    zero_variables(space, known)

The variables that are zero: those given as zero, and those the equations then force to zero. When
the zeros make one side of an equation vanish, the other side must vanish too: if exactly one of its
factors is an unknown that may be zero, that unknown is zero, and if none is, the input is
contradictory and this throws.
"""
function zero_variables(space::ParameterSpace, known)
    zeros = Set(n for (n, v) in known if iszero(v))
    while true
        added = false
        subs = Dict{Num, Any}(space.vars[space.index[n]] => 0 for n in zeros)
        for eq in space.equations
            sides = (expand(substitute(Num(eq.lhs), subs)), expand(substitute(Num(eq.rhs), subs)))
            vanishing = is_zero.(sides)
            xor(vanishing...) || continue   # neither side vanishes, or both do
            other = vanishing[1] ? sides[2] : sides[1]
            factors = [n for n in names_in(space, other)
                       if is_zero(substitute(other, Dict(space.vars[space.index[n]] => 0)))]
            candidates = [n for n in factors if !haskey(known, n) && space.domains[space.index[n]] !== :positive]
            isempty(candidates) &&
                throw(ArgumentError("contradictory input: with $(join(sort!(collect(zeros)), ", ")) = 0 one side of $eq " *
                                    "vanishes, so the other must too, but none of its factors " *
                                    "($(join(factors, ", "))) can be zero"))
            length(candidates) == 1 || continue   # several candidates: no conclusion
            push!(zeros, only(candidates))
            added = true
        end
        added || return zeros
    end
end
#---

#+++ A generic point on the equations
"""
    point_on_equations(space, zeros)

A point that satisfies every equation, with the variables named in `zeros` at zero and the others at
generic values: a seeded random start projected onto the equations by Gauss–Newton steps in the
logarithms of the variables. Cached for each set of zeros.
"""
function point_on_equations(space::ParameterSpace, zeros)
    key = sort!(collect(Symbol, zeros))
    return get!(space.points, key) do
        x = generic_values(space)
        fixed = [n in key for n in space.names]
        x[fixed] .= 0
        free = findall(.!fixed)
        for _ in 1:200
            a = arguments_at(space, x)
            l, r = space.lhs(a), space.rhs(a)
            w = 1 ./ max.(abs.(l) .+ abs.(r), floatmin())
            residual = (l .- r) .* w
            maximum(abs, residual; init = 0.0) < 1e-13 && return x
            J = space.jacobian(a)[:, free] .* w .* transpose(x[free])
            x[free] .*= exp.(clamp.(-(pinv(J) * residual), -1, 1))
        end
        throw(ArgumentError("found no point that satisfies the equations" *
                            (isempty(key) ? "" : " with $(join(key, ", ")) = 0")))
    end
end
#---

#+++ Structure: which unknowns the given values determine
"""
    structure(space, known)

For the given values `known` (a Dict from names to numbers), returns a NamedTuple with:
- `unknowns`: the indices of the other variables, except those forced to zero (`forced`);
- `independent`: a maximal set of independent equations among those that do not vanish identically
  once the known values and the zeros are substituted, and `redundant`, the others;
- `determined` and `undetermined`: the variables that the known values fix (zeros included), and the
  unknowns they leave free;
- `nullity`: the number of free directions, and `missing_inputs`: how many more values the required
  (non-optional) unknowns need.

It uses the rank and null space of the Jacobian with respect to the unknowns, in relative terms, at a
generic point that satisfies the equations with the same zeros.
"""
function structure(space::ParameterSpace, known::AbstractDict; rtol = 1e-9)
    zeros = zero_variables(space, known)
    forced = [n for n in space.names if n in zeros && !haskey(known, n)]
    unknowns = [j for (j, n) in enumerate(space.names) if !haskey(known, n) && n ∉ zeros]
    subs = substitutions(space, merge(known, Dict(n => 0.0 for n in forced)))
    active = [i for (i, r) in enumerate(space.residuals) if !is_zero(substitute(r, subs))]

    #++++ Jacobian with respect to the unknowns, scaled by their values, rows of comparable size
    x = point_on_equations(space, zeros)
    J = space.jacobian(arguments_at(space, x))[active, unknowns] .* transpose(x[unknowns])
    for row in eachrow(J)
        s = maximum(abs, row; init = 0.0)
        s > 0 && (row ./= s)
    end
    #----

    r = isempty(J) ? 0 : rank(J; rtol)
    independent = r == 0 ? Int[] : sort(active[qr(permutedims(J), ColumnNorm()).p[1:r]])
    N = nullspace(J; rtol)   # free directions of the unknowns
    free = [maximum(abs, N[k, :]; init = 0.0) > 1e-8 for k in eachindex(unknowns)]
    required = [!space.optional[j] for j in unknowns]
    missing_inputs = size(N, 2) == 0 ? 0 : rank(N[required, :]; atol = 1e-8)
    return (; unknowns, forced, independent, redundant = setdiff(active, independent),
              determined = sort(vcat(unknowns[.!free], [space.index[n] for n in forced])),
              undetermined = unknowns[free], nullity = size(N, 2), missing_inputs)
end
#---

#+++ Solving the equations, one unknown or one smallest block at a time
"""
    smallest_block(unknowns_of)

The indices of a smallest set of equations whose unknowns are exactly as many as the equations;
`unknowns_of[k]` lists the unknowns of equation `k`.
"""
function smallest_block(unknowns_of)
    n = length(unknowns_of)
    for k in 2:n, mask in 1:(2^n - 1)
        count_ones(mask) == k || continue
        block = [i for i in 1:n if isodd(mask >> (i - 1))]
        length(union(unknowns_of[block]...)) == k && return block
    end
    return collect(1:n)
end

"""
    two_terms(eq, u)

`(a, n, b, m)` when `eq` has the form a uⁿ + b uᵐ with n > m ≥ 0 and a and b free of `u`, and `nothing`
otherwise (including when `eq` is not a polynomial in `u`).
"""
function two_terms(eq, u)
    try
        coefficients, residual = Symbolics.polynomial_coeffs(expand(eq), [u])
        is_zero(residual) || return nothing
        terms = [(Symbolics.degree(Num(monomial), u), c) for (monomial, c) in coefficients if !is_zero(c)]
        length(terms) == 2 || return nothing
        (n, a), (m, b) = sort(terms; by = first, rev = true)
        n > m || return nothing
        return (a, n, b, m)
    catch
        return nothing
    end
end

"-b/a, written as b/(-a) when a carries the minus sign, so that it prints with fewer signs."
negative_ratio(b, a) = startswith(string(a), "-") ? b / expand(-a) : -b / a

"""
    power_roots(eq, u)

The roots of `eq = 0` in `u` that can be real, when the equation has the form a uⁿ + b uᵐ = 0 with a
and b free of `u` (see [`two_terms`](@ref)): zero if m > 0, and with k = n - m and r = -b/a, r for
k = 1, ±√r for k = 2, ∛r for k = 3 and ±r^(1/k) above. These give cleaner formulas than the general
solver. Returns `nothing` for other equations.
"""
function power_roots(eq, u)
    terms = two_terms(eq, u)
    terms === nothing && return nothing
    a, n, b, m = terms
    k, r = n - m, negative_ratio(b, a)
    zero_root = m > 0 ? Any[0] : Any[]
    try
        k == 1 && return vcat(zero_root, [r])
        k == 2 && return vcat(zero_root, [sqrt(r), -sqrt(r)])
        k == 3 && return vcat(zero_root, [cbrt(r)])
        return vcat(zero_root, [r^(1//k), -r^(1//k)])
    catch
        return nothing
    end
end

"""
    linear_elimination(space, equations, unknowns_of, unknowns_in)

An equation among `equations` that is linear in one of its unknowns with a coefficient that involves no
unknown, a uⱼ + b = 0: returns `(i, j, -b/a)`, the index of the equation, the unknown, and the
expression of that unknown in terms of the others, which stays polynomial. Returns `nothing` if no
equation qualifies.
"""
function linear_elimination(space::ParameterSpace, equations, unknowns_of, unknowns_in)
    for (i, eq) in enumerate(equations), j in unknowns_of[i]
        terms = two_terms(eq, space.vars[j])
        terms === nothing && continue
        a, n, b, m = terms
        (n == 1 && m == 0 && isempty(unknowns_in(a))) && return (i, j, negative_ratio(b, a))
    end
    return nothing
end

"""
    propagate(space, equations, unknowns; value, accept, symbolic, partial = false)

Solve `equations` (expressions equal to zero) for the variables with indices `unknowns`. An equation
with a single unknown is solved first. When none is left, an equation that is linear in one unknown
with a known coefficient expresses that unknown through the others, which are then solved first
([`linear_elimination`](@ref)); otherwise a smallest coupled block is solved as a system.
`value(root)` turns a root into a complex number, and the root kept for variable `j` is the only
distinct one with `accept(j, value(root))`. With `symbolic = true` the kept root itself is
substituted into the remaining equations (a formula); otherwise its exact numeric value is. With
`partial = true` a step that fails is skipped, and its unknowns stay unsolved.

Returns a Dict from variable index to `(root, value(root))`.
"""
function propagate(space::ParameterSpace, equations, unknowns; value, accept, symbolic, partial = false)
    index_of = Dict(Symbolics.unwrap(space.vars[j]) => j for j in unknowns)
    unknowns_in(eq) = sort!(unique([index_of[v] for v in Symbolics.get_variables(eq) if haskey(index_of, v)]))
    equations = collect(Num, equations)
    chosen = Dict{Int, Any}()
    eliminated = Pair{Int, Num}[]   # unknowns expressed through others, evaluated at the end
    while !isempty(equations)
        unknowns_of = map(unknowns_in, equations)
        k = findfirst(isempty, unknowns_of)
        k === nothing || (deleteat!(equations, k); continue)   # nothing left to solve in this one
        k = findfirst(u -> length(u) == 1, unknowns_of)

        #++++ No single unknown left: eliminate one through an equation that is linear in it
        elimination = k === nothing ? linear_elimination(space, equations, unknowns_of, unknowns_in) : nothing
        if elimination !== nothing
            i, j, expression = elimination
            push!(eliminated, j => expression)
            deleteat!(equations, i)
            delete!(index_of, Symbolics.unwrap(space.vars[j]))
            subs = Dict(space.vars[j] => expression)
            equations = [symbolic ? substitute(eq, subs) : expand(substitute(eq, subs)) for eq in equations]
            continue
        end
        #----

        block = k === nothing ? smallest_block(unknowns_of) : [k]
        vars = sort!(unique(reduce(vcat, unknowns_of[block])))
        solved = try
            solve_block(space, equations[block], vars; value, accept)
        catch e
            (partial && e isa ArgumentError) || rethrow()
            nothing
        end
        deleteat!(equations, block)
        if solved === nothing   # skipped: these unknowns stay symbolic, and formulas that need them are dropped
            foreach(j -> delete!(index_of, Symbolics.unwrap(space.vars[j])), vars)
            continue
        end
        merge!(chosen, solved)
        subs = Dict(space.vars[j] => (symbolic ? chosen[j][1] : exact(real(chosen[j][2]))) for j in vars)
        equations = [symbolic ? substitute(eq, subs) : expand(substitute(eq, subs)) for eq in equations]
    end

    #++++ The eliminated unknowns, from the last one back (each depends only on unknowns solved after it)
    for (j, expression) in Iterators.reverse(eliminated)
        subs = Dict(space.vars[k] => (symbolic ? root : exact(real(z))) for (k, (root, z)) in chosen)
        root = substitute(expression, subs)
        z = try value(root) catch; ComplexF64(NaN, NaN) end
        if accept(j, z)
            chosen[j] = (root, z)
        elseif !partial
            throw(ArgumentError("no admissible value for $(space.names[j]): $(expression) gives $z"))
        end
    end
    #----
    return chosen
end

"Solve one equation for one unknown, or a block of equations for as many unknowns, and keep the accepted root."
function solve_block(space::ParameterSpace, equations, vars; value, accept)
    names = join(space.names[vars], ", ")
    safe_value(root) = try value(root) catch; ComplexF64(NaN, NaN) end
    if length(vars) == 1
        j = only(vars)
        roots = power_roots(only(equations), space.vars[j])
        roots === nothing && (roots = symbolic_solve(only(equations), space.vars[j]; warns = false))
        roots === nothing && throw(ArgumentError("symbolic_solve could not solve $(only(equations)) for $names"))
        candidates = [Dict(j => (root, safe_value(root))) for root in roots]
    else
        solutions = symbolic_solve(equations, space.vars[vars]; warns = false)
        solutions === nothing && throw(ArgumentError("symbolic_solve could not solve the system $equations for $names"))
        candidates = [Dict(j => (s[space.vars[j]], safe_value(s[space.vars[j]])) for j in vars) for s in solutions]
    end
    accepted = filter(c -> all(accept(j, z) for (j, (_, z)) in c), candidates)
    distinct = unique(c -> [round(real(c[j][2]); sigdigits = 10) for j in vars], accepted)
    length(distinct) == 1 && return only(distinct)
    found = [Tuple(c[j][2] for j in vars) for c in candidates]
    isempty(distinct) && throw(ArgumentError("no admissible value for $names among the roots $found"))
    throw(ArgumentError("several admissible values for $names: $([Tuple(real(c[j][2]) for j in vars) for c in distinct]); " *
                        "give it, or a variable that tells them apart"))
end
#---

#+++ Solutions
"""
    ParameterSolution

Result of [`solve_parameters`](@ref). `sol[:x]` is the value of variable `x`, given or solved, and
`sol.values` holds them all. `sol.undetermined` lists the variables that the given values leave free,
and `sol.missing_inputs` counts the values that the required (non-optional) ones still need;
[`closing_inputs`](@ref) suggests which. [`formula`](@ref) gives closed-form expressions.
"""
mutable struct ParameterSolution
    space          :: ParameterSpace
    given          :: Dict{Symbol, Float64}
    values         :: Dict{Symbol, Float64}               # given and determined
    undetermined   :: Vector{Symbol}
    missing_inputs :: Int
    closing        :: Dict{Symbol, Float64}               # arbitrary values for the undetermined variables
    all_values     :: Dict{Symbol, Float64}               # every variable at the solve's point
    independent    :: Vector{Int}                         # the equations the solve used
    forced         :: Vector{Symbol}                      # unknowns that the given zeros force to zero
    formulas       :: Union{Nothing, Dict{Symbol, Any}}   # filled by the first call to formula
end

function Base.getindex(sol::ParameterSolution, name::Symbol)
    haskey(sol.values, name) && return sol.values[name]
    haskey(sol.space.index, name) || throw(ArgumentError("$name is not one of the variables"))
    picks, alternatives = closing_inputs(sol; optional = sol.space.optional[sol.space.index[name]])
    others = filter(!=(name), picks)
    hint = isempty(picks)  ? "" :
           name ∉ picks    ? "; also giving $(join(picks, ", ")) would determine it" :
           "; give it" * (isempty(alternatives[name]) ? "" : " (or one of $(join(alternatives[name], ", ")))") *
                         (isempty(others) ? "" : " together with $(join(others, ", "))")
    throw(ArgumentError("$name is not determined by the given values$hint"))
end

function Base.show(io::IO, ::MIME"text/plain", sol::ParameterSolution)
    space = sol.space
    nsolved = length(sol.values) - length(sol.given)
    print(io, "ParameterSolution: $(length(sol.given)) given, $nsolved solved, $(length(sol.undetermined)) undetermined")
    for (j, n) in enumerate(space.names)
        haskey(sol.values, n) || continue
        print(io, "\n  ", haskey(sol.given, n) ? "given " : "solved", "  ", rpad(n, 18),
                  round(sol.values[n]; sigdigits = 6), " ", format_units(space.units[j]))
    end
    isempty(sol.undetermined) && return
    print(io, "\n  undetermined: ", join(sol.undetermined, ", "))
    k = sol.missing_inputs
    print(io, k > 0 ? " ($k more input$(k == 1 ? "" : "s") needed; see closing_inputs)" : " (all optional)")
end
#---

#+++ solve_parameters
"""
    solve_parameters(space; rtol = 1e-6, given...)

The variables of `space` that follow from the `given` values (keyword arguments named after the
variables), as a [`ParameterSolution`](@ref).

Given values must lie in their variables' domains. Over-specified input is checked: every equation
left out of the solve must hold at the solution to within `rtol`, relative to the larger side, or the
solve throws an `ArgumentError` that names the equation. The default tolerance accepts rounded inputs
such as Ri = 6.37105.
"""
function solve_parameters(space::ParameterSpace; rtol = 1e-6, given...)
    #++++ Check the given values
    known = Dict{Symbol, Float64}()
    for (n, v) in pairs(given)
        haskey(space.index, n) || throw(ArgumentError("$n is not one of the variables: $(join(space.names, ", "))"))
        (v isa Real && isfinite(v)) || throw(ArgumentError("$n = $v is not a finite real number"))
        domain = space.domains[space.index[n]]
        domain === :positive && v <= 0 && throw(ArgumentError("$n = $v must be positive"))
        domain === :nonnegative && v < 0 && throw(ArgumentError("$n = $v must not be negative"))
        known[n] = v
    end
    #----

    #++++ Fix the undetermined variables at arbitrary values; the determined ones do not depend on them
    s0 = structure(space, known)
    s, closing, rng = s0, Dict{Symbol, Float64}(), Xoshiro(2)
    while !isempty(s.undetermined)
        closing[space.names[first(s.undetermined)]] = exp(randn(rng) / 2)
        s = structure(space, merge(known, closing))
    end
    #----

    #++++ Solve, then check the redundant equations
    work = merge(known, closing, Dict(n => 0.0 for n in s.forced))
    subs = substitutions(space, work)
    equations = [expand(substitute(space.residuals[i], subs)) for i in s.independent]
    chosen = propagate(space, equations, s.unknowns; value = to_number, symbolic = false,
                       accept = (j, z) -> admissible(space, j, z))
    all_values = merge(work, Dict(space.names[j] => real(z) for (j, (_, z)) in chosen))
    check_consistency(space, s.redundant, all_values, rtol)
    #----

    values = merge(known, Dict(space.names[j] => all_values[space.names[j]] for j in s0.determined))
    return ParameterSolution(space, known, values, space.names[s0.undetermined], s0.missing_inputs,
                             closing, all_values, s.independent, s0.forced, nothing)
end

"Throw if one of the equations `redundant` does not hold, to within `rtol`, at `values` (a Dict from names to numbers)."
function check_consistency(space::ParameterSpace, redundant, values, rtol)
    isempty(redundant) && return nothing
    a = arguments_at(space, [values[n] for n in space.names])
    l, r = space.lhs(a), space.rhs(a)
    for i in redundant
        scale = max(abs(l[i]), abs(r[i]))
        abs(l[i] - r[i]) <= rtol * scale && continue
        throw(ArgumentError("inconsistent input: $(space.equations[i]) gives $(round(l[i]; sigdigits = 6)) on the " *
                            "left and $(round(r[i]; sigdigits = 6)) on the right (relative difference " *
                            "$(round(abs(l[i] - r[i]) / scale; sigdigits = 2)), tolerance $rtol)"))
    end
    return nothing
end
#---

#+++ Which extra inputs would close the system
"""
    closing_inputs(sol; optional = false)

A smallest set of extra inputs that would determine every required variable (every variable, with
`optional = true`), picked in the space's priority order, and for each pick the variables that could
replace it. Returns `(picks, alternatives)`, with `alternatives` a Dict from each pick to a vector.
"""
function closing_inputs(sol::ParameterSolution; optional = false)
    space, rng = sol.space, Xoshiro(3)
    open_directions(known) = (s = structure(space, known); optional ? s.nullity : s.missing_inputs)
    with(known, n) = merge(known, Dict(n => exp(randn(rng) / 2)))

    known, picks = Dict{Symbol, Float64}(sol.given), Symbol[]
    remaining = open_directions(known)
    for n in sol.undetermined
        remaining == 0 && break
        trial = with(known, n)
        open_directions(trial) < remaining || continue
        known, remaining = trial, open_directions(trial)
        push!(picks, n)
    end
    alternatives = Dict(p => [n for n in sol.undetermined if n ∉ picks &&
                              open_directions(with(delete!(copy(known), p), n)) == 0] for p in picks)
    return picks, alternatives
end
#---

#+++ Formulas
"""
    formula(sol, name)

Closed-form expression of variable `name` in terms of the given variables (and constants such as π),
or `nothing` if the variable is undetermined or no formula was found. A given variable is its own
formula.
"""
function formula(sol::ParameterSolution, name::Symbol)
    haskey(sol.space.index, name) || throw(ArgumentError("$name is not one of the variables"))
    sol.formulas === nothing && (sol.formulas = solve_symbolically(sol))
    return get(sol.formulas, name, nothing)
end

"""
    solve_symbolically(sol)

Formulas for the variables of `sol`: the propagation of [`solve_parameters`](@ref) repeated with the
given variables kept symbolic, choosing at each step the root whose value matches the numeric solution.
"""
function solve_symbolically(sol::ParameterSolution)
    space = sol.space
    knowns = merge(sol.given, sol.closing)
    unknowns = [j for (j, n) in enumerate(space.names) if !haskey(knowns, n) && n ∉ sol.forced]
    x = [sol.all_values[n] for n in space.names]
    point = merge(Dict{Num, Float64}(zip(space.vars, x)), space.constants)
    zeros = Dict{Num, Any}(space.vars[space.index[n]] => 0 for n in sol.forced)
    equations = [substitute(space.residuals[i], zeros) for i in sol.independent]
    chosen = propagate(space, equations, unknowns; symbolic = true, partial = true,
                       value = root -> to_number(substitute(root, point)),
                       accept = (j, z) -> abs(z - x[j]) <= 1e-4 * abs(x[j]))   # loose: some roots lose digits

    #++++ Keep the formulas of determined variables that contain only given variables and constants
    allowed = Set(vcat(collect(keys(sol.given)), [Symbolics.tosymbol(c) for c in keys(space.constants)]))
    formulas = Dict{Symbol, Any}(n => space.vars[space.index[n]] for n in keys(sol.given))
    for n in sol.forced
        formulas[n] = Num(0)
    end
    for (j, (root, _)) in chosen
        n = space.names[j]
        haskey(sol.values, n) || continue
        f = simplify(Num(root))
        all(v -> Symbolics.tosymbol(v) in allowed, Symbolics.get_variables(f)) && (formulas[n] = f)
    end
    #----
    return formulas
end
#---
