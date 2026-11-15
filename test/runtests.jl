# Tests of ParameterSpaces on small canonical parameter spaces. Run from the repository root with
#
#     julia --project -e 'using Pkg; Pkg.test()'

#+++ Imports and helpers
using Test
using ParameterSpaces
using Symbolics

"The value of the expression `ex` of `space` with the named `values` substituted (cbrt and sqrt are evaluated too)."
function evaluate_at(space, ex; values...)
    subs = merge(Dict{Num, Any}(space.vars[space.index[n]] => v for (n, v) in values), space.constants)
    return real(ParameterSpaces.to_number(substitute(ex, subs)))
end

"The message of the exception that `f()` throws."
error_message(f) = try f(); "" catch e; sprint(showerror, e) end
#---

#+++ Canonical spaces
"Pipe flow: the Reynolds number, with the density ρ and the dynamic viscosity μ = ρν as optional variables."
function pipe_flow()
    @variables U D ν Re ρ μ
    return ParameterSpace([U, D, ν, Re, ρ, μ], [Re * ν ~ U * D, μ ~ ρ * ν];
                          units = (U = "m/s", D = "m", ν = "m^2/s", ρ = "kg/m^3", μ = "kg/(m*s)"),
                          optional = (:ρ, :μ))
end

"Ship-model similitude: matching both the Reynolds and the Froude number couples the speed U and length L."
function ship_model()
    @variables U L ν g Re Fr
    return ParameterSpace([U, L, ν, g, Re, Fr], [Re * ν ~ U * L, Fr^2 * g * L ~ U^2];
                          units = (U = "m/s", L = "m", ν = "m^2/s", g = "m/s^2"))
end

"Simple pendulum: T = 2π √(L/g), with π kept symbolic."
function pendulum()
    @variables T L g π
    return ParameterSpace([T, L, g], [T^2 * g ~ 4π^2 * L];
                          units = (T = "s", L = "m", g = "m/s^2"), constants = Dict(π => Base.π))
end

"Pitot tube (Bernoulli): the total pressure p₀ = p + ρU²/2, a relation that is not a power law."
function pitot()
    @variables p₀ p ρ U
    return ParameterSpace([p₀, p, ρ, U], [2p₀ ~ 2p + ρ * U^2];   # p₀ = p + ρU²/2, written without the fraction
                          units = (p₀ = "kg/(m*s^2)", p = "kg/(m*s^2)", ρ = "kg/m^3", U = "m/s"))
end

"Ohm's law V = I R, with a voltage and a current of either sign and a positive resistance."
function ohm()
    @variables V I R
    return ParameterSpace([V, I, R], [V ~ I * R];
                          units = (V = "kg*m^2/(s^3*A)", I = "A", R = "kg*m^2/(s^3*A^2)"),
                          domains = (V = :real, I = :real))
end
#---

@testset "ParameterSpaces" begin

    #+++ Construction
    @testset "construction" begin
        @test degrees_of_freedom(pipe_flow()) == 4    # 6 variables, 2 equations
        @test degrees_of_freedom(ship_model()) == 4
        @test degrees_of_freedom(pendulum()) == 2
        @test degrees_of_freedom(pitot()) == 3
        @test ParameterSpaces.parse_units("kg*m^2/(s^3*A)") == (1//1, 2//1, -3//1, -1//1, 0//1, 0//1, 0//1)
        @test ParameterSpaces.parse_units("1/s") == (0//1, 0//1, -1//1, 0//1, 0//1, 0//1, 0//1)
        @test ParameterSpaces.parse_units("1") == ParameterSpaces.NO_UNITS
        @test ParameterSpaces.format_units(ParameterSpaces.parse_units("kg/(m*s^2)")) == "kg m⁻¹ s⁻²"

        @variables a b c
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b^2]; units = (a = "m", b = "m"))   # m against m²
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b]; units = (a = "kg", b = "m"))    # kg against m
        @test_throws ArgumentError ParameterSpace([a], [a ~ b])                                     # b is not declared
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b]; units = (c = "m",))              # c is not a variable
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b]; domains = (a = :imaginary,))
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b]; units = (a = "ft", b = "ft"))    # not an SI base unit
    end
    #---

    #+++ One equation, solved for any of its variables
    @testset "Reynolds number" begin
        space = pipe_flow()
        sol = solve_parameters(space; U = 2.0, D = 0.05, ν = 1e-6)
        @test sol[:Re] ≈ 1e5
        @test sol.missing_inputs == 0
        @test Set(sol.undetermined) == Set([:ρ, :μ])   # optional, so not missing
        @test_throws ArgumentError sol[:μ]

        inverse = solve_parameters(space; Re = 1e5, D = 0.05, ν = 1e-6)
        @test inverse[:U] ≈ 2.0
        @test solve_parameters(space; Re = 1e5, U = 2.0, ν = 1e-6)[:D] ≈ 0.05
        @test evaluate_at(space, formula(inverse, :U); Re = 2.0, D = 3.0, ν = 5.0) ≈ 2 * 5 / 3
        @test isequal(formula(inverse, :Re), space.vars[space.index[:Re]])   # a given variable is its own formula
        @test formula(inverse, :μ) === nothing                                # undetermined

        water = solve_parameters(space; Re = 1e5, D = 0.05, ρ = 1000.0, μ = 1e-3)   # ν from μ and ρ
        @test water[:ν] ≈ 1e-6
        @test water[:U] ≈ 2.0
        @test isempty(water.undetermined)
    end
    #---

    #+++ Unknowns that appear together in several equations
    @testset "coupled unknowns (ship model)" begin
        space = ship_model()
        sol = solve_parameters(space; Re = 1e7, Fr = 0.3, ν = 1e-6, g = 9.81)
        U = cbrt(1e7 * 1e-6 * 0.3^2 * 9.81)
        @test sol[:U] ≈ U
        @test sol[:L] ≈ U^2 / (0.3^2 * 9.81)
        x = (Re = 2.0, Fr = 0.5, ν = 1.5, g = 3.0)
        Ux = cbrt(x.Re * x.ν * x.Fr^2 * x.g)
        @test evaluate_at(space, formula(sol, :U); x...) ≈ Ux
        @test evaluate_at(space, formula(sol, :L); x...) ≈ Ux^2 / (x.Fr^2 * x.g)
    end
    #---

    #+++ Constants, too few values and too many
    @testset "constants, missing and conflicting inputs (pendulum)" begin
        space = pendulum()
        sol = solve_parameters(space; L = 1.0, g = 9.81)
        @test sol[:T] ≈ 2π * sqrt(1 / 9.81)
        @test occursin("π", string(formula(sol, :T)))
        @test evaluate_at(space, formula(sol, :T); L = 2.0, g = 3.0) ≈ 2π * sqrt(2 / 3)
        @test solve_parameters(space; T = 2.0, g = 9.81)[:L] ≈ 9.81 * (2 / 2π)^2

        only_T = solve_parameters(space; T = 2.0)
        @test only_T.missing_inputs == 1
        picks, alternatives = closing_inputs(only_T)
        @test picks == [:L]
        @test alternatives[:L] == [:g]
        @test occursin("give it (or one of g)", error_message(() -> only_T[:L]))

        @test solve_parameters(space; T = 2π * sqrt(1 / 9.81), L = 1.0, g = 9.81) isa ParameterSolution
        @test occursin("inconsistent", error_message(() -> solve_parameters(space; T = 2.1, L = 1.0, g = 9.81)))
        @test_throws ArgumentError solve_parameters(space; T = 2.1, L = 1.0, g = 9.81)
        @test_throws ArgumentError solve_parameters(space; L = -1.0)   # outside its domain
    end
    #---

    #+++ A relation that is not a power law
    @testset "sums (Pitot tube)" begin
        space = pitot()
        sol = solve_parameters(space; p₀ = 101825.0, p = 101325.0, ρ = 1.2)
        @test sol[:U] ≈ sqrt(2 * 500 / 1.2)   # the positive root
        @test solve_parameters(space; p₀ = 101825.0, U = sol[:U], ρ = 1.2)[:p] ≈ 101325.0
        @test evaluate_at(space, formula(sol, :U); p₀ = 3.0, p = 1.0, ρ = 2.0) ≈ sqrt(2 * 2 / 2)
    end
    #---

    #+++ Signs and zeros
    @testset "signs and zeros (Ohm's law)" begin
        space = ohm()
        @test solve_parameters(space; V = -1.5, R = 3.0)[:I] ≈ -0.5
        @test solve_parameters(space; I = -0.5, R = 3.0)[:V] ≈ -1.5
        @test_throws ArgumentError solve_parameters(space; V = 1.5, I = -0.5)   # R would be negative
        @test_throws ArgumentError solve_parameters(space; R = -3.0)

        open = solve_parameters(space; V = 0.0)   # with R > 0, V = 0 forces I = 0
        @test open[:I] == 0
        @test open.undetermined == [:R]
        @test isequal(formula(open, :I), Num(0))
        @test_throws ArgumentError solve_parameters(space; V = 0.0, I = 0.5)   # I R = 0 is impossible
    end
    #---

    #+++ Printing
    @testset "printing" begin
        sol = solve_parameters(pitot(); p₀ = 101825.0, p = 101325.0, ρ = 1.2)
        text = sprint(show, MIME"text/plain"(), sol)
        @test occursin("3 given, 1 solved", text)
        @test occursin("kg m⁻¹ s⁻²", text)
        @test occursin("ParameterSpace with 3 variables, 1 equation and 2 degrees of freedom", sprint(show, pendulum()))
    end
    #---
end
