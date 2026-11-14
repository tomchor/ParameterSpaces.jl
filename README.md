# ParameterSpaces.jl

ParameterSpaces holds a set of variables linked by equations and computes the other variables from any sufficient set of values, on either side of an equation. It reports when the given values are not enough and which extra inputs would close the system, checks values that over-specify it, and gives closed-form formulas for what it solves. It is built on [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

## Installation

The package is not registered yet. Add it from this repository (it needs Julia 1.13 or later):

```julia
using Pkg
Pkg.add(url = "https://github.com/tomchor/ParameterSpaces.jl")
```

## Quick start

A space is a list of variables and the equations that link them. The Reynolds number of a pipe flow links the velocity U, the diameter D and the kinematic viscosity ν:

```julia
using ParameterSpaces

@variables U D ν Re
space = ParameterSpace([U, D, ν, Re], [Re * ν ~ U * D];   # Re = U D / ν
                       units = (U = "m/s", D = "m", ν = "m^2/s"))

solve_parameters(space; U = 2, D = 0.05, ν = 1e-6)[:Re]   # 100000.0

sol = solve_parameters(space; Re = 1e5, D = 0.05, ν = 1e-6)
sol[:U]                      # 2.0
formula(sol, :U)             # (Re*ν) / D
latexify(formula(sol, :U))   # the same formula in LaTeX
```

Any variable can be given or solved for. The equations must be polynomial in the variables, without division by a variable or radicals, because Symbolics' `symbolic_solve` requires it, so Re = U D / ν is written `Re * ν ~ U * D`.

## Coupled unknowns

A ship model that matches both the Reynolds and the Froude number of the full-size ship needs one particular speed U and length L. Both equations contain both unknowns, and the solver handles them together:

```julia
@variables U L ν g Re Fr
ship = ParameterSpace([U, L, ν, g, Re, Fr],
                      [Re * ν ~ U * L,        # Re = U L / ν
                       Fr^2 * g * L ~ U^2];   # Fr = U / √(g L)
                      units = (U = "m/s", L = "m", ν = "m^2/s", g = "m/s^2"))

sol = solve_parameters(ship; Re = 1e7, Fr = 0.3, ν = 1e-6, g = 9.81)
sol[:U], sol[:L]    # (2.0668, 4.8383)
formula(sol, :U)    # cbrt((Fr^2)*Re*g*ν)
formula(sol, :L)    # (cbrt((Fr^2)*Re*g*ν)^2) / ((Fr^2)*g)
```

## Missing and conflicting inputs

The period of a simple pendulum, T = 2π √(L/g), links three variables. Declaring π as a symbol, with its value in `constants`, keeps it in the formulas. Building the space inside a function keeps that symbol local, so it doesn't replace `π` elsewhere:

```julia
function pendulum_space()
    @variables T L g π
    return ParameterSpace([T, L, g], [T^2 * g ~ 4π^2 * L];   # T = 2π √(L/g)
                          units = (T = "s", L = "m", g = "m/s^2"), constants = Dict(π => Base.π))
end
pendulum = pendulum_space()

sol = solve_parameters(pendulum; L = 1, g = 9.81)
sol[:T]             # 2.0061
formula(sol, :T)    # sqrt((4L*(π^2)) / g)
```

Values that leave some variables free still give everything they determine. `closing_inputs` returns the extra inputs that would close the system, each with the variables that could replace it:

```julia
only_T = solve_parameters(pendulum; T = 2)
only_T.missing_inputs     # 1
closing_inputs(only_T)    # ([:L], Dict(:L => [:g]))
only_T[:L]                # ArgumentError: L is not determined by the given values; give it (or one of g)
```

Values that disagree are reported with the equation they break:

```julia
solve_parameters(pendulum; T = 2.1, L = 1, g = 9.81)
# ArgumentError: inconsistent input: (T^2)*g ~ 4L*(π^2) gives 43.2621 on the left and 39.4784 on the right (relative difference 0.087, tolerance 1.0e-6)
```

## Relations that are not power laws

Equations can have sums. A Pitot tube measures the total pressure p₀ = p + ρU²/2, from which Bernoulli's equation gives the speed. The velocity's domain (positive by default) picks the positive root:

```julia
@variables p₀ p ρ U
pitot = ParameterSpace([p₀, p, ρ, U], [2p₀ ~ 2p + ρ * U^2];   # p₀ = p + ρU²/2
                       units = (p₀ = "kg/(m*s^2)", p = "kg/(m*s^2)", ρ = "kg/m^3", U = "m/s"))

sol = solve_parameters(pitot; p₀ = 101825, p = 101325, ρ = 1.2)
sol[:U]             # 28.868
formula(sol, :U)    # sqrt((-2p + 2p₀) / ρ)
```

Printing a solution shows each value with its units:

```
ParameterSolution: 3 given, 1 solved, 0 undetermined
  given   p₀                101825.0 kg m⁻¹ s⁻²
  given   p                 101325.0 kg m⁻¹ s⁻²
  given   ρ                 1.2 kg m⁻³
  solved  U                 28.8675 m s⁻¹
```

## Signs and zeros

In Ohm's law, V = I R, the voltage and the current can have either sign, while the resistance stays positive. A zero voltage then forces a zero current:

```julia
@variables V I R
ohm = ParameterSpace([V, I, R], [V ~ I * R];
                     units = (V = "kg*m^2/(s^3*A)", I = "A", R = "kg*m^2/(s^3*A^2)"),
                     domains = (V = :real, I = :real))   # R keeps the default, :positive

solve_parameters(ohm; V = -1.5, R = 3)[:I]   # -0.5
open_circuit = solve_parameters(ohm; V = 0)
open_circuit[:I]                             # 0.0
open_circuit.undetermined                    # [:R]
```

## Defining a space

`ParameterSpace(vars, equations; units, domains, optional, constants)` takes:

- `vars`: the variables, in priority order. This is the order in which `closing_inputs` suggests missing inputs, so put the inputs you would usually give first.
- `equations`: polynomial relations written with `~`. Write a root as a power: w = ∛(Q H) becomes `w^3 ~ Q * H`, and the solver keeps the single real root.
- `units`: unit strings over the SI base units (kg, m, s, A, K, mol, cd), such as `"kg/(m*s^2)"`. Variables that are not listed are dimensionless. Every equation is checked for dimensional consistency when the space is built.
- `domains`: `:positive` (the default), `:nonnegative` or `:real`. They pick the physical root among those `symbolic_solve` returns and reject impossible inputs.
- `optional`: variables that may stay undetermined without counting as missing inputs. For example, adding the dynamic viscosity μ = ρν to the pipe flow with ρ and μ optional lets you give either ν or both ρ and μ.
- `constants`: symbols with fixed values, such as π, which stay symbolic in formulas. Constants must be dimensionless; give a dimensional constant, such as g, as a variable.

## How it works

- Structure: the rank and null space of the Jacobian with respect to the unknowns tell which unknowns the given values determine, how many inputs are missing, and which equations are redundant. They are evaluated at a generic point that satisfies the equations, found by Gauss–Newton projection of a random start, because at an arbitrary point they give wrong answers.
- Solving: equations with one unknown are solved one at a time. When the remaining unknowns are coupled, an equation that is linear in one of them, with a known coefficient, expresses it through the others; otherwise a smallest block of equations is solved as a system with `symbolic_solve`. Equations of the form a uⁿ + b uᵐ = 0 are solved directly, which gives cleaner formulas. Each variable's domain picks the root.
- Consistency: the redundant equations are checked at the solution with `rtol = 1e-6`, which accepts rounded inputs.
- Formulas: the same propagation runs with the given variables kept symbolic, and at each step it keeps the root whose value matches the numeric solution.

## Limitations

- Equations must be polynomial, so roots go in as powers, and divisions by numbers give less tidy formulas (the Pitot example is written without the 1/2 for that reason).
- Formulas are correct but not always fully simplified: the pendulum's period prints as `sqrt((4L*(π^2)) / g)` rather than 2π √(L/g).
- Constants must be dimensionless.
- Loading is slow. The first precompilation of Symbolics, Nemo and Groebner takes a few minutes, and each session then spends about a minute loading and compiling.

## Tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
