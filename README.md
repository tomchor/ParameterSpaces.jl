# ParameterSpaces.jl

ParameterSpaces holds a set of variables linked by equations and computes the other variables from any sufficient set of values, dimensional or nondimensional, on either side of an equation. It reports when the given values are not enough and which extra inputs would close the system, checks values that over-specify it, and gives closed-form formulas for what it solves. It is built on [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

It was written for large-eddy simulations of submesoscale ocean fronts, whose parameter space it defines as `front_space()`, but the machinery is generic.

## Installation

The package is not registered yet. Add it from this repository (it needs Julia 1.13 or later):

```julia
using Pkg
Pkg.add(url = "https://github.com/tomchor/ParameterSpaces.jl")
```

## Quick start

Define the variables with `@variables` (re-exported from Symbolics), write the equations with `~`, and declare units over m and s. The equations must be polynomial, without division or radicals, because Symbolics' `symbolic_solve` requires it. In this example H and Lf are solved from the variables on the other side of both equations:

```julia
using ParameterSpaces

@variables f H Lf Δb Ro δ
space = ParameterSpace([f, H, Lf, Δb, Ro, δ],
                       [Ro * f^2 * Lf^2 ~ Δb * H,   # Ro = Δb H / (f² Lf²)
                        δ * Lf ~ H];                # δ = H / Lf
                       units = (f = "1/s", H = "m", Lf = "m", Δb = "m/s^2"))

sol = solve_parameters(space; f = 1e-4, Δb = 1.962e-4, Ro = 3.924, δ = 0.1)
sol[:H], sol[:Lf]            # (50.0, 500.0)
formula(sol, :H)             # (Δb*(δ^2)) / (Ro*(f^2))
latexify(formula(sol, :H))   # the same formula in LaTeX
```

## Missing and conflicting inputs

Values that leave some variables free still give everything they determine, and `closing_inputs` suggests what else to give. Values that disagree are reported with the equation they break. With the built-in front space:

```julia
front = front_space()

groups = solve_parameters(front; Ro₀ = 3.924, Roₘₑₛₒ = 0, Ri = 6.37105, δ = 0.1, R = 2.598e-3, Pr = 1)
groups.missing_inputs        # 2: the groups fix neither a time scale nor a length scale
groups[:frontal_strength]    # ≈ 39.24, which is Ro₀/δ
picks, alternatives = closing_inputs(groups)
picks                        # [:f₀, :H]
alternatives[:H]             # [:Lf, :Δb, :Qb, :U_g, :w_conv, :w_subm, :ℓ_s]
groups[:Δb]                  # ArgumentError: Δb is not determined by the given values; also giving f₀, H would determine it

solve_parameters(front; f₀ = 1e-4, H = 50, Lf = 500, Δb = 1.962e-4, Ro₀ = 4)
# ArgumentError: inconsistent input: (Lf^2)*Ro₀*(f₀^2) ~ H*Δb gives 0.01 on the left and 0.00981 on the right (relative difference 0.019, tolerance 1.0e-6)
```

## The front parameter space

`front_space()` links 23 variables with 14 equations, which leaves 9 degrees of freedom:

- Dimensional inputs: the Coriolis frequency `f₀`, the mixed-layer depth `H`, the front width `Lf`, the buoyancy difference across the front `Δb`, the stratification below the mixed layer `N²∞`, the surface buoyancy flux `Qb` (positive for cooling), the mesoscale strain rate `α`, and the eddy viscosity `ν` and diffusivity `κ`.
- Groups: the initial frontal Rossby number `Ro₀ = Δb H / (f₀² Lf²)`, the mesoscale Rossby number `Roₘₑₛₒ = α / f₀`, the balanced Richardson number `Ri = N²∞ f₀² Lf² / Δb²`, the aspect ratio `δ = H / Lf`, the restratification ratio `R = Qb f₀ Lf² / (Δb² H²)`, the Prandtl number `Pr = ν / κ` and the Ekman number `Ek = ν / (f₀ H²)`.
- Derived scales: the peak cross-front gradient `M² = Δb / Lf`, the frontal strength `frontal_strength = M² / f₀²`, the thermal-wind velocity `U_g = M² H / f₀`, the convective velocity `w_conv = ∛(Qb H)`, the submesoscale velocity `w_subm = ∛(H³ M⁴ / f₀)`, the deformation radius `ℓ_s = N∞ H / f₀` and the inertial period `T_inertial = 2π / f₀`.

Solved from the inputs of the unstrained simulations:

```julia
julia> solve_parameters(front_space(); f₀ = 1e-4, H = 50, Lf = 500, Ro₀ = 3.924, Ri = 6.37105, Roₘₑₛₒ = 0, Qb = 1e-8, Pr = 1)
ParameterSolution: 8 given, 12 solved, 3 undetermined
  given   f₀                0.0001 s⁻¹
  given   H                 50.0 m
  given   Lf                500.0 m
  solved  Δb                0.0001962 m s⁻²
  solved  N²∞               9.81e-5 s⁻²
  given   Qb                1.0e-8 m² s⁻³
  solved  α                 0.0 s⁻¹
  given   Ro₀               3.924
  given   Roₘₑₛₒ            0.0
  given   Ri                6.37105
  solved  δ                 0.1
  solved  R                 0.00259778
  given   Pr                1.0
  solved  M²                3.924e-7 s⁻²
  solved  frontal_strength  39.24
  solved  U_g               0.1962 m s⁻¹
  solved  w_conv            0.00793701 m s⁻¹
  solved  w_subm            0.0577372 m s⁻¹
  solved  ℓ_s               4952.27 m
  solved  T_inertial        62831.9 s
  undetermined: ν, κ, Ek (all optional)
```

## Defining a space

[`src/front_space.jl`](src/front_space.jl) is a complete definition. The parts of `ParameterSpace(vars, equations; units, domains, optional, constants)`:

- `vars`: the variables, in priority order. This is the order in which `closing_inputs` suggests missing inputs.
- `equations`: polynomial relations written with `~`. For example, w = ∛(Q H) is written `w^3 ~ Q * H`, and the solver keeps its single real root.
- `units`: unit strings over m and s, such as `"m/s^2"`. Variables that are not listed are dimensionless. Every equation is checked for dimensional consistency when the space is built.
- `domains`: `:positive` (the default), `:nonnegative` or `:real`. They pick the physical root among those `symbolic_solve` returns, and reject impossible inputs. A given zero can force others to zero: in the front space, Qb = 0 gives R = w_conv = 0.
- `optional`: variables that may stay undetermined without counting as missing inputs, such as the eddy viscosity of an LES closure.
- `constants`: symbols with fixed values, such as a symbolic π, which stay symbolic in formulas (`T_inertial = 2π / f₀`).

## How it works

- Structure: the rank and null space of the Jacobian with respect to the unknowns tell which unknowns the given values determine, how many inputs are missing, and which equations are redundant. They are evaluated at a generic point that satisfies the equations, found by Gauss–Newton projection of a random start, because at an arbitrary point they give wrong answers.
- Solving: equations with one unknown are solved one at a time, and coupled unknowns as a smallest block, with `symbolic_solve`. Equations of the form a uⁿ + b = 0 are solved directly, which gives cleaner formulas. Each variable's domain picks the root.
- Consistency: the redundant equations are checked at the solution with `rtol = 1e-6`, which accepts rounded inputs such as Ri = 6.37105.
- Formulas: the same propagation runs with the given variables kept symbolic, and at each step it keeps the root whose value matches the numeric solution.

## Limitations

- Equations must be polynomial, so cube roots go in as `w_conv^3 ~ Qb*H`.
- Some formulas are correct but not fully simplified: in the front space, ℓ_s comes out as `(1//2)*√(4(Lf^2)*Ri*(Ro₀^2))`, which is Lf Ro₀ √Ri.
- Loading is slow. The first precompilation of Symbolics, Nemo and Groebner takes a few minutes, and each session then spends about 50 s loading and compiling.

## Tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
