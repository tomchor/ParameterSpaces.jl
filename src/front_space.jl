# front_space.jl
#
# The parameter space of large-eddy simulations of a submesoscale ocean front, either strained by
# a mesoscale flow or released after a nudged spin-up: nine dimensional inputs, the seven
# nondimensional groups of their Buckingham pi analysis, and derived scales. H is the depth of the
# mixed layer, |hᵢ| in the simulations' notation, where the mixed-layer base hᵢ is negative.

"""
    front_space()

The parameter space of the front simulations: 23 variables linked by 14 equations, with 9 degrees of
freedom. The equations are written in polynomial form, as `symbolic_solve` requires; their usual
form is in the comments of `front_space.jl`.

Dimensional inputs, in SI units: Coriolis frequency `f₀`, mixed-layer depth `H` (= |hᵢ|), front
width `Lf`, buoyancy difference across the front `Δb`, pycnocline stratification `N²∞`, surface
buoyancy flux `Qb` (positive for cooling), mesoscale strain rate `α`, and the eddy viscosity `ν` and
diffusivity `κ` (optional: they may stay undetermined).

Groups: initial frontal Rossby number `Ro₀`, mesoscale Rossby number `Roₘₑₛₒ`, balanced Richardson
number `Ri`, aspect ratio `δ`, restratification ratio `R`, Prandtl number `Pr` and Ekman number `Ek`
(optional).

Derived scales: peak cross-front gradient `M²`, frontal strength `frontal_strength` (M²/f₀²),
thermal-wind velocity `U_g`, convective velocity `w_conv`, submesoscale velocity `w_subm`,
deformation radius `ℓ_s` and inertial period `T_inertial`.
"""
function front_space()
    #+++ Variables, in priority order: dimensional inputs, groups, derived scales
    @variables f₀ H Lf Δb N²∞ Qb α ν κ
    @variables Ro₀ Roₘₑₛₒ Ri δ R Pr Ek
    @variables M² frontal_strength U_g w_conv w_subm ℓ_s T_inertial
    π_ = Symbolics.variable(:π)   # π, kept symbolic in formulas
    #---

    #+++ Equations, in polynomial form
    equations = [
        Ro₀ * f₀^2 * Lf^2 ~ Δb * H,           # Ro₀ = Δb H / (f₀² Lf²)      initial frontal Rossby number
        Roₘₑₛₒ * f₀ ~ α,                       # Roₘₑₛₒ = α / f₀             mesoscale Rossby number
        Ri * Δb^2 ~ N²∞ * f₀^2 * Lf^2,         # Ri = N²∞ f₀² Lf² / Δb²      balanced Richardson number
        δ * Lf ~ H,                            # δ = H / Lf                  aspect ratio of the front
        R * Δb^2 * H^2 ~ Qb * f₀ * Lf^2,       # R = Qb f₀ Lf² / (Δb² H²)    restratification ratio
        Pr * κ ~ ν,                            # Pr = ν / κ                  Prandtl number
        Ek * f₀ * H^2 ~ ν,                     # Ek = ν / (f₀ H²)            Ekman number
        M² * Lf ~ Δb,                          # M² = Δb / Lf                peak cross-front gradient
        frontal_strength * f₀^2 ~ M²,          # frontal_strength = M² / f₀²
        U_g * f₀ ~ M² * H,                     # U_g = M² H / f₀             thermal-wind velocity
        w_conv^3 ~ Qb * H,                     # w_conv = ∛(Qb H)            convective velocity
        w_subm^3 * f₀ ~ H^3 * M²^2,            # w_subm = ∛(H³ M⁴ / f₀)      submesoscale velocity
        ℓ_s^2 * f₀^2 ~ N²∞ * H^2,              # ℓ_s = N∞ H / f₀             deformation radius
        T_inertial * f₀ ~ 2π_,                 # T_inertial = 2π / f₀        inertial period
    ]
    #---

    #+++ The space
    return ParameterSpace([f₀, H, Lf, Δb, N²∞, Qb, α, ν, κ, Ro₀, Roₘₑₛₒ, Ri, δ, R, Pr, Ek,
                           M², frontal_strength, U_g, w_conv, w_subm, ℓ_s, T_inertial], equations;
                          units = (f₀ = "1/s", H = "m", Lf = "m", Δb = "m/s^2", N²∞ = "1/s^2",
                                   Qb = "m^2/s^3", α = "1/s", ν = "m^2/s", κ = "m^2/s", M² = "1/s^2",
                                   U_g = "m/s", w_conv = "m/s", w_subm = "m/s", ℓ_s = "m", T_inertial = "s"),
                          domains = (N²∞ = :nonnegative, Qb = :real, α = :nonnegative, Roₘₑₛₒ = :nonnegative,
                                     Ri = :nonnegative, R = :real, w_conv = :real, ℓ_s = :nonnegative),
                          optional = (:ν, :κ, :Ek),
                          constants = Dict(π_ => π))
    #---
end
