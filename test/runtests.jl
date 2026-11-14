# Tests of ParameterSpaces: the generic solver, on the front parameter space and on small spaces.
# Run from the repository root with
#
#     julia --project -e 'using Pkg; Pkg.test()'

#+++ Imports and helpers
using Test
using ParameterSpaces
using Symbolics

const space = front_space()
var(name) = space.vars[space.index[name]]

"The value of the expression `ex` with the named values in `point` substituted (cbrt is evaluated too)."
evaluate_at(ex, point) = real(ParameterSpaces.to_number(substitute(ex, Dict(var(n) => v for (n, v) in point))))

"A random point with values of order one for the variables named in `names`."
random_point(names) = Dict(n => exp(randn() / 2) for n in names)
#---

@testset "ParameterSpaces" begin

    #+++ Construction
    @testset "construction" begin
        @test length(space.vars) == 23
        @test length(space.equations) == 14
        @test degrees_of_freedom(space) == 9
        @test ParameterSpaces.parse_units("m^2/s^3") == (2//1, -3//1)
        @test ParameterSpaces.parse_units("1/s") == (0//1, -1//1)
        @test ParameterSpaces.parse_units("1") == (0//1, 0//1)

        @variables a b c
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b^2]; units = (a = "m", b = "m"))   # m against m²
        @test_throws ArgumentError ParameterSpace([a], [a ~ b])                                     # b is not declared
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b]; units = (c = "m",))              # c is not a variable
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b]; domains = (a = :imaginary,))
        @test_throws ArgumentError ParameterSpace([a, b], [a ~ b]; units = (a = "kg", b = "kg"))    # only m and s
    end
    #---

    #+++ Forward: the dimensional inputs of the simulations' defaults give their groups
    @testset "simulation defaults, forward" begin
        u = solve_parameters(space; f₀ = 1e-4, H = 50.0, Lf = 500.0, Δb = 1.962e-4, N²∞ = 9.81e-5,
                                    Qb = 1e-8, α = 0.0, Pr = 1.0)
        @test u[:Ro₀] ≈ 3.924
        @test u[:Ri] ≈ 6.37105 rtol = 1e-6
        @test u[:R] ≈ 2.598e-3 rtol = 1e-3
        @test u[:frontal_strength] ≈ 39.24
        @test u[:δ] ≈ 0.1
        @test u[:Roₘₑₛₒ] == 0
        @test u[:T_inertial] ≈ 2π / 1e-4
        @test u.missing_inputs == 0
        @test Set(u.undetermined) == Set([:ν, :κ, :Ek])   # optional: no value of ν or κ was given
        @test_throws ArgumentError u[:Ek]
        @test occursin("π", string(formula(u, :T_inertial)))

        s = solve_parameters(space; f₀ = 1e-4, H = 50.0, Lf = 8000.0, Δb = 3.2e-4, N²∞ = 9.81e-5,
                                    Qb = 1e-8, α = 2e-6, Pr = 1.0)
        @test s[:Ro₀] ≈ 0.025
        @test s[:Ri] ≈ 613.125
        @test s[:R] ≈ 0.25
        @test s[:frontal_strength] ≈ 4
        @test s[:Roₘₑₛₒ] ≈ 0.02
    end
    #---

    #+++ Inverse: the inversions that the simulations' setup code does by hand
    @testset "flag inversions" begin
        p = (f₀ = 1e-4, H = 50.0, Lf = 8000.0, Ro₀ = 0.025, Ri = 613.125, Roₘₑₛₒ = 0.02, Qb = 1e-8, Pr = 1.0)
        sol = solve_parameters(space; p...)
        Δb = p.Ro₀ * p.f₀^2 * p.Lf^2 / p.H
        @test sol[:Δb] ≈ Δb rtol = 1e-12
        @test sol[:N²∞] ≈ p.Ri * Δb^2 / (p.f₀ * p.Lf)^2 rtol = 1e-12
        @test sol[:α] ≈ p.Roₘₑₛₒ * p.f₀ rtol = 1e-12
        @test sol[:R] ≈ p.Qb / (p.Ro₀^2 * p.f₀^3 * p.Lf^2) rtol = 1e-12

        #++++ Closed forms, compared with the expected ones at a random point of the given variables
        x = random_point(keys(p))
        @test evaluate_at(formula(sol, :Δb), x) ≈ x[:Ro₀] * x[:f₀]^2 * x[:Lf]^2 / x[:H]
        @test evaluate_at(formula(sol, :N²∞), x) ≈ x[:Ri] * x[:Ro₀]^2 * x[:f₀]^2 * x[:Lf]^2 / x[:H]^2
        @test evaluate_at(formula(sol, :R), x) ≈ x[:Qb] / (x[:Ro₀]^2 * x[:f₀]^3 * x[:Lf]^2)
        @test evaluate_at(formula(sol, :w_conv), x) ≈ cbrt(x[:Qb] * x[:H])
        @test isequal(formula(sol, :Ro₀), var(:Ro₀))     # a given variable is its own formula
        @test formula(sol, :ν) === nothing               # undetermined
        @test occursin("frac", string(latexify(formula(sol, :Δb))))
        #----
    end
    #---

    #+++ The reference values of the Buckingham pi analysis of the front parameters, and back from the groups
    @testset "Buckingham pi reference values and round trip" begin
        inputs = (f₀ = 1e-4, H = 60.0, Lf = 1000.0, Δb = 6e-5, N²∞ = 3.6e-7, Qb = 4.79e-10, α = 1e-5,
                  ν = 1e-2, κ = 1e-2)
        groups = solve_parameters(space; inputs...)
        @test groups[:Ro₀] ≈ 0.36
        @test groups[:Ri] ≈ 1
        @test groups[:Roₘₑₛₒ] ≈ 0.1
        @test groups[:δ] ≈ 0.06
        @test groups[:R] ≈ 3.696e-3 rtol = 1e-3
        @test groups[:Pr] ≈ 1
        @test groups[:Ek] ≈ 1e-2 / (1e-4 * 60^2)
        @test isempty(groups.undetermined)

        back = solve_parameters(space; (n => groups[n] for n in (:Ro₀, :Roₘₑₛₒ, :Ri, :δ, :R, :Pr, :Ek))...,
                                       f₀ = 1e-4, H = 60.0)
        for (n, v) in pairs(inputs)
            @test back[n] ≈ v rtol = 1e-10
        end
    end
    #---

    #+++ A coupled block: H and Lf from f₀, Δb, Ro₀ and δ
    @testset "coupled block" begin
        p = (f₀ = 1e-4, Δb = 1.962e-4, Ro₀ = 3.924, δ = 0.1)
        sol = solve_parameters(space; p...)
        @test sol[:Lf] ≈ p.Δb * p.δ / (p.f₀^2 * p.Ro₀)      # 500 m
        @test sol[:H] ≈ p.δ^2 * p.Δb / (p.f₀^2 * p.Ro₀)     # 50 m
        x = random_point(keys(p))
        @test evaluate_at(formula(sol, :H), x) ≈ x[:δ]^2 * x[:Δb] / (x[:f₀]^2 * x[:Ro₀])
        @test evaluate_at(formula(sol, :Lf), x) ≈ x[:δ] * x[:Δb] / (x[:f₀]^2 * x[:Ro₀])
    end
    #---

    #+++ Too few values: what is missing, and what would close the system
    @testset "too few values" begin
        sol = solve_parameters(space; Ro₀ = 3.924, Roₘₑₛₒ = 0.0, Ri = 6.37105, δ = 0.1, R = 2.598e-3, Pr = 1.0)
        @test sol.missing_inputs == 2                        # a time and a length scale
        @test sol[:frontal_strength] ≈ 3.924 / 0.1           # = Ro₀/δ, fixed by the groups alone
        @test sol[:α] == 0                                   # = Roₘₑₛₒ f₀ for any f₀
        picks, alternatives = closing_inputs(sol)
        @test picks == [:f₀, :H]
        @test :Lf in alternatives[:H]
        @test :T_inertial in alternatives[:f₀]
        @test_throws ArgumentError sol[:Δb]

        with_f₀ = solve_parameters(space; Ro₀ = 3.924, Roₘₑₛₒ = 0.0, Ri = 6.37105, δ = 0.1, R = 2.598e-3,
                                          Pr = 1.0, f₀ = 1e-4)
        @test with_f₀.missing_inputs == 1
        @test with_f₀[:T_inertial] ≈ 2π / 1e-4
    end
    #---

    #+++ Too many values: consistent ones pass, others name the equation they break
    @testset "too many values" begin
        full = solve_parameters(space; f₀ = 1e-4, H = 50.0, Lf = 8000.0, Δb = 3.2e-4, N²∞ = 9.81e-5,
                                       Qb = 1e-8, α = 2e-6, Pr = 1.0, ν = 1e-3)
        @test length(full.values) == 23
        @test solve_parameters(space; full.values...) isa ParameterSolution

        wrong = merge(full.values, Dict(:R => 1.02 * full.values[:R]))
        @test_throws ArgumentError solve_parameters(space; wrong...)
        message = try solve_parameters(space; wrong...); "" catch e; sprint(showerror, e) end
        @test occursin("inconsistent", message) && occursin("R", message)

        # The simulations' rounded nondimensional inputs agree with their dimensional values within the default tolerance
        @test solve_parameters(space; f₀ = 1e-4, H = 50.0, Lf = 500.0, Ro₀ = 3.924, Ri = 6.37105,
                                      Δb = 1.962e-4, N²∞ = 9.81e-5, Qb = 1e-8, Pr = 1.0) isa ParameterSolution
    end
    #---

    #+++ Zeros and signs
    @testset "zeros" begin
        base = (f₀ = 1e-4, H = 50.0, Lf = 500.0, Ro₀ = 3.924, Ri = 6.37105, Pr = 1.0)
        sol = solve_parameters(space; base..., Roₘₑₛₒ = 0.0, Qb = 0.0)
        @test sol[:α] == 0
        @test sol[:R] == 0
        @test sol[:w_conv] == 0
        @test isequal(formula(sol, :R), Num(0))
        @test solve_parameters(space; base..., R = 0.0, α = 0.0)[:Qb] == 0
        @test solve_parameters(space; base..., Roₘₑₛₒ = 0.0, α = 0.0, Qb = 1e-8)[:R] ≈ 2.598e-3 rtol = 1e-3
        @test_throws ArgumentError solve_parameters(space; f₀ = 1e-4, H = 0.0)   # H must be positive
        @test_throws ArgumentError solve_parameters(space; base..., Qb = 0.0, R = 1e-3, Roₘₑₛₒ = 0.0)
    end

    @testset "signs" begin
        base = (f₀ = 1e-4, H = 50.0, Lf = 500.0, Ro₀ = 3.924, Ri = 6.37105, Roₘₑₛₒ = 0.0, Pr = 1.0)
        heating = solve_parameters(space; base..., Qb = -1e-8)
        @test heating[:R] ≈ -2.598e-3 rtol = 1e-3
        @test heating[:w_conv] ≈ cbrt(-1e-8 * 50)
        @test solve_parameters(space; base..., R = heating[:R])[:Qb] ≈ -1e-8
        @test_throws ArgumentError solve_parameters(space; H = -50.0)
        @test_throws ArgumentError solve_parameters(space; Roₘₑₛₒ = -0.02)
    end
    #---

    #+++ Relations other than power laws
    @testset "sums" begin
        @variables w_scale w_conv w_subm
        mixing = ParameterSpace([w_scale, w_conv, w_subm], [w_scale ~ w_conv + w_subm];
                                units = (w_scale = "m/s", w_conv = "m/s", w_subm = "m/s"))
        @test solve_parameters(mixing; w_scale = 0.3, w_conv = 0.1)[:w_subm] ≈ 0.2
        @test_throws ArgumentError ParameterSpace([w_scale, w_conv, w_subm], [w_scale ~ w_conv + w_subm^2];
                                                  units = (w_scale = "m/s", w_conv = "m/s", w_subm = "m/s"))
    end
    #---
end
