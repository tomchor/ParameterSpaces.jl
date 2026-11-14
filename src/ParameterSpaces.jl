"""
    ParameterSpaces

Sets of variables linked by equations. Give values for enough of the variables (dimensional or
nondimensional, on either side of an equation) and [`solve_parameters`](@ref) computes the rest,
reports which inputs are still missing, and checks over-specified input for consistency. Closed-form
expressions for the solved variables come from [`formula`](@ref).

Built on Symbolics.jl: `symbolic_solve` solves the equations (Nemo for polynomial roots, Groebner
for systems), and the Jacobian of the equations, at a generic point that satisfies them, tells
which variables a set of inputs determines. [`front_space`](@ref) is the parameter space of the
submesoscale-front simulations the package was written for.
"""
module ParameterSpaces

#+++ Imports
using Symbolics
using Nemo: Nemo          # loading it enables symbolic_solve for polynomial roots
using Groebner: Groebner  # loading it enables symbolic_solve for systems of equations
using Latexify: latexify
using LinearAlgebra: rank, nullspace, qr, ColumnNorm, pinv
using Random: Xoshiro
#---

#+++ Exports
export ParameterSpace, ParameterSolution
export solve_parameters, closing_inputs, formula, degrees_of_freedom, latexify
export front_space
export @variables   # from Symbolics, so that scripts can define a space with `using ParameterSpaces` alone
#---

#+++ Source files
include("parameter_space.jl")  # the generic machinery
include("front_space.jl")      # the parameter space of the front simulations
#---

end # module
