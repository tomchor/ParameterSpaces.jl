"""
    ParameterSpaces

Sets of variables linked by equations. Give values for enough of the variables (dimensional or
nondimensional, on either side of an equation) and [`solve_parameters`](@ref) computes the rest,
reports which inputs are still missing, and checks over-specified input for consistency. Closed-form
expressions for the solved variables come from [`formula`](@ref).

Built on Symbolics.jl: `symbolic_solve` solves the equations (Nemo for polynomial roots, Groebner
for systems), and the Jacobian of the equations, at a generic point that satisfies them, tells
which variables a set of inputs determines.
"""
module ParameterSpaces

#+++ Imports
using Symbolics
using Nemo: Nemo          # loading it enables symbolic_solve for polynomial roots
using Groebner: Groebner  # loading it enables symbolic_solve for systems of equations
using LinearAlgebra: rank, nullspace, qr, ColumnNorm, pinv
using Random: Xoshiro
#---

#+++ Exports
export ParameterSpace, ParameterSolution
export solve_parameters, closing_inputs, formula, degrees_of_freedom
export @variables   # from Symbolics, so that scripts can define a space with `using ParameterSpaces` alone
#---

#+++ Source files
include("parameter_space.jl")
#---

end # module
