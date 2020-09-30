using Documenter
using Pkg
Pkg.activate("..")
#push!(LOAD_PATH,"../src/")
using SymmetricPolynomials

makedocs(sitename="SymmetricPolynomials")
