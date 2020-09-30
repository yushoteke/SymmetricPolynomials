# SymmetricPolynomials

### Installation
To install, type in the repl
```julia
] add SymmetricPolynomials
```

### Usage
```julia
using SymmetricPolynomials
```
Factor x^3+y^3+z^3
```julia
x = semi_elementary_monomial((0,0,3),(0,0,0))
decompose(x)
```

Factor w^3+x^3+y^3+z^3
```julia
x = semi_elementary_monomial((0,0,0,3),(0,0,0,0))
decompose(x)
```

Factor x^3(y+z)+y^3(x+z)+z^3(x+y)
```julia
x = semi_elementary_monomial((0,1,3),(0,0,0))
decompose(x)
```
