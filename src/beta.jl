"""
    dirichlet_beta()

 Calculates Dirichlet beta function,
     [https://en.wikipedia.org/wiki/Dirichlet_beta_function](https://en.wikipedia.org/wiki/Dirichlet_beta_function)
  
## Arguments
* ``s`` `::Number`: it should work for any type of number, but mainly tested for `Complex{Float64}`

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> dirichlet_beta(1.5)
0.8645026534612017
```
"""
function dirichlet_beta(s::Number)
    β = 4.0^(-s) * ( SpecialFunctions.zeta(s,0.25) - SpecialFunctions.zeta(s,0.75) )
end
