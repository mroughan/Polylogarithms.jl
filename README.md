# Polylogarithms

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mroughan.github.io/Polylogarithms.jl/stable)
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://matthew.roughan@adelaide.edu.au.github.io/Polylogarithms.jl/dev) -->
[![Build Status](https://travis-ci.com/mroughan/Polylogarithms.jl.svg?branch=master)](https://travis-ci.com/mroughan/Polylogarithms.jl)
[![Coverage](https://codecov.io/gh/mroughan/Polylogarithms.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mroughan/Polylogarithms.jl)

This implements the
[Polylogarithm](https://en.wikipedia.org/wiki/Polylogarithm#Relationship_to_other_functions)
and some related functions that were needed (Harmonic numbers,
Stieltjes constants, and Bernoulli numbers and polynomials, and Euler numbers).

The code is aimed at calculating Li_s(z) for all (complex) s and z. 

This is still a little experimental, but there is a fairly large test
set that all works nicely.

Note that the aimed for accuracy is 1.0e-12 relative error, but that
occasional errors as large as 1.0e-11 have been seen. 

# Functions

 + `polylog(s, z)` the polylogarithm function
 
 + `bernoulli(n)`  Provides the first 59 Bernoulli numbers as exact rationals
 + `bernoulli(n,x)`  Provides the Bernoulli polynomials
 
 + `euler(n)`  Provides the first 61 Euler numbers as BigInt 

 + `harmonic(n)` Provides the Harmonic numbers
 + `harmonic(n,r)` Provides the generalised Harmonic numbers
 
 + `stieltjes(n)` Provides the first 10 [Stieltjes
 constants](https://en.wikipedia.org/wiki/Stieltjes_constants) (see
 Abramowitz and Stegun, 23.2.5), also known as the generalized
 Euler-Mascheroni constants.
 
 + `dirichlet_beta(z)` Provides the Dirichlet beta function
 

# Examples

```
julia> using Polylogarithms
julia> polylog(2.0, 1.0)
1.6449340668482273
```

# Details and Accuracy

Extended details of the algorithms being used here are provided at
https://arxiv.org/abs/2010.09860.

Accuracy has been tested over a wide range of scenarios. It starts to
falter for large $z$ and $\Re(s) > 8$. 

# Note

According to naming conventions, this package should have been called
Polylogarithm, but there is already an older package doing
polylogarithms, https://github.com/Expander/polylogarithm but it's
using C/CPP/Fortran bindings, and only appears to do s=2,3,4,5,6.

There is a new package
[PolyLog](https://github.com/Expander/PolyLog.jl), which does some
newer stuff for dilogarithms, so that might be for you if you care
about that case. 
