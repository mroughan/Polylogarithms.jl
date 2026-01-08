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
 
 + `dirichlet_beta(z)` Provides the Dirichlet beta function (the eta
    function is included in SpecialFunctions)

 + `rogers(z)` Rogers L-function

 + `spence(z)` and `dilog` -- aliases for `polylog(2, z)` 

 + `trilog(z)` and `tetralog` -- aliaes for  `polylog(3, z)` and `polylog(3, z)` 

## Recent additions

   Currently in the process of adding functions to calculate the
   derivatives of the polylogarithm WRT to both z and s. Doing this
   for a wide range of arguments requires derivatives of the Riemann
   zeta function so this is also in the process of being added.
   
## Notes

   There are few extra functions exported that are used as components
   here. These are currently being used in testing, but ultimately
   will taken out of the export list, so don't use these unless you
   are clear about what you are doing.
   

# Examples

```
julia> using Polylogarithms
julia> polylog(2.0, 1.0)
1.6449340668482273
```

# Details and Accuracy

Extended details of the algorithms being used here are provided at
https://arxiv.org/abs/2010.09860. However note that this reflected
performance of v0.1, which has been substantially improved for v0.2.

Relative accuracy has been tested to be within 1.0e-12 over a wide
range of inputs for 

  + s in the 16x16 square in the complex plane bounded by +- 8
  + z in the square in the complex plane bounded by +- 1.0e20

It fails to achieve this accuracy in <<1% of cases.

Many of the tests are in the tests directory, but the larger body of
tests, and plots visualising the results are the `extended_tests`
directory using the precomputed reference data in the `data`
directory. Running presumes you can  write files (eg plots)

# Note

According to naming conventions, this package should have been called
Polylogarithm, but there is already an older package doing
polylogarithms, https://github.com/Expander/polylogarithm but it's
using C/CPP/Fortran bindings, and only appears to do s=2,3,4,5,6.

There is a new package
[PolyLog](https://github.com/Expander/PolyLog.jl), which does some
newer stuff for dilogarithms, so that might be for you if you only
care about that case.

# Changes

Changes from v0.1 to v0.2 are in `Changes.txt`, but most of the
changes are improvements. You should see few changes to any interfaces
unless you are using the diagnostic output or using non-exported
functions.

