# Standard tests

This directory contains the standard tests of the functions and series
developed in this package.

The main tests are all run by `runtests.jl` including

+ `bernoulli_test.jl` -- test the Bernoulli numbers and polynomials
+ `gamma_test.jl` -- most gamma functionality is provided through SpecialFunctions, so this just tests a small number of bits
+ `harmonic_test.jl` -- harmonic and generalized harmonic series
+ `polylog_test1.jl' -- test functionality and identities for the polylogarithm function
+ `polylog_test1.jl' -- test the polylogarithm on data (not all tests will pass perfectly)
+ `stieltjes_test.jl` -- test the Stieltjes constants
+ `polylog_derivative_test.jl` -- test the derivatives of the polylog function

Many of these will use reference data from `../data/`


