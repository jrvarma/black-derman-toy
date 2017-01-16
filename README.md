Overview
========

This software uses the Black-Derman-Toy (BDT) model to value Options on Bonds (Interest Rate Options) or  bonds with embedded interest rate options (put/call options). A single factor binomial interest rate tree is built calibrated to the specified yield curve and volatility curve and this is used to value the options.

The C++ code was written in the late 1990s and has not been actively used and maintained for many years now. It has been revised to use the open source non linear optimization library [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt), and it now compiles using recent versions of `gcc` (6.3.1).

Installation
============

1. Install [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt)

2. Compile `black-derman-toy.cpp` telling the linker to use the `nlopt` library and the `C math` library (the `gcc` switches are `-lm -lnlopt -lm`). Alternatively, you can simply run `make` using the included `Makefile`.

Usage
=====

`black-derman-toy` takes one mandatory argument &ndash; the name of the input file. The two sample input files `bond-option-inputs.dat` and `embedded-option-inputs.dat` illustrate the format of these input files. See the comments in these files for an explanation. 

`bond-option-results.out` and `embedded-option-results.out` contain the outputs that would be obtained by running `black-derman-toy bond-option-inputs.dat` and `black-derman-toy embedded-option-inputs.dat` respectively.

COPYRIGHT
---------

The program is copyrighted and distributed under GNU GPL. Copyright (C) 2001 Prof. Jayanth R. Varma, jrvarma@iima.ac.in, Indian Institute of Management, Ahmedabad 380 015, INDIA. 




