# NLO DIS in dipole picture

This code implements the DIS  structure functions at NLO in dipole picture.

## Reference
References
* H. Hänninen, H. Mäntysaari, R. Paatelainen, J. Penttala, [Phys. Rev. Lett. 130 (2023) 19, 19](https://doi.org/10.1103/PhysRevLett.130.192301), [arXiv:2211.03504](https://arxiv.org/abs/2211.03504)

This is a cleaned-up version of the code originally written by H. Hänninen [see nlobkfitter branch](https://github.com/hejajama/nlobkfitter), non-zero masses implemented by H. Hänninen and J. Penttala.


## Compile
1. Edit CMakeLists.txt: paths for `nlodisfit_bayesian` should point to correct directories. Fetch that library from [here](https://github.com/hejajama/bayesian-nlodisfit-dipole)
2. Do `mkdir build; cd build; cmake ..; make

## Usage
See `src/main.cpp`


