# NLO DIS in dipole picture

This code implements the computation of the DIS structure functions at NLO in the dipole picture. Computation is implemented in C and C++.

**IMPORTANT NOTE:** The source code in this repository is considered *under development*, and it is not guaranteed to always work or to reproduce previously published results. For the highest reliability and accuracy in scientific computation, it is **recommended** to use the source code available on the [Releases page](https://github.com/hejajama/nlodis/releases/latest). This can also be done by selecting the latest release tag after cloning the repository:

```
git clone https://github.com/hejajama/nlodis
cd nlodis
git checkout v1.0.0
```


## Citation
When you use this code, please cite the following publications:

- [1] H. Hänninen, H. Mäntysaari and J. Penttala, "A numerical implementation of the NLO DIS structure functions in the dipole picture", Submitted to SciPost Physics Codebases, [arXiv:2604.09071](https://arxiv.org/abs/2604.09071)
- [2] H. Hänninen, H. Mäntysaari and J. Penttala, "NLO DIS in dipole picture" (software), Zenodo [doi:10.5281/zenodo.19367635](https://doi.org/10.5281/zenodo.19367635)


## Reference
This is a cleaned-up version of the [code](https://github.com/hejajama/nlobkfitter) originally written by H. Hänninen, with non-zero masses implemented by H. Hänninen and J. Penttala.

This code has been developed and used in the following publications:
* H. Hänninen, H. Mäntysaari, R. Paatelainen, J. Penttala, [Phys. Rev. Lett. 130 (2023) 19, 19](https://doi.org/10.1103/PhysRevLett.130.192301), [arXiv:2211.03504](https://arxiv.org/abs/2211.03504)
* G. Beuf, H. Hänninen, T. Lappi and H. Mäntysaari, [Phys. Rev. D102 (2020) 074028 ](), [arXiv:2007.01645](https://arxiv.org/abs/2007.01645)
* B. Ducloué, H. Hänninen, T. Lappi, Y. Zhu, [Phys. Rev.D 96 (2017) 9, 094017](https://doi.org/10.1103/PhysRevD.96.094017), [arXiv:1708.07328](https://arxiv.org/abs/1708.07328) 

The NLO DIS calculation in the dipole picture with massive quarks, implemented in this code, has been published in
* G. Beuf, T. Lappi, R. Paatelainen, [Phys.Rev.D 106 (2022) 3, 034013](https://doi.org/10.1103/PhysRevD.106.034013), [arXiv:2204.02486](https://arxiv.org/abs/2204.02486)
* G. Beuf, T. Lappi, R. Paatelainen, [Phys. Rev. Lett. 129 (2022) 7, 072001](https://doi.org/10.1103/PhysRevLett.129.072001), [arXiv:2112.03158](https://arxiv.org/abs/2112.03158)
* G. Beuf, T. Lappi, R. Paatelainen, [Phys.Rev.D 104 (2021) 5, 056032](https://doi.org/10.1103/PhysRevD.104.056032), [arXiv:2103.14549]


## Compilation
1. Dependencies: CMake and [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). The code also uses the [Cuba library](https://feynarts.de/cuba/) for multidimensional numerical integration. Cuba-4.2.2 is included in this code package and compiled automatically.
1. Do `mkdir build; cd build; cmake ..; make`


## Usage
This code requires a dipole-proton scattering amplitude (that should satisfy the NLO BK evolution equation). For example, one can use dipole amplitudes obtained in C. Casuga, H. Hänninen, H. Mäntysaari, [Phys. Rev. D112 (2025) 3, 034003](https://doi.org/10.1103/54zd-hyvg), [arXiv:2506.00487](https://arxiv.org/abs/2506.00487). Data files compatible with this code can be found in the [Zenodo repository DOI:10.5281/zenodo.15552940](https://doi.org/10.5281/zenodo.15552940)

For further instructions, see Ref. [1].



### Example program
An example program is provided and described in the published article [1], and can be executed after compilation with:

```./build/bin/nlodis```

See [./src/main.cpp](./src/main.cpp) for the source code. 

A simple example program is implemented in [./src/simple_example.cpp](./src/simple_example.cpp) and automatically compiled. It computes $F_2$ at $x=10^{-3}, Q^2=8.5$ GeV$^2$:
```bash
./build/bin/simple_example
=== NLODIS Configuration Summary ===
Order: NLO
Subtraction Scheme: Unsubtracted (UNSUB)
Nc Scheme: Finite Nc
Running Coupling Scale: Smallest dipole
Heavy quark x scheme: MassIndependentX, x = x_Bj = Q^2/(Q^2 + W^2)
IR Freezing Scheme: Smooth
Maximum dipole size (maxr): 30 GeV^-1
Running coupling C^2 factor: 1.74
Non-perturbative scale (Q0^2): 1 GeV^2
Proton transverse area: 23.3192 GeV^-2 = 9.08067 mb
Quark flavors and masses:
  light, m=0.005 GeV (e=0.816497)
  c, m=1.24 GeV (e=0.666667)
Active flavors for running coupling: determined from quark list
Dipole: Data read from file nlobkdatafiles/zenodo.15552940/balsd/bk_map.dat, minr: 1e-06 GeV^(-1), maxr: 27.5255 GeV^(-1), rpoints: 200 initial rapidity 0, maximum rapidity 15.44 Q_{s,0}^2(initial rapidity) = 0.223376 GeV^2, Q_{s,0}^2(Y=ln 1/0.01) = 0.195612 GeV^2 [ N(r^2=2/Q_s^2) = 0.393469]
Cuba integration method: vegas, maxeval 2000000, relaccuracy 0.001 cores not set (using default)
===================================

F_2(Q^2=8.5 GeV^2, xbj=0.001) = 0.875895
``` 
Using default Monte Carlo integration settings, this takes $\mathcal{O}(2)$ minutes on a typical laptop (tested  using MacBook Pro from 2022).

Before using this program, the user is encourated to verify that the current version passes unit tests by running

```./bin/nlodis_tests```

in the `build` directory. This runs unit tests implemented in [src/tests](./src/tests) folder (see different .cpp files). In order to compute $F_L$, one should change 
```c++
double F2 = dis.F2(Q2,xbj);
```
to
```c++
double F2 = dis.F2(Q2,xbj);
```
The different options that the user can control (e.g. running coupling prescription, dipole amplitude etc) are documented in [arXiv:2604.09071](https://arxiv.org/abs/2604.09071).


### Note about numerical accuracy
The aligned jet contribution gives a large logarithm $\sim \ln Q^2/m_f^2$ for the transverse cross section. This (the dipole contribution) is numerically challenging to evaluate if $Q^2/m_f^2$ is large. As such, in
the large $Q^2$ region one should always check numerical stability. This can be done, e.g., by running the code with 10 times larger `--mcintpoints` (default value is `2e6`).

### Documentation 

The API documentation can be generated as follows
```bash
cd docs
doxygen
```

## License
This project is licensed under the [MIT License](https://opensource.org/license/mit). You are free to use, modify, and distribute the code, including for commercial purposes, provided that the original copyright notice and license are included.

If you use this code in scientific work, please cite the associated publication (see above).

The unmodified version of the Cuba library, which is distributed as part of this package (`src/Cuba-4.2.2`), is available under the GNU LGPL license. See `src/Cuba-4.2.2/COPYING`.


## Future developments
* Include support for $b$-dependent dipoles
* "Undo" $z_2$ integral in $\sigma_{\mathrm{dip}}$, allowing the user to use a $z_2$-dependent evolution rapidity also in $\sigma_\mathrm{dip}$ similarly as in $\sigma_{qg}$. See [arXiv:2112.08818](https://arxiv.org/abs/2112.08818) Sec. 3.3.3 (for the $m_q=0$ case)
* Include zero quark mass limit results
