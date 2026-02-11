# NLO DIS in dipole picture

This code implements the DIS  structure functions at NLO in dipole picture.

**IMPORTANT NOTE** This package is currently under development, and it is not quaranteed to work or to reproduce previously published results.

## Reference
This is a cleaned-up version of the [code](https://github.com/hejajama/nlobkfitter) originally written by H. Hänninen, with non-zero masses implemented by H. Hänninen and J. Penttala.

This code has been developed and used in the following publications
* H. Hänninen, H. Mäntysaari, R. Paatelainen, J. Penttala, [Phys. Rev. Lett. 130 (2023) 19, 19](https://doi.org/10.1103/PhysRevLett.130.192301), [arXiv:2211.03504](https://arxiv.org/abs/2211.03504)
* G. Beuf, H. Hänninen, T. Lappi and H. Mäntysaari, [Phys. Rev. D102 (2020) 074028 ](), [arXiv:2007.01645](https://arxiv.org/abs/2007.01645)
* B. Ducloué, H. Hänninen, T. Lappi, Y. Zhu, [Phys. Rev.D 96 (2017) 9, 094017](https://doi.org/10.1103/PhysRevD.96.094017), [arXiv:1708.07328](https://arxiv.org/abs/1708.07328) 

The NLO DIS calculation in dipole picture with massive quarks, implemented in this code, has been published in
* G. Beuf, T. Lappi, R. Paatelainen, [Phys. Rev. Lett. 129 (2022) 7, 072001](https://doi.org/10.1103/PhysRevLett.129.072001), [arXIv:2112.03158](https://arxiv.org/abs/2112.03158)


## Compile
1. Dependencies: CMake and GSL. The code also uses the [Cuba library](https://feynarts.de/cuba/) for multi dimensional numericla integration, Cuba-4.2.2 is included in this code package and compiled automatically.
1. Do `mkdir build; cd build; cmake ..; make`

## Usage
This code requires dipole-proton scattering amplitude (that should satisfy the NLO BK evolution equation). For example, one can use dipole amplitudes obtained in C: Casuga, H. Hänninen, H. Mäntysaari, [Phys. Rev. D112 (2025) 3, 034003](https://doi.org/10.1103/54zd-hyvg), [arXiv:2506.00487](https://arxiv.org/abs/2506.00487). Datafiles compatible with this code can be found from the [Zenodo repository DOI:10.5281/zenodo.15552940](https://doi.org/10.5281/zenodo.15552940)

Run the example program: `./build/bin/nlodis`

See `src/main.cpp`

## License
This code is available under the TODO license.

Unmodified version of the Cuba library, that is distributed as a part of this package (`src/Cuba-4.2.2`), is available under the GNU LGPL license. See `src/Cuba-4.2.2/COPYING`



## Future developments
* Include support for $b$-dependent dipoles
* "Undo" $z_2$ integral in $\sigma_{\mathrm{dip}}$, allowing the user to use a $z_2$-dependent evolution rapidity also in $\sigma_\mathrm{dip}$ similarly as in $\sigma_{qg}$. See [arXiv:2112.08818](https://arxiv.org/abs/2112.08818) Sec. 3.3.3 (for the $m_q=0$ case)
* Include zero quark mass limit results