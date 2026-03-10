
#include "dipole/bkdipole/bkdipole.hpp"
#include "nlodis.hpp"
#include "qcd.hpp"
#include <iostream>

int main()
{
  NLODIS dis;
  // Read in dipole amplitude
  // This datafile can be downloaded from 
  // https://doi.org/10.5281/zenodo.15552940.
 dis.SetDipole(std::make_unique<BKDipole>("data/balsd/bk_map.dat"));
 // Perform NLO calculation
  dis.SetOrder(Order::NLO);

 // Set other parameters according to the applied fit, see arXiv:2506.00487 table 1

  // Running coupling scale
  dis.SetRunningCouplingC2(1.74); 
  // The distance scale is set by the smallest dipole size
  dis.SetRunningCouplingScheme(RunningCouplingScheme::SMALLEST);
  
  // Charm quark mass
  dis.SetQuarkMass(Quark::Type::C, 1.24);
  // Light quark mass is exactly 0 in arXiv:2506.00487
  // Use a finite value here for numerical stability
  dis.SetQuarkMass(Quark::Type::LIGHT, 0.005);
  
  // Proton transverse area
  dis.SetProtonTransverseArea(9.08, Unit::MB);
  
  // Print parameters to stdout
  dis.PrintConfiguration(); 

  // Compute F2
  double Q2=8.5;
  double xbj=1e-3;
  double F2 = dis.F2(Q2,xbj);

  std::cout << "F_2(Q^2=" << Q2 << " GeV^2, xbj=" << xbj << ") = " << F2 << std::endl;

  return 0;
}
