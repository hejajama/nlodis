#include <iostream>
#include <string>
#include <vector>
#include "nlodis.hpp"

int main(int argc, char* argv[]) {
    
    NLODIS dis("gbw.dat");
    cout << "sigma(gamma+p;Q^2=10,xbj=1e-3,pol=L)=" << dis.Photon_proton_cross_section(10, 1e-3,L) << endl;

    return 0;
}