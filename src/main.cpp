#include <iostream>
#include <string>
#include <vector>
#include "nlodis.hpp"
#include <gsl/gsl_errno.h>

int main(int argc, char* argv[]) {
    
    // Suppress GSL error handler for underflow errors during integration
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

    NLODIS dis(argv[1]);
    dis.SetOrder(NLO);
    double Q2=10;
    double xbj=0.1;
    cout << "#sigma(gamma+p;Q^2="<< Q2 << ",xbj="<< xbj << ",pol=L)=" << endl;
    double res = dis.Photon_proton_cross_section(Q2, xbj,L);
    cout << res << endl;

    return 0;
}