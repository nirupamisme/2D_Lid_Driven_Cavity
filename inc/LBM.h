#ifndef __LBM_H
#define __LBM_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "UserTypes.h"

// Class for the LBM solver
class LBM {
private:
    int rows, cols;
    double uBC, omega_bgk;
    
    // Primitive variables
    array2D ux, uy, rho;

    // Pre-collision populations
    array2D f00, fp0, fpm, f0m, fmm, fm0, fmp, f0p, fpp;

    // Post-collision populations
    array2D f00S, fp0S, fpmS, f0mS, fmmS, fm0S, fmpS, f0pS, fppS;

public:
    LBM(int, int, double, double);

    ~LBM();

    // Run function for lid driven cavity (BGK-Raw)
    void runBGKRaw(int);

    // Run function for lid driven cavity (MRT-Central)
    void runMRTCentral(int, double, double, double);

    // Run function for lid driven cavity (BGK-Central)
    void runBGKCentral(int);

    // Write the U velocity
    void writeUx();

    // Write the V velocity
    void writeUy();

    // Residual calculation function
    double residuals(array2D&, array2D&);
};

#endif      // __LBM_H
