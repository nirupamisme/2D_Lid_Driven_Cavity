#include "LBM.h"
#include "IOobject.h"
#include "UserTypes.h"

using namespace std;

LBM::LBM(int r, int c, double u, double omg): rows(r), cols(c), uBC(u),
    omega_bgk(omg),
    ux(c, array1D(r, 0.0)), uy(c, array1D(r, 0.0)),
    rho(c, array1D(r, 1.0)), f00(c, array1D(r, 4.0 / 9.0)),
    fp0(c, array1D(r, 1.0 / 9.0)), fpm(c, array1D(r, 1.0 / 36.0)),
    f0m(c, array1D(r, 1.0 / 9.0)), fmm(c, array1D(r, 1.0 / 36.0)),
    fm0(c, array1D(r, 1.0 / 9.0)), fmp(c, array1D(r, 1.0 / 36.0)),
    f0p(c, array1D(r, 1.0 / 9.0)), fpp(c, array1D(r, 1.0 / 36.0)) {
    f00S = f00;
    fp0S = fp0;
    fpmS = fpm;
    f0mS = f0m;
    fmmS = fmm;
    fm0S = fm0;
    fmpS = fmp;
    f0pS = f0p;
    fppS = fpp;

    // Initalize the top boundary
    for (int i = 1; i < cols-1; i++)
        ux[i][rows-1] = uBC;
} 

LBM::~LBM() {}

// Run function for lid driven cavity (BGK-Raw)
void LBM::runBGKRaw(int timeStep) {
    double feq00, feqp0, feqpm, feq0m, feqmm, feqm0, feqmp, feq0p, feqpp;

    for (int i = 0; i <= timeStep; i++) {
        for (int j = 1; j < cols-1; j++) {
            for (int k = 1; k < rows-1; k++) {
                rho[j][k] = f00[j][k] + (((fp0[j][k] + fm0[j][k]) +
                    (f0p[j][k] + f0m[j][k])) + ((fpp[j][k] +
                    fmm[j][k]) + (fmp[j][k] + fpm[j][k])));
                ux[j][k] = ((((fp0[j][k] - fm0[j][k])) +
                    ((fpp[j][k] - fmm[j][k]) + (-fmp[j][k] +
                    fpm[j][k])))) / rho[j][k];
                uy[j][k] = (((f0p[j][k] - f0m[j][k])) +
                    ((fpp[j][k] - fmm[j][k]) + (fmp[j][k] -
                    fpm[j][k]))) / rho[j][k];

                feq00 = (rho[j][k] * (-2.0 + 3.0 * (ux[j][k] * ux[j][k])) * 
                    (-2.0 + 3.0 * (uy[j][k] * uy[j][k]))) / 9.0;
                feqm0 = -(rho[j][k] * (1.0 + 3.0 * (-1.0 + ux[j][k]) *
                    ux[j][k]) * (-2.0 + 3.0 * 
                    (uy[j][k] * uy[j][k]))) / 18.0;
                feqp0 = -(rho[j][k] * (1.0 + 3.0 * ux[j][k] *
                    (1.0 + ux[j][k])) * (-2.0 + 3.0 * (uy[j][k] *
                    uy[j][k]))) / 18.0;
                feq0m = -(rho[j][k] * (1.0 + 3.0 * (-1.0 + uy[j][k]) *
                    uy[j][k]) * (-2.0 + 3.0 *
                    (ux[j][k] * ux[j][k]))) / 18.0;
                feq0p = -(rho[j][k] * (1.0 + 3.0 * uy[j][k] *
                    (1.0 + uy[j][k])) * (-2.0 + 3.0 * (ux[j][k] *
                    ux[j][k]))) / 18.0;
                feqmm = (rho[j][k] * (1.0 + 3.0 * (-1.0 + ux[j][k]) *
                    ux[j][k]) * (1.0 + 3.0 * (-1.0 + uy[j][k]) *
                    uy[j][k])) / 36.0;
                feqpp = (rho[j][k] * (1.0 + 3.0 * ux[j][k] *
                    (1.0 + ux[j][k])) * (1.0 + 3.0 * uy[j][k] *
                    (1.0 + uy[j][k]))) / 36.0;
                feqpm = (rho[j][k] * (1.0 + 3.0 * ux[j][k] *
                    (1.0 + ux[j][k])) * (1.0 + 3.0 *
                    (-1.0 + uy[j][k]) * uy[j][k])) / 36.0;
                feqmp = (rho[j][k] * (1.0 + 3.0 *
                    (-1.0 + ux[j][k]) * ux[j][k]) * (1.0 + 3.0 *
                    uy[j][k] * (1.0 + uy[j][k]))) / 36.0;

                // Collision and streaming
                f00S[j][k] = f00[j][k] + omega_bgk * (feq00 - f00[j][k]);
                fm0S[j-1][k] = fm0[j][k] + omega_bgk * (feqm0 - fm0[j][k]);
                fp0S[j+1][k] = fp0[j][k] + omega_bgk * (feqp0 - fp0[j][k]);
                f0mS[j][k-1] = f0m[j][k] + omega_bgk * (feq0m - f0m[j][k]);
                f0pS[j][k+1] = f0p[j][k] + omega_bgk * (feq0p - f0p[j][k]);
                fmmS[j-1][k-1] = fmm[j][k] + omega_bgk * (feqmm - fmm[j][k]);
                fppS[j+1][k+1] = fpp[j][k] + omega_bgk * (feqpp - fpp[j][k]);
                fpmS[j+1][k-1] = fpm[j][k] + omega_bgk * (feqpm - fpm[j][k]);
                fmpS[j-1][k+1] = fmp[j][k] + omega_bgk * (feqmp - fmp[j][k]);
            }
        }

        for (int i = 1; i < rows-1; i++) {
            // Left
            // Delayed
            fp0S[1][i] = fm0[0][i];
            fppS[1][i+1] = fmm[0][i];
            fpmS[1][i-1] = fmp[0][i];

            // Right
            fm0S[cols-2][i] = fp0[cols-1][i];
            fmpS[cols-2][i+1] = fpm[cols-1][i];
            fmmS[cols-2][i-1] = fpp[cols-1][i];
        }

        for (int i = 1; i < cols-1; i++) {
            // Top boundary
            f0mS[i][rows-2] = f0p[i][rows-1];
            fmmS[i-1][rows-2] = fpp[i][rows-1] - 6.0 / 36.0 * uBC;
            fpmS[i+1][rows-2] = fmp[i][rows-1] + 6.0 / 36.0 * uBC;

            // Bottom boundary modified corner
            f0pS[i][1] = f0m[i][0];
            fmpS[i-1][1] = fpm[i][0];
            fppS[i+1][1] = fmm[i][0];
        }

        // Corners
        fppS[1][1] = fmm[0][0];
        fpmS[1][rows-2] = fmp[0][rows-1] - 6.0 / 36.0 * uBC;
        fmpS[cols-2][1] = fpm[cols-1][0];
        fmmS[cols-2][rows-2] = fpp[cols-1][rows-1] + 6.0 / 36.0 * uBC;

        double res = residuals(f00S, f00);

        f00S.swap(f00);
        fp0S.swap(fp0);
        fm0S.swap(fm0);
        f0pS.swap(f0p);
        f0mS.swap(f0m);
        fmmS.swap(fmm);
        fppS.swap(fpp);
        fmpS.swap(fmp);
        fpmS.swap(fpm);

        if (i % 100 == 0 && i != 0)
            cout << "Time = " << i << ", " << "Residual = " << res << endl;
    }

    // Scale the U and V velocities
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            ux[i][j] = ux[i][j] / uBC;
            uy[i][j] = uy[i][j] / uBC;
        }
    }
}

// Run function for lid driven cavity (MRT-Central)
void LBM::runMRTCentral(int timeStep, double omega_p, double omega_3, double omega_4) {
    double k00, k11, k02, k10, k01, k12, k20, k21, k22, e, p;

    for (int i = 0; i <= timeStep; i++) {
        for (int j = 1; j < cols-1; j++) {
            for (int k = 1; k < rows-1; k++) {
                rho[j][k] = f00[j][k] + (((fp0[j][k] + fm0[j][k]) +
                    (f0p[j][k] + f0m[j][k])) + ((fpp[j][k] +
                    fmm[j][k]) + (fmp[j][k] + fpm[j][k])));
                ux[j][k] = ((((fp0[j][k] - fm0[j][k])) +
                    ((fpp[j][k] - fmm[j][k]) + (-fmp[j][k] +
                    fpm[j][k])))) / rho[j][k];
                uy[j][k] = (((f0p[j][k] - f0m[j][k])) +
                    ((fpp[j][k] - fmm[j][k]) + (fmp[j][k] -
                    fpm[j][k]))) / rho[j][k];
                k00 = f00[j][k] + f0m[j][k] + f0p[j][k] + fm0[j][k] + fmm[j][k] + fmp[j][k] + fp0[j][k] + fpm[j][k] + fpp[j][k];
                k01 = -f00[j][k]*uy[j][k] + f0m[j][k]*(-uy[j][k] - 1) + f0p[j][k]*(1 - uy[j][k]) - fm0[j][k]*uy[j][k] + fmm[j][k]*(-uy[j][k] - 1) + fmp[j][k]*(1 - uy[j][k]) - fp0[j][k]*uy[j][k] + fpm[j][k]*(-uy[j][k] - 1) + fpp[j][k]*(1 - uy[j][k]);
                k02 = f00[j][k]*(uy[j][k]) * (uy[j][k]) + f0m[j][k]*(-uy[j][k] - 1) * (-uy[j][k] - 1) + f0p[j][k]*(1 - uy[j][k]) * (1 - uy[j][k]) + fm0[j][k]*(uy[j][k]) * (uy[j][k]) + fmm[j][k]*(-uy[j][k] - 1) * (-uy[j][k] - 1) + fmp[j][k]*(1 - uy[j][k]) * (1 - uy[j][k]) + fp0[j][k]*(uy[j][k]) * (uy[j][k]) + fpm[j][k]*(-uy[j][k] - 1) * (-uy[j][k] - 1) + fpp[j][k]*(1 - uy[j][k]) * (1 - uy[j][k]);
                k10 = -f00[j][k]*ux[j][k] - f0m[j][k]*ux[j][k] - f0p[j][k]*ux[j][k] + fm0[j][k]*(-ux[j][k] - 1) + fmm[j][k]*(-ux[j][k] - 1) + fmp[j][k]*(-ux[j][k] - 1) + fp0[j][k]*(1 - ux[j][k]) + fpm[j][k]*(1 - ux[j][k]) + fpp[j][k]*(1 - ux[j][k]);
                k11 = f00[j][k]*ux[j][k]*uy[j][k] - f0m[j][k]*ux[j][k]*(-uy[j][k] - 1) - f0p[j][k]*ux[j][k]*(1 - uy[j][k]) - fm0[j][k]*uy[j][k]*(-ux[j][k] - 1) + fmm[j][k]*(-ux[j][k] - 1)*(-uy[j][k] - 1) + fmp[j][k]*(1 - uy[j][k])*(-ux[j][k] - 1) - fp0[j][k]*uy[j][k]*(1 - ux[j][k]) + fpm[j][k]*(1 - ux[j][k])*(-uy[j][k] - 1) + fpp[j][k]*(1 - ux[j][k])*(1 - uy[j][k]);
                k12 = -f00[j][k]*ux[j][k]*(uy[j][k]) * (uy[j][k]) - f0m[j][k]*ux[j][k]*(-uy[j][k] - 1) * (-uy[j][k] - 1) - f0p[j][k]*ux[j][k]*(1 - uy[j][k]) * (1 - uy[j][k]) + fm0[j][k]*(uy[j][k]) * (uy[j][k])*(-ux[j][k] - 1) + fmm[j][k]*(-ux[j][k] - 1)*(-uy[j][k] - 1) * (-uy[j][k] - 1) + fmp[j][k]*(1 - uy[j][k]) * (1 - uy[j][k])*(-ux[j][k] - 1) + fp0[j][k]*(uy[j][k]) * (uy[j][k])*(1 - ux[j][k]) + fpm[j][k]*(1 - ux[j][k])*(-uy[j][k] - 1) * (-uy[j][k] - 1) + fpp[j][k]*(1 - ux[j][k])*(1 - uy[j][k]) * (1 - uy[j][k]);
                k20 = f00[j][k]*(ux[j][k]) * (ux[j][k]) + f0m[j][k]*(ux[j][k]) * (ux[j][k]) + f0p[j][k]*(ux[j][k]) * (ux[j][k]) + fm0[j][k]*(-ux[j][k] - 1) * (-ux[j][k] - 1) + fmm[j][k]*(-ux[j][k] - 1) * (-ux[j][k] - 1) + fmp[j][k]*(-ux[j][k] - 1) * (-ux[j][k] - 1) + fp0[j][k]*(1 - ux[j][k]) * (1 - ux[j][k]) + fpm[j][k]*(1 - ux[j][k]) * (1 - ux[j][k]) + fpp[j][k]*(1 - ux[j][k]) * (1 - ux[j][k]);
                k21 = -f00[j][k]*(ux[j][k]) * (ux[j][k])*uy[j][k] + f0m[j][k]*(ux[j][k]) * (ux[j][k])*(-uy[j][k] - 1) + f0p[j][k]*(ux[j][k]) * (ux[j][k])*(1 - uy[j][k]) - fm0[j][k]*uy[j][k]*(-ux[j][k] - 1) * (-ux[j][k] - 1) + fmm[j][k]*(-ux[j][k] - 1) * (-ux[j][k] - 1)*(-uy[j][k] - 1) + fmp[j][k]*(1 - uy[j][k])*(-ux[j][k] - 1) * (-ux[j][k] - 1) - fp0[j][k]*uy[j][k]*(1 - ux[j][k]) * (1 - ux[j][k]) + fpm[j][k]*(1 - ux[j][k]) * (1 - ux[j][k])*(-uy[j][k] - 1) + fpp[j][k]*(1 - ux[j][k]) * (1 - ux[j][k])*(1 - uy[j][k]);
                k22 = f00[j][k]*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) + f0m[j][k]*(ux[j][k]) * (ux[j][k])*(-uy[j][k] - 1) * (-uy[j][k] - 1) + f0p[j][k]*(ux[j][k]) * (ux[j][k])*(1 - uy[j][k]) * (1 - uy[j][k]) + fm0[j][k]*(uy[j][k]) * (uy[j][k])*(-ux[j][k] - 1) * (-ux[j][k] - 1) + fmm[j][k]*(-ux[j][k] - 1) * (-ux[j][k] - 1)*(-uy[j][k] - 1) * (-uy[j][k] - 1) + fmp[j][k]*(1 - uy[j][k]) * (1 - uy[j][k])*(-ux[j][k] - 1) * (-ux[j][k] - 1) + fp0[j][k]*(uy[j][k]) * (uy[j][k])*(1 - ux[j][k]) * (1 - ux[j][k]) + fpm[j][k]*(1 - ux[j][k]) * (1 - ux[j][k])*(-uy[j][k] - 1) * (-uy[j][k] - 1) + fpp[j][k]*(1 - ux[j][k]) * (1 - ux[j][k])*(1 - uy[j][k]) * (1 - uy[j][k]);
                e = k20 - k02;
                p = k20 + k02;

                // Collision
                k11 = k11 + omega_bgk * (-k11);
                e = e + omega_bgk * (-e);
                p = p + omega_p * (2.0 * rho[j][k] / 3.0 - p);
                k21 = k21 + omega_3 * (-k21);
                k12 = k12 + omega_3 * (-k12);
                k22 = k22 + omega_4 * (rho[j][k] / 9.0 - k22);

                k20 = (e + p) / 2.0;
                k02 = (p - e) / 2.0;

                // Transform and stream
                f00S[j][k] = k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) - k00*(ux[j][k]) * (ux[j][k]) - k00*(uy[j][k]) * (uy[j][k]) + k00 + 2*k01*(ux[j][k]) * (ux[j][k])*uy[j][k] - 2*k01*uy[j][k] + k02*(ux[j][k]) * (ux[j][k]) - k02 + 2*k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) - 2*k10*ux[j][k] + 4*k11*ux[j][k]*uy[j][k] + 2*k12*ux[j][k] + k20*(uy[j][k]) * (uy[j][k]) - k20 + 2*k21*uy[j][k] + k22;
                f0mS[j][k-1] = -1.0/2.0*k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k00*(ux[j][k]) * (ux[j][k])*uy[j][k] + (1.0/2.0)*k00*(uy[j][k]) * (uy[j][k]) - 1.0/2.0*k00*uy[j][k] - k01*(ux[j][k]) * (ux[j][k])*uy[j][k] + (1.0/2.0)*k01*(ux[j][k]) * (ux[j][k]) + k01*uy[j][k] - 1.0/2.0*k01 - 1.0/2.0*k02*(ux[j][k]) * (ux[j][k]) + (1.0/2.0)*k02 - k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) + k10*ux[j][k]*uy[j][k] - 2*k11*ux[j][k]*uy[j][k] + k11*ux[j][k] - k12*ux[j][k] - 1.0/2.0*k20*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k20*uy[j][k] - k21*uy[j][k] + (1.0/2.0)*k21 - 1.0/2.0*k22;
                f0pS[j][k+1] = -1.0/2.0*k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) - 1.0/2.0*k00*(ux[j][k]) * (ux[j][k])*uy[j][k] + (1.0/2.0)*k00*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k00*uy[j][k] - k01*(ux[j][k]) * (ux[j][k])*uy[j][k] - 1.0/2.0*k01*(ux[j][k]) * (ux[j][k]) + k01*uy[j][k] + (1.0/2.0)*k01 - 1.0/2.0*k02*(ux[j][k]) * (ux[j][k]) + (1.0/2.0)*k02 - k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) - k10*ux[j][k]*uy[j][k] - 2*k11*ux[j][k]*uy[j][k] - k11*ux[j][k] - k12*ux[j][k] - 1.0/2.0*k20*(uy[j][k]) * (uy[j][k]) - 1.0/2.0*k20*uy[j][k] - k21*uy[j][k] - 1.0/2.0*k21 - 1.0/2.0*k22;
                fm0S[j-1][k] = -1.0/2.0*k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k00*(ux[j][k]) * (ux[j][k]) + (1.0/2.0)*k00*ux[j][k]*(uy[j][k]) * (uy[j][k]) - 1.0/2.0*k00*ux[j][k] - k01*(ux[j][k]) * (ux[j][k])*uy[j][k] + k01*ux[j][k]*uy[j][k] - 1.0/2.0*k02*(ux[j][k]) * (ux[j][k]) + (1.0/2.0)*k02*ux[j][k] - k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) + k10*ux[j][k] + (1.0/2.0)*k10*(uy[j][k]) * (uy[j][k]) - 1.0/2.0*k10 - 2*k11*ux[j][k]*uy[j][k] + k11*uy[j][k] - k12*ux[j][k] + (1.0/2.0)*k12 - 1.0/2.0*k20*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k20 - k21*uy[j][k] - 1.0/2.0*k22;
                fmmS[j-1][k-1] = (1.0/4.0)*k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) - 1.0/4.0*k00*(ux[j][k]) * (ux[j][k])*uy[j][k] - 1.0/4.0*k00*ux[j][k]*(uy[j][k]) * (uy[j][k]) + (1.0/4.0)*k00*ux[j][k]*uy[j][k] + (1.0/2.0)*k01*(ux[j][k]) * (ux[j][k])*uy[j][k] - 1.0/4.0*k01*(ux[j][k]) * (ux[j][k]) - 1.0/2.0*k01*ux[j][k]*uy[j][k] + (1.0/4.0)*k01*ux[j][k] + (1.0/4.0)*k02*(ux[j][k]) * (ux[j][k]) - 1.0/4.0*k02*ux[j][k] + (1.0/2.0)*k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) - 1.0/2.0*k10*ux[j][k]*uy[j][k] - 1.0/4.0*k10*(uy[j][k]) * (uy[j][k]) + (1.0/4.0)*k10*uy[j][k] + k11*ux[j][k]*uy[j][k] - 1.0/2.0*k11*ux[j][k] - 1.0/2.0*k11*uy[j][k] + (1.0/4.0)*k11 + (1.0/2.0)*k12*ux[j][k] - 1.0/4.0*k12 + (1.0/4.0)*k20*(uy[j][k]) * (uy[j][k]) - 1.0/4.0*k20*uy[j][k] + (1.0/2.0)*k21*uy[j][k] - 1.0/4.0*k21 + (1.0/4.0)*k22;
                fmpS[j-1][k+1] = (1.0/4.0)*k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) + (1.0/4.0)*k00*(ux[j][k]) * (ux[j][k])*uy[j][k] - 1.0/4.0*k00*ux[j][k]*(uy[j][k]) * (uy[j][k]) - 1.0/4.0*k00*ux[j][k]*uy[j][k] + (1.0/2.0)*k01*(ux[j][k]) * (ux[j][k])*uy[j][k] + (1.0/4.0)*k01*(ux[j][k]) * (ux[j][k]) - 1.0/2.0*k01*ux[j][k]*uy[j][k] - 1.0/4.0*k01*ux[j][k] + (1.0/4.0)*k02*(ux[j][k]) * (ux[j][k]) - 1.0/4.0*k02*ux[j][k] + (1.0/2.0)*k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k10*ux[j][k]*uy[j][k] - 1.0/4.0*k10*(uy[j][k]) * (uy[j][k]) - 1.0/4.0*k10*uy[j][k] + k11*ux[j][k]*uy[j][k] + (1.0/2.0)*k11*ux[j][k] - 1.0/2.0*k11*uy[j][k] - 1.0/4.0*k11 + (1.0/2.0)*k12*ux[j][k] - 1.0/4.0*k12 + (1.0/4.0)*k20*(uy[j][k]) * (uy[j][k]) + (1.0/4.0)*k20*uy[j][k] + (1.0/2.0)*k21*uy[j][k] + (1.0/4.0)*k21 + (1.0/4.0)*k22;
                fp0S[j+1][k] = -1.0/2.0*k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k00*(ux[j][k]) * (ux[j][k]) - 1.0/2.0*k00*ux[j][k]*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k00*ux[j][k] - k01*(ux[j][k]) * (ux[j][k])*uy[j][k] - k01*ux[j][k]*uy[j][k] - 1.0/2.0*k02*(ux[j][k]) * (ux[j][k]) - 1.0/2.0*k02*ux[j][k] - k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) + k10*ux[j][k] - 1.0/2.0*k10*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k10 - 2*k11*ux[j][k]*uy[j][k] - k11*uy[j][k] - k12*ux[j][k] - 1.0/2.0*k12 - 1.0/2.0*k20*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k20 - k21*uy[j][k] - 1.0/2.0*k22;
                fpmS[j+1][k-1] = (1.0/4.0)*k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) - 1.0/4.0*k00*(ux[j][k]) * (ux[j][k])*uy[j][k] + (1.0/4.0)*k00*ux[j][k]*(uy[j][k]) * (uy[j][k]) - 1.0/4.0*k00*ux[j][k]*uy[j][k] + (1.0/2.0)*k01*(ux[j][k]) * (ux[j][k])*uy[j][k] - 1.0/4.0*k01*(ux[j][k]) * (ux[j][k]) + (1.0/2.0)*k01*ux[j][k]*uy[j][k] - 1.0/4.0*k01*ux[j][k] + (1.0/4.0)*k02*(ux[j][k]) * (ux[j][k]) + (1.0/4.0)*k02*ux[j][k] + (1.0/2.0)*k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) - 1.0/2.0*k10*ux[j][k]*uy[j][k] + (1.0/4.0)*k10*(uy[j][k]) * (uy[j][k]) - 1.0/4.0*k10*uy[j][k] + k11*ux[j][k]*uy[j][k] - 1.0/2.0*k11*ux[j][k] + (1.0/2.0)*k11*uy[j][k] - 1.0/4.0*k11 + (1.0/2.0)*k12*ux[j][k] + (1.0/4.0)*k12 + (1.0/4.0)*k20*(uy[j][k]) * (uy[j][k]) - 1.0/4.0*k20*uy[j][k] + (1.0/2.0)*k21*uy[j][k] - 1.0/4.0*k21 + (1.0/4.0)*k22;
                fppS[j+1][k+1] = (1.0/4.0)*k00*(ux[j][k]) * (ux[j][k])*(uy[j][k]) * (uy[j][k]) + (1.0/4.0)*k00*(ux[j][k]) * (ux[j][k])*uy[j][k] + (1.0/4.0)*k00*ux[j][k]*(uy[j][k]) * (uy[j][k]) + (1.0/4.0)*k00*ux[j][k]*uy[j][k] + (1.0/2.0)*k01*(ux[j][k]) * (ux[j][k])*uy[j][k] + (1.0/4.0)*k01*(ux[j][k]) * (ux[j][k]) + (1.0/2.0)*k01*ux[j][k]*uy[j][k] + (1.0/4.0)*k01*ux[j][k] + (1.0/4.0)*k02*(ux[j][k]) * (ux[j][k]) + (1.0/4.0)*k02*ux[j][k] + (1.0/2.0)*k10*ux[j][k]*(uy[j][k]) * (uy[j][k]) + (1.0/2.0)*k10*ux[j][k]*uy[j][k] + (1.0/4.0)*k10*(uy[j][k]) * (uy[j][k]) + (1.0/4.0)*k10*uy[j][k] + k11*ux[j][k]*uy[j][k] + (1.0/2.0)*k11*ux[j][k] + (1.0/2.0)*k11*uy[j][k] + (1.0/4.0)*k11 + (1.0/2.0)*k12*ux[j][k] + (1.0/4.0)*k12 + (1.0/4.0)*k20*(uy[j][k]) * (uy[j][k]) + (1.0/4.0)*k20*uy[j][k] + (1.0/2.0)*k21*uy[j][k] + (1.0/4.0)*k21 + (1.0/4.0)*k22;
            }
        }

        for (int i = 1; i < rows-1; i++) {
            // Left
            // Delayed
            fp0S[1][i] = fm0[0][i];
            fppS[1][i+1] = fmm[0][i];
            fpmS[1][i-1] = fmp[0][i];

            // Right
            fm0S[cols-2][i] = fp0[cols-1][i];
            fmpS[cols-2][i+1] = fpm[cols-1][i];
            fmmS[cols-2][i-1] = fpp[cols-1][i];
        }

        for (int i = 1; i < cols-1; i++) {
            // Top boundary
            f0mS[i][rows-2] = f0p[i][rows-1];
            fmmS[i-1][rows-2] = fpp[i][rows-1] - 6.0 / 36.0 * uBC;
            fpmS[i+1][rows-2] = fmp[i][rows-1] + 6.0 / 36.0 * uBC;

            // Bottom boundary modified corner
            f0pS[i][1] = f0m[i][0];
            fmpS[i-1][1] = fpm[i][0];
            fppS[i+1][1] = fmm[i][0];
        }

        // Corners
        fppS[1][1] = fmm[0][0];
        fpmS[1][rows-2] = fmp[0][rows-1] - 6.0 / 36.0 * uBC;
        fmpS[cols-2][1] = fpm[cols-1][0];
        fmmS[cols-2][rows-2] = fpp[cols-1][rows-1] + 6.0 / 36.0 * uBC;

        double res = residuals(f00S, f00);

        f00S.swap(f00);
        fp0S.swap(fp0);
        fm0S.swap(fm0);
        f0pS.swap(f0p);
        f0mS.swap(f0m);
        fmmS.swap(fmm);
        fppS.swap(fpp);
        fmpS.swap(fmp);
        fpmS.swap(fpm);

        if (i % 100 == 0 && i != 0)
            cout << "Time = " << i << ", " << "Residual = " << res << endl;
    }

    // Scale the U and V velocities
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            ux[i][j] = ux[i][j] / uBC;
            uy[i][j] = uy[i][j] / uBC;
        }
    }
}

// Run function for lid driven cavity (BGK-Central)
void LBM::runBGKCentral(int timeStep) {
    double feq00, feqp0, feqpm, feq0m, feqmm, feqm0, feqmp, feq0p, feqpp;

    for (int i = 0; i <= timeStep; i++) {
        for (int j = 1; j < cols-1; j++) {
            for (int k = 1; k < rows-1; k++) {
                rho[j][k] = f00[j][k] + (((fp0[j][k] + fm0[j][k]) +
                    (f0p[j][k] + f0m[j][k])) + ((fpp[j][k] +
                    fmm[j][k]) + (fmp[j][k] + fpm[j][k])));
                ux[j][k] = ((((fp0[j][k] - fm0[j][k])) +
                    ((fpp[j][k] - fmm[j][k]) + (-fmp[j][k] +
                    fpm[j][k])))) / rho[j][k];
                uy[j][k] = (((f0p[j][k] - f0m[j][k])) +
                    ((fpp[j][k] - fmm[j][k]) + (fmp[j][k] -
                    fpm[j][k]))) / rho[j][k];

                // Calculation of equilibrium function
                feq00 = rho[j][k] * ((ux[j][k]) * (ux[j][k]) * 
                    (uy[j][k]) * (uy[j][k]) - 0.66666666666666674 *
                    (ux[j][k]) * (ux[j][k]) - 0.66666666666666674 * 
                    (uy[j][k]) * (uy[j][k]) + 0.44444444444444453);
                feq0m = (1.0 / 2.0) * rho[j][k] * (-(ux[j][k]) * 
                    (ux[j][k]) * (uy[j][k]) * (uy[j][k]) + (ux[j][k]) * 
                    (ux[j][k]) * uy[j][k] - 0.33333333333333331 * (ux[j][k]) *
                    (ux[j][k]) + 0.66666666666666674 * (uy[j][k]) * 
                    (uy[j][k]) - 0.66666666666666674 * uy[j][k] + 
                    0.22222222222222221);
                feq0p = (1.0 / 2.0) * rho[j][k] * (-(ux[j][k]) *
                    (ux[j][k]) * (uy[j][k]) * (uy[j][k]) - (ux[j][k]) *
                    (ux[j][k]) * uy[j][k] - 0.33333333333333331 * 
                    (ux[j][k]) * (ux[j][k]) + 0.66666666666666674 * 
                    (uy[j][k]) * (uy[j][k]) + 0.66666666666666674 * 
                    uy[j][k] + 0.22222222222222221);
                feqm0 = (1.0 / 2.0) * rho[j][k] * (-(ux[j][k]) *
                    (ux[j][k]) * (uy[j][k]) * (uy[j][k]) + 
                    0.66666666666666674 * (ux[j][k]) * (ux[j][k]) + ux[j][k] *
                    (uy[j][k]) * (uy[j][k]) - 0.66666666666666674 * ux[j][k] - 
                    0.33333333333333331 * (uy[j][k]) * (uy[j][k]) +
                    0.22222222222222221);
                feqmm = (1.0 / 4.0) * rho[j][k] * ((ux[j][k]) *
                    (ux[j][k]) * (uy[j][k]) * (uy[j][k]) - (ux[j][k]) *
                    (ux[j][k]) * uy[j][k] + 0.33333333333333331 * (ux[j][k]) *
                    (ux[j][k]) - ux[j][k] * (uy[j][k]) * (uy[j][k]) + 
                    ux[j][k] * uy[j][k] - 0.33333333333333331 * ux[j][k] +
                    0.33333333333333331 * (uy[j][k]) * (uy[j][k]) - 
                    0.33333333333333331 * uy[j][k] + 0.1111111111111111);
                feqmp = (1.0 / 4.0) * rho[j][k] * ((ux[j][k]) *
                    (ux[j][k]) * (uy[j][k]) * (uy[j][k]) + (ux[j][k]) *
                    (ux[j][k]) * uy[j][k] + 0.33333333333333331 * (ux[j][k]) *
                    (ux[j][k]) - ux[j][k]*(uy[j][k]) * (uy[j][k]) - ux[j][k] *
                    uy[j][k] - 0.33333333333333331 * ux[j][k] +
                    0.33333333333333331 * (uy[j][k]) * (uy[j][k]) +
                    0.33333333333333331 * uy[j][k] + 0.1111111111111111);
                feqp0 = (1.0 / 2.0) * rho[j][k] * (-(ux[j][k]) *
                    (ux[j][k]) * (uy[j][k]) * (uy[j][k]) +
                    0.66666666666666674 * (ux[j][k]) * (ux[j][k]) - ux[j][k] *
                    (uy[j][k]) * (uy[j][k]) + 0.66666666666666674 * ux[j][k] -
                    0.33333333333333331 * (uy[j][k]) * (uy[j][k]) +
                    0.22222222222222221);
                feqpm = (1.0 / 4.0) * rho[j][k] * ((ux[j][k]) *
                    (ux[j][k]) * (uy[j][k]) * (uy[j][k]) - (ux[j][k]) *
                    (ux[j][k]) * uy[j][k] + 0.33333333333333331 * (ux[j][k]) *
                    (ux[j][k]) + ux[j][k] * (uy[j][k]) * (uy[j][k]) -
                    ux[j][k] * uy[j][k] + 0.33333333333333331 * ux[j][k] +
                    0.33333333333333331 * (uy[j][k]) * (uy[j][k]) -
                    0.33333333333333331 * uy[j][k] + 0.1111111111111111);
                feqpp = (1.0 / 4.0) * rho[j][k] * ((ux[j][k]) *
                    (ux[j][k]) * (uy[j][k]) * (uy[j][k]) + (ux[j][k]) *
                    (ux[j][k]) * uy[j][k] + 0.33333333333333331 * (ux[j][k]) *
                    (ux[j][k]) + ux[j][k] * (uy[j][k]) * (uy[j][k]) +
                    ux[j][k] * uy[j][k] + 0.33333333333333331 * ux[j][k] +
                    0.33333333333333331 * (uy[j][k]) * (uy[j][k]) +
                    0.33333333333333331 * uy[j][k] + 0.1111111111111111);

                // Collision and streaming
                f00S[j][k] = f00[j][k] + omega_bgk * (feq00 - f00[j][k]);
                fm0S[j-1][k] = fm0[j][k] + omega_bgk * (feqm0 - fm0[j][k]);
                fp0S[j+1][k] = fp0[j][k] + omega_bgk * (feqp0 - fp0[j][k]);
                f0mS[j][k-1] = f0m[j][k] + omega_bgk * (feq0m - f0m[j][k]);
                f0pS[j][k+1] = f0p[j][k] + omega_bgk * (feq0p - f0p[j][k]);
                fmmS[j-1][k-1] = fmm[j][k] + omega_bgk * (feqmm - fmm[j][k]);
                fppS[j+1][k+1] = fpp[j][k] + omega_bgk * (feqpp - fpp[j][k]);
                fpmS[j+1][k-1] = fpm[j][k] + omega_bgk * (feqpm - fpm[j][k]);
                fmpS[j-1][k+1] = fmp[j][k] + omega_bgk * (feqmp - fmp[j][k]);
            }
        }

        for (int i = 1; i < rows-1; i++) {
            // Left
            // Delayed
            fp0S[1][i] = fm0[0][i];
            fppS[1][i+1] = fmm[0][i];
            fpmS[1][i-1] = fmp[0][i];

            // Right
            fm0S[cols-2][i] = fp0[cols-1][i];
            fmpS[cols-2][i+1] = fpm[cols-1][i];
            fmmS[cols-2][i-1] = fpp[cols-1][i];
        }

        for (int i = 1; i < cols-1; i++) {
            // Top boundary
            f0mS[i][rows-2] = f0p[i][rows-1];
            fmmS[i-1][rows-2] = fpp[i][rows-1] - 6.0 / 36.0 * uBC;
            fpmS[i+1][rows-2] = fmp[i][rows-1] + 6.0 / 36.0 * uBC;

            // Bottom boundary modified corner
            f0pS[i][1] = f0m[i][0];
            fmpS[i-1][1] = fpm[i][0];
            fppS[i+1][1] = fmm[i][0];
        }

        // Corners
        fppS[1][1] = fmm[0][0];
        fpmS[1][rows-2] = fmp[0][rows-1] - 6.0 / 36.0 * uBC;
        fmpS[cols-2][1] = fpm[cols-1][0];
        fmmS[cols-2][rows-2] = fpp[cols-1][rows-1] + 6.0 / 36.0 * uBC;

        double res = residuals(f00S, f00);

        f00S.swap(f00);
        fp0S.swap(fp0);
        fm0S.swap(fm0);
        f0pS.swap(f0p);
        f0mS.swap(f0m);
        fmmS.swap(fmm);
        fppS.swap(fpp);
        fmpS.swap(fmp);
        fpmS.swap(fpm);

        if (i % 100 == 0 && i != 0)
            cout << "Time = " << i << ", " << "Residual = " << res << endl;
    }

    // Scale the U and V velocities
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            ux[i][j] = ux[i][j] / uBC;
            uy[i][j] = uy[i][j] / uBC;
        }
    }
}

// Write the U velocity
void LBM::writeUx() {
    WriteData w("./ux.csv");

    w.write(ux);
}

// Write the V velocity
void LBM::writeUy() {
    WriteData w("./uy.csv");

    w.write(uy);
}

// Residual calculation function
double LBM::residuals(array2D& f1, array2D& f2) {
    double res = 0;
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            res += pow(f1[i][j] - f2[i][j], 2);
        }
    }
    res = sqrt(res);

    return res;
}
