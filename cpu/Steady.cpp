#include <math.h>

extern "C" void ScaLBL_D3Q19_AAeven_SteadyPhase(double *dist, int start, int finish,
                                        int Np, double Fx, double Fy, double Fz,
                                        double *ColorGrad, int *Phi, double tau_A, double tau_B, double *Velocity) {
    double rho, jx, jy, jz;
    double m1, m2, m4, m6, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18;

    constexpr double mrt_V1 = 0.05263157894736842;
    constexpr double mrt_V2 = 0.012531328320802;
    constexpr double mrt_V3 = 0.04761904761904762;
    constexpr double mrt_V4 = 0.004594820384294068;
    constexpr double mrt_V5 = 0.01587301587301587;
    constexpr double mrt_V6 = 0.0555555555555555555555555;
    constexpr double mrt_V7 = 0.02777777777777778;
    constexpr double mrt_V8 = 0.08333333333333333;
    constexpr double mrt_V9 = 0.003341687552213868;
    constexpr double mrt_V10 = 0.003968253968253968;
    constexpr double mrt_V11 = 0.01388888888888889;
    constexpr double mrt_V12 = 0.04166666666666666;

    for (int n = start; n < finish; n++) {

        double phi = (double)Phi[n];
        double tau = (tau_A + tau_B) * 0.5 + tanh(phi * 60.0) * (tau_A - tau_B) * 0.5;
        double rlx_setA = 1.0 / tau;
        double rlx_setB = 8.f * (2.f - rlx_setA) / (8.f - rlx_setA);

        // q=0
        double fq = dist[n];
        rho = fq;
        m1 = -30.0 * fq;
        m2 = 12.0 * fq;

        // q=1
        fq = dist[2 * Np + n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jx = fq;
        m4 = -4.0 * fq;
        m9 = 2.0 * fq;
        m10 = -4.0 * fq;

        // f2 = dist[10*Np+n];
        fq = dist[1 * Np + n];
        rho += fq;
        m1 -= 11.0 * (fq);
        m2 -= 4.0 * (fq);
        jx -= fq;
        m4 += 4.0 * (fq);
        m9 += 2.0 * (fq);
        m10 -= 4.0 * (fq);

        // q=3
        fq = dist[4 * Np + n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jy = fq;
        m6 = -4.0 * fq;
        m9 -= fq;
        m10 += 2.0 * fq;
        m11 = fq;
        m12 = -2.0 * fq;

        // q = 4
        fq = dist[3 * Np + n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jy -= fq;
        m6 += 4.0 * fq;
        m9 -= fq;
        m10 += 2.0 * fq;
        m11 += fq;
        m12 -= 2.0 * fq;

        // q=5
        fq = dist[6 * Np + n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jz = fq;
        m8 = -4.0 * fq;
        m9 -= fq;
        m10 += 2.0 * fq;
        m11 -= fq;
        m12 += 2.0 * fq;

        // q = 6
        fq = dist[5 * Np + n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jz -= fq;
        m8 += 4.0 * fq;
        m9 -= fq;
        m10 += 2.0 * fq;
        m11 -= fq;
        m12 += 2.0 * fq;

        // q=7
        fq = dist[8 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx += fq;
        m4 += fq;
        jy += fq;
        m6 += fq;
        m9 += fq;
        m10 += fq;
        m11 += fq;
        m12 += fq;
        m13 = fq;
        m16 = fq;
        m17 = -fq;

        // q = 8
        fq = dist[7 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx -= fq;
        m4 -= fq;
        jy -= fq;
        m6 -= fq;
        m9 += fq;
        m10 += fq;
        m11 += fq;
        m12 += fq;
        m13 += fq;
        m16 -= fq;
        m17 += fq;

        // q=9
        fq = dist[10 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx += fq;
        m4 += fq;
        jy -= fq;
        m6 -= fq;
        m9 += fq;
        m10 += fq;
        m11 += fq;
        m12 += fq;
        m13 -= fq;
        m16 += fq;
        m17 += fq;

        // q = 10
        fq = dist[9 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx -= fq;
        m4 -= fq;
        jy += fq;
        m6 += fq;
        m9 += fq;
        m10 += fq;
        m11 += fq;
        m12 += fq;
        m13 -= fq;
        m16 -= fq;
        m17 -= fq;

        // q=11
        fq = dist[12 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx += fq;
        m4 += fq;
        jz += fq;
        m8 += fq;
        m9 += fq;
        m10 += fq;
        m11 -= fq;
        m12 -= fq;
        m15 = fq;
        m16 -= fq;
        m18 = fq;

        // q=12
        fq = dist[11 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx -= fq;
        m4 -= fq;
        jz -= fq;
        m8 -= fq;
        m9 += fq;
        m10 += fq;
        m11 -= fq;
        m12 -= fq;
        m15 += fq;
        m16 += fq;
        m18 -= fq;

        // q=13
        fq = dist[14 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx += fq;
        m4 += fq;
        jz -= fq;
        m8 -= fq;
        m9 += fq;
        m10 += fq;
        m11 -= fq;
        m12 -= fq;
        m15 -= fq;
        m16 -= fq;
        m18 -= fq;

        // q=14
        fq = dist[13 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx -= fq;
        m4 -= fq;
        jz += fq;
        m8 += fq;
        m9 += fq;
        m10 += fq;
        m11 -= fq;
        m12 -= fq;
        m15 -= fq;
        m16 += fq;
        m18 += fq;

        // q=15
        fq = dist[16 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jy += fq;
        m6 += fq;
        jz += fq;
        m8 += fq;
        m9 -= 2.0 * fq;
        m10 -= 2.0 * fq;
        m14 = fq;
        m17 += fq;
        m18 -= fq;

        // q=16
        fq = dist[15 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jy -= fq;
        m6 -= fq;
        jz -= fq;
        m8 -= fq;
        m9 -= 2.0 * fq;
        m10 -= 2.0 * fq;
        m14 += fq;
        m17 -= fq;
        m18 += fq;

        // q=17
        fq = dist[18 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jy += fq;
        m6 += fq;
        jz -= fq;
        m8 -= fq;
        m9 -= 2.0 * fq;
        m10 -= 2.0 * fq;
        m14 -= fq;
        m17 += fq;
        m18 += fq;

        // q=18
        fq = dist[17 * Np + n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jy -= fq;
        m6 -= fq;
        jz += fq;
        m8 += fq;
        m9 -= 2.0 * fq;
        m10 -= 2.0 * fq;
        m14 -= fq;
        m17 -= fq;
        m18 -= fq;

        //.............................interface force + velocity correction........................................
        double nx = ColorGrad[n];
        double ny = ColorGrad[Np + n];
        double nz = ColorGrad[2 * Np + n];

        double jx_eq = jx + 0.5 * Fx;
        double jy_eq = jy + 0.5 * Fy;
        double jz_eq = jz + 0.5 * Fz;

        double dot = (jx_eq * nx + jy_eq * ny + jz_eq * nz) / rho;

        double gx = -dot * nx;
        double gy = -dot * ny;
        double gz = -dot * nz;

        double Fx_run = Fx + gx;
        double Fy_run = Fy + gy;
        double Fz_run = Fz + gz;

        jx_eq = jx + 0.5 * Fx_run;
        jy_eq = jy + 0.5 * Fy_run;
        jz_eq = jz + 0.5 * Fz_run;
        //......................................................................     

        //........................................................................
        //					READ THE DISTRIBUTIONS
        //		(read from opposite array due to previous swap operation)
        //........................................................................

        //..............incorporate external force................................................
        //..............carry out relaxation process...............................................
        m1 = m1 +
             rlx_setA *
                 ((19 * (jx_eq * jx_eq + jy_eq * jy_eq + jz_eq * jz_eq) / rho - 11 * rho) - m1);
        m2 = m2 +
             rlx_setA *
                 ((3 * rho - 5.5 * (jx_eq * jx_eq + jy_eq * jy_eq + jz_eq * jz_eq) / rho) - m2);
        m4 = m4 + rlx_setB * ((-0.6666666666666666 * jx_eq) - m4);
        m6 = m6 + rlx_setB * ((-0.6666666666666666 * jy_eq) - m6);
        m8 = m8 + rlx_setB * ((-0.6666666666666666 * jz_eq) - m8);
        m9 = m9 + rlx_setA * (((2 * jx_eq * jx_eq - jy_eq * jy_eq - jz_eq * jz_eq) / rho) - m9);
        m10 =
            m10 +
            rlx_setA * (-0.5 * ((2 * jx_eq * jx_eq - jy_eq * jy_eq - jz_eq * jz_eq) / rho) - m10);
        m11 = m11 + rlx_setA * (((jy_eq * jy_eq - jz_eq * jz_eq) / rho) - m11);
        m12 = m12 + rlx_setA * (-0.5 * ((jy_eq * jy_eq - jz_eq * jz_eq) / rho) - m12);
        m13 = m13 + rlx_setA * ((jx_eq * jy_eq / rho) - m13);
        m14 = m14 + rlx_setA * ((jy_eq * jz_eq / rho) - m14);
        m15 = m15 + rlx_setA * ((jx_eq * jz_eq / rho) - m15);
        m16 = m16 + rlx_setB * (-m16);
        m17 = m17 + rlx_setB * (-m17);
        m18 = m18 + rlx_setB * (-m18);
        //.......................................................................................................

        //.................inverse transformation......................................................

        // q=0
        fq = mrt_V1 * rho - mrt_V2 * m1 + mrt_V3 * m2;
        dist[n] = fq;

        // q = 1
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jx - m4) +
             mrt_V6 * (m9 - m10) + 0.16666666 * Fx_run;
        dist[1 * Np + n] = fq;

        // q=2
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m4 - jx) +
             mrt_V6 * (m9 - m10) - 0.16666666 * Fx_run;
        dist[2 * Np + n] = fq;

        // q = 3
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jy - m6) +
             mrt_V7 * (m10 - m9) + mrt_V8 * (m11 - m12) + 0.16666666 * Fy_run;
        dist[3 * Np + n] = fq;

        // q = 4
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m6 - jy) +
             mrt_V7 * (m10 - m9) + mrt_V8 * (m11 - m12) - 0.16666666 * Fy_run;
        dist[4 * Np + n] = fq;

        // q = 5
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jz - m8) +
             mrt_V7 * (m10 - m9) + mrt_V8 * (m12 - m11) + 0.16666666 * Fz_run;
        dist[5 * Np + n] = fq;

        // q = 6
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m8 - jz) +
             mrt_V7 * (m10 - m9) + mrt_V8 * (m12 - m11) - 0.16666666 * Fz_run;
        dist[6 * Np + n] = fq;

        // q = 7
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx + jy) +
             0.025 * (m4 + m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
             mrt_V12 * m12 + 0.25 * m13 + 0.125 * (m16 - m17) +
             0.08333333333 * (Fx_run + Fy_run);
        dist[7 * Np + n] = fq;

        // q = 8
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jx + jy) -
             0.025 * (m4 + m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
             mrt_V12 * m12 + 0.25 * m13 + 0.125 * (m17 - m16) -
             0.08333333333 * (Fx_run + Fy_run);
        dist[8 * Np + n] = fq;

        // q = 9
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx - jy) +
             0.025 * (m4 - m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
             mrt_V12 * m12 - 0.25 * m13 + 0.125 * (m16 + m17) +
             0.08333333333 * (Fx_run - Fy_run);
        dist[9 * Np + n] = fq;

        // q = 10
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy - jx) +
             0.025 * (m6 - m4) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
             mrt_V12 * m12 - 0.25 * m13 - 0.125 * (m16 + m17) -
             0.08333333333 * (Fx_run - Fy_run);
        dist[10 * Np + n] = fq;

        // q = 11
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx + jz) +
             0.025 * (m4 + m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
             mrt_V12 * m12 + 0.25 * m15 + 0.125 * (m18 - m16) +
             0.08333333333 * (Fx_run + Fz_run);
        dist[11 * Np + n] = fq;

        // q = 12
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jx + jz) -
             0.025 * (m4 + m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
             mrt_V12 * m12 + 0.25 * m15 + 0.125 * (m16 - m18) -
             0.08333333333 * (Fx_run + Fz_run);
        dist[12 * Np + n] = fq;

        // q = 13
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx - jz) +
             0.025 * (m4 - m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
             mrt_V12 * m12 - 0.25 * m15 - 0.125 * (m16 + m18) +
             0.08333333333 * (Fx_run - Fz_run);
        dist[13 * Np + n] = fq;

        // q= 14
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jz - jx) +
             0.025 * (m8 - m4) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
             mrt_V12 * m12 - 0.25 * m15 + 0.125 * (m16 + m18) -
             0.08333333333 * (Fx_run - Fz_run);

        dist[14 * Np + n] = fq;

        // q = 15
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy + jz) +
             0.025 * (m6 + m8) - mrt_V6 * m9 - mrt_V7 * m10 + 0.25 * m14 +
             0.125 * (m17 - m18) + 0.08333333333 * (Fy_run + Fz_run);
        dist[15 * Np + n] = fq;

        // q = 16
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jy + jz) -
             0.025 * (m6 + m8) - mrt_V6 * m9 - mrt_V7 * m10 + 0.25 * m14 +
             0.125 * (m18 - m17) - 0.08333333333 * (Fy_run + Fz_run);
        dist[16 * Np + n] = fq;

        // q = 17
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy - jz) +
             0.025 * (m6 - m8) - mrt_V6 * m9 - mrt_V7 * m10 - 0.25 * m14 +
             0.125 * (m17 + m18) + 0.08333333333 * (Fy_run - Fz_run);
        dist[17 * Np + n] = fq;

        // q = 18
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jz - jy) +
             0.025 * (m8 - m6) - mrt_V6 * m9 - mrt_V7 * m10 - 0.25 * m14 -
             0.125 * (m17 + m18) - 0.08333333333 * (Fy_run - Fz_run);
        dist[18 * Np + n] = fq;

        //........................................................................
    }
}

extern "C" void ScaLBL_D3Q19_AAodd_SteadyPhase(int *neighborList, double *dist,
                                       int start, int finish, int Np,
                                       double Fx, double Fy, double Fz,
                                       double *ColorGrad, int *Phi, double tauA, double tauB, double *Velocity) {
    double rho, jx, jy, jz;
    double m1, m2, m4, m6, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18;
    constexpr double mrt_V1 = 0.05263157894736842;
    constexpr double mrt_V2 = 0.012531328320802;
    constexpr double mrt_V3 = 0.04761904761904762;
    constexpr double mrt_V4 = 0.004594820384294068;
    constexpr double mrt_V5 = 0.01587301587301587;
    constexpr double mrt_V6 = 0.0555555555555555555555555;
    constexpr double mrt_V7 = 0.02777777777777778;
    constexpr double mrt_V8 = 0.08333333333333333;
    constexpr double mrt_V9 = 0.003341687552213868;
    constexpr double mrt_V10 = 0.003968253968253968;
    constexpr double mrt_V11 = 0.01388888888888889;
    constexpr double mrt_V12 = 0.04166666666666666;

    int nread;
    for (int n = start; n < finish; n++) {
        double phi = (double)Phi[n];
        double tau = (tauA + tauB) * 0.5 + tanh(phi * 60.0) * (tauA - tauB) * 0.5;
        double rlx_setA = 1.0 / tau;
        double rlx_setB = 8.f * (2.f - rlx_setA) / (8.f - rlx_setA);
        // q=0
        double fq = dist[n];
        rho = fq;
        m1 = -30.0 * fq;
        m2 = 12.0 * fq;

        // q=1
        nread = neighborList[n];
        fq = dist[nread]; 
        //fp = dist[10*Np+n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jx = fq;
        m4 = -4.0 * fq;
        m9 = 2.0 * fq;
        m10 = -4.0 * fq;

        // f2 = dist[10*Np+n];
        nread =
            neighborList[n + Np];
        fq = dist[nread];  
        //fq = dist[Np+n];
        rho += fq;
        m1 -= 11.0 * (fq);
        m2 -= 4.0 * (fq);
        jx -= fq;
        m4 += 4.0 * (fq);
        m9 += 2.0 * (fq);
        m10 -= 4.0 * (fq);

        // q=3
        nread = neighborList[n + 2 * Np];
        fq = dist[nread];
        //fq = dist[11*Np+n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jy = fq;
        m6 = -4.0 * fq;
        m9 -= fq;
        m10 += 2.0 * fq;
        m11 = fq;
        m12 = -2.0 * fq;

        // q = 4
        nread = neighborList[n + 3 * Np]; 
        fq = dist[nread];
        //fq = dist[2*Np+n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jy -= fq;
        m6 += 4.0 * fq;
        m9 -= fq;
        m10 += 2.0 * fq;
        m11 += fq;
        m12 -= 2.0 * fq;

        // q=5
        nread = neighborList[n + 4 * Np];
        fq = dist[nread];
        //fq = dist[12*Np+n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jz = fq;
        m8 = -4.0 * fq;
        m9 -= fq;
        m10 += 2.0 * fq;
        m11 -= fq;
        m12 += 2.0 * fq;

        // q = 6
        nread = neighborList[n + 5 * Np];
        fq = dist[nread];
        //fq = dist[3*Np+n];
        rho += fq;
        m1 -= 11.0 * fq;
        m2 -= 4.0 * fq;
        jz -= fq;
        m8 += 4.0 * fq;
        m9 -= fq;
        m10 += 2.0 * fq;
        m11 -= fq;
        m12 += 2.0 * fq;

        // q=7
        nread = neighborList[n + 6 * Np];
        fq = dist[nread];
        //fq = dist[13*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx += fq;
        m4 += fq;
        jy += fq;
        m6 += fq;
        m9 += fq;
        m10 += fq;
        m11 += fq;
        m12 += fq;
        m13 = fq;
        m16 = fq;
        m17 = -fq;

        // q = 8
        nread = neighborList[n + 7 * Np];
        fq = dist[nread];
        //fq = dist[4*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx -= fq;
        m4 -= fq;
        jy -= fq;
        m6 -= fq;
        m9 += fq;
        m10 += fq;
        m11 += fq;
        m12 += fq;
        m13 += fq;
        m16 -= fq;
        m17 += fq;

        // q=9
        nread = neighborList[n + 8 * Np];
        fq = dist[nread];
        //fq = dist[14*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx += fq;
        m4 += fq;
        jy -= fq;
        m6 -= fq;
        m9 += fq;
        m10 += fq;
        m11 += fq;
        m12 += fq;
        m13 -= fq;
        m16 += fq;
        m17 += fq;

        // q = 10
        nread = neighborList[n + 9 * Np];
        fq = dist[nread];
        //fq = dist[5*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx -= fq;
        m4 -= fq;
        jy += fq;
        m6 += fq;
        m9 += fq;
        m10 += fq;
        m11 += fq;
        m12 += fq;
        m13 -= fq;
        m16 -= fq;
        m17 -= fq;

        // q=11
        nread = neighborList[n + 10 * Np];
        fq = dist[nread];
        //fq = dist[15*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx += fq;
        m4 += fq;
        jz += fq;
        m8 += fq;
        m9 += fq;
        m10 += fq;
        m11 -= fq;
        m12 -= fq;
        m15 = fq;
        m16 -= fq;
        m18 = fq;

        // q=12
        nread = neighborList[n + 11 * Np];
        fq = dist[nread];
        //fq = dist[6*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx -= fq;
        m4 -= fq;
        jz -= fq;
        m8 -= fq;
        m9 += fq;
        m10 += fq;
        m11 -= fq;
        m12 -= fq;
        m15 += fq;
        m16 += fq;
        m18 -= fq;

        // q=13
        nread = neighborList[n + 12 * Np];
        fq = dist[nread];
        //fq = dist[16*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx += fq;
        m4 += fq;
        jz -= fq;
        m8 -= fq;
        m9 += fq;
        m10 += fq;
        m11 -= fq;
        m12 -= fq;
        m15 -= fq;
        m16 -= fq;
        m18 -= fq;

        // q=14
        nread = neighborList[n + 13 * Np];
        fq = dist[nread];
        //fq = dist[7*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jx -= fq;
        m4 -= fq;
        jz += fq;
        m8 += fq;
        m9 += fq;
        m10 += fq;
        m11 -= fq;
        m12 -= fq;
        m15 -= fq;
        m16 += fq;
        m18 += fq;

        // q=15
        nread = neighborList[n + 14 * Np];
        fq = dist[nread];
        //fq = dist[17*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jy += fq;
        m6 += fq;
        jz += fq;
        m8 += fq;
        m9 -= 2.0 * fq;
        m10 -= 2.0 * fq;
        m14 = fq;
        m17 += fq;
        m18 -= fq;

        // q=16
        nread = neighborList[n + 15 * Np];
        fq = dist[nread];
        //fq = dist[8*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jy -= fq;
        m6 -= fq;
        jz -= fq;
        m8 -= fq;
        m9 -= 2.0 * fq;
        m10 -= 2.0 * fq;
        m14 += fq;
        m17 -= fq;
        m18 += fq;

        // q=17
        //fq = dist[18*Np+n];
        nread = neighborList[n + 16 * Np];
        fq = dist[nread];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jy += fq;
        m6 += fq;
        jz -= fq;
        m8 -= fq;
        m9 -= 2.0 * fq;
        m10 -= 2.0 * fq;
        m14 -= fq;
        m17 += fq;
        m18 += fq;

        // q=18
        nread = neighborList[n + 17 * Np];
        fq = dist[nread];
        //fq = dist[9*Np+n];
        rho += fq;
        m1 += 8.0 * fq;
        m2 += fq;
        jy -= fq;
        m6 -= fq;
        jz += fq;
        m8 += fq;
        m9 -= 2.0 * fq;
        m10 -= 2.0 * fq;
        m14 -= fq;
        m17 -= fq;
        m18 -= fq;

        //.............................interface force + velocity correction........................................
        double nx = ColorGrad[n];
        double ny = ColorGrad[Np + n];
        double nz = ColorGrad[2 * Np + n];

        double jx_eq = jx + 0.5 * Fx;
        double jy_eq = jy + 0.5 * Fy;
        double jz_eq = jz + 0.5 * Fz;

        double dot = (jx_eq * nx + jy_eq * ny + jz_eq * nz) / rho;

        double gx = -dot * nx;
        double gy = -dot * ny;
        double gz = -dot * nz;

        double Fx_run = Fx + gx;
        double Fy_run = Fy + gy;
        double Fz_run = Fz + gz;

        jx_eq = jx + 0.5 * Fx_run;
        jy_eq = jy + 0.5 * Fy_run;
        jz_eq = jz + 0.5 * Fz_run;
        //......................................................................     

        //..............incorporate external force................................................
        //..............carry out relaxation process...............................................
        m1 = m1 +
             rlx_setA *
                 ((19 * (jx_eq * jx_eq + jy_eq * jy_eq + jz_eq * jz_eq) / rho - 11 * rho) - m1);
        m2 = m2 +
             rlx_setA *
                 ((3 * rho - 5.5 * (jx_eq * jx_eq + jy_eq * jy_eq + jz_eq * jz_eq) / rho) - m2);
        m4 = m4 + rlx_setB * ((-0.6666666666666666 * jx_eq) - m4);
        m6 = m6 + rlx_setB * ((-0.6666666666666666 * jy_eq) - m6);
        m8 = m8 + rlx_setB * ((-0.6666666666666666 * jz_eq) - m8);
        m9 = m9 + rlx_setA * (((2 * jx_eq * jx_eq - jy_eq * jy_eq - jz_eq * jz_eq) / rho) - m9);
        m10 =
            m10 +
            rlx_setA * (-0.5 * ((2 * jx_eq * jx_eq - jy_eq * jy_eq - jz_eq * jz_eq) / rho) - m10);
        m11 = m11 + rlx_setA * (((jy_eq * jy_eq - jz_eq * jz_eq) / rho) - m11);
        m12 = m12 + rlx_setA * (-0.5 * ((jy_eq * jy_eq - jz_eq * jz_eq) / rho) - m12);
        m13 = m13 + rlx_setA * ((jx_eq * jy_eq / rho) - m13);
        m14 = m14 + rlx_setA * ((jy_eq * jz_eq / rho) - m14);
        m15 = m15 + rlx_setA * ((jx_eq * jz_eq / rho) - m15);
        m16 = m16 + rlx_setB * (-m16);
        m17 = m17 + rlx_setB * (-m17);
        m18 = m18 + rlx_setB * (-m18);
        //.......................................................................................................

        //.................inverse transformation......................................................

        // q=0
        fq = mrt_V1 * rho - mrt_V2 * m1 + mrt_V3 * m2;
        dist[n] = fq;

        // q = 1
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jx - m4) +
             mrt_V6 * (m9 - m10) + 0.16666666 * Fx_run;
        nread = neighborList[n + Np];
        dist[nread] = fq;

        // q=2
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m4 - jx) +
             mrt_V6 * (m9 - m10) - 0.16666666 * Fx_run;
        nread = neighborList[n];
        dist[nread] = fq;

        // q = 3
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jy - m6) +
             mrt_V7 * (m10 - m9) + mrt_V8 * (m11 - m12) + 0.16666666 * Fy_run;
        nread = neighborList[n + 3 * Np];
        dist[nread] = fq;

        // q = 4
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m6 - jy) +
             mrt_V7 * (m10 - m9) + mrt_V8 * (m11 - m12) - 0.16666666 * Fy_run;
        nread = neighborList[n + 2 * Np];
        dist[nread] = fq;

        // q = 5
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (jz - m8) +
             mrt_V7 * (m10 - m9) + mrt_V8 * (m12 - m11) + 0.16666666 * Fz_run;
        nread = neighborList[n + 5 * Np];
        dist[nread] = fq;

        // q = 6
        fq = mrt_V1 * rho - mrt_V4 * m1 - mrt_V5 * m2 + 0.1 * (m8 - jz) +
             mrt_V7 * (m10 - m9) + mrt_V8 * (m12 - m11) - 0.16666666 * Fz_run;
        nread = neighborList[n + 4 * Np];
        dist[nread] = fq;

        // q = 7
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx + jy) +
             0.025 * (m4 + m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
             mrt_V12 * m12 + 0.25 * m13 + 0.125 * (m16 - m17) +
             0.08333333333 * (Fx_run + Fy_run);
        nread = neighborList[n + 7 * Np];
        dist[nread] = fq;

        // q = 8
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jx + jy) -
             0.025 * (m4 + m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
             mrt_V12 * m12 + 0.25 * m13 + 0.125 * (m17 - m16) -
             0.08333333333 * (Fx_run + Fy_run);
        nread = neighborList[n + 6 * Np];
        dist[nread] = fq;

        // q = 9
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx - jy) +
             0.025 * (m4 - m6) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
             mrt_V12 * m12 - 0.25 * m13 + 0.125 * (m16 + m17) +
             0.08333333333 * (Fx_run - Fy_run);
        nread = neighborList[n + 9 * Np];
        dist[nread] = fq;

        // q = 10
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy - jx) +
             0.025 * (m6 - m4) + mrt_V7 * m9 + mrt_V11 * m10 + mrt_V8 * m11 +
             mrt_V12 * m12 - 0.25 * m13 - 0.125 * (m16 + m17) -
             0.08333333333 * (Fx_run - Fy_run);
        nread = neighborList[n + 8 * Np];
        dist[nread] = fq;

        // q = 11
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx + jz) +
             0.025 * (m4 + m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
             mrt_V12 * m12 + 0.25 * m15 + 0.125 * (m18 - m16) +
             0.08333333333 * (Fx_run + Fz_run);
        nread = neighborList[n + 11 * Np];
        dist[nread] = fq;

        // q = 12
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jx + jz) -
             0.025 * (m4 + m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
             mrt_V12 * m12 + 0.25 * m15 + 0.125 * (m16 - m18) -
             0.08333333333 * (Fx_run + Fz_run);
        nread = neighborList[n + 10 * Np];
        dist[nread] = fq;

        // q = 13
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jx - jz) +
             0.025 * (m4 - m8) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
             mrt_V12 * m12 - 0.25 * m15 - 0.125 * (m16 + m18) +
             0.08333333333 * (Fx_run - Fz_run);
        nread = neighborList[n + 13 * Np];
        dist[nread] = fq;

        // q= 14
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jz - jx) +
             0.025 * (m8 - m4) + mrt_V7 * m9 + mrt_V11 * m10 - mrt_V8 * m11 -
             mrt_V12 * m12 - 0.25 * m15 + 0.125 * (m16 + m18) -
             0.08333333333 * (Fx_run - Fz_run);
        nread = neighborList[n + 12 * Np];
        dist[nread] = fq;

        // q = 15
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy + jz) +
             0.025 * (m6 + m8) - mrt_V6 * m9 - mrt_V7 * m10 + 0.25 * m14 +
             0.125 * (m17 - m18) + 0.08333333333 * (Fy_run + Fz_run);
        nread = neighborList[n + 15 * Np];
        dist[nread] = fq;

        // q = 16
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 - 0.1 * (jy + jz) -
             0.025 * (m6 + m8) - mrt_V6 * m9 - mrt_V7 * m10 + 0.25 * m14 +
             0.125 * (m18 - m17) - 0.08333333333 * (Fy_run + Fz_run);
        nread = neighborList[n + 14 * Np];
        dist[nread] = fq;

        // q = 17
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jy - jz) +
             0.025 * (m6 - m8) - mrt_V6 * m9 - mrt_V7 * m10 - 0.25 * m14 +
             0.125 * (m17 + m18) + 0.08333333333 * (Fy_run - Fz_run);
        nread = neighborList[n + 17 * Np];
        dist[nread] = fq;

        // q = 18
        fq = mrt_V1 * rho + mrt_V9 * m1 + mrt_V10 * m2 + 0.1 * (jz - jy) +
             0.025 * (m8 - m6) - mrt_V6 * m9 - mrt_V7 * m10 - 0.25 * m14 -
             0.125 * (m17 + m18) - 0.08333333333 * (Fy_run - Fz_run);
        nread = neighborList[n + 16 * Np];
        dist[nread] = fq;
    }
}

extern "C" void SteadyComputeVelocity(double *dist, double *vel, double *pressure_out, int Np, double Fx, double Fy, double Fz, double *ColorGrad) {
    int n;
    int N = Np;
    double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9;
    double f10, f11, f12, f13, f14, f15, f16, f17, f18;
    double rho, vx, vy, vz;

    for (n = 0; n < N; n++) {
        f0 = dist[n];
        f1 = dist[N + n];
        f2 = dist[2 * N + n];
        f3 = dist[3 * N + n];
        f4 = dist[4 * N + n];
        f5 = dist[5 * N + n];
        f6 = dist[6 * N + n];
        f7 = dist[7 * N + n];
        f8 = dist[8 * N + n];
        f9 = dist[9 * N + n];
        f10 = dist[10 * N + n];
        f11 = dist[11 * N + n];
        f12 = dist[12 * N + n];
        f13 = dist[13 * N + n];
        f14 = dist[14 * N + n];
        f15 = dist[15 * N + n];
        f16 = dist[16 * N + n];
        f17 = dist[17 * N + n];
        f18 = dist[18 * N + n];

        rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
        
        vx = f1 - f2 + f7 - f8 + f9 - f10 + f11 - f12 + f13 - f14;
        vy = f3 - f4 + f7 - f8 - f9 + f10 + f15 - f16 + f17 - f18;
        vz = f5 - f6 + f11 - f12 - f13 + f14 + f15 - f16 - f17 + f18;

        double nx = ColorGrad[n];
        double ny = ColorGrad[N + n];
        double nz = ColorGrad[2 * N + n];

        double jx_star = vx + 0.5 * Fx;
        double jy_star = vy + 0.5 * Fy;
        double jz_star = vz + 0.5 * Fz;

        double dot = (jx_star * nx + jy_star * ny + jz_star * nz) / rho;

        double gx = -dot * nx;
        double gy = -dot * ny;
        double gz = -dot * nz;

        vel[n] = (jx_star + 0.5 * gx) / rho;
        vel[N + n] = (jy_star + 0.5 * gy) / rho;
        vel[2 * N + n] = (jz_star + 0.5 * gz) / rho;
        
        pressure_out[n] = rho / 3.0;
    }
}