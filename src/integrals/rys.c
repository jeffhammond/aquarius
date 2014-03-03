/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "internal.h"

/**
 * generate the roots and weights of the quadrature from the eigenvalues and first element of each eigenvector
 */
void quadrature(int n, double* restrict a, double* restrict b, double mu0, double* restrict rt, double* restrict wt)
{
    double Z[n * n] __attribute__((aligned(64)));

    int info = dstev('V', n, a, b, Z, n);
    assert(info == 0);

    for (int i = 0;i < n;i++) {
        rt[i] = a[i];
        wt[i] = Z[i * n] * Z[i * n] * mu0;
    }
}

/**
 * generate the roots and weights of the Rys quadrature
 *
 * see Golub, G. H.; Welsch, J. H. Math. Comput. 23, 221-230 (1969)
 *     K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 */
void rysquad(double T, int n, double* restrict rt, double* restrict wt)
{
    double a[n] __attribute__((aligned(64)));
    double b[n - 1] __attribute__((aligned(64)));
    double R[n+1][n+1] __attribute__((aligned(64)));
    const int n1 = n+1;

    for (int i = 0;i < n1;i++) {
        for (int j = 0;j < n1;j++) {
            R[i][j] = fm(T, i+j);
        }
    }

    for (int i = 0;i < n1;i++) {
        for (int j = i;j < n1;j++) {
            double rij = R[i][j];
            for (int k = 0;k <= i - 1;k++) {
                double rki = R[k][i];
                double rkj = R[k][j];
                double rkk = R[k][k];
                rij -= rki * rkj / rkk;
            }
            R[i][j] = rij;
        }
    }

    for (int i = 0;i < n1;i++) {
        double abs = fabs(R[i][i])
        double tmp = sqrt(1 / abs);
        for (int j = 0;j < n1;j++) {
            R[i][j] *= tmp;
        }
    }

    a[0] = R[0][1] / R[0][0];

    for (int i = 0;i < n - 1;i++) {
        double r00  = R[i][i];
        double r01  = R[i][i+1];
        double ir00 = 1/r00;
        double r11  = R[i+1][i+1];
        double r12  = R[i+1][i+2];
        double ir11 = 1/r11;
        b[i]   = r11 * ir00;
        a[i+1] = r12 * ir11 - r01 * ir00;
    }

    {
        double fm0 = fm(T, 0);
        quadrature(n, a, b, fm0, rt, wt);
    }
}

