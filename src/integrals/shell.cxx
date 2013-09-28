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

#include "input/molecule.hpp"

#include "shell.hpp"
#include "internal.h"

/**
 * Compute the index of a function in cartesian angular momentum.
 */
#define FUNC_CART(x,y,z) ((((x)*(3+(x)+2*((y)+(z))))/2) + (y))
/**
 * Compute the index of a function in spherical harmonic angular momentum.
 *
 * Regular spherical harmonics are referenced by n=l, l>=m>=-l. Contaminants may also be referenced by
 * n>l>=0, n-l even.
 */
#define FUNC_SPHER(n,l,m) ((((n)-(l))*((n)+(l)-1))/2 + 2*(n) + ((m) > 0 ? -2*(m) : 2*(m)+1))

using namespace std;
using namespace aquarius;
using namespace aquarius::integrals;
using namespace aquarius::input;
using namespace aquarius::symmetry;

Shell::Shell(const Shell& other)
: center(other.center), L(other.L), nprim(other.nprim), ncontr(other.ncontr), nfunc(other.nfunc),
  ndegen(other.ndegen), nfunc_per_irrep(other.nfunc_per_irrep), spherical(other.spherical),
  keep_contaminants(other.keep_contaminants), func_irrep(other.func_irrep),
  irrep_pos(other.irrep_pos), irreps(other.irreps), parity(other.parity), cart2spher(NULL)
{
    exponents = SAFE_MALLOC(double, nprim);
    coefficients = SAFE_MALLOC(double, nprim*ncontr);

    copy(other.exponents, other.exponents+nprim, exponents);
    copy(other.coefficients, other.coefficients+nprim*ncontr, coefficients);

    if (spherical)
    {
        cart2spher = SAFE_MALLOC(double, nfunc*(L+1)*(L+2)/2);
        copy(other.cart2spher, other.cart2spher+nfunc*(L+1)*(L+2)/2, cart2spher);
    }
}

Shell::Shell(const Shell& other, bool spherical, bool keep_contaminants)
: center(other.center), L(other.L), nprim(other.nprim), ncontr(other.ncontr), nfunc(other.nfunc),
  ndegen(other.ndegen), nfunc_per_irrep(other.nfunc_per_irrep), spherical(spherical),
  keep_contaminants(keep_contaminants), func_irrep(other.func_irrep),
  irrep_pos(other.irrep_pos), irreps(other.irreps), parity(other.parity), cart2spher(NULL)
{
    exponents = SAFE_MALLOC(double, nprim);
    coefficients = SAFE_MALLOC(double, nprim*ncontr);

    copy(other.exponents, other.exponents+nprim, exponents);
    copy(other.coefficients, other.coefficients+nprim*ncontr, coefficients);

    if (spherical)
    {
        cart2spher = SAFE_MALLOC(double, nfunc*(L+1)*(L+2)/2);
        copy(other.cart2spher, other.cart2spher+nfunc*(L+1)*(L+2)/2, cart2spher);
    }
}

Shell::Shell(const Center& pos, int L, int nprim, int ncontr, bool spherical, bool keep_contaminants,
             const double* exponents, const double* coefficients)
: center(pos), L(L), nprim(nprim), ncontr(ncontr), nfunc(spherical && !keep_contaminants ? 2*L+1 : (L+1)*(L+2)/2),
  ndegen(pos.getCenters().size()), spherical(spherical), keep_contaminants(keep_contaminants), cart2spher(NULL)
{
    this->exponents = SAFE_MALLOC(double, nprim);
    this->coefficients = SAFE_MALLOC(double, nprim*ncontr);

    copy(exponents, exponents+nprim, this->exponents);
    copy(coefficients, coefficients+nprim*ncontr, this->coefficients);

    const PointGroup& group = pos.getPointGroup();
    int nirrep = group.getNumIrreps();
    int order = group.getOrder();

    vector<int> proj(order);

    irreps.resize(nfunc, vector<int>(nirrep));
    func_irrep.resize(nfunc, vector<int>(nirrep));
    irrep_pos.resize(nfunc, vector<int>(nirrep));

    nfunc_per_irrep.resize(nirrep, 0);

    parity.resize(nfunc, vector<int>(order));

    for (int op = 0;op < order;op++)
    {
        if (spherical)
        {
            int f = 0;
            for (int l = L;l >= (keep_contaminants ? 0 : L);l -= 2)
            {
                for (int m = l;m > 0;m--)
                {
                    parity[f++][op] = (group.sphericalParity(l,  m, op) < 0 ? -1 : 1);
                    parity[f++][op] = (group.sphericalParity(l, -m, op) < 0 ? -1 : 1);
                }
                parity[f++][op] = (group.sphericalParity(l, 0, op) < 0 ? -1 : 1);
            }
        }
        else
        {
            int f = 0;
            for (int x = 0;x <= L;x++)
            {
                for (int y = 0;y <= L-x;y++)
                {
                    int z = L-x-y;
                    parity[f++][op] = (group.cartesianParity(x, y, z, op) < 0 ? -1 : 1);
                }
            }
        }
    }

    /*
     * the tricky part: determine the irrep of each final SO function
     * each AO function will give n SO functions, where n is the number of symmetry-equivalent atoms associated to this shell
     *
     * loop through the functions and for each irrep, determine the projection of this function onto the degenerate centers
     * if this is non-zero, then this is one of the irreps for this function
     *
     * this could potentially be simplified by using the DCR of the atom's stabilizer, but this works
     */
    for (int func = 0;func < nfunc;func++)
    {
        int i = 0;

        for (int irrep = 0;irrep < order;irrep++)
        {
            for (int j = 0;j < ndegen;j++)
            {
                proj[j] = 0;
            }

            /*
             * do the projection, using the characters of the irrep and the parity of either the cartesian or spherical functions
             */
            for (int op = 0;op < order;op++)
            {
                proj[pos.getCenterAfterOp(op)] += (group.character(irrep, op)*parity[func][op] < 0 ? -1 : 1);
            }

            for (int j = 0;j < ndegen;j++)
            {
                proj[j] /= (order/ndegen);
            }

            /*
             * check if the projection is non-zero
             */
            int nonzero = 0;
            for (int j = 0;j < ndegen;j++)
            {
                nonzero += abs(proj[j]);
            }
            if (nonzero > 0)
            {
                irrep_pos[func][irrep] = i;
                func_irrep[func][i] = nfunc_per_irrep[irrep];
                irreps[func][i] = irrep;
                nfunc_per_irrep[irrep]++;
                i++;
            }
        }
    }

    /*
     * Normalize the shell
     */
    const double PI2_N34 = 0.25197943553838073034791409490358;

    for (int i = 0;i < ncontr;i++)
    {
        double norm = 0.0;
        for (int j = 0;j < nprim;j++)
        {
            for (int k = 0;k < nprim;k++)
            {
                double zeta = sqrt(exponents[j]*exponents[k])/(exponents[j]+exponents[k]);
                norm += coefficients[i*nprim+j]*coefficients[i*nprim+k]*pow(2*zeta,(double)L+1.5);
            }
        }

        for (int j = 0;j < nprim;j++)
        {
            this->coefficients[i*nprim+j] *= PI2_N34*pow(4*exponents[j],((double)L+1.5)/2)/sqrt(norm);
        }
    }

    /*
     * Generate the cartesian -> spherical harmonic transformation
     */
    if (spherical)
    {
        cart2spher = SAFE_MALLOC(double, nfunc*(L+1)*(L+2)/2);
        int k = 0;
        for (int l = L;l >= (keep_contaminants ? 0 : L);l -= 2)
        {
            for (int m = l;m > 0;m--)
            {
                for (int x = 0;x <= L;x++)
                    for (int y = 0;y <= L-x;y++)
                        cart2spher[k++] = cartcoef(l, m, x, y, L-x-y);
                for (int x = 0;x <= L;x++)
                    for (int y = 0;y <= L-x;y++)
                        cart2spher[k++] = cartcoef(l, -m, x, y, L-x-y);
            }
            for (int x = 0;x <= L;x++)
                for (int y = 0;y <= L-x;y++)
                    cart2spher[k++] = cartcoef(l, 0, x, y, L-x-y);
        }
        assert(k == nfunc*(L+1)*(L+2)/2);
    }
}

Shell::~Shell()
{
    FREE(exponents);
    FREE(coefficients);
    if (cart2spher != NULL) FREE(cart2spher);
}

vector<vector<int> > Shell::setupIndices(const Context& ctx, const Molecule& m)
{
    int nirrep = m.getGroup().getNumIrreps();

    vector<vector<int> > idx;
    vector<int> nfunc(nirrep, (int)0);

    for (Molecule::const_shell_iterator s = m.getShellsBegin();s != m.getShellsEnd();++s)
    {
        idx.push_back(vector<int>(nirrep));
        vector<int>& index = idx.back();

        if (ctx.getOrdering() == Context::ISCF || ctx.getOrdering() == Context::ISFC)
        {
            for (int i = 0;i < nirrep;i++) index[i] = nfunc[i];
        }
        else
        {
            if (ctx.getOrdering() == Context::SICF || ctx.getOrdering() == Context::SIFC)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i] + s->getNFuncInIrrep(i)*s->getNContr();
                }
            }
            else if (ctx.getOrdering() == Context::SFIC ||
                     ctx.getOrdering() == Context::SFCI ||
                     ctx.getOrdering() == Context::SCFI)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i];
                }
            }
            else if (ctx.getOrdering() == Context::SCIF)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i] + s->getNFuncInIrrep(i);
                }
            }

            int nfunctot = 0;
            for (int i = 0;i < nirrep;i++) nfunctot += nfunc[i];
            for (int i = 0;i < nirrep;i++) index[i] += nfunctot;
        }

        for (int i = 0;i < nirrep;i++) nfunc[i] += s->getNFuncInIrrep(i)*s->getNContr();
    }

    for (int i = 1;i < nirrep;i++) nfunc[i] += nfunc[i-1];

    if (ctx.getOrdering() == Context::ISCF || ctx.getOrdering() == Context::ISFC)
    {
        for (vector<vector<int> >::iterator s = idx.begin();s != idx.end();++s)
        {
            vector<int>& index = *s;
            for (int i = 1;i < nirrep;i++) index[i] += nfunc[i-1];
        }
    }

    return idx;
}

int Shell::getIndex(const Context& ctx, vector<int> idx, int func, int contr, int degen) const
{
    int irrep = irreps[func][degen];

    int ffunc;
    if (spherical)
    {
        ffunc = ctx.getSphericalOrdering(L)[func];
    }
    else
    {
        ffunc = ctx.getCartesianOrdering(L)[func];
    }

    switch (ctx.getOrdering())
    {
        case Context::ISCF:
        case Context::SICF:
            return idx[irrep] + contr*nfunc_per_irrep[irrep] + func_irrep[func][degen];
        case Context::ISFC:
        case Context::SIFC:
            return idx[irrep] + contr + func_irrep[func][degen]*ncontr;
        case Context::SCIF:
            return idx[irrep] + contr*nfunc*ndegen + func_irrep[func][degen];
        case Context::SCFI:
            return idx[irrep] + (contr*nfunc + ffunc)*ndegen + irrep_pos[func][irrep];
        case Context::SFIC:
            return idx[irrep] + contr + (ffunc*ndegen + irrep_pos[func][irrep])*ncontr;
        case Context::SFCI:
            return idx[irrep] + (contr + ffunc*ncontr)*ndegen + irrep_pos[func][irrep];
    }

    return 0;
}

double Shell::cartcoef(int l, int m, int lx, int ly, int lz)
{
    int am = abs(m);
    int j = lx+ly-am;
    if ((j&1) == 1) return 0.0;
    j /= 2;

    double c = sqrt((double)(binom(2*lx,lx)*binom(2*ly,ly)*binom(2*lz,lz)*binom(l+am,am))/
                    (double)(binom(2*l,l)*binom(l,am)*binom(lx+ly+lz,lx)*binom(ly+lz,ly))/
                    (double)(dfact(2*lx-1)*dfact(2*ly-1)*dfact(2*lz-1)))/pow(2.0,l);
    if (m != 0) c *= sqrt(2.0);

    if (m >= 0)
    {
        if (((am-lx)&1) == 1) return 0.0;
        if (((am-lx)&3) == 2) c = -c;
    }
    else
    {
        if (((am-lx)&1) == 0) return 0.0;
        if (((am-lx)&3) == 3) c = -c;
    }

    double sum = 0.0;
    for (int i = 0;i <= (l-am)/2;i++)
    {
        for (int k = 0;k <= j;k++)
        {
            double tmp = binom(2*l-2*i,l+am)*binom(l,i)*binom(i,j)*binom(j,k)*binom(am,lx-2*k);
            if (((i+k)&1) == 1) tmp = -tmp;
            sum += tmp;
        }
    }

    return sum*c;
}
