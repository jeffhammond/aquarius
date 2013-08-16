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

#include "limits.h"

#include "aomoints.hpp"

#include "util/util.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::scf;
using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::slide;

/*
static bool sortIntsByRS(const idx4_t& i1, const idx4_t& i2)
{
    if (i1.l < i2.l)
    {
        return true;
    }
    else if (i1.l > i2.l)
    {
        return false;
    }
    if (i1.k < i2.k)
    {
        return true;
    }
    else
    {
        return false;
    }
}
*/

template <typename T>
AOMOIntegrals<T>::AOMOIntegrals(const AOUHF<T>& uhf)
: MOIntegrals<T>(uhf)
{
    doTransformation(uhf.ints);
}

template <typename T>
AOMOIntegrals<T>::pqrs_integrals::pqrs_integrals(const AOIntegrals<T>& aoints)
: Distributed<T>(aoints.ctf)
{
    ns = nr = nq = np = aoints.molecule.getNumOrbitals();

    size_t noldints = aoints.getNumInts();
    const T *oldints = aoints.getInts();
    const idx4_t *oldidxs = aoints.getIndices();
    nints = noldints;

    for (size_t i = 0;i < noldints;i++)
    {
        idx4_t idx = oldidxs[i];

        if (!((idx.i == idx.k && idx.j == idx.l) ||
              (idx.i == idx.l && idx.j == idx.k)))
        {
            nints++;
        }
    }

    ints = SAFE_MALLOC(T, nints);
    idxs = SAFE_MALLOC(idx4_t, nints);

    size_t j = 0;
    for (size_t i = 0;i < noldints;i++)
    {
        T val = oldints[i];
        idx4_t idx = oldidxs[i];

        if (idx.i > idx.j) swap(idx.i, idx.j);
        if (idx.k > idx.l) swap(idx.k, idx.l);

        assert(idx.i >= 0 && idx.i < np);
        assert(idx.j >= 0 && idx.j < nq);
        assert(idx.k >= 0 && idx.k < nr);
        assert(idx.l >= 0 && idx.l < ns);

        assert(j < nints);
        ints[j] = val;
        idxs[j] = idx;
        j++;

        if (idx.i != idx.k || idx.j != idx.l)
        {
            swap(idx.i, idx.k);
            swap(idx.j, idx.l);
            assert(j < nints);
            ints[j] = val;
            idxs[j] = idx;
            j++;
        }
    }
    assert(j == nints);
}

template <typename T>
AOMOIntegrals<T>::pqrs_integrals::pqrs_integrals(abrs_integrals& abrs)
: Distributed<T>(abrs.ctf)
{
    np = abrs.nr;
    nq = abrs.ns;
    nr = abrs.na;
    ns = abrs.nb;

    nints = abrs.nints;
    ints = abrs.ints;
    idxs = SAFE_MALLOC(idx4_t, nints);
    size_t ipqrs, irs;
    for (ipqrs = 0, irs = 0;irs < abrs.nrs;irs++)
    {
        int p = abrs.rs[irs].i;
        int q = abrs.rs[irs].j;
        assert(p >= 0 && p < np);
        assert(q >= 0 && q < nq);

        for (int s = 0;s < ns;s++)
        {
            for (int r = 0;r < nr;r++)
            {
                assert(ipqrs < nints);
                idxs[ipqrs].i = p;
                idxs[ipqrs].j = q;
                idxs[ipqrs].k = r;
                idxs[ipqrs].l = s;
                ipqrs++;
            }
        }
    }
    assert(ipqrs == nints);

    FREE(abrs.rs);
}

template <typename T>
void AOMOIntegrals<T>::pqrs_integrals::free()
{
    FREE(ints);
    FREE(idxs);
}

template <typename T>
void AOMOIntegrals<T>::pqrs_integrals::sortInts(bool rles, size_t& nrs, size_t*& rscount)
{
    if (rles)
    {
        nrs = nr*(nr+1)/2;
    }
    else
    {
        nrs = nr*ns;
    }

    rscount = SAFE_MALLOC(size_t, nrs);
    fill(rscount, rscount+nrs, 0);

    for (size_t i = 0;i < nints;i++)
    {
        if (rles)
        {
            assert(idxs[i].k+idxs[i].l*(idxs[i].l+1)/2 < nrs);
            rscount[idxs[i].k+idxs[i].l*(idxs[i].l+1)/2]++;
        }
        else
        {
            assert(idxs[i].k+idxs[i].l*nr < nrs);
            rscount[idxs[i].k+idxs[i].l*nr]++;
        }
    }

    size_t *rsoff = SAFE_MALLOC(size_t, nrs);

    rsoff[0] = 0;
    for (size_t i = 1;i < nrs;i++) rsoff[i] = rsoff[i-1]+rscount[i-1];
    assert(rsoff[nrs-1]+rscount[nrs-1] == nints);

    T* newints = SAFE_MALLOC(T, nints);
    idx4_t* newidxs = SAFE_MALLOC(idx4_t, nints);

    for (size_t i = 0;i < nints;i++)
    {
        size_t rs = (rles ? idxs[i].k+idxs[i].l*(idxs[i].l+1)/2
                          : idxs[i].k+idxs[i].l*nr);

        newints[rsoff[rs]] = ints[i];
        newidxs[rsoff[rs]] = idxs[i];

        rsoff[rs]++;
    }

    swap(ints, newints);
    swap(idxs, newidxs);
    FREE(newints);
    FREE(newidxs);
    FREE(rsoff);
}

/*
 * Redistribute integrals such that each node has all pq for each rs pair
 */
template <typename T>
void AOMOIntegrals<T>::pqrs_integrals::collect(bool rles)
{
    assert(sizeof(short) == sizeof(int16_t));
    MPI::Datatype IDX4_T_TYPE = MPI::Datatype(MPI_SHORT).Create_contiguous(4);
    IDX4_T_TYPE.Commit();

    size_t nrs;
    size_t *rscount;
    sortInts(rles, nrs, rscount);

    size_t nnewints = 0;
    size_t *sendcount = SAFE_MALLOC(size_t, nproc);
    size_t *sendoff = SAFE_MALLOC(size_t, nproc);
    size_t *recvcount = SAFE_MALLOC(size_t, nproc);
    size_t *recvoff = SAFE_MALLOC(size_t, nproc);
    fill(sendcount, sendcount+nproc, 0);
    sendoff[0] = 0;
    fill(recvcount, recvcount+nproc, 0);
    recvoff[0] = 0;

    for (int i = 0;i < nproc;i++)
    {
        for (size_t rs = (nrs*i)/nproc;rs < (nrs*(i+1))/nproc;rs++)
        {
            assert(rs >= 0 && rs < nrs);
            sendcount[i] += rscount[rs];
        }
        if (i > 0) sendoff[i] = sendoff[i-1]+sendcount[i-1];
    }

    FREE(rscount);
    this->comm.Alltoall(sendcount, 1, MPI_TYPE_<size_t>::value(),
                        recvcount, 1, MPI_TYPE_<size_t>::value());

    for (int i = 1;i < nproc;i++)
    {
        recvoff[i] = recvoff[i-1]+recvcount[i-1];
    }
    nnewints = recvoff[nproc-1]+recvcount[nproc-1];

    assert(recvoff[nproc-1]+recvcount[nproc-1] == nnewints);
    assert(sendoff[nproc-1]+sendcount[nproc-1] == nints);

    assert(allsum(nints) == allsum(nnewints));

    T* newints = SAFE_MALLOC(T, nnewints);
    idx4_t* newidxs = SAFE_MALLOC(idx4_t, nnewints);

    int *realsendcount = SAFE_MALLOC(int, nproc);
    int *realrecvcount = SAFE_MALLOC(int, nproc);
    int *realsendoff = SAFE_MALLOC(int, nproc);
    int *realrecvoff = SAFE_MALLOC(int, nproc);

    for (int i = 0;i < nproc;i++)
    {
        assert(sendcount[i] <= INT_MAX);
        assert(recvcount[i] <= INT_MAX);
        assert(sendoff[i] <= INT_MAX);
        assert(recvoff[i] <= INT_MAX);
        realsendcount[i] = (int)sendcount[i];
        realrecvcount[i] = (int)recvcount[i];
        realsendoff[i] = (int)sendoff[i];
        realrecvoff[i] = (int)recvoff[i];
    }

    this->comm.Alltoallv(   ints, realsendcount, realsendoff, MPI_TYPE_<T>::value(),
                         newints, realrecvcount, realrecvoff, MPI_TYPE_<T>::value());

    this->comm.Alltoallv(   idxs, realsendcount, realsendoff, IDX4_T_TYPE,
                         newidxs, realrecvcount, realrecvoff, IDX4_T_TYPE);

    FREE(realsendcount);
    FREE(realsendoff);
    FREE(realrecvcount);
    FREE(realrecvoff);
    FREE(sendcount);
    FREE(sendoff);
    FREE(recvcount);
    FREE(recvoff);

    nints = nnewints;
    swap(ints, newints);
    swap(idxs, newidxs);
    FREE(newints);
    FREE(newidxs);

    sortInts(rles, nrs, rscount);
    FREE(rscount);
}

template <typename T>
AOMOIntegrals<T>::abrs_integrals::abrs_integrals(pqrs_integrals& pqrs, const bool pleq)
: Distributed<T>(pqrs.ctf)
{
    na = pqrs.np;
    nb = pqrs.nq;
    nr = pqrs.nr;
    ns = pqrs.ns;

    nrs = (pqrs.nints > 0 ? 1 : 0);
    for (size_t ipqrs = 1;ipqrs < pqrs.nints;ipqrs++)
    {
        assert( pqrs.idxs[ipqrs].l  > pqrs.idxs[ipqrs-1].l ||
               (pqrs.idxs[ipqrs].l == pqrs.idxs[ipqrs-1].l &&
                pqrs.idxs[ipqrs].k >= pqrs.idxs[ipqrs-1].k));
        if (pqrs.idxs[ipqrs].k != pqrs.idxs[ipqrs-1].k ||
            pqrs.idxs[ipqrs].l != pqrs.idxs[ipqrs-1].l) nrs++;
    }

    assert(allsum((long)nrs) <= nr*ns);

    rs = SAFE_MALLOC(idx2_t, nrs);
    nints = nrs*na*nb;
    ints = SAFE_MALLOC(T, nints);
    fill(ints, ints+nints, (T)0);

    size_t ipqrs = 0, irs = 0, iabrs = 0;
    if (nrs > 0)
    {
        rs[0].i = pqrs.idxs[0].k;
        rs[0].j = pqrs.idxs[0].l;
        assert(rs[0].i >= 0 && rs[0].i < nr);
        assert(rs[0].j >= 0 && rs[0].j < ns);
    }
    for (;ipqrs < pqrs.nints;ipqrs++)
    {
        assert(ipqrs >= 0 && ipqrs < pqrs.nints);
        if (ipqrs > 0 &&
            (pqrs.idxs[ipqrs].k != pqrs.idxs[ipqrs-1].k ||
             pqrs.idxs[ipqrs].l != pqrs.idxs[ipqrs-1].l))
        {
            irs++;
            iabrs += na*nb;
            assert(irs >= 0 && irs < nrs);
            rs[irs].i = pqrs.idxs[ipqrs].k;
            rs[irs].j = pqrs.idxs[ipqrs].l;
            assert(rs[irs].i >= 0 && rs[irs].i < nr);
            assert(rs[irs].j >= 0 && rs[irs].j < ns);
        }

        int p = pqrs.idxs[ipqrs].i;
        int q = pqrs.idxs[ipqrs].j;

        assert(iabrs >= 0 && iabrs+na*nb <= nints);
        assert(p >= 0 && p < na);
        assert(q >= 0 && q < nb);
        ints[iabrs+p+q*na] = pqrs.ints[ipqrs];
        if (p != q && pleq)
            ints[iabrs+q+p*na] = pqrs.ints[ipqrs];
    }
    assert(irs == nrs-1 || nrs == 0);
    assert(iabrs == nints-na*nb || nrs == 0);
    assert(ipqrs == pqrs.nints);

    pqrs.free();
}

template <typename T>
typename AOMOIntegrals<T>::abrs_integrals AOMOIntegrals<T>::abrs_integrals::transform(Index index, const char trans, const int nc, const double* C, const int ldc)
{
    abrs_integrals out = *this;

    out.rs = SAFE_MALLOC(idx2_t, nrs);
    copy(rs, rs+nrs, out.rs);

    if (index == A)
    {
        out.na = nc;
        out.nints = nrs*nc*nb;
        out.ints = SAFE_MALLOC(T, out.nints);

        if (nc == 0) return out;

        if (trans == 'N')
        {
            gemm('T', 'N', nc, nb*nrs, na, 1.0,        C, ldc,
                                                    ints,  na,
                                           0.0, out.ints,  nc);
        }
        else
        {
            gemm('N', 'N', nc, nb*nrs, na, 1.0,        C, ldc,
                                                    ints,  na,
                                           0.0, out.ints,  nc);
        }
    }
    else
    {
        out.nb = nc;
        out.nints = nrs*na*nc;
        out.ints = SAFE_MALLOC(T, out.nints);

        if (nc == 0) return out;

        T* tmp = SAFE_MALLOC(T, max(nints, out.nints));
        transpose(na, nb*nrs, 1.0, ints, na, 0.0, tmp, nb*nrs);

        if (trans == 'N')
        {
            gemm('T', 'N', nc, na*nrs, nb, 1.0,        C, ldc,
                                                     tmp,  nb,
                                           0.0, out.ints,  nc);
        }
        else
        {
            gemm('N', 'N', nc, na*nrs, nb, 1.0,        C, ldc,
                                                     tmp,  nb,
                                           0.0, out.ints,  nc);
        }

        FREE(tmp);

        tmp = SAFE_MALLOC(T, out.nints);
        transpose(nc*nrs, na, 1.0, out.ints, nc*nrs, 0.0, tmp, na);
        swap(tmp, out.ints);
        FREE(tmp);

        /*
        size_t iin = 0;
        size_t iout = 0;
        for (size_t irs = 0;irs < nrs;irs++)
        {
            if (trans == 'N')
            {
                gemm('N', 'N', na, nc, nb, 1.0,      ints+iin,  na,
                                                            C, ldc,
                                           0.0, out.ints+iout,  na);
            }
            else
            {
                gemm('N', 'T', na, nc, nb, 1.0,      ints+iin,  na,
                                                            C, ldc,
                                           0.0, out.ints+iout,  na);
            }

            iin += na*nb;
            iout += na*nc;
        }
        */
    }

    return out;
}

template <typename T>
void AOMOIntegrals<T>::abrs_integrals::transcribe(DistTensor<T>& tensor, bool assymij, bool assymkl, bool reverse)
{
    assert(nr*ns == allsum((long)nrs));

    if (reverse)
    {
        if (assymij) assert(ns == nb);
        if (assymkl) assert(nr == na);

        vector<tkv_pair<T> > pairs;

        /*
         * (ab|rs) -> <sb|ra>
         */
        for (size_t idx = 0,irs = 0;irs < nrs;irs++)
        {
            int r = rs[irs].i;
            int s = rs[irs].j;

            for (int b = 0;b < nb;b++)
            {
                for (int a = 0;a < na;a++)
                {
                    if ((!assymij || s < b) && (!assymkl || r < a))
                        pairs.push_back(tkv_pair<T>(((a*nr+r)*nb+b)*ns+s, ints[idx]));
                    idx++;
                }
            }
        }

        if (assymij) assert(tensor.getSymmetry()[0] == AS);
        if (assymkl) assert(tensor.getSymmetry()[2] == AS);
        assert(tensor.getLengths()[0] == ns);
        assert(tensor.getLengths()[1] == nb);
        assert(tensor.getLengths()[2] == nr);
        assert(tensor.getLengths()[3] == na);
        tensor.writeRemoteData(1, 0, pairs.size(), pairs.data());
        pairs.clear();

        if (assymij)
        {
            /*
             * -(ab|rs) -> <bs|ra>
             */
            for (size_t idx = 0,irs = 0;irs < nrs;irs++)
            {
                int r = rs[irs].i;
                int s = rs[irs].j;

                for (int b = 0;b < nb;b++)
                {
                    for (int a = 0;a < na;a++)
                    {
                        if (b < s && (!assymkl || r < a))
                            pairs.push_back(tkv_pair<T>(((((int64_t)a)*nr+r)*ns+s)*nb+b, ints[idx]));
                        idx++;
                    }
                }
            }

            tensor.writeRemoteData(-1, 1, pairs.size(), pairs.data());
        }
        else if (assymkl)
        {
            /*
             * -(ab|rs) -> <sb|ar>
             */
            for (size_t idx = 0,irs = 0;irs < nrs;irs++)
            {
                int r = rs[irs].i;
                int s = rs[irs].j;

                for (int b = 0;b < nb;b++)
                {
                    for (int a = 0;a < na;a++)
                    {
                        if (a < r)
                            pairs.push_back(tkv_pair<T>(((((int64_t)r)*na+a)*nb+b)*ns+s, ints[idx]));
                        idx++;
                    }
                }
            }

            tensor.writeRemoteData(-1, 1, pairs.size(), pairs.data());
        }
    }
    else
    {
        if (assymij) assert(na == nr);
        if (assymkl) assert(nb == ns);

        vector<tkv_pair<T> > pairs;

        /*
         * (ab|rs) -> <ar|bs>
         */
        for (size_t idx = 0,irs = 0;irs < nrs;irs++)
        {
            int r = rs[irs].i;
            int s = rs[irs].j;

            for (int b = 0;b < nb;b++)
            {
                for (int a = 0;a < na;a++)
                {
                    if ((!assymij || a < r) && (!assymkl || b < s))
                    {
                        pairs.push_back(tkv_pair<T>(((((int64_t)s)*nb+b)*nr+r)*na+a, ints[idx]));
                    }
                    idx++;
                }
            }
        }

        if (assymij) assert(tensor.getSymmetry()[0] == AS);
        if (assymkl) assert(tensor.getSymmetry()[2] == AS);
        assert(tensor.getLengths()[0] == na);
        assert(tensor.getLengths()[1] == nr);
        assert(tensor.getLengths()[2] == nb);
        assert(tensor.getLengths()[3] == ns);
        tensor.writeRemoteData(1, 0, pairs.size(), pairs.data());
        pairs.clear();

        if (assymij)
        {
            /*
             * -(ab|rs) -> <ra|bs>
             */
            for (size_t idx = 0,irs = 0;irs < nrs;irs++)
            {
                int r = rs[irs].i;
                int s = rs[irs].j;

                for (int b = 0;b < nb;b++)
                {
                    for (int a = 0;a < na;a++)
                    {
                        if (r < a && (!assymkl || b < s))
                            pairs.push_back(tkv_pair<T>(((((int64_t)s)*nb+b)*na+a)*nr+r, ints[idx]));
                        idx++;
                    }
                }
            }

            tensor.writeRemoteData(-1, 1, pairs.size(), pairs.data());
        }
        else if (assymkl)
        {
            /*
             * -(ab|rs) -> <ar|sb>
             */
            for (size_t idx = 0,irs = 0;irs < nrs;irs++)
            {
                int r = rs[irs].i;
                int s = rs[irs].j;

                for (int b = 0;b < nb;b++)
                {
                    for (int a = 0;a < na;a++)
                    {
                        if (s < b)
                            pairs.push_back(tkv_pair<T>(((((int64_t)b)*ns+s)*nr+r)*na+a, ints[idx]));
                        idx++;
                    }
                }
            }

            tensor.writeRemoteData(-1, 1, pairs.size(), pairs.data());
        }
    }
}

template <typename T>
void AOMOIntegrals<T>::abrs_integrals::free()
{
    FREE(ints);
    FREE(rs);
}

template <typename T>
void AOMOIntegrals<T>::doTransformation(const AOIntegrals<T>& ints)
{
    int N = this->uhf.getMolecule().getNumOrbitals();
    int nI = this->uhf.getMolecule().getNumAlphaElectrons();
    int ni = this->uhf.getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    int sizeAAII[] = {nA, nA, nI, nI};
    int sizeaaii[] = {na, na, ni, ni};
    int shapeNNNN[] = {NS, NS, NS, NS};

    DistTensor<T> ABIJ__(this->ctf, 4, sizeAAII, shapeNNNN, false);
    DistTensor<T> abij__(this->ctf, 4, sizeaaii, shapeNNNN, false);

    int64_t npair;
    T *cA, *ca, *cI, *ci;

    /*
     * Read transformation coefficients
     */
    this->uhf.getCA().getAllData(npair, cA);
    assert(npair == N*nA);
    this->uhf.getCa().getAllData(npair, ca);
    assert(npair == N*na);
    this->uhf.getCI().getAllData(npair, cI);
    assert(npair == N*nI);
    this->uhf.getCi().getAllData(npair, ci);
    assert(npair == N*ni);

    /*
     * Resort integrals so that each node has (pq|r_k s_l) where pq
     * are dense blocks for each sparse rs pair
     */
    pqrs_integrals pqrs(ints);
    pqrs.collect(true);
    abrs_integrals PQrs(pqrs, true);

    /*
     * First quarter-transformation
     */
    abrs_integrals PArs = PQrs.transform(B, 'N', nA, cA, N);
    abrs_integrals Pars = PQrs.transform(B, 'N', na, ca, N);
    abrs_integrals PIrs = PQrs.transform(B, 'N', nI, cI, N);
    abrs_integrals Pirs = PQrs.transform(B, 'N', ni, ci, N);
    PQrs.free();

    /*
     * Second quarter-transformation
     */
    abrs_integrals ABrs = PArs.transform(A, 'N', nA, cA, N);
    PArs.free();
    abrs_integrals abrs = Pars.transform(A, 'N', na, ca, N);
    Pars.free();
    abrs_integrals AIrs = PIrs.transform(A, 'N', nA, cA, N);
    abrs_integrals IJrs = PIrs.transform(A, 'N', nI, cI, N);
    PIrs.free();
    abrs_integrals airs = Pirs.transform(A, 'N', na, ca, N);
    abrs_integrals ijrs = Pirs.transform(A, 'N', ni, ci, N);
    Pirs.free();

    /*
     * Make <AB||CD>
     */
    pqrs_integrals rsAB(ABrs);
    rsAB.collect(false);

    abrs_integrals RSAB(rsAB, true);
    abrs_integrals RDAB = RSAB.transform(B, 'N', nA, cA, N);
    RSAB.free();

    abrs_integrals CDAB = RDAB.transform(A, 'N', nA, cA, N);
    RDAB.free();
    CDAB.transcribe(*this->ABCD_, true, true, false);
    CDAB.free();

    /*
     * Make <Ab|Cd> and <ab||cd>
     */
    pqrs_integrals rsab(abrs);
    rsab.collect(false);

    abrs_integrals RSab(rsab, true);
    abrs_integrals RDab = RSab.transform(B, 'N', nA, cA, N);
    abrs_integrals Rdab = RSab.transform(B, 'N', na, ca, N);
    RSab.free();

    abrs_integrals CDab = RDab.transform(A, 'N', nA, cA, N);
    RDab.free();
    CDab.transcribe(*this->AbCd_, false, false, false);
    CDab.free();

    abrs_integrals cdab = Rdab.transform(A, 'N', na, ca, N);
    Rdab.free();
    cdab.transcribe(*this->abcd_, true, true, false);
    cdab.free();

    /*
     * Make <AB||CI>, <aB|cI>, and <AB|IJ>
     */
    pqrs_integrals rsAI(AIrs);
    rsAI.collect(false);

    abrs_integrals RSAI(rsAI, true);
    abrs_integrals RCAI = RSAI.transform(B, 'N', nA, cA, N);
    abrs_integrals RcAI = RSAI.transform(B, 'N', na, ca, N);
    abrs_integrals RJAI = RSAI.transform(B, 'N', nI, cI, N);
    RSAI.free();

    abrs_integrals BCAI = RCAI.transform(A, 'N', nA, cA, N);
    RCAI.free();
    BCAI.transcribe(*this->ABCI_, true, false, false);
    BCAI.free();

    abrs_integrals bcAI = RcAI.transform(A, 'N', na, ca, N);
    RcAI.free();
    bcAI.transcribe(*this->aBcI_, false, false, false);
    bcAI.free();

    abrs_integrals BJAI = RJAI.transform(A, 'N', nA, cA, N);
    RJAI.free();
    BJAI.transcribe(ABIJ__, false, false, false);
    BJAI.free();

    /*
     * Make <Ab|Ci>, <ab||ci>, <Ab|Ij>, and <ab|ij>
     */
    pqrs_integrals rsai(airs);
    rsai.collect(false);

    abrs_integrals RSai(rsai, true);
    abrs_integrals RCai = RSai.transform(B, 'N', nA, cA, N);
    abrs_integrals Rcai = RSai.transform(B, 'N', na, ca, N);
    abrs_integrals RJai = RSai.transform(B, 'N', nI, cI, N);
    abrs_integrals Rjai = RSai.transform(B, 'N', ni, ci, N);
    RSai.free();

    abrs_integrals BCai = RCai.transform(A, 'N', nA, cA, N);
    RCai.free();
    BCai.transcribe(*this->AbCi_, false, false, false);
    BCai.free();

    abrs_integrals bcai = Rcai.transform(A, 'N', na, ca, N);
    Rcai.free();
    bcai.transcribe(*this->abci_, true, false, false);
    bcai.free();

    abrs_integrals BJai = RJai.transform(A, 'N', nA, cA, N);
    RJai.free();
    BJai.transcribe(*this->AbIj_, false, false, false);
    BJai.free();

    abrs_integrals bjai = Rjai.transform(A, 'N', na, ca, N);
    Rjai.free();
    bjai.transcribe(abij__, false, false, false);
    bjai.free();

    /*
     * Make <IJ||KL>, <IJ||KA>, <Ij|Ka>, <aI|bJ>, and <AI|BJ>
     */
    pqrs_integrals rsIJ(IJrs);
    rsIJ.collect(false);

    abrs_integrals RSIJ(rsIJ, true);
    abrs_integrals RBIJ = RSIJ.transform(B, 'N', nA, cA, N);
    abrs_integrals RbIJ = RSIJ.transform(B, 'N', na, ca, N);
    abrs_integrals RLIJ = RSIJ.transform(B, 'N', nI, cI, N);
    abrs_integrals RlIJ = RSIJ.transform(B, 'N', ni, ci, N);
    RSIJ.free();

    abrs_integrals ABIJ = RBIJ.transform(A, 'N', nA, cA, N);
    RBIJ.free();
    ABIJ.transcribe(*this->AIBJ_, false, false, false);
    ABIJ.free();

    abrs_integrals abIJ = RbIJ.transform(A, 'N', na, ca, N);
    RbIJ.free();
    abIJ.transcribe(*this->aIbJ_, false, false, false);
    abIJ.free();

    abrs_integrals akIJ = RlIJ.transform(A, 'N', na, ca, N);
    RlIJ.free();
    akIJ.transcribe(*this->IjKa_, false, false, true);
    akIJ.free();

    abrs_integrals AKIJ = RLIJ.transform(A, 'N', nA, cA, N);
    abrs_integrals KLIJ = RLIJ.transform(A, 'N', nI, cI, N);
    RLIJ.free();
    AKIJ.transcribe(*this->IJKA_, true, false, true);
    AKIJ.free();
    KLIJ.transcribe(*this->IJKL_, true, true, false);
    KLIJ.free();

    /*
     * Make <Ij|Kl>, <ij||kl>, <iJ|kA>, <ij||ka>, <Ai|Bj>, and <ai|bj>
     */
    pqrs_integrals rsij(ijrs);
    rsij.collect(false);

    abrs_integrals RSij(rsij, true);
    abrs_integrals RBij = RSij.transform(B, 'N', nA, cA, N);
    abrs_integrals Rbij = RSij.transform(B, 'N', na, ca, N);
    abrs_integrals RLij = RSij.transform(B, 'N', nI, cI, N);
    abrs_integrals Rlij = RSij.transform(B, 'N', ni, ci, N);
    RSij.free();

    abrs_integrals ABij = RBij.transform(A, 'N', nA, cA, N);
    RBij.free();
    ABij.transcribe(*this->AiBj_, false, false, false);
    ABij.free();

    abrs_integrals abij = Rbij.transform(A, 'N', na, ca, N);
    Rbij.free();
    abij.transcribe(*this->aibj_, false, false, false);
    abij.free();

    abrs_integrals AKij = RLij.transform(A, 'N', nA, cA, N);
    abrs_integrals KLij = RLij.transform(A, 'N', nI, cI, N);
    RLij.free();
    AKij.transcribe(*this->iJkA_, false, false, true);
    AKij.free();
    KLij.transcribe(*this->IjKl_, false, false, false);
    KLij.free();

    abrs_integrals akij = Rlij.transform(A, 'N', na, ca, N);
    abrs_integrals klij = Rlij.transform(A, 'N', ni, ci, N);
    Rlij.free();
    akij.transcribe(*this->ijka_, true, false, true);
    akij.free();
    klij.transcribe(*this->ijkl_, true, true, false);
    klij.free();

    /*
     * Make <AI||BJ> and <ai||bj>
     */
    (*this->AIBJ_)["AIBJ"] -= ABIJ__["ABJI"];
    (*this->aibj_)["aibj"] -= abij__["abji"];

    /*
     * Make <AB||IJ> and <ab||ij>
     */
    //(*this->ABIJ_)["ABIJ"] = 0.5*ABIJ__["ABIJ"];
    //(*this->abij_)["abij"] = 0.5*abij__["abij"];
    (*this->ABIJ_)["ABIJ"]  = ABIJ__["ABIJ"];
    (*this->ABIJ_)["ABIJ"] -= ABIJ__["ABJI"];
    (*this->abij_)["abij"]  = abij__["abij"];
    (*this->abij_)["abij"] -= abij__["abji"];

    /*
     * Make <Ai|bJ> = -<Ab|Ji> and <aI|Bj> = -<Ba|Ij>
     */
    (*this->AibJ_)["AbJi"] -= (*this->AbIj_)["AbJi"];
    (*this->aIBj_)["BaIj"] -= (*this->AbIj_)["BaIj"];

    if (nA > 0) FREE(cA);
    if (na > 0) FREE(ca);
    if (nI > 0) FREE(cI);
    if (ni > 0) FREE(ci);
}

INSTANTIATE_SPECIALIZATIONS(AOMOIntegrals);