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

#include "time/time.hpp"
#include "util/util.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::scf;
using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;

template <typename T>
AOMOIntegrals<T>::AOMOIntegrals(const string& name, const Config& config)
: MOIntegrals<T>("aomoints", name, config)
{
    this->getProduct("H").addRequirement(Requirement("eri","I"));
}

template <typename T>
AOMOIntegrals<T>::pqrs_integrals::pqrs_integrals(int norb, const ERI& aoints)
: Distributed(aoints.arena)
{
    PROFILE_FUNCTION

    ns = nr = nq = np = norb;

    const vector<T>& oldints = aoints.ints;
    const vector<idx4_t>& oldidxs = aoints.idxs;
    size_t noldints = oldints.size();
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

        assert(idx.i < np);
        assert(idx.j < nq);
        assert(idx.k < nr);
        assert(idx.l < ns);

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

    PROFILE_STOP
}

template <typename T>
AOMOIntegrals<T>::pqrs_integrals::pqrs_integrals(abrs_integrals& abrs)
: Distributed(abrs.arena)
{
    PROFILE_FUNCTION

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

    PROFILE_STOP
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
    PROFILE_FUNCTION

    assert(sizeof(short) == sizeof(int16_t));
    MPI::Datatype IDX4_T_TYPE = MPI::Datatype(MPI_SHORT).Create_contiguous(4);
    IDX4_T_TYPE.Commit();

    size_t nrs;
    size_t *rscount;
    sortInts(rles, nrs, rscount);

    size_t *sendcount = SAFE_MALLOC(size_t, nproc);
    size_t *recvcount = SAFE_MALLOC(size_t, nproc);
    fill(sendcount, sendcount+nproc, 0);
    fill(recvcount, recvcount+nproc, 0);

    for (int i = 0;i < nproc;i++)
    {
        for (size_t rs = (nrs*i)/nproc;rs < (nrs*(i+1))/nproc;rs++)
        {
            assert(rs >= 0 && rs < nrs);
            sendcount[i] += rscount[rs];
        }
    }

    FREE(rscount);

    PROFILE_SECTION(collect_comm)
    this->arena.Alltoall(sendcount, recvcount, 1);
    PROFILE_STOP

    size_t nnewints = 0;
    for (int i = 0;i < nproc;i++) nnewints += recvcount[i];

//    assert(allsum(nints) == allsum(nnewints));

    T* newints = SAFE_MALLOC(T, nnewints);
    idx4_t* newidxs = SAFE_MALLOC(idx4_t, nnewints);

    int *realsendcount = SAFE_MALLOC(int, nproc);
    int *realrecvcount = SAFE_MALLOC(int, nproc);

    for (int i = 0;i < nproc;i++)
    {
        assert(sendcount[i] <= INT_MAX);
        assert(recvcount[i] <= INT_MAX);
        realsendcount[i] = (int)sendcount[i];
        realrecvcount[i] = (int)recvcount[i];
    }

    PROFILE_SECTION(collect_comm)
    this->arena.Alltoallv(ints, realsendcount, newints, realrecvcount);

    this->arena.Alltoallv(idxs, realsendcount, newidxs, realrecvcount, IDX4_T_TYPE);
    PROFILE_STOP

    FREE(realsendcount);
    FREE(realrecvcount);
    FREE(sendcount);
    FREE(recvcount);

    nints = nnewints;
    swap(ints, newints);
    swap(idxs, newidxs);
    FREE(newints);
    FREE(newidxs);

    sortInts(rles, nrs, rscount);
    FREE(rscount);

    PROFILE_STOP
}

template <typename T>
AOMOIntegrals<T>::abrs_integrals::abrs_integrals(pqrs_integrals& pqrs, const bool pleq)
: Distributed(pqrs.arena)
{
    PROFILE_FUNCTION

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

    rs = SAFE_MALLOC(idx2_t, nrs);
    nints = nrs*na*nb;
    ints = SAFE_MALLOC(T, nints);
    fill(ints, ints+nints, (T)0);

    size_t ipqrs = 0, irs = 0, iabrs = 0;
    if (nrs > 0)
    {
        rs[0].i = pqrs.idxs[0].k;
        rs[0].j = pqrs.idxs[0].l;
        assert(rs[0].i < nr);
        assert(rs[0].j < ns);
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
            assert(rs[irs].i < nr);
            assert(rs[irs].j < ns);
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

    PROFILE_STOP
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

        PROFILE_SECTION(transform)
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
        PROFILE_STOP
    }
    else
    {
        out.nb = nc;
        out.nints = nrs*na*nc;
        out.ints = SAFE_MALLOC(T, out.nints);

        if (nc == 0) return out;

        PROFILE_SECTION(transform)
        /*
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
        */

        ///*
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
        //*/
        PROFILE_STOP
    }

    return out;
}

template <typename T>
void AOMOIntegrals<T>::abrs_integrals::transcribe(DistTensor<T>& tensor, bool assymij, bool assymkl, bool reverse)
{
    PROFILE_FUNCTION

    if (reverse)
    {
        if (assymij) assert(ns == nb);
        if (assymkl) assert(nr == na);

        vector<tkv_pair<T> > pairs;

        /*
         * (ab|rs) -> <sb|ra>
         */
        PROFILE_SECTION(1)
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
        PROFILE_STOP

        PROFILE_SECTION(2)
        if (assymij) assert(tensor.getSymmetry()[0] == AS);
        if (assymkl) assert(tensor.getSymmetry()[2] == AS);
        assert(tensor.getLengths()[0] == ns);
        assert(tensor.getLengths()[1] == nb);
        assert(tensor.getLengths()[2] == nr);
        assert(tensor.getLengths()[3] == na);
        tensor.writeRemoteData(1, 0, pairs);
        pairs.clear();
        PROFILE_STOP

        if (assymij)
        {
            /*
             * -(ab|rs) -> <bs|ra>
             */
            PROFILE_SECTION(3)
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
            PROFILE_STOP

            PROFILE_SECTION(4)
            tensor.writeRemoteData(-1, 1, pairs);
            PROFILE_STOP
        }
        else if (assymkl)
        {
            /*
             * -(ab|rs) -> <sb|ar>
             */
            PROFILE_SECTION(5)
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
            PROFILE_STOP

            PROFILE_SECTION(6)
            tensor.writeRemoteData(-1, 1, pairs);
            PROFILE_STOP
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
        PROFILE_SECTION(7)
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
        PROFILE_STOP

        PROFILE_SECTION(8)
        if (assymij) assert(tensor.getSymmetry()[0] == AS);
        if (assymkl) assert(tensor.getSymmetry()[2] == AS);
        assert(tensor.getLengths()[0] == na);
        assert(tensor.getLengths()[1] == nr);
        assert(tensor.getLengths()[2] == nb);
        assert(tensor.getLengths()[3] == ns);
        tensor.writeRemoteData(1, 0, pairs);
        pairs.clear();
        PROFILE_STOP

        if (assymij)
        {
            /*
             * -(ab|rs) -> <ra|bs>
             */
            PROFILE_SECTION(9)
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
            PROFILE_STOP

            PROFILE_SECTION(10)
            tensor.writeRemoteData(-1, 1, pairs);
            PROFILE_STOP
        }
        else if (assymkl)
        {
            /*
             * -(ab|rs) -> <ar|sb>
             */
            PROFILE_SECTION(11)
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
            PROFILE_STOP

            PROFILE_SECTION(12)
            tensor.writeRemoteData(-1, 1, pairs);
            PROFILE_STOP
        }
    }

    PROFILE_STOP
}

template <typename T>
void AOMOIntegrals<T>::abrs_integrals::free()
{
    FREE(ints);
    FREE(rs);
}

template <typename T>
void AOMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    const MOSpace<T>& occ = this->template get<MOSpace<T> >("occ");
    const MOSpace<T>& vrt = this->template get<MOSpace<T> >("vrt");

    const DistTensor<T>& Fa = this->template get<DistTensor<T> >("Fa");
    const DistTensor<T>& Fb = this->template get<DistTensor<T> >("Fb");

    this->put("H", new TwoElectronOperator<T>(OneElectronOperator<T>(occ, vrt, Fa, Fb, true)));
    TwoElectronOperator<T>& H = this->template get<TwoElectronOperator<T> >("H");

    const ERI& ints = this->template get<ERI>("I");

    const DistTensor<T>& cA_ = vrt.Calpha;
    const DistTensor<T>& ca_ = vrt.Cbeta;
    const DistTensor<T>& cI_ = occ.Calpha;
    const DistTensor<T>& ci_ = occ.Cbeta;

    int N = occ.nao;
    int nI = occ.nalpha;
    int ni = occ.nbeta;
    int nA = vrt.nalpha;
    int na = vrt.nbeta;

    DistTensor<T> ABIJ__(arena, 4, vec(nA,nA,nI,nI), vec(NS,NS,NS,NS), false);
    DistTensor<T> abij__(arena, 4, vec(na,na,ni,ni), vec(NS,NS,NS,NS), false);

    vector<T> cA, ca, cI, ci;

    /*
     * Read transformation coefficients
     */
    cA_.getAllData(cA);
    assert(cA.size() == N*nA);
    ca_.getAllData(ca);
    assert(ca.size() == N*na);
    cI_.getAllData(cI);
    assert(cI.size() == N*nI);
    ci_.getAllData(ci);
    assert(ci.size() == N*ni);

    /*
     * Resort integrals so that each node has (pq|r_k s_l) where pq
     * are dense blocks for each sparse rs pair
     */
    pqrs_integrals pqrs(N, ints);
    pqrs.collect(true);
    abrs_integrals PQrs(pqrs, true);

    /*
     * First quarter-transformation
     */
    abrs_integrals PArs = PQrs.transform(B, 'N', nA, cA.data(), N);
    abrs_integrals Pars = PQrs.transform(B, 'N', na, ca.data(), N);
    abrs_integrals PIrs = PQrs.transform(B, 'N', nI, cI.data(), N);
    abrs_integrals Pirs = PQrs.transform(B, 'N', ni, ci.data(), N);
    PQrs.free();

    /*
     * Second quarter-transformation
     */
    abrs_integrals ABrs = PArs.transform(A, 'N', nA, cA.data(), N);
    PArs.free();
    abrs_integrals abrs = Pars.transform(A, 'N', na, ca.data(), N);
    Pars.free();
    abrs_integrals AIrs = PIrs.transform(A, 'N', nA, cA.data(), N);
    abrs_integrals IJrs = PIrs.transform(A, 'N', nI, cI.data(), N);
    PIrs.free();
    abrs_integrals airs = Pirs.transform(A, 'N', na, ca.data(), N);
    abrs_integrals ijrs = Pirs.transform(A, 'N', ni, ci.data(), N);
    Pirs.free();

    /*
     * Make <AB||CD>
     */
    pqrs_integrals rsAB(ABrs);
    rsAB.collect(false);

    abrs_integrals RSAB(rsAB, true);
    abrs_integrals RDAB = RSAB.transform(B, 'N', nA, cA.data(), N);
    RSAB.free();

    abrs_integrals CDAB = RDAB.transform(A, 'N', nA, cA.data(), N);
    RDAB.free();
    PROFILE_SECTION(j)
    CDAB.transcribe(H.getABCD()(2,0,2,0), true, true, false);
    PROFILE_STOP
    CDAB.free();

    /*
     * Make <Ab|Cd> and <ab||cd>
     */
    pqrs_integrals rsab(abrs);
    rsab.collect(false);

    abrs_integrals RSab(rsab, true);
    abrs_integrals RDab = RSab.transform(B, 'N', nA, cA.data(), N);
    abrs_integrals Rdab = RSab.transform(B, 'N', na, ca.data(), N);
    RSab.free();

    abrs_integrals CDab = RDab.transform(A, 'N', nA, cA.data(), N);
    RDab.free();
    PROFILE_SECTION(k)
    CDab.transcribe(H.getABCD()(1,0,1,0), false, false, false);
    PROFILE_STOP
    CDab.free();

    abrs_integrals cdab = Rdab.transform(A, 'N', na, ca.data(), N);
    Rdab.free();
    PROFILE_SECTION(l)
    cdab.transcribe(H.getABCD()(0,0,0,0), true, true, false);
    PROFILE_STOP
    cdab.free();

    /*
     * Make <AB||CI>, <aB|cI>, and <AB|IJ>
     */
    pqrs_integrals rsAI(AIrs);
    rsAI.collect(false);

    abrs_integrals RSAI(rsAI, true);
    abrs_integrals RCAI = RSAI.transform(B, 'N', nA, cA.data(), N);
    abrs_integrals RcAI = RSAI.transform(B, 'N', na, ca.data(), N);
    abrs_integrals RJAI = RSAI.transform(B, 'N', nI, cI.data(), N);
    RSAI.free();

    abrs_integrals BCAI = RCAI.transform(A, 'N', nA, cA.data(), N);
    RCAI.free();
    PROFILE_SECTION(m)
    BCAI.transcribe(H.getABCI()(2,0,1,1), true, false, false);
    PROFILE_STOP
    BCAI.free();

    abrs_integrals bcAI = RcAI.transform(A, 'N', na, ca.data(), N);
    RcAI.free();
    PROFILE_SECTION(n)
    bcAI.transcribe(H.getABCI()(1,0,0,1), false, false, false);
    PROFILE_STOP
    bcAI.free();

    abrs_integrals BJAI = RJAI.transform(A, 'N', nA, cA.data(), N);
    RJAI.free();
    PROFILE_SECTION(o)
    BJAI.transcribe(ABIJ__, false, false, false);
    PROFILE_STOP
    BJAI.free();

    /*
     * Make <Ab|Ci>, <ab||ci>, <Ab|Ij>, and <ab|ij>
     */
    pqrs_integrals rsai(airs);
    rsai.collect(false);

    abrs_integrals RSai(rsai, true);
    abrs_integrals RCai = RSai.transform(B, 'N', nA, cA.data(), N);
    abrs_integrals Rcai = RSai.transform(B, 'N', na, ca.data(), N);
    abrs_integrals RJai = RSai.transform(B, 'N', nI, cI.data(), N);
    abrs_integrals Rjai = RSai.transform(B, 'N', ni, ci.data(), N);
    RSai.free();

    abrs_integrals BCai = RCai.transform(A, 'N', nA, cA.data(), N);
    RCai.free();
    PROFILE_SECTION(p)
    BCai.transcribe(H.getABCI()(1,0,1,0), false, false, false);
    PROFILE_STOP
    BCai.free();

    abrs_integrals bcai = Rcai.transform(A, 'N', na, ca.data(), N);
    Rcai.free();
    PROFILE_SECTION(q)
    bcai.transcribe(H.getABCI()(0,0,0,0), true, false, false);
    PROFILE_STOP
    bcai.free();

    abrs_integrals BJai = RJai.transform(A, 'N', nA, cA.data(), N);
    RJai.free();
    PROFILE_SECTION(r)
    BJai.transcribe(H.getABIJ()(1,0,0,1), false, false, false);
    PROFILE_STOP
    BJai.free();

    abrs_integrals bjai = Rjai.transform(A, 'N', na, ca.data(), N);
    Rjai.free();
    PROFILE_SECTION(s)
    bjai.transcribe(abij__, false, false, false);
    PROFILE_STOP
    bjai.free();

    /*
     * Make <IJ||KL>, <IJ||KA>, <Ij|Ka>, <aI|bJ>, and <AI|BJ>
     */
    pqrs_integrals rsIJ(IJrs);
    rsIJ.collect(false);

    abrs_integrals RSIJ(rsIJ, true);
    abrs_integrals RBIJ = RSIJ.transform(B, 'N', nA, cA.data(), N);
    abrs_integrals RbIJ = RSIJ.transform(B, 'N', na, ca.data(), N);
    abrs_integrals RLIJ = RSIJ.transform(B, 'N', nI, cI.data(), N);
    abrs_integrals RlIJ = RSIJ.transform(B, 'N', ni, ci.data(), N);
    RSIJ.free();

    abrs_integrals ABIJ = RBIJ.transform(A, 'N', nA, cA.data(), N);
    RBIJ.free();
    PROFILE_SECTION(t)
    ABIJ.transcribe(H.getAIBJ()(1,1,1,1), false, false, false);
    PROFILE_STOP
    ABIJ.free();

    abrs_integrals abIJ = RbIJ.transform(A, 'N', na, ca.data(), N);
    RbIJ.free();
    PROFILE_SECTION(u)
    abIJ.transcribe(H.getAIBJ()(0,1,0,1), false, false, false);
    PROFILE_STOP
    abIJ.free();

    abrs_integrals akIJ = RlIJ.transform(A, 'N', na, ca.data(), N);
    RlIJ.free();
    PROFILE_SECTION(v)
    akIJ.transcribe(H.getIJKA()(0,1,0,1), false, false, true);
    PROFILE_STOP
    akIJ.free();

    abrs_integrals AKIJ = RLIJ.transform(A, 'N', nA, cA.data(), N);
    abrs_integrals KLIJ = RLIJ.transform(A, 'N', nI, cI.data(), N);
    RLIJ.free();
    PROFILE_SECTION(w)
    AKIJ.transcribe(H.getIJKA()(0,2,1,1), true, false, true);
    PROFILE_STOP
    AKIJ.free();
    PROFILE_SECTION(x)
    KLIJ.transcribe(H.getIJKL()(0,2,0,2), true, true, false);
    PROFILE_STOP
    KLIJ.free();

    /*
     * Make <Ij|Kl>, <ij||kl>, <iJ|kA>, <ij||ka>, <Ai|Bj>, and <ai|bj>
     */
    pqrs_integrals rsij(ijrs);
    rsij.collect(false);

    abrs_integrals RSij(rsij, true);
    abrs_integrals RBij = RSij.transform(B, 'N', nA, cA.data(), N);
    abrs_integrals Rbij = RSij.transform(B, 'N', na, ca.data(), N);
    abrs_integrals RLij = RSij.transform(B, 'N', nI, cI.data(), N);
    abrs_integrals Rlij = RSij.transform(B, 'N', ni, ci.data(), N);
    RSij.free();

    abrs_integrals ABij = RBij.transform(A, 'N', nA, cA.data(), N);
    RBij.free();
    PROFILE_SECTION(a)
    ABij.transcribe(H.getAIBJ()(1,0,1,0), false, false, false);
    PROFILE_STOP
    ABij.free();

    abrs_integrals abij = Rbij.transform(A, 'N', na, ca.data(), N);
    Rbij.free();
    PROFILE_SECTION(b)
    abij.transcribe(H.getAIBJ()(0,0,0,0), false, false, false);
    PROFILE_STOP
    abij.free();

    abrs_integrals AKij = RLij.transform(A, 'N', nA, cA.data(), N);
    abrs_integrals KLij = RLij.transform(A, 'N', nI, cI.data(), N);
    RLij.free();
    PROFILE_SECTION(c)
    AKij.transcribe(H.getIJKA()(0,1,1,0), false, false, true);
    PROFILE_STOP
    AKij.free();
    PROFILE_SECTION(d)
    KLij.transcribe(H.getIJKL()(0,1,0,1), false, false, false);
    PROFILE_STOP
    KLij.free();

    abrs_integrals akij = Rlij.transform(A, 'N', na, ca.data(), N);
    abrs_integrals klij = Rlij.transform(A, 'N', ni, ci.data(), N);
    Rlij.free();
    PROFILE_SECTION(e)
    akij.transcribe(H.getIJKA()(0,0,0,0), true, false, true);
    PROFILE_STOP
    akij.free();
    PROFILE_SECTION(f)
    klij.transcribe(H.getIJKL()(0,0,0,0), true, true, false);
    PROFILE_STOP
    klij.free();

    /*
     * Make <AI||BJ> and <ai||bj>
     */
    PROFILE_SECTION(g)
    H.getAIBJ()(1,1,1,1)["AIBJ"] -= ABIJ__["ABJI"];
    H.getAIBJ()(0,0,0,0)["aibj"] -= abij__["abji"];
    PROFILE_STOP

    /*
     * Make <AB||IJ> and <ab||ij>
     */
    //(*this->ABIJ_)["ABIJ"] = 0.5*ABIJ__["ABIJ"];
    //(*H.getABIJ()_)["abij"] = 0.5*abij__["abij"];
    PROFILE_SECTION(h)
    H.getABIJ()(2,0,0,2)["ABIJ"]  = ABIJ__["ABIJ"];
    H.getABIJ()(2,0,0,2)["ABIJ"] -= ABIJ__["ABJI"];
    H.getABIJ()(0,0,0,0)["abij"]  = abij__["abij"];
    H.getABIJ()(0,0,0,0)["abij"] -= abij__["abji"];
    PROFILE_STOP

    /*
     * Make <Ai|bJ> = -<Ab|Ji> and <aI|Bj> = -<Ba|Ij>
     */
    PROFILE_SECTION(i)
    H.getAIBJ()(1,0,0,1)["AbJi"] = -H.getABIJ()(1,0,0,1)["AbJi"];
    H.getAIBJ()(0,1,1,0)["BaIj"] = -H.getABIJ()(1,0,0,1)["BaIj"];
    PROFILE_STOP
}

INSTANTIATE_SPECIALIZATIONS(AOMOIntegrals);
REGISTER_TASK(AOMOIntegrals<double>,"aomoints");
