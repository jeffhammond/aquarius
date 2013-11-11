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

#include "ctf_tensor.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::tensor;

template <typename T>
map<const tCTF_World<T>*,pair<int,CTFTensor<T>*> > CTFTensor<T>::scalars;

/*
 * Create a scalar (0-dimensional tensor)
 */
template <typename T>
CTFTensor<T>::CTFTensor(const Arena& arena, T scalar)
: IndexableTensor< CTFTensor<T>,T >(), Resource(arena), len(0), sym(0)
{
    allocate();
    *dt = scalar;
    register_scalar();
}

/*
 * Create a scalar (0-dimensional tensor) on the same arena as A
 */
template <typename T>
CTFTensor<T>::CTFTensor(const CTFTensor<T>& A, T scalar)
: IndexableTensor< CTFTensor<T>,T >(), Resource(A.arena),
  len(0), sym(0)
{
    allocate();
    *dt = scalar;
    register_scalar();
}

/*
 * Create a tensor of the same size and shape as A, optionally copying or zeroing the data
 */
template <typename T>
CTFTensor<T>::CTFTensor(const CTFTensor<T>& A, bool copy, bool zero)
: IndexableTensor< CTFTensor<T>,T >(A.ndim), Resource(A.arena),
  len(A.len), sym(A.sym)
{
    allocate();

    if (copy)
    {
        *this = A;
    }
    else if (zero)
    {
        *dt = (T)0;
    }

    register_scalar();
}

template <typename T>
CTFTensor<T>::CTFTensor(CTFTensor<T>* A)
: IndexableTensor< CTFTensor<T>,T >(A->ndim), Resource(A->arena),
  len(A->len), sym(A->sym)
{
    dt = A->dt;
    delete A;
    register_scalar();
}

template <typename T>
CTFTensor<T>::CTFTensor(const CTFTensor<T>& A, const vector<int>& start_A, const vector<int>& len_A)
: IndexableTensor< CTFTensor<T>,T >(A.ndim), Resource(A.arena),
  len(len_A), sym(A.sym)
{
    allocate();
    slice((T)1, false, A, start_A, (T)0);
    register_scalar();
}

/*
 * Create a tensor of the specified size and shape, optionally zeroing the data
 */
template <typename T>
CTFTensor<T>::CTFTensor(const Arena& arena, int ndim, const vector<int>& len, const vector<int>& sym,
                          bool zero)
: IndexableTensor< CTFTensor<T>,T >(ndim), Resource(arena),
  len(len), sym(sym)
{
    assert(len.size() == ndim);
    assert(sym.size() == ndim);

    #ifdef VALIDATE_INPUTS
    validate_tensor(ndim, len.data(), NULL, sym.data());
    #endif //VALIDATE_INPUTS

    allocate();
    if (zero) *dt = (T)0;

    register_scalar();
}

template <typename T>
CTFTensor<T>::~CTFTensor()
{
    unregister_scalar();
    free();
}

template <typename T>
void CTFTensor<T>::allocate()
{
    dt = new tCTF_Tensor<T>(ndim, len.data(), sym.data(), arena.ctf<T>(), "FIXME", 1);
}

template <typename T>
void CTFTensor<T>::free()
{
    delete dt;
}

template <typename T>
void CTFTensor<T>::register_scalar()
{
    if (scalars.find(&arena.ctf<T>()) == scalars.end())
    {
        /*
         * If we are the first, make a new entry and put a new
         * scalar in it. The entry in scalars must be made FIRST,
         * since the new scalar will call this constructor too.
         */
        scalars[&arena.ctf<T>()].first = -1;
        scalars[&arena.ctf<T>()].second = new CTFTensor<T>(arena);
    }

    scalars[&arena.ctf<T>()].first++;
    //cout << "creating: " << scalars[&arena.ctf<T>()].first << endl;
}

template <typename T>
void CTFTensor<T>::unregister_scalar()
{
    /*
     * The last tensor (besides the scalar in scalars)
     * will delete the entry, so if it does not exist
     * then we must be that scalar and nothing needs to be done.
     */
    if (scalars.find(&arena.ctf<T>()) == scalars.end()) return;

    //cout << "deleting: " << scalars[&arena.ctf<T>()].first << endl;
    if (--scalars[&arena.ctf<T>()].first == 0)
    {
        /*
         * The entry must be deleted FIRST, so that the scalar
         * knows to do nothing.
         */
        CTFTensor<T>* scalar = scalars[&arena.ctf<T>()].second;
        scalars.erase(&arena.ctf<T>());
        delete scalar;
    }
}

template <typename T>
CTFTensor<T>& CTFTensor<T>::scalar() const
{
    return *scalars[&arena.ctf<T>()].second;
}

template <typename T>
void CTFTensor<T>::resize(int ndim, const vector<int>& len, const vector<int>& sym, bool zero)
{
    assert(len.size() == ndim);
    assert(sym.size() == ndim);

    this->ndim = ndim;
    this->len = len;
    this->sym = sym;

    free();
    allocate();
    if (zero) *dt = (T)0;
}

template <typename T>
T* CTFTensor<T>::getRawData(int64_t& size)
{
    return const_cast<T*>(const_cast<const CTFTensor<T>&>(*this).getRawData(size));
}

template <typename T>
const T* CTFTensor<T>::getRawData(int64_t& size) const
{
    long_int size_;
    T* data = dt->get_raw_data(&size_);
    size = size_;
    return data;
}

template <typename T>
void CTFTensor<T>::slice(T alpha, bool conja, const CTFTensor<T>& A,
                          const vector<int>& start_A, T beta)
{
    slice(alpha, conja, A, start_A, beta, vector<int>(this->ndim, 0), len);
}

template <typename T>
void CTFTensor<T>::slice(T alpha, bool conja, const CTFTensor<T>& A,
                          T beta, const vector<int>& start_B)
{
    slice(alpha, conja, A, vector<int>(this->ndim, 0), beta, start_B, A.len);
}

template <typename T>
void CTFTensor<T>::slice(T alpha, bool conja, const CTFTensor<T>& A, const vector<int>& start_A,
                          T  beta,                                     const vector<int>& start_B,
                                                                       const vector<int>& len)
{
    assert(this->ndim == A.ndim);

    vector<int> end_A(this->ndim);
    vector<int> end_B(this->ndim);

    for (int i = 0;i < this->ndim;i++)
    {
        end_A[i] = start_A[i]+len[i];
        end_B[i] = start_B[i]+len[i];
        assert(sym[i] == A.sym[i] && sym[i] == NS);
        assert(start_A[i] >= 0);
        assert(start_B[i] >= 0);
        assert(end_A[i] <= A.len[i]);
        assert(end_B[i] <= this->len[i]);
    }

    dt->sum_slice(start_B.data(), end_B.data(), beta, *A.dt, start_A.data(), end_A.data(), alpha);
}

template <typename T>
void CTFTensor<T>::div(T alpha, bool conja, const CTFTensor<T>& A,
                                 bool conjb, const CTFTensor<T>& B, T beta)
{
    const_cast<tCTF_Tensor<T>*>(A.dt)->align(*dt);
    const_cast<tCTF_Tensor<T>*>(B.dt)->align(*dt);
    int64_t size, size_A, size_B;
    T* raw_data = getRawData(size);
    const T* raw_data_A = A.getRawData(size_A);
    const T* raw_data_B = B.getRawData(size_B);
    assert(size == size_A);
    assert(size == size_B);
    if (conja)
    {
        if (conjb)
        {
            for (int64_t i = 0;i < size;i++)
            {
                if (abs(raw_data_B[i]) > DBL_MIN)
                {
                    raw_data[i] = beta*raw_data[i] + alpha*conj(raw_data_A[i])/conj(raw_data_B[i]);
                }
            }
        }
        else
        {
            for (int64_t i = 0;i < size;i++)
            {
                if (abs(raw_data_B[i]) > DBL_MIN)
                {
                    raw_data[i] = beta*raw_data[i] + alpha*conj(raw_data_A[i])/raw_data_B[i];
                }
            }
        }
    }
    else
    {
        if (conjb)
        {
            for (int64_t i = 0;i < size;i++)
            {
                if (abs(raw_data_B[i]) > DBL_MIN)
                {
                    raw_data[i] = beta*raw_data[i] + alpha*raw_data_A[i]/conj(raw_data_B[i]);
                }
            }
        }
        else
        {
            for (int64_t i = 0;i < size;i++)
            {
                if (abs(raw_data_B[i]) > DBL_MIN)
                {
                    raw_data[i] = beta*raw_data[i] + alpha*raw_data_A[i]/raw_data_B[i];
                }
            }
        }
    }
}

template <typename T>
void CTFTensor<T>::invert(T alpha, bool conja, const CTFTensor<T>& A, T beta)
{
    dt->align(*A.dt);
    int64_t size, size_A;
    T* raw_data = getRawData(size);
    const T* raw_data_A = A.getRawData(size_A);
    assert(size == size_A);
    if (conja)
    {
        for (int64_t i = 0;i < size;i++)
        {
            if (abs(raw_data_A[i]) > DBL_MIN)
            {
                raw_data[i] = beta*raw_data[i] + alpha/conj(raw_data_A[i]);
            }
        }
    }
    else
    {
        for (int64_t i = 0;i < size;i++)
        {
            if (abs(raw_data_A[i]) > DBL_MIN)
            {
                raw_data[i] = beta*raw_data[i] + alpha/raw_data_A[i];
            }
        }
    }
}

template <typename T>
void CTFTensor<T>::print(FILE* fp, double cutoff) const
{
    dt->print(fp, cutoff);
}

template <typename T>
void CTFTensor<T>::compare(FILE* fp, const CTFTensor<T>& other, double cutoff) const
{
    dt->compare(*other.dt, fp, cutoff);
}

template <typename T>
typename real_type<T>::type CTFTensor<T>::norm(int p) const
{
    T ans = (T)0;
    if (p == 00)
    {
        ans = dt->reduce(CTF_OP_NORM_INFTY);
    }
    else if (p == 1)
    {
        ans = dt->reduce(CTF_OP_NORM1);
    }
    else if (p == 2)
    {
        ans = dt->reduce(CTF_OP_NORM2);
    }
    return abs(ans);
}

template <typename T>
void CTFTensor<T>::mult(T alpha, bool conja, const CTFTensor<T>& A, const string& idx_A,
                                  bool conjb, const CTFTensor<T>& B, const string& idx_B,
                         T  beta,                                     const string& idx_C)
{
    dt->contract(alpha, *A.dt, idx_A.c_str(),
                        *B.dt, idx_B.c_str(),
                  beta,        idx_C.c_str());
}

template <typename T>
void CTFTensor<T>::sum(T alpha, bool conja, const CTFTensor<T>& A, const string& idx_A,
                        T  beta,                                     const string& idx_B)
{
    dt->sum(alpha, *A.dt, idx_A.c_str(),
             beta,        idx_B.c_str());
}

template <typename T>
void CTFTensor<T>::scale(T alpha, const string& idx_A)
{
    dt->scale(alpha, idx_A.c_str());
}

template <typename T>
T CTFTensor<T>::dot(bool conja, const CTFTensor<T>& A, const string& idx_A,
                    bool conjb,                        const string& idx_B) const
{
    CTFTensor<T>& dt = scalar();
    vector<T> val;
    dt.mult(1, conja,     A, idx_A,
               conjb, *this, idx_B,
            0,                  "");
    dt.getAllData(val);
    assert(val.size()==1);
    return val[0];
}

template <typename T>
void CTFTensor<T>::weight(const vector<const vector<T>*>& d)
{
    if (this->ndim == 0) return;

    assert(d.size() == this->ndim);
    for (int i = 0;i < d.size();i++) assert(d[i]->size() == len[i]);

    vector<tkv_pair<T> > pairs;
    getLocalData(pairs);

    for (int i = 0;i < pairs.size();i++)
    {
        int64_t k = pairs[i].k;

        T den = 0;
        for (int j = 0;j < this->ndim;j++)
        {
            int o = k%len[j];
            k = k/len[j];
            den += (*d[j])[o];
        }

        pairs[i].d /= den;
    }

    writeRemoteData(pairs);
}

INSTANTIATE_SPECIALIZATIONS(CTFTensor);