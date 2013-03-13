/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice,This list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice,This list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_TENSOR_HPP_
#define _AQUARIUS_TENSOR_HPP_

#include <stdexcept>

#include "stl_ext/stl_ext.hpp"
#include "util/util.h"
#include "util/fortran.h"

namespace aquarius
{
namespace tensor
{

template <class Derived, class T> class Tensor;
template <class Derived, class T> class ScaledTensor;
template <class Derived, class T> class InvertedTensor;
template <class Derived, class T> class TensorMult;
template <class Derived, class T> class TensorDiv;

class TensorError;
class OutOfBoundsError;
class LengthMismatchError;
class IndexMismatchError;
class InvalidNdimError;
class InvalidLengthError;
class InvalidLdError;
class LdTooSmallError;
class SymmetryMismatchError;
class InvalidSymmetryError;
class InvalidStartError;

#define INHERIT_FROM_TENSOR(Derived,T) \
    protected: \
        using aquarius::tensor::Tensor< Derived,T >::derived; \
    public: \
        using aquarius::tensor::Tensor< Derived,T >::operator=; \
        using aquarius::tensor::Tensor< Derived,T >::operator+=; \
        using aquarius::tensor::Tensor< Derived,T >::operator-=; \
        using aquarius::tensor::Tensor< Derived,T >::operator*=; \
        using aquarius::tensor::Tensor< Derived,T >::operator/=; \
        using aquarius::tensor::Tensor< Derived,T >::operator*; \
        using aquarius::tensor::Tensor< Derived,T >::operator/; \
        Derived & operator=(const Derived & other) \
        { \
            static_cast< aquarius::tensor::Tensor< Derived,T >& >(*this).operator=(other); \
            return *this; \
        } \
    private:

template <class Derived, typename T>
class Tensor
{
    protected:
        Derived& derived;

    public:
        typedef T dtype;

        Tensor(Derived& derived) : derived(derived) {}

        virtual ~Tensor() {}

        const Derived& getDerived() const { return derived; }

        /**********************************************************************
         *
         * Operators with scalars
         *
         *********************************************************************/
        Tensor<Derived,T>& operator=(const T val)
        {
            sum(val, (T)0);
            return *this;
        }

        Tensor<Derived,T>& operator+=(const T val)
        {
            sum(val, (T)1);
            return *this;
        }

        Tensor<Derived,T>& operator-=(const T val)
        {
            sum(-val, (T)1);
            return *this;
        }

        Tensor<Derived,T>& operator*=(const T val)
        {
            mult(val);
            return *this;
        }

        Tensor<Derived,T>& operator/=(const T val)
        {
            mult(1.0/val);
            return *this;
        }

        /**********************************************************************
         *
         * Binary operations (multiplication and division)
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator=(const TensorMult<cvDerived,T>& other)
        {
            mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator+=(const TensorMult<cvDerived,T>& other)
        {
            mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)1);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator-=(const TensorMult<cvDerived,T>& other)
        {
            mult(-other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)1);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator=(const TensorDiv<cvDerived,T>& other)
        {
            div(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator+=(const TensorDiv<cvDerived,T>& other)
        {
            div(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)1);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator-=(const TensorDiv<cvDerived,T>& other)
        {
            div(-other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)1);
            return *this;
        }

        /**********************************************************************
         *
         * Unary operations (assignment, summation, multiplication, and division)
         *
         *********************************************************************/
        Tensor<Derived,T>& operator=(const Tensor<Derived,T>& other)
        {
            sum((T)1, false, other, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator=(const Tensor<cvDerived,T>& other)
        {
            sum((T)1, false, other, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator+=(const Tensor<cvDerived,T>& other)
        {
            sum((T)1, false, other, (T)1);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator-=(const Tensor<cvDerived,T>& other)
        {
            sum((T)(-1), false, other, (T)1);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator*=(const Tensor<cvDerived,T>& other)
        {
            mult((T)1, false, *this, false, other, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator/=(const Tensor<cvDerived,T>& other)
        {
            div((T)1, false, *this, false, other, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator=(const ScaledTensor<cvDerived,T>& other)
        {
            sum(other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator+=(const ScaledTensor<cvDerived,T>& other)
        {
            sum(other.factor_, other.conj_, other.tensor_, (T)1);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator-=(const ScaledTensor<cvDerived,T>& other)
        {
            sum(-other.factor_, other.conj_, other.tensor_, (T)1);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator*=(const ScaledTensor<cvDerived,T>& other)
        {
            mult(other.factor_, false, *this, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator/=(const ScaledTensor<cvDerived,T>& other)
        {
            div((T)1/other.factor_, false, *this, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator=(const InvertedTensor<cvDerived,T>& other)
        {
            invert(other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator+=(const InvertedTensor<cvDerived,T>& other)
        {
            invert(other.factor_, other.conj_, other.tensor_, (T)1);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator-=(const InvertedTensor<cvDerived,T>& other)
        {
            invert(-other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator*=(const InvertedTensor<cvDerived,T>& other)
        {
            div(other.factor_, false, *this, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(Tensor<Derived,T>&))
        operator/=(const InvertedTensor<cvDerived,T>& other)
        {
            mult((T)1/other.factor_, false, *this, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        /**********************************************************************
         *
         * Intermediate operations
         *
         *********************************************************************/
        friend ScaledTensor<Derived,T> operator*(const double factor, Tensor<Derived,T>& other)
        {
            return ScaledTensor<Derived,T>(other.derived, factor);
        }

        friend ScaledTensor<const Derived,T> operator*(const double factor, const Tensor<Derived,T>& other)
        {
            return ScaledTensor<const Derived,T>(other.derived, factor);
        }

        ScaledTensor<Derived,T> operator*(const T factor)
        {
            return ScaledTensor<Derived,T>(derived, factor);
        }

        ScaledTensor<const Derived,T> operator*(const T factor) const
        {
            return ScaledTensor<const Derived,T>(derived, factor);
        }

        friend InvertedTensor<Derived,T> operator/(const double factor, const Tensor<Derived,T>& other)
        {
            return InvertedTensor<Derived,T>(other.derived, factor);
        }

        ScaledTensor<Derived,T> operator/(const T factor)
        {
            return ScaledTensor<Derived,T>(derived, (T)1/factor);
        }

        ScaledTensor<const Derived,T> operator/(const T factor) const
        {
            return ScaledTensor<const Derived,T>(derived, (T)1/factor);
        }

        ScaledTensor<Derived,T> operator-()
        {
            return ScaledTensor<Derived,T>(derived, (T)(-1));
        }

        ScaledTensor<const Derived,T> operator-() const
        {
            return ScaledTensor<const Derived,T>(derived, (T)(-1));
        }

        friend ScaledTensor<const Derived,T> conj(const Tensor<Derived,T>& t)
        {
            return ScaledTensor<const Derived,T>(t.derived, (T)1, true);
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorMult<Derived,T>))
        operator*(const Tensor<cvDerived,T>& other) const
        {
            return TensorMult<Derived,T>(ScaledTensor<const Derived,T>(derived, (T)1),
                                         ScaledTensor<const Derived,T>(other.derived, (T)1));
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorMult<Derived,T>))
        operator/(const Tensor<cvDerived,T>& other) const
        {
            return TensorDiv<Derived,T>(ScaledTensor<const Derived,T>(derived, (T)1),
                                        ScaledTensor<const Derived,T>(other.derived, (T)1));
        }

        /**********************************************************************
         *
         * Stubs
         *
         *********************************************************************/

        /*
         * this = alpha*this + beta*A*B
         */
        virtual void mult(const T alpha, bool conja, const Tensor<Derived,T>& A,
                                         bool conjb, const Tensor<Derived,T>& B, const T beta) = 0;

        /*
         * this = alpha*this
         */
        virtual void mult(const T alpha) = 0;

        /*
         * this = alpha*this + beta*A/B
         */
        virtual void div(const T alpha, bool conja, const Tensor<Derived,T>& A,
                                        bool conjb, const Tensor<Derived,T>& B, const T beta) = 0;

        /*
         * this = alpha*this + beta*A
         */
        virtual void sum(const T alpha, bool conja, const Tensor<Derived,T>& A, const T beta) = 0;

        /*
         * this = alpha*this + beta
         */
        virtual void sum(const T alpha, const T beta) = 0;

        /*
         * this = alpha*this + beta/A
         */
        virtual void invert(const T alpha, bool conja, const Tensor<Derived,T>& A, const T beta) = 0;
};

template <class Derived, typename T>
class ScaledTensor
{
    public:
        Derived& tensor_;
        T factor_;
        bool conj_;

        template <typename cvDerived>
        ScaledTensor(const ScaledTensor<cvDerived,T>& other)
        : tensor_(other.tensor_), factor_(other.factor_), conj_(other.conj_) {}

        ScaledTensor(Derived& tensor, const T factor, const bool conj=false)
        : tensor_(tensor), factor_(factor), conj_(conj) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        ScaledTensor<Derived,T> operator-() const
        {
            ScaledTensor<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return ret;
        }

        friend ScaledTensor<const Derived,T> conj(const ScaledTensor<Derived,T>& st)
        {
            ScaledTensor<Derived,T> ret(st);
            ret.conj = !ret.conj;
            return ret;
        }

        /**********************************************************************
         *
         * Unary tensor operations
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const Tensor<cvDerived,T>& other)
        {
            tensor_.sum((T)1, false, other, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const Tensor<cvDerived,T>& other)
        {
            tensor_.sum((T)1, false, other, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const Tensor<cvDerived,T>& other)
        {
            tensor_.sum((T)(-1), false, other, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator*=(const Tensor<cvDerived,T>& other)
        {
            tensor_.mult(factor_, false, tensor_, false, other, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator/=(const Tensor<cvDerived,T>& other)
        {
            tensor_.div(factor_, false, tensor_, false, other, (T)0);
            return *this;
        }

        ScaledTensor<Derived,T>& operator=(const ScaledTensor<Derived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.sum(-other.factor_, other.conj_, other.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator*=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.mult(factor_*other.factor_, false, tensor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator/=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.div(factor_/other.factor_, false, tensor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.invert(other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.invert(other.factor_, other.conj_, other.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.invert(-other.factor_, other.conj_, other.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator*=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.div(factor_*other.factor_, false, tensor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator/=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.mult(factor_/other.factor_, false, tensor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        /**********************************************************************
         *
         * Binary tensor operations
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const TensorMult<cvDerived,T>& other)
        {
            tensor_.mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const TensorMult<cvDerived,T>& other)
        {
            tensor_.mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const TensorMult<cvDerived,T>& other)
        {
            tensor_.mult(-other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const TensorDiv<cvDerived,T>& other)
        {
            tensor_.div(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const TensorDiv<cvDerived,T>& other)
        {
            tensor_.div(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const TensorDiv<cvDerived,T>& other)
        {
            tensor_.div(-other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorMult<Derived,T>))
        operator*(const ScaledTensor<cvDerived,T>& other) const
        {
            return TensorMult<Derived,T>(*this, other);
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorMult<Derived,T>))
        operator*(const Tensor<cvDerived,T>& other) const
        {
            return TensorMult<Derived,T>(*this, ScaledTensor<const Derived,T>(other.getDerived(), (T)1));
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorDiv<Derived,T>))
        operator/(const ScaledTensor<cvDerived,T>& other) const
        {
            return TensorDiv<Derived,T>(*this, other);
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorDiv<Derived,T>))
        operator/(const Tensor<cvDerived,T>& other) const
        {
            return TensorDiv<Derived,T>(*this, ScaledTensor<const Derived,T>(other.getDerived(), (T)1));
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        ScaledTensor<Derived,T> operator*(const T factor) const
        {
            ScaledTensor<Derived,T> it(*this);
            it.factor_ *= factor;
            return it;
        }

        friend ScaledTensor<Derived,T> operator*(const T factor, const ScaledTensor<Derived,T>& other)
        {
            return other*factor;
        }

        ScaledTensor<Derived,T> operator/(const T factor) const
        {
            ScaledTensor<Derived,T> it(*this);
            it.factor_ /= factor;
            return it;
        }

        friend InvertedTensor<Derived,T> operator/(const T factor, const ScaledTensor<Derived,T>& other)
        {
            return InvertedTensor<Derived,T>(other.tensor_, factor/other.factor_);
        }

        ScaledTensor<Derived,T>& operator=(const T val)
        {
            tensor_.sum(val, (T)0);
            return *this;
        }

        ScaledTensor<Derived,T>& operator+=(const T val)
        {
            tensor_.sum(val, factor_);
            return *this;
        }

        ScaledTensor<Derived,T>& operator-=(const T val)
        {
            tensor_.sum(-val, factor_);
            return *this;
        }

        ScaledTensor<Derived,T>& operator*=(const T val)
        {
            tensor_.mult(val);
            return *this;
        }

        ScaledTensor<Derived,T>& operator/=(const T val)
        {
            tensor_.mult((T)1/val);
            return *this;
        }
};

template <class Derived1, class Derived2, class T>
typename std::enable_if<std::is_same<const Derived1, const Derived2>::value,TensorMult<Derived1,T> >::type
operator*(const Tensor<Derived1,T>& t1, const ScaledTensor<Derived2,T>& t2)
{
    return TensorMult<Derived1,T>(ScaledTensor<const Derived1,T>(t1.getDerived(), (T)1), t2);
}

template <class Derived1, class Derived2, class T>
typename std::enable_if<std::is_same<const Derived1, const Derived2>::value,TensorDiv<Derived1,T> >::type
operator/(const Tensor<Derived1,T>& t1, const ScaledTensor<Derived2,T>& t2)
{
    return TensorDiv<Derived1,T>(ScaledTensor<const Derived1,T>(t1.getDerived(), (T)1), t2);
}

template <class Derived, typename T>
class InvertedTensor
{
    private:
        const InvertedTensor& operator=(const InvertedTensor<Derived,T>& other);

    public:
        Derived& tensor_;
        T factor_;
        bool conj_;

        InvertedTensor(Derived& tensor, const T factor, const bool conj=false)
        : tensor_(tensor), factor_(factor), conj_(conj) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        InvertedTensor<Derived,T> operator-() const
        {
            InvertedTensor<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return *this;
        }

        friend InvertedTensor<Derived,T> conj(const InvertedTensor<Derived,T>& tm)
        {
            InvertedTensor<Derived,T> ret(tm);
            ret.conj_ = !ret.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        InvertedTensor<Derived,T> operator*(const T factor) const
        {
            InvertedTensor<Derived,T> ret(*this);
            ret.factor_ *= factor;
            return ret;
        }

        InvertedTensor<Derived,T> operator/(const T factor) const
        {
            InvertedTensor<Derived,T> ret(*this);
            ret.factor_ /= factor;
            return ret;
        }

        friend InvertedTensor<Derived,T> operator*(const T factor, const InvertedTensor<Derived,T>& other)
        {
            return other*factor;
        }
};

template <class Derived, typename T>
class TensorMult
{
    private:
        const TensorMult& operator=(const TensorMult<Derived,T>& other);

    public:
        ScaledTensor<const Derived,T> A_;
        ScaledTensor<const Derived,T> B_;
        T factor_;

        template <class Derived1, class Derived2>
        TensorMult(const ScaledTensor<Derived1,T>& A, const ScaledTensor<Derived2,T>& B)
        : A_(A), B_(B), factor_(A_.factor_*B_.factor_) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        TensorMult<Derived,T> operator-() const
        {
            TensorMult<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return *this;
        }

        friend TensorMult<Derived,T> conj(const TensorMult<Derived,T>& tm)
        {
            TensorMult<Derived,T> ret(tm);
            ret.A_.conj_ = !ret.A_.conj_;
            ret.B_.conj_ = !ret.B_.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        TensorMult<Derived,T> operator*(const T factor) const
        {
            TensorMult<Derived,T> ret(*this);
            ret.factor_ *= factor;
            return ret;
        }

        TensorMult<Derived,T> operator/(const T factor) const
        {
            TensorMult<Derived,T> ret(*this);
            ret.factor_ /= factor;
            return ret;
        }

        friend TensorMult<Derived,T> operator*(const T factor, const TensorMult<Derived,T>& other)
        {
            return other*factor;
        }
};

template <class Derived, typename T>
class TensorDiv
{
    private:
        const TensorDiv& operator=(const TensorDiv<Derived,T>& other);

    public:
        ScaledTensor<const Derived,T> A_;
        ScaledTensor<const Derived,T> B_;
        T factor_;

        template <class Derived1, class Derived2>
        TensorDiv(const ScaledTensor<Derived1,T>& A, const ScaledTensor<Derived2,T>& B)
        : A_(A), B_(B), factor_(A_.factor_/B_.factor_) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        TensorDiv<Derived,T> operator-() const
        {
            TensorDiv<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return *this;
        }

        friend TensorDiv<Derived,T> conj(const TensorDiv<Derived,T>& tm)
        {
            TensorDiv<Derived,T> ret(tm);
            ret.A_.conj_ = !ret.A_.conj_;
            ret.B_.conj_ = !ret.B_.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        TensorDiv<Derived,T> operator*(const T factor) const
        {
            TensorDiv<Derived,T> ret(*this);
            ret.factor_ *= factor;
            return ret;
        }

        TensorDiv<Derived,T> operator/(const T factor) const
        {
            TensorDiv<Derived,T> ret(*this);
            ret.factor_ /= factor;
            return ret;
        }

        friend TensorDiv<Derived,T> operator*(const T factor, const TensorDiv<Derived,T>& other)
        {
            return other*factor;
        }
};

class TensorError : public std::exception
{
    public:
        virtual const char* what() const throw() = 0;
};

class OutOfBoundsError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "out-of-bounds read or write"; }
};

class LengthMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "length mismatch error"; }
};

class IndexMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "index mismatch error"; }
};

class InvalidNdimError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid number of dimensions"; }
};

class InvalidLengthError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid length"; }
};

class InvalidLdError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid leading dimension"; }
};

class LdTooSmallError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "leading dimension is too small"; }
};

class SymmetryMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "symmetry mismatch error"; }
};

class InvalidSymmetryError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid symmetry value"; }
};

class InvalidStartError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid start value"; }
};

}
}

#endif
