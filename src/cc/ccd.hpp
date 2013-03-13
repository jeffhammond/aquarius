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

#ifndef _AQUARIUS_CC_CCD_HPP_
#define _AQUARIUS_CC_CCD_HPP_

#include "ccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCD : public CCSD<U>
{
    public:
        CCD(const input::Config& config, op::TwoElectronOperator<U>& moints)
        : tensor::Tensor<op::ExcitationOperator<U,2>,U>(*this),
          CCSD<U>(config, moints)
        {
            this->T(1) = 0;

            energy = 0.25*scalar(moints.getABIJ()*T(2));

            conv =          conv,T(2).getSpinCase(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,T(2).getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,T(2).getSpinCase(2).reduce(CTF_OP_MAXABS));
        }

        void _iterate()
        {
            using CCSD<U>::T;
            using CCSD<U>::D;
            using CCSD<U>::Z;

            op::TwoElectronOperator H(moints, op::TwoElectronOperator::AB|
                                              op::TwoElectronOperator::IJ|
                                              op::TwoElectronOperator::IJKL|
                                              op::TwoElectronOperator::AIBJ);

            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FAE = H.getAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FMI = H.getIJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNEF = H.getIJAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEF = H.getABCD();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIJ = H.getIJKL();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEI = H.getAIBJ();

            FAE["aa"] = 0.0;
            FMI["ii"] = 0.0;

            /**************************************************************************
             *
             * Intermediates
             */
            PROFILE_SECTION(calc_FMI)
            FMI["mi"] += 0.5*WMNEF["mnef"]*T(2)["efin"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WIJKL)
            WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*T(2)["efij"];
            PROFILE_STOP

            PROFILE_SECTION(calc_FAE)
            FAE["ae"] -= 0.5*WMNEF["mnef"]*T(2)["afmn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WMBEJ)
            WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
            PROFILE_STOP
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * T(1)->T(2) and T(2)->T(2)
             */
            PROFILE_SECTION(calc_WMNEF)
            Z(2)["abij"] = WMNEF["ijab"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_FAE)
            Z(2)["abij"] += FAE["af"]*T(2)["fbij"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_FMI)
            Z(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_ABCD)
            Z(2)["abij"] += 0.5*WABEF["abef"]*T(2)["efij"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_IJKL)
            Z(2)["abij"] += 0.5*WMNIJ["mnij"]*T(2)["abmn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_RING)
            Z(2)["abij"] -= WAMEI["amei"]*T(2)["ebmj"];
            PROFILE_STOP
            /*
             *************************************************************************/

            PROFILE_SECTION(calc_EN)

            Z *= D;
            Z -= T;
            T += Z;

            energy = 0.25*scalar(moints.getVABIJ()*T(2));

            conv =               Z(1).getSpinCase(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,Z(1).getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2).getSpinCase(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2).getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2).getSpinCase(2).reduce(CTF_OP_MAXABS));

            this->diis.extrapolate(T, Z);

            PROFILE_STOP
        }
};

}
}

#endif
