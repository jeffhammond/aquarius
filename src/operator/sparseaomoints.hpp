#ifndef _AQUARIUS_OPERATOR_SPARSEAOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_SPARSEAOMOINTS_HPP_

#include "util/global.hpp"

#include "scf/aouhf.hpp"
#include "integrals/2eints.hpp"

#include "moints.hpp"
#include "aomoints.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class SparseAOMOIntegrals : public MOIntegrals<T>
{
    public:
        SparseAOMOIntegrals(const string& name, input::Config& config);

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
