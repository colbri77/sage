#ifndef __MATRIXFACTORY_FN_H__
#define __MATRIXFACTORY_FN_H__

#include "i_matrixFactory.h"

namespace matrix {

class MatrixFactory_Fn : public IMatrixFactory<uint64_t>
{
    public:
        MatrixFactory_Fn(uint64_t minPolynomial);
        virtual ~MatrixFactory_Fn();

        virtual TAKE_OWN IMatrix<uint64_t>* create(uint r,uint c) const;

    protected:

    private:
        uint64_t polynomial;
};

} // namespace matrix

#endif // __MATRIXFACTORY_FN_H__
