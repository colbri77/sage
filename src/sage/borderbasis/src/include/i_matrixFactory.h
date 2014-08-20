#ifndef __I_MATRIXFACTORY_H__
#define __I_MATRIXFACTORY_H__

#include "matrix.h"

namespace matrix {

template<typename T>
class IMatrixFactory
{
    public:
        IMatrixFactory() {}
        virtual ~IMatrixFactory() {}

        virtual TAKE_OWN IMatrix<T>* create(uint r,uint c) const = 0;
};

} // namespace matrix

#endif // __I_MATRIXFACTORY_H__
