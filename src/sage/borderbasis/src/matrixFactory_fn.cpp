#include "include/matrixFactory_fn.h"

namespace matrix {

MatrixFactory_Fn::MatrixFactory_Fn(uint64_t minPolynomial)
: IMatrixFactory<uint64_t>(),
polynomial(minPolynomial)
{
    ENSURE(minPolynomial>1, "MatrixFactory_Fn(minPolynomial): invalid polynomial");
}

MatrixFactory_Fn::~MatrixFactory_Fn()
{

}

TAKE_OWN IMatrix<uint64_t>* MatrixFactory_Fn::create(uint r,uint c) const
{
    IMatrix<uint64_t>* result = NULL;

    if(polynomial<4) {
        result = new MatrixF2(r,c);
    }
    else if(polynomial<0x800) {
        result = new MatrixFn(polynomial,r,c);
    }
    else {
        ASSERT_NOT_REACHED;
    }

    return result;
}

} // namespace matrix
