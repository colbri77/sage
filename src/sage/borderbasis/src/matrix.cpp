#include "include/matrix.h"

namespace matrix {

//########## IMatrix ####################################

template<typename T>
IMatrix<T>::IMatrix(uint r,uint c)
: rows(r),
columns(c)
{

}

template<typename T>
IMatrix<T>::~IMatrix()
{

}

template<typename T>
uint IMatrix<T>::getRows() const
{
    return rows;
}

template<typename T>
uint IMatrix<T>::getColumns() const
{
    return columns;
}

template class IMatrix<uint64_t>;
template class IMatrix<int64_t>;


//########## MatrixF2 ###################################

MatrixF2::MatrixF2(uint r,uint c)
    : IMatrix<uint64_t>(r,c),
    matrix(mzd_init(r,c)),
    inRowEchelonForm(false)
{

}

MatrixF2::~MatrixF2()
{
    mzd_free(matrix);
}

uint64_t MatrixF2::get(uint r,uint c) const
{
    return mzd_read_bit(matrix,r,c);
}

void MatrixF2::set(uint r,uint c,uint64_t value)
{
    inRowEchelonForm = false;
    mzd_write_bit(matrix,r,c,value);
}

void MatrixF2::toRowEchelon(bool full)
{
    mzd_echelonize_m4ri(matrix,full ? 1 : 0,0);
    inRowEchelonForm = true;
}


//########## MatrixFn ##################################

MatrixFn::MatrixFn(uint64_t polynomial,uint r,uint c)
    : IMatrix<uint64_t>(r,c),
    inRowEchelonForm(false)
{
    field = gf2e_init(polynomial);
    matrix = mzed_init(field,r,c);
}

MatrixFn::~MatrixFn()
{
    mzed_free(matrix);
    gf2e_free(field);
}

uint64_t MatrixFn::get(uint r,uint c) const
{
    return mzed_read_elem(matrix,r,c);
}

void MatrixFn::set(uint r,uint c,uint64_t value)
{
    inRowEchelonForm = false;
    mzed_write_elem(matrix,r,c,value);
}

void MatrixFn::toRowEchelon(bool full)
{
    mzed_echelonize(matrix,full ? 1 : 0);
    inRowEchelonForm = true;
}


} // namespace matrix
