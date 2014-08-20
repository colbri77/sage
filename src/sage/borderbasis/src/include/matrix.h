#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "m4ri/m4ri.h"
#include "m4rie/m4rie.h"
#include "definitions.h"

namespace matrix {

template<typename T>
class IMatrix
{
public:
    IMatrix(uint r,uint c);
    virtual ~IMatrix();

    virtual T get(uint r,uint c) const = 0;
    virtual void set(uint r,uint c,T value) = 0;
    virtual void toRowEchelon(bool full) = 0;

    uint getRows() const;
    uint getColumns() const;

private:
    uint rows;
    uint columns;
};

class MatrixF2 : public IMatrix<uint64_t>
{
    public:
        MatrixF2(uint r,uint c);
        virtual ~MatrixF2();

        virtual uint64_t get(uint r,uint c) const OVERRIDE;
        virtual void set(uint r,uint c,uint64_t value) OVERRIDE;
        virtual void toRowEchelon(bool full) OVERRIDE;

    private:
        mzd_t* matrix;
        bool inRowEchelonForm;
};

class MatrixFn : public IMatrix<uint64_t>
{
    public:
        MatrixFn(uint64_t polynomial,uint r,uint c);
        virtual ~MatrixFn();

        virtual uint64_t get(uint r,uint c) const OVERRIDE;
        virtual void set(uint r,uint c,uint64_t value) OVERRIDE;
        virtual void toRowEchelon(bool full) OVERRIDE;

    private:
        gf2e* field;
        mzed_t* matrix;
        bool inRowEchelonForm;
};

} // namespace matrix

#endif // __MATRIX_H__
