#ifndef __FIELD_H__
#define __FIELD_H__

#include "definitions.h"
#include "../galois/gf.h"

namespace math {

template<typename T>
class IField
{
public:
    IField(){};
    virtual ~IField(){};

    virtual T add(T x,T y) const = 0;
    virtual T subtract(T x,T y) const = 0;
    virtual T multiply(T x,T y) const = 0;
    virtual T divide(T x,T y) const = 0;
    virtual T pow(T x,T y) const = 0;
};

class FieldFn : public IField<uint64_t>
{
public:
    FieldFn(uint64_t pol);
    virtual ~FieldFn();

    virtual uint64_t add(uint64_t x,uint64_t y) const OVERRIDE;
    virtual uint64_t subtract(uint64_t x,uint64_t y) const OVERRIDE;
    virtual uint64_t multiply(uint64_t x,uint64_t y) const OVERRIDE;
    virtual uint64_t divide(uint64_t x,uint64_t y) const OVERRIDE;
    virtual uint64_t pow(uint64_t x,uint64_t y) const OVERRIDE;

private:
    state* s;
};

} // namespace math

#endif // __FIELD_H__

