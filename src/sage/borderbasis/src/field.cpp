#include "include/field.h"

namespace math {

//----- FieldFN -----------------------------------

FieldFn::FieldFn(uint64_t pol)
{
    uint exp = 1;
    for(uint64_t tmp=3;tmp<pol;tmp=tmp*2+1)
        exp++;
    s = gf_init(exp,(uint)pol);
}

FieldFn::~FieldFn()
{
    gf_uninit(s);
}

uint64_t FieldFn::add(uint64_t x,uint64_t y) const
{
    return x^y;
}

uint64_t FieldFn::subtract(uint64_t x,uint64_t y) const
{
    return x^y;
}

uint64_t FieldFn::multiply(uint64_t x,uint64_t y) const
{
    return (uint64_t)gf_mul(s,(uint)x,(uint)y);
}

uint64_t FieldFn::divide(uint64_t x,uint64_t y) const
{
    return (uint64_t)gf_div(s,(uint)x,(uint)y);
}

uint64_t FieldFn::pow(uint64_t x,uint64_t y) const
{
    uint64_t result = 1;
    for(uint64_t i=0;i<y;i++)
        result = (uint64_t)gf_mul(s,(uint)result,(uint)x);
    return result;
}

} // namespace math
