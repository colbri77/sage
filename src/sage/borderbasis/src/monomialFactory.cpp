#include "include/monomialFactory.h"

#include "include/degLexMonomial.h"

namespace polynomial {

//-----MonomialFactoryNoOrderPos-------------------------------

MonomialFactoryNoOrderPos::MonomialFactoryNoOrderPos(uint indet)
:IMonomialFactory(),
indet(indet)
{

}

MonomialFactoryNoOrderPos::~MonomialFactoryNoOrderPos()
{

}

bool MonomialFactoryNoOrderPos::supportsGetPos() const
{
    return false;
}

TAKE_OWN IMonomial* MonomialFactoryNoOrderPos::create() const
{
    return (IMonomial*)new DegLexMonomialNoOrderPos(indet);
}

TAKE_OWN IMonomial* MonomialFactoryNoOrderPos::create(uint64_t pos) const
{
    NOT_IMPLEMENTED;
}

//----------MonomialFactoryDegLex---------------------------------

MonomialFactoryDegLex::MonomialFactoryDegLex(uint indet)
: IMonomialFactory(),
monomials(new FastFlexibleArray()),
indet(indet)
{
    DegLexMonomial* one = new DegLexMonomial(indet,monomials);
    monomials->add(0,one);
}

MonomialFactoryDegLex::~MonomialFactoryDegLex()
{
    delete monomials;
}

bool MonomialFactoryDegLex::supportsGetPos() const
{
    return true;
}

TAKE_OWN IMonomial* MonomialFactoryDegLex::create() const
{
    return (IMonomial*)(monomials->get(0));
}

TAKE_OWN IMonomial* MonomialFactoryDegLex::create(uint64_t pos) const
{
    IMonomial* result = (IMonomial*)(monomials->get(pos));

    if(result == NULL) {
        DegLexMonomial* elem = new DegLexMonomial(pos,indet,monomials);
        monomials->add(pos,elem);
        result = elem;
    }

    return result;
}

} // namespace polynomial
