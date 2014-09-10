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
    create(0);
}

MonomialFactoryDegLex::MonomialFactoryDegLex()
: IMonomialFactory(),
monomials(new FastFlexibleArray()),
indet(0)
{

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


//-----MonomialFactoryNoOrderPosGF2-------------------------------

MonomialFactoryNoOrderPosGF2::MonomialFactoryNoOrderPosGF2(uint indet)
:MonomialFactoryNoOrderPos(indet)
{

}

MonomialFactoryNoOrderPosGF2::~MonomialFactoryNoOrderPosGF2()
{

}

TAKE_OWN IMonomial* MonomialFactoryNoOrderPosGF2::create() const
{
    return (IMonomial*)new DegLexMonomialNoOrderPosGF2(indet);
}

TAKE_OWN IMonomial* MonomialFactoryNoOrderPosGF2::create(uint64_t pos) const
{
    NOT_IMPLEMENTED;
}

//----------MonomialFactoryDegLexGF2---------------------------------

MonomialFactoryDegLexGF2::MonomialFactoryDegLexGF2(uint indet)
: MonomialFactoryDegLex()
{
    MonomialFactoryDegLex::indet = indet;
    create(0);
}

MonomialFactoryDegLexGF2::~MonomialFactoryDegLexGF2()
{
}

TAKE_OWN IMonomial* MonomialFactoryDegLexGF2::create() const
{
    return (IMonomial*)(monomials->get(0));
}

TAKE_OWN IMonomial* MonomialFactoryDegLexGF2::create(uint64_t pos) const
{
    IMonomial* result = (IMonomial*)(monomials->get(pos));

    if(result == NULL) {
        DegLexMonomialGF2* elem = new DegLexMonomialGF2(pos,indet,monomials);
        monomials->add(pos,elem);
        result = elem;
    }
    return result;
}



} // namespace polynomial
