#include "include/monomialFactory.h"

#include "include/degLexMonomial.h"
#include "include/degRevLexMonomial.h"

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
    return NULL;
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
:MonomialFactoryNoOrderPos(indet),
excludedIndets(new bool[indet])
{
    for(uint i=0;i<indet;i++)
        excludedIndets[i] = false;
}

MonomialFactoryNoOrderPosGF2::~MonomialFactoryNoOrderPosGF2()
{
    delete excludedIndets;
}

TAKE_OWN IMonomial* MonomialFactoryNoOrderPosGF2::create() const
{
    return (IMonomial*)new DegLexMonomialNoOrderPosGF2(indet,excludedIndets);
}

TAKE_OWN IMonomial* MonomialFactoryNoOrderPosGF2::create(uint64_t pos) const
{
    NOT_IMPLEMENTED;
}

void MonomialFactoryNoOrderPosGF2::excludeIndet(uint index)
{
    excludedIndets[index] = true;
}

//----------MonomialFactoryDegLexGF2---------------------------------

MonomialFactoryDegLexGF2::MonomialFactoryDegLexGF2(uint indet)
: MonomialFactoryDegLex(),
excludedIndets(new bool[indet])
{
    for(uint i=0;i<indet;i++)
        excludedIndets[i] = false;
    MonomialFactoryDegLex::indet = indet;
    create(0);
}

MonomialFactoryDegLexGF2::MonomialFactoryDegLexGF2()
: MonomialFactoryDegLex(),
excludedIndets(NULL)
{
    MonomialFactoryDegLex::indet = 0;
}

MonomialFactoryDegLexGF2::~MonomialFactoryDegLexGF2()
{
    delete excludedIndets;
}

TAKE_OWN IMonomial* MonomialFactoryDegLexGF2::create(uint64_t pos) const
{
    IMonomial* result = (IMonomial*)(monomials->get(pos));
    if(result == NULL) {
        DegLexMonomialGF2* elem = new DegLexMonomialGF2(pos,indet,monomials,excludedIndets);
        monomials->add(pos,elem);
        result = elem;
    }
    return result;
}

void MonomialFactoryDegLexGF2::excludeIndet(uint index)
{
    excludedIndets[index] = true;
}

//-----MonomialFactoryDegRevLex--------------------------------
MonomialFactoryDegRevLex::MonomialFactoryDegRevLex(uint indet)
: MonomialFactoryDegLex()
{
    MonomialFactoryDegLex::indet = indet;
    create(0);
}

MonomialFactoryDegRevLex::~MonomialFactoryDegRevLex()
{

}

TAKE_OWN IMonomial* MonomialFactoryDegRevLex::create(uint64_t pos) const
{
    IMonomial* result = (IMonomial*)(MonomialFactoryDegLex::monomials->get(pos));

    if(result == NULL) {
        DegRevLexMonomial* elem = new DegRevLexMonomial(pos,MonomialFactoryDegLex::indet,MonomialFactoryDegLex::monomials);
        MonomialFactoryDegLex::monomials->add(pos,elem);
        result = elem;
    }
    return result;
}

//-----MonomialFactoryDegRevLexGF2-----------------------------
MonomialFactoryDegRevLexGF2::MonomialFactoryDegRevLexGF2(uint indet)
: MonomialFactoryDegLexGF2()
{
    MonomialFactoryDegLexGF2::excludedIndets = new bool[indet];
    for(uint i=0;i<indet;i++)
        MonomialFactoryDegLexGF2::excludedIndets[i] = false;
    MonomialFactoryDegLex::indet = indet;
    create(0);
}

MonomialFactoryDegRevLexGF2::~MonomialFactoryDegRevLexGF2()
{

}

TAKE_OWN IMonomial* MonomialFactoryDegRevLexGF2::create(uint64_t pos) const
{
    IMonomial* result = (IMonomial*)(MonomialFactoryDegLex::monomials->get(pos));

    if(result == NULL) {
        DegRevLexMonomialGF2* elem = new DegRevLexMonomialGF2(pos,MonomialFactoryDegLex::indet,MonomialFactoryDegLex::monomials,MonomialFactoryDegLexGF2::excludedIndets);
        MonomialFactoryDegLex::monomials->add(pos,elem);
        result = elem;
    }

    return result;
}

} // namespace polynomial
