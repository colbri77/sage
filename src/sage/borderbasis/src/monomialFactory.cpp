#include "include/monomialFactory.h"

#include "include/degLexMonomial.h"

namespace polynomial {

MonomialFactory::MonomialFactory(MonomialType type)
: type(type)
{

}

MonomialFactory::~MonomialFactory()
{

}

bool MonomialFactory::supportsGetPos() const
{
    IMonomial* tmp = create(1);
    bool result = tmp->supportsGetPos();
    delete tmp;
    return result;
}

TAKE_OWN IMonomial* MonomialFactory::create(uint indet) const
{
    if(type==MONOMIALTYPE_DEGLEX)
        return (IMonomial*)new DegLexMonomial(indet);
    else
        return (IMonomial*)new DegLexMonomialNoOrderPos(indet);
}

TAKE_OWN IMonomial* MonomialFactory::create(uint64_t pos, uint indet) const
{
    if(type==MONOMIALTYPE_DEGLEX) {
        return (IMonomial*)new DegLexMonomial(pos, indet);
    } else {
        NOT_IMPLEMENTED;    // we can't do this with DEGLEX_NO_ORDER_POS
        return NULL;
    }
}

} // namespace polynomial
