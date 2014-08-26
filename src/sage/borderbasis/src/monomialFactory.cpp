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
    // at the moment, there is only one type...
    return (IMonomial*)new DegLexMonomial(indet);
}

TAKE_OWN IMonomial* MonomialFactory::create(uint64_t pos, uint indet) const
{
    // at the moment, there is only one type...
    return (IMonomial*)new DegLexMonomial(pos, indet);
}

} // namespace polynomial
