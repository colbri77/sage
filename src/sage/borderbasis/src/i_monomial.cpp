#include "include/i_monomial.h"

#include <cstring>

namespace polynomial {

IMonomial::IMonomial(IMonomial::termOrder termOrdering)
: termOrdering(termOrdering)
{

}

IMonomial::~IMonomial()
{

}

int IMonomial::compare(const IMonomial* other) const
{
    ENSURE(termOrdering == other->termOrdering, "IMonomial::compare(other): monomials have different term ordering");
    ENSURE(getIndet() == other->getIndet(), "IMonomial::compare(other): monomials have a different amount of variables");
    // pos is uint64_t => simple subtraction might cause an overflow
    return getPos()>other->getPos() ? 1 :
            (getPos()<other->getPos() ? -1 : 0);
}

void IMonomial::del()
{
    delete this;
}


} // namespace polynomial
