#include "include/i_monomial.h"

#include <cstring>

namespace polynomial {

template<typename T>
IMonomial<T>::IMonomial(IMonomial::termOrder termOrdering)
: termOrdering(termOrdering)
{

}

template<typename T>
IMonomial<T>::~IMonomial()
{

}

template<typename T>
int IMonomial<T>::compare(const IMonomial<T>* other) const
{
    ENSURE(termOrdering == other->termOrdering, "IMonomial::compare(other): monomials have different term ordering");
    ENSURE(getIndet() == other->getIndet(), "IMonomial::compare(other): monomials have a different amount of variables");
    // pos is uint64_t => simple subtraction might cause an overflow
    return getPos()>other->getPos() ? 1 :
            (getPos()<other->getPos() ? -1 : 0);
}

template class IMonomial<uint64_t>;
template class IMonomial<int64_t>;

} // namespace polynomial
