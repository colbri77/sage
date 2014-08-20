#include "include/monomialFactory.h"

#include "include/degLexMonomial.h"

namespace polynomial {

template<typename T>
MonomialFactory<T>::MonomialFactory(MonomialType type)
: type(type)
{

}

template<typename T>
MonomialFactory<T>::~MonomialFactory()
{

}

template<typename T>
TAKE_OWN IMonomial<T>* MonomialFactory<T>::create(uint indet) const
{
    // at the moment, there is only one type...
    return (IMonomial<T>*)new DegLexMonomial(indet);
}

template<typename T>
TAKE_OWN IMonomial<T>* MonomialFactory<T>::create(uint64_t pos, uint indet) const
{
    // at the moment, there is only one type...
    return (IMonomial<T>*)new DegLexMonomial(pos, indet);
}

template class MonomialFactory<uint64_t>;
template class MonomialFactory<int64_t>;

} // namespace polynomial
