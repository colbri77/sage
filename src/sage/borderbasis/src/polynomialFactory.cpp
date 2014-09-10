#include "include/polynomialFactory.h"

#include "include/polynomial.h"

namespace polynomial {

template<typename T>
PolynomialFactory<T>::PolynomialFactory(PolType type)
: type(type)
{

}

template<typename T>
PolynomialFactory<T>::~PolynomialFactory()
{

}

template<typename T>
TAKE_OWN IPolynomial<T>* PolynomialFactory<T>::create(uint indet) const
{
    if(type==POLTYPE_VECTOR)
        return new Polynomial<T>(indet);
    else
        return new PolynomialGF2<T>(indet);
}

template class PolynomialFactory<uint64_t>;
template class PolynomialFactory<int64_t>;

} // namespace polynomial
