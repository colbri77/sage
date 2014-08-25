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
    // at the moment, there is only one type...
    return new Polynomial<T>(indet);
}

template class PolynomialFactory<uint64_t>;
template class PolynomialFactory<int64_t>;

} // namespace polynomial
