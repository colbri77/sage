#ifndef __POLYNOMIALFACTORY_H__
#define __POLYNOMIALFACTORY_H__

#include "i_polynomial.h"

namespace polynomial {

enum PolType {
	POLTYPE_VECTOR,
	POLTYPE_VECTOR_GF2
};

template<typename T>
class PolynomialFactory
{
public:
    PolynomialFactory(PolType type);
    virtual ~PolynomialFactory();

    virtual TAKE_OWN IPolynomial<T>* create(uint indet) const;

private:
    PolType type;
};

} // namespace polynomial

#endif // __POLYNOMIALFACTORY_H__
