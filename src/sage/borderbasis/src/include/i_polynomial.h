#ifndef __I_POLYNOMIAL_H__
#define __I_POLYNOMIAL_H__

#include "term.h"
#include "field.h"

namespace polynomial {

template<typename T>
class IPolynomial
{
    public:
        IPolynomial();
        virtual ~IPolynomial();

        virtual TAKE_OWN IPolynomial* copy() const = 0;

        virtual void push(TAKE_OWN Term<T>* term) = 0;
        virtual void push_back(TAKE_OWN Term<T>* term) = 0;
        virtual void remove(uint index) = 0;

        virtual uint getIndet() const = 0;
        virtual uint size() const = 0;
        virtual bool isZero() const = 0;
        virtual void clear() = 0;

        virtual Term<T>* at(uint index) = 0;
        virtual const Term<T>* at(uint index) const = 0;

        virtual void incrementAtIndet(uint index) = 0;

        virtual int compare(const IPolynomial* other) const;

        virtual void hash(uint64_t* out /*size 2*/) const;
        virtual void subtract(IPolynomial<T>* other,IField<T>* f);
};

typedef IPolynomial<uint64_t> IPolynomial_uint64;

} // namespace polynomial

#endif // __I_POLYNOMIAL_H__
