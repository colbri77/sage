#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include "owningVector.h"
#include "i_polynomial.h"

namespace polynomial {

template<typename T>
class Polynomial : public IPolynomial<T>
{
    public:
        Polynomial(uint indet);
        virtual ~Polynomial();

        virtual TAKE_OWN Polynomial<T>* copy() const OVERRIDE;

        virtual void push(TAKE_OWN Term<T>* term) OVERRIDE;
        virtual void push_back(TAKE_OWN Term<T>* term) OVERRIDE;

        virtual uint getIndet() const OVERRIDE;
        virtual uint size() const OVERRIDE;
        virtual bool isZero() const OVERRIDE;
        virtual void clear() OVERRIDE;

        virtual Term<T>* at(uint index) OVERRIDE;
        virtual const Term<T>* at(uint index) const OVERRIDE;

        virtual void incrementAtIndet(uint index) OVERRIDE;

    private:
        uint indet;
        OwningVector<Term<T>*>* rep;
};

} // namespace polynomial

#endif // __POLYNOMIAL_H__
