#ifndef __TERM_H__
#define __TERM_H__

#include "i_monomial.h"

namespace polynomial {

template<typename T>
class Term
{
public:
    Term(T coef, TAKE_OWN IMonomial<T>* monomial);
    Term();
    ~Term();

    T getCoef() const;
    void setCoef(T coef);

    IMonomial<T>* getMonomial();
    const IMonomial<T>* getMonomial() const;
    void setMonomial(TAKE_OWN IMonomial<T>* monomial);

    int compare(const Term* other) const;

private:
    IMonomial<T>* monomial;
    T coef;
};

} // namespace polynomial

#endif // __TERM_H__
