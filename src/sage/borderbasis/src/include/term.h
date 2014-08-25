#ifndef __TERM_H__
#define __TERM_H__

#include "i_monomial.h"

namespace polynomial {

template<typename T>
class Term
{
public:
    Term(T coef, TAKE_OWN IMonomial* monomial);
    Term();
    ~Term();

    T getCoef() const;
    void setCoef(T coef);

    IMonomial* getMonomial();
    const IMonomial* getMonomial() const;
    void setMonomial(TAKE_OWN IMonomial* monomial);

    int compare(const Term* other) const;

private:
    IMonomial* monomial;
    T coef;
};

} // namespace polynomial

#endif // __TERM_H__
