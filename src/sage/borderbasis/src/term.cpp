#include "include/term.h"

namespace polynomial {

template<typename T>
Term<T>::Term(T coef, TAKE_OWN IMonomial* monomial)
: monomial(monomial),
coef(coef)
{

}

template<typename T>
Term<T>::Term()
: monomial(NULL),
coef(0)
{

}

template<typename T>
Term<T>::~Term()
{
    if(monomial!=NULL)
        monomial->del();
}

template<typename T>
T Term<T>::getCoef() const
{
    return coef;
}

template<typename T>
void Term<T>::setCoef(T coef)
{
    this->coef = coef;
}

template<typename T>
IMonomial* Term<T>::getMonomial()
{
    return monomial;
}

template<typename T>
const IMonomial* Term<T>::getMonomial() const
{
    return monomial;
}

template<typename T>
void Term<T>::setMonomial(TAKE_OWN IMonomial* monomial)
{
    this->monomial = monomial;
}

template<typename T>
int Term<T>::compare(const Term* other) const
{
    int result = 0;
    result = monomial->compare(other->monomial);
    if(result!=0)
        return result;
    result = coef-other->coef;
    return result;
}

template class Term<uint64_t>;
template class Term<int64_t>;
template class Term<uint8_t>;

} // namespace polynomial
