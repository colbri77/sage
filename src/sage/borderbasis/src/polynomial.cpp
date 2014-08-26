#include "include/polynomial.h"

#include "include/owningVector.h"

namespace polynomial {

template<typename T>
Polynomial<T>::Polynomial(uint indet)
: indet(indet),
rep(new OwningVector<Term<T>*>())
{

}

template<typename T>
Polynomial<T>::~Polynomial()
{
    delete rep;
}

template<typename T>
TAKE_OWN Polynomial<T>* Polynomial<T>::copy() const
{
    Polynomial<T>* result = new Polynomial<T>(indet);
    for(uint i=0;i<rep->size();i++) {
        Term<T>* term = new Term<T>();
        term->setCoef(rep->at(i)->getCoef());
        term->setMonomial(rep->at(i)->getMonomial()->copy());
        result->rep->push_back(term);
    }
    return result;
}

template<typename T>
void Polynomial<T>::push(TAKE_OWN Term<T>* term)
{
    if(term->getCoef()==0)
        return;

    for(typename vector<Term<T>*>::iterator it=rep->begin();it!=rep->end();it++) {
        if((*it)->compare(term)<=0) {
            rep->insert(it,term);
            return;
        }
    }
    rep->insert(rep->end(),term);
}

template<typename T>
void Polynomial<T>::push_back(TAKE_OWN Term<T>* term)
{
    if(term->getCoef()==0)
        return;

    rep->push_back(term);
}

template<typename T>
uint Polynomial<T>::getIndet() const
{
    return indet;
}

template<typename T>
uint Polynomial<T>::size() const
{
    return rep->size();
}

template<typename T>
bool Polynomial<T>::isZero() const
{
    for(uint i=0;i<rep->size();i++) {
        if(rep->at(i)->getCoef()>0)
            return false;
    }

    return true;
}

template<typename T>
void Polynomial<T>::clear()
{
    rep->clear();
}

template<typename T>
Term<T>* Polynomial<T>::at(uint index)
{
    return rep->at(index);
}

template<typename T>
const Term<T>* Polynomial<T>::at(uint index) const
{
    return (const Term<T>*)(rep->at(index));
}

template<typename T>
void Polynomial<T>::incrementAtIndet(uint index)
{
    for(uint i=0;i<rep->size();i++) {
        rep->at(i)->getMonomial()->extend(index,1);
    }
}

template<typename T>
void Polynomial<T>::remove(uint index)
{
    typename vector<Term<T>*>::iterator it=rep->begin();
    it += index;
    delete rep->at(index);
    rep->erase(it);
}

template class Polynomial<uint64_t>;
template class Polynomial<int64_t>;

} // namespace polynomial
