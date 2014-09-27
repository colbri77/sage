#include "include/polynomial.h"

#include "include/owningVector.h"

#include <stack>

namespace polynomial {

//-----Polynomial<T>-----------------------------------

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
TAKE_OWN IPolynomial<T>* Polynomial<T>::copy() const
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
    if(term->getCoef()==0) {
        delete term;
        return;
    }

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
void Polynomial<T>::remove(uint index)
{
    typename vector<Term<T>*>::iterator it=rep->begin();
    it += index;
    delete rep->at(index);
    rep->erase(it);
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
        IMonomial* m = rep->at(i)->getMonomial()->extend(index,1);
        rep->at(i)->setMonomial(m);
    }
}

template class Polynomial<uint64_t>;
template class Polynomial<int64_t>;

//-----PolynomialGF2------------------------------------
template<typename T>
PolynomialGF2<T>::PolynomialGF2(uint indet)
: Polynomial<T>(indet)
{

}

template<typename T>
PolynomialGF2<T>::~PolynomialGF2()
{

}

template<typename T>
void PolynomialGF2<T>::incrementAtIndet(uint index)
{
    Polynomial<T>::incrementAtIndet(index);

    for(uint i=0;i<Polynomial<T>::rep->size();i++) {
        for(uint k=i+1;k<Polynomial<T>::rep->size();k++) {
            int cmpResult = Polynomial<T>::rep->at(i)->getMonomial()->compare(Polynomial<T>::rep->at(k)->getMonomial());
            if(cmpResult==0) {
                Polynomial<T>::rep->remove(i);
                Polynomial<T>::rep->remove(k-1);
                i--;
                break;
            }
            else if(cmpResult<0) {
                IMonomial* monTmp = Polynomial<T>::rep->at(i)->getMonomial();
                Polynomial<T>::rep->at(i)->setMonomial(Polynomial<T>::rep->at(k)->getMonomial());
                Polynomial<T>::rep->at(k)->setMonomial(monTmp);
                i--;
                break;
            }
        }
    }
}

template<typename T>
TAKE_OWN IPolynomial<T>* PolynomialGF2<T>::copy() const
{
    PolynomialGF2<T>* result = new PolynomialGF2<T>(Polynomial<T>::indet);
    for(uint i=0;i<Polynomial<T>::rep->size();i++) {
        Term<T>* term = new Term<T>();
        term->setCoef(Polynomial<T>::rep->at(i)->getCoef());
        term->setMonomial(Polynomial<T>::rep->at(i)->getMonomial()->copy());
        result->rep->push_back(term);
    }
    return result;
}

template<typename T>
void PolynomialGF2<T>::push(TAKE_OWN Term<T>* term)
{
    if(term->getCoef()==0) {
        delete term;
        return;
    }

    int posCtr = 0;
    for(typename vector<Term<T>*>::iterator it=Polynomial<T>::rep->begin();it!=Polynomial<T>::rep->end();it++,posCtr++) {
        int cmp = (*it)->compare(term);
        if(cmp<0) {
            Polynomial<T>::rep->insert(it,term);
            return;
        }
        else if(cmp==0) {
            Polynomial<T>::rep->remove(posCtr);
            return;
        }
    }
    Polynomial<T>::rep->insert(Polynomial<T>::rep->end(),term);
}

template<typename T>
TAKE_OWN IPolynomial<T>* PolynomialGF2<T>::getLinearReducible(int* index,const bool* indexMap,const IField<T>* f) const
{
    *index = -1;
    Term<T>* t = NULL;
    uint i = 0;
    for(uint i_end=Polynomial<T>::rep->size();i<i_end;i++) {
        t = Polynomial<T>::rep->at(i);
        IMonomial* m = t->getMonomial();
        if(m->getDegree()>1)
            return NULL;
        int varIndex = m->getDegree()==1 ? (int)m->getLV() : -1;
        if(varIndex!=-1 && !indexMap[varIndex]) {
            *index = varIndex;
            break;
        }
    }
    if((*index) == -1)
        return NULL;

    IPolynomial<T>* result = copy();
    result->remove(i);

    return result;
}

template<typename T>
void PolynomialGF2<T>::substitute(int indet,const IPolynomial<T>* replacement,const IField<T>* f)
{
    stack<Term<T>*> _stack = stack<Term<T>*>();
    for(uint i=0;i<Polynomial<T>::rep->size();i++) {
        Term<T>* t = Polynomial<T>::rep->at(i);
        if(t->getMonomial()->at(indet)>0) {
            Term<T>* tNew = new Term<T>(1,t->getMonomial()->copy()->set(indet,0));

            // actual substitution
            for(uint k=0,k_end=replacement->size();k<k_end;k++) {
                Term<T>* r = new Term<T>();
                const IMonomial* rMonom = replacement->at(k)->getMonomial();
                r->setCoef(1);
                IMonomial* m = tNew->getMonomial()->copy();
                IMonomial* mTmp = NULL;
                for(uint d=0;d<this->indet;d++) {
                    if(rMonom->at(d)>0) {
                        mTmp = m;
                        m = m->extend(d,1);
                        if(m != mTmp)
                            mTmp->del();
                    }
                }
                r->setMonomial(m);
                _stack.push(r);
            }

            delete tNew;

            Polynomial<T>::remove(i);
            i--;
        }
    }

    // now append new Terms to the polynomial
    while(!_stack.empty()) {
        Term<T>* t = _stack.top();
        _stack.pop();
        push(t);
    }
}

template class PolynomialGF2<uint64_t>;
template class PolynomialGF2<int64_t>;


} // namespace polynomial
