#include "include/compUniverse.h"

#include <stack>
#include <cstring>
#include "include/owningVector.h"
#include "include/degLexMonomial.h"
#include "include/term.h"

namespace borderbasis {

//########## ICompUniverse ################################

template<typename T>
ICompUniverse<T>::ICompUniverse(uint indet)
: indet(indet)
{

}

template<typename T>
ICompUniverse<T>::~ICompUniverse()
{

}

template<typename T>
bool ICompUniverse<T>::contains(IPolynomial<T>* pol) const
{
    for(uint i=0,end_i=pol->size(); i<end_i; i++) {
        if(!contains(pol->at(i)->getMonomial()))
            return false;
    }
    return true;
}

template<typename T>
bool ICompUniverse<T>::contains(IMonomial* monomial) const
{
    return contains(monomial->getPos());
}

template<typename T>
void ICompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions)
{
    add(additions,0);
}

template<typename T>
void ICompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions, uint start)
{
    for(uint i=start,end_i=additions->size();i<end_i;i++) {
        IPolynomial<T>* p = additions->at(i);
        for(uint k=0,end_k=p->size();k<end_k;k++) {
            add(p->at(k)->getMonomial());
        }
    }
}

template<typename T>
void ICompUniverse<T>::add(IMonomial* monomial)
{
    add(monomial->getPos());
}

//########## LinearCompUniverse ############################################

template<typename T>
LinearCompUniverse<T>::LinearCompUniverse(uint indet)
: ICompUniverse<T>(indet),
limit(0)
{

}

template<typename T>
LinearCompUniverse<T>::~LinearCompUniverse()
{

}

template<typename T>
void LinearCompUniverse<T>::clear()
{
    limit = 0;
}

template<typename T>
void LinearCompUniverse<T>::extend(uint64_t limitDegree)
{
    limit = limitDegree;
}

template<typename T>
uint LinearCompUniverse<T>::getMaxDegree() const
{
    return limit;
}

template<typename T>
uint64_t LinearCompUniverse<T>::getMaxPos() const
{
    DegLexMonomial* monomial = new DegLexMonomial(ICompUniverse<T>::indet);
    monomial->set(0,limit);
    uint64_t result = monomial->getPos();
    delete monomial;

    return result;
}

template<typename T>
bool LinearCompUniverse<T>::contains(uint64_t pos) const
{
    DegLexMonomial* monomial = new DegLexMonomial(pos,ICompUniverse<T>::indet);
    uint degree = monomial->getDegree();
    delete monomial;

    return degree<=limit;
}

template<typename T>
bool LinearCompUniverse<T>::contains(IMonomial* monomial) const
{
    return monomial->getDegree()<=limit;
}

template<typename T>
void LinearCompUniverse<T>::add(uint64_t pos)
{
    // since we're only accepting degree hard borders, processing this request might cause confusion.
    NOT_IMPLEMENTED;
}

template<typename T>
void LinearCompUniverse<T>::addBorder()
{
    limit++;
}

template<typename T>
void LinearCompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions,uint start)
{
    // this actually expands the universe at least the degree of the biggest monomial in the set
    for(uint i=start,end_i=additions->size();i<end_i;i++) {
        uint degNew = additions->at(i)->at(0)->getMonomial()->getDegree();
        if(degNew>limit)
            limit = degNew;
    }
}

template class LinearCompUniverse<uint64_t>;
template class LinearCompUniverse<int64_t>;


//########## SpecificCompUniverse ######################################

template<typename T>
SpecificCompUniverse<T>::SpecificCompUniverse(uint indet)
: ICompUniverse<T>(indet),
U(new Polynomial<T>(indet))
{

}

template<typename T>
SpecificCompUniverse<T>::~SpecificCompUniverse()
{
    delete U;
}

template<typename T>
void SpecificCompUniverse<T>::clear()
{
    U->clear();
}

template<typename T>
void SpecificCompUniverse<T>::extend(uint64_t limitDegree)
{
    NOT_IMPLEMENTED;
}

template<typename T>
uint SpecificCompUniverse<T>::getMaxDegree() const
{
    if(U->size()==0)
        return 0;
    return U->at(0)->getMonomial()->getDegree();
}

template<typename T>
uint64_t SpecificCompUniverse<T>::getMaxPos() const
{
    if(U->size()==0)
        return 0;
    return U->at(0)->getMonomial()->getPos();
}

template<typename T>
bool SpecificCompUniverse<T>::contains(uint64_t pos) const
{
    NOT_IMPLEMENTED;
    return NULL;
}

template<typename T>
bool SpecificCompUniverse<T>::contains(IMonomial* monomial) const
{
    for(uint i=0,end=U->size();i<end;i++) {
        if(U->at(i)->getMonomial()->divides(monomial))
            return true;
    }
    return false;
}

template<typename T>
void SpecificCompUniverse<T>::add(uint64_t pos)
{
    NOT_IMPLEMENTED;
}

template<typename T>
void SpecificCompUniverse<T>::addBorder()
{
    OwningVector<IPolynomial<T>*>* newElements = new OwningVector<IPolynomial<T>*>();
    for(uint i=0,end_i=U->size();i<end_i;i++) {
        Polynomial<T>* p = new Polynomial<T>(ICompUniverse<T>::indet);
        for(uint k=0,end_k=U->getIndet();k<end_k;k++) {
            IMonomial* m = U->at(i)->getMonomial()->copy();
            m->extend(k,1);
            p->push(new Term<T>(1,m));
        }
        newElements->push_back(p);
    }
    ICompUniverse<T>::add(newElements);
    delete newElements;
}

template<typename T>
void SpecificCompUniverse<T>::add(IMonomial* monomial)
{
    for(int i=0;i<U->size();i++) {
        if(monomial->divides(U->at(i)->getMonomial())) {
            return;
        }
        else if(U->at(i)->getMonomial()->divides(monomial)) {
            U->remove(i);
            i--;
        }
    }
    U->push(new Term<T>(1,monomial->copy()));
}

template class SpecificCompUniverse<uint64_t>;
template class SpecificCompUniverse<int64_t>;

} // namespace borderbasis
