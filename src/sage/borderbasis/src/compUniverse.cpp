#include "include/compUniverse.h"

#include <stack>
#include <cstring>
#include "include/owningVector.h"
#include "include/degLexMonomial.h"

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
bool ICompUniverse<T>::contains(IMonomial<T>* monomial) const
{
    return contains(monomial->getPos());
}

template<typename T>
void ICompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions,uint start)
{
    for(uint i=0,end_i=additions->size(); i<end_i; i++) {
        IPolynomial<T>* pol = additions->at(i);
        for(uint k=0,end_k=pol->size(); k<end_k; k++) {
            add(pol->at(k)->getMonomial()->getPos());
        }
    }
}

template<typename T>
void ICompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions)
{
    add(additions,0);
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
bool LinearCompUniverse<T>::contains(IMonomial<T>* monomial) const
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
U(new uint8_t[0]),
uLen(0),
uBlocks(0),
lastUBorderCandidates(new OwningVector<IMonomial<T>*>())
{

}

template<typename T>
SpecificCompUniverse<T>::~SpecificCompUniverse()
{
    delete U;
    delete lastUBorderCandidates;
}

template<typename T>
void SpecificCompUniverse<T>::clear()
{
    delete U;
    U = new uint8_t[0];
    uLen = 0;
    uBlocks = 0;
    lastUBorderCandidates->clear();
}

template<typename T>
void SpecificCompUniverse<T>::extend(uint64_t limitDegree)
{
    // construct a monomial at the edge
    DegLexMonomial* monomial = new DegLexMonomial(ICompUniverse<T>::indet);
    monomial->set(0,limitDegree);
    uint64_t pos = monomial->getPos();
    delete monomial;

    for(; !contains(pos) && pos>0; pos--) {
        add(pos);
    }
}

template<typename T>
uint SpecificCompUniverse<T>::getMaxDegree() const
{
    uint64_t pos = uLen;
    for(; !contains(pos) && pos>0; pos--);
    DegLexMonomial* monomial = new DegLexMonomial(pos,ICompUniverse<T>::indet);
    uint result = monomial->getDegree();
    delete monomial;

    return result;
}

template<typename T>
uint64_t SpecificCompUniverse<T>::getMaxPos() const
{
    return uLen-1;
}

template<typename T>
bool SpecificCompUniverse<T>::contains(uint64_t pos) const
{
    return (uLen>pos) && ((U[(pos)>>3]&(1<<((pos)&7)))!=0);
}

template<typename T>
void SpecificCompUniverse<T>::add(uint64_t pos)
{
    if(pos>(uLen+1)) {
        // plus 2 is not specifically calculated - better safe than sorry.
        uint newUBlocks = (pos+9)/8+2;
        uint8_t* uNew = new uint8_t[newUBlocks]();
        memcpy(uNew,U,uBlocks);
        uBlocks = newUBlocks;
        uLen = pos+1;
        delete U;
        U = uNew;
    }
    (U[(pos)>>3]|=(1<<((pos)&7)));
}

template<typename T>
void SpecificCompUniverse<T>::addBorder()
{
    stack<IMonomial<T>*> _stack = stack<IMonomial<T>*>();
    uint64_t maxPos = 0;

    // fill the stack with border candidates
    while(lastUBorderCandidates->size()>0) {
        IMonomial<T>* t = lastUBorderCandidates->pop_lift();
        for(uint i=0;i<ICompUniverse<T>::indet;i++) {
            IMonomial<T>* tNew = t->copy();
            tNew->extend(i,1);
            if(tNew->getPos()>maxPos)
                maxPos = tNew->getPos();
            _stack.push(tNew);
        }
        delete t;
    }

    // extend by the maxima first, otherwise we might have inefficient many memory reallocations
    add(maxPos);

    // finally, add the remaining monomials
    while(!_stack.empty()) {
        IMonomial<T>* t = _stack.top();
        uint64_t pos = t->getPos();
        _stack.pop();
        if(!contains(pos)) {
            add(pos);
            lastUBorderCandidates->push_back(t);
        } else {
            delete t;
        }
    }
}

template<typename T>
void SpecificCompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions,uint start)
{
    stack<IMonomial<T>*> _stack = stack<IMonomial<T>*>();

    // fill stack with extension candidates
    uint64_t maxPos = 0;
    for(uint iPol=start,ts=additions->size();iPol<ts;iPol++) {
        IPolynomial<T>* p = additions->at(iPol);
        for(uint iMonomial=0,ts2=p->size();iMonomial<ts2;iMonomial++) {
            IMonomial<T>* t = p->at(iMonomial)->getMonomial();
            uint64_t pos = t->getPos();
            if(pos>maxPos)
                maxPos = pos;
            _stack.push(t->copy());
        }
    }

    // extend by the maxima first, otherwise we might have inefficient many memory reallocations
    add(maxPos);

    // finally, add the remaining monomials and their divisors
    while(!_stack.empty()) {
        IMonomial<T>* t = _stack.top();
        uint64_t pos = t->getPos();
        _stack.pop();

        if(!contains(pos)) {
            add(pos);
            lastUBorderCandidates->push_back(t->copy());
            for(uint i=0;i<ICompUniverse<T>::indet;i++) {
                if(t->at(i)>0) {
                    IMonomial<T>* tNew = t->copy();
                    tNew->extend(i,-1);
                    _stack.push(tNew);
                }
            }
        }

        delete t;
    }
}

template class SpecificCompUniverse<uint64_t>;
template class SpecificCompUniverse<int64_t>;


//########## SpecificCompUniverseNoBorderLog ###############################

template<typename T>
SpecificCompUniverseNoBorderLog<T>::SpecificCompUniverseNoBorderLog(uint indet)
: SpecificCompUniverse<T>(indet)
{

}

template<typename T>
SpecificCompUniverseNoBorderLog<T>::~SpecificCompUniverseNoBorderLog()
{

}

template<typename T>
void SpecificCompUniverseNoBorderLog<T>::addBorder()
{
    // this class is more efficient when adding polynomials by ignoring the need to take care of the border.
    // This means, of course, that the border can't be used => Exception
    NOT_IMPLEMENTED;
}

template<typename T>
void SpecificCompUniverseNoBorderLog<T>::add(IOwningList<IPolynomial<T>*>* additions,uint start)
{
    stack<IMonomial<T>*> _stack = stack<IMonomial<T>*>();

    // fill stack with extension candidates
    uint64_t maxPos = 0;
    for(uint iPol=start,ts=additions->size();iPol<ts;iPol++) {
        IPolynomial<T>* p = additions->at(iPol);
        for(uint iMonomial=0,ts2=p->size();iMonomial<ts2;iMonomial++) {
            IMonomial<T>* t = p->at(iMonomial)->getMonomial();
            uint64_t pos = t->getPos();
            if(pos>maxPos)
                maxPos = pos;
            _stack.push(t->copy());
        }
    }

    // extend by the maxima first, otherwise we might have inefficient many memory reallocations
    SpecificCompUniverse<T>::add(maxPos);

    // finally, add the remaining monomials and their divisors
    while(!_stack.empty()) {
        IMonomial<T>* t = _stack.top();
        uint64_t pos = t->getPos();
        _stack.pop();

        if(!SpecificCompUniverse<T>::contains(pos)) {
            SpecificCompUniverse<T>::add(pos);
            for(uint i=0;i<ICompUniverse<T>::indet;i++) {
                if(t->at(i)>0) {
                    IMonomial<T>* tNew = t->copy();
                    tNew->extend(i,-1);
                    _stack.push(tNew);
                }
            }
        }

        delete t;
    }
}

template class SpecificCompUniverseNoBorderLog<uint64_t>;
template class SpecificCompUniverseNoBorderLog<int64_t>;

} // namespace borderbasis
