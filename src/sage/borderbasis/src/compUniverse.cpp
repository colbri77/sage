#include "include/compUniverse.h"

#include <stack>
#include <cstring>
#include "include/owningVector.h"
#include "include/degLexMonomial.h"

namespace borderbasis {

//########## ICompUniverse ################################

template<typename T>
ICompUniverse<T>::ICompUniverse(uint indet)
: indet(indet),
exclusions(new vector<IMonomial*>())
{

}

template<typename T>
ICompUniverse<T>::~ICompUniverse()
{
    for(uint i=0;i<exclusions->size();i++) {
        exclusions->at(i)->del();
    }
    delete exclusions;
}

template<typename T>
bool ICompUniverse<T>::exclude(IMonomial* monomial)
{
    for(uint i=0;i<exclusions->size();i++) {
        if(exclusions->at(i)->divides(monomial))
            return false;
        if(monomial->divides(exclusions->at(i))) {
            exclusions->at(i)->del();
            exclusions->at(i) = monomial;
            return true;
        }
    }
    exclusions->push_back(monomial->copy());
    return true;
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
    bool result = contains(monomial->getPos());
    if(result) {
        for(uint i=0,i_end=exclusions->size();i<i_end;i++) {
            if(exclusions->at(i)->divides(monomial) && exclusions->at(i)->compare(monomial)!=0) {
                return false;
            }
        }
    }
    return result;
}

template<typename T>
void ICompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions,uint start)
{
    for(uint i=0,end_i=additions->size(); i<end_i; i++) {
        IPolynomial<T>* pol = additions->at(i);
        for(uint k=0,end_k=pol->size(); k<end_k; k++) {
            add(pol->at(k)->getMonomial());
        }
    }
}

template<typename T>
void ICompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions)
{
    add(additions,0);
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
bool LinearCompUniverse<T>::contains(uint64_t pos) const
{
    NOT_IMPLEMENTED;
    return false;
}

template<typename T>
bool LinearCompUniverse<T>::contains(IMonomial* monomial) const
{
     bool result = monomial->getDegree()<=limit;
     if(result) {
        for(uint i=0,i_end=ICompUniverse<T>::exclusions->size();i<i_end;i++) {
            if(ICompUniverse<T>::exclusions->at(i)->divides(monomial) && ICompUniverse<T>::exclusions->at(i)->compare(monomial)!=0) {
                return false;
            }
        }
    }
    return result;
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

template<typename T>
bool LinearCompUniverse<T>::beyondLastElement(IMonomial* monomial) const
{
    return monomial->getDegree()>limit;
}

template class LinearCompUniverse<uint64_t>;
template class LinearCompUniverse<int64_t>;


//########## SpecificCompUniverse ######################################

template<typename T>
SpecificCompUniverse<T>::SpecificCompUniverse(uint indet)
: ICompUniverse<T>(indet),
U(new uint8_t[1]),
uLen(0),
uBlocks(0),
maxDegree(0),
lastUBorderCandidates(new vector<IMonomial*>())
{

}

template<typename T>
SpecificCompUniverse<T>::~SpecificCompUniverse()
{
    delete U;
    for(uint i=0;i<lastUBorderCandidates->size();i++)
        lastUBorderCandidates->at(i)->del();
    delete lastUBorderCandidates;
}

template<typename T>
void SpecificCompUniverse<T>::clear()
{
    delete U;
    U = new uint8_t[1];
    uLen = 0;
    uBlocks = 0;
    maxDegree = 0;
    lastUBorderCandidates->clear();
}

template<typename T>
void SpecificCompUniverse<T>::extend(uint64_t limitDegree)
{
    NOT_IMPLEMENTED;
}

template<typename T>
uint SpecificCompUniverse<T>::getMaxDegree() const
{
    return maxDegree;
}

template<typename T>
bool SpecificCompUniverse<T>::contains(uint64_t pos) const
{
    return (uLen>pos) && ((U[(pos)>>3]&(1<<((pos)&7)))!=0);
}

template<typename T>
void SpecificCompUniverse<T>::add(uint64_t pos)
{
    if(pos>=uLen) {
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
    stack<IMonomial*> _stack = stack<IMonomial*>();
    uint64_t maxPos = 0;

    // fill the stack with border candidates
    while(lastUBorderCandidates->size()>0) {
        IMonomial* t = lastUBorderCandidates->back();
        lastUBorderCandidates->pop_back();
        for(uint i=0;i<ICompUniverse<T>::indet;i++) {
            IMonomial* tNew = t->copy();
            tNew = tNew->extend(i,1);

            bool excluded = false;
            for(uint i=0;i<ICompUniverse<T>::exclusions->size() && !excluded;i++) {
                excluded |= ICompUniverse<T>::exclusions->at(i)->divides(tNew);
                excluded |= (ICompUniverse<T>::exclusions->at(i)->compare(tNew)==0);
            }

            if(!excluded) {
                if(tNew->getPos()>maxPos) {
                    maxPos = tNew->getPos();
                    maxDegree = tNew->getDegree();
                }
                _stack.push(tNew);
            } else {
                tNew->del();
            }
        }
        t->del();
    }

    // extend by the maxima first, otherwise we might have inefficient many memory reallocations
    add(maxPos);

    // finally, add the remaining monomials
    while(!_stack.empty()) {
        IMonomial* t = _stack.top();
        uint64_t pos = t->getPos();
        _stack.pop();
        if(!contains(pos) || pos==maxPos) {
            // only allow children of maxPos once
            if(pos==maxPos) maxPos++;
            add(pos);
            lastUBorderCandidates->push_back(t);
        } else {
            t->del();
        }
    }
}

template<typename T>
void SpecificCompUniverse<T>::add(IOwningList<IPolynomial<T>*>* additions,uint start)
{
    stack<IMonomial*> _stack = stack<IMonomial*>();

    // fill stack with extension candidates
    uint64_t maxPos = 0;
    for(uint iPol=start,ts=additions->size();iPol<ts;iPol++) {
        IPolynomial<T>* p = additions->at(iPol);
        for(uint iMonomial=0,ts2=p->size();iMonomial<ts2;iMonomial++) {
            IMonomial* t = p->at(iMonomial)->getMonomial();
            uint64_t pos = t->getPos();
            if(pos>maxPos) {
                maxPos = pos;
                if(maxDegree < t->getDegree())
                    maxDegree = t->getDegree();
            }
            _stack.push(t->copy());
        }
    }

    // extend by the maxima first, otherwise we might have inefficient many memory reallocations
    add(maxPos);

    // finally, add the remaining monomials and their divisors
    while(!_stack.empty()) {
        IMonomial* t = _stack.top();
        uint64_t pos = t->getPos();
        _stack.pop();

        if(!contains(pos) || pos==maxPos) {
            // only allow children of maxPos once
            if(pos==maxPos) maxPos++;
            add(pos);
            lastUBorderCandidates->push_back(t->copy());
            for(uint i=0;i<ICompUniverse<T>::indet;i++) {
                if(t->at(i)>0) {
                    IMonomial* tNew = t->copy();
                    tNew = tNew->extend(i,-1);
                    _stack.push(tNew);
                }
            }
        }

        t->del();
    }
}

template<typename T>
void SpecificCompUniverse<T>::add(IMonomial* monomial)
{
    if(monomial->getDegree()>maxDegree)
        maxDegree = monomial->getDegree();
    add(monomial->getPos());
}

template<typename T>
bool SpecificCompUniverse<T>::beyondLastElement(IMonomial* monomial) const
{
    return monomial->getPos()>uLen;
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
    stack<IMonomial*> _stack = stack<IMonomial*>();
    IMonomial* t = NULL;

    // fill stack with extension candidates
    uint64_t maxPos = 0;
    for(uint iPol=start,ts=additions->size();iPol<ts;iPol++) {
        IPolynomial<T>* p = additions->at(iPol);
        for(uint iMonomial=0,ts2=p->size();iMonomial<ts2;iMonomial++) {
            t = p->at(iMonomial)->getMonomial();
            uint64_t pos = t->getPos();
            if(pos>maxPos) {
                maxPos = pos;
                if(SpecificCompUniverse<T>::maxDegree < t->getDegree())
                    SpecificCompUniverse<T>::maxDegree = t->getDegree();
            }
            _stack.push(t->copy());
        }
    }

    // extend by the maxima first, otherwise we might have inefficient many memory reallocations
    SpecificCompUniverse<T>::add(maxPos);

    // finally, add the remaining monomials and their divisors
    while(!_stack.empty()) {
        t = _stack.top();
        uint64_t pos = t->getPos();
        _stack.pop();

        if(!SpecificCompUniverse<T>::contains(pos) || pos==maxPos) {
             // only allow children of maxPos once
            if(pos==maxPos) maxPos++;
            SpecificCompUniverse<T>::add(pos);
            for(uint i=0;i<ICompUniverse<T>::indet;i++) {
                if(t->at(i)>0) {
                    IMonomial* tNew = t->copy();
                    tNew = tNew->extend(i,-1);
                    _stack.push(tNew);
                }
            }
        }

        t->del();
    }
}

template class SpecificCompUniverseNoBorderLog<uint64_t>;
template class SpecificCompUniverseNoBorderLog<int64_t>;


//########## SpecificCompUniverseNoOrderPos ######################################

template<typename T>
SpecificCompUniverseNoOrderPos<T>::SpecificCompUniverseNoOrderPos(uint indet)
: ICompUniverse<T>(indet),
U(new Polynomial<T>(indet))
{

}

template<typename T>
SpecificCompUniverseNoOrderPos<T>::~SpecificCompUniverseNoOrderPos()
{
    delete U;
}

template<typename T>
void SpecificCompUniverseNoOrderPos<T>::clear()
{
    U->clear();
}

template<typename T>
void SpecificCompUniverseNoOrderPos<T>::extend(uint64_t limitDegree)
{
    NOT_IMPLEMENTED;
}

template<typename T>
uint SpecificCompUniverseNoOrderPos<T>::getMaxDegree() const
{
    if(U->size()==0)
        return 0;
    return U->at(0)->getMonomial()->getDegree();
}

template<typename T>
bool SpecificCompUniverseNoOrderPos<T>::contains(uint64_t pos) const
{
    NOT_IMPLEMENTED;
    return NULL;
}

template<typename T>
bool SpecificCompUniverseNoOrderPos<T>::contains(IMonomial* monomial) const
{
    bool result = false;
    for(uint i=0,end=U->size();i<end;i++) {
        if(monomial->divides(U->at(i)->getMonomial()))
            result = true;
            break;
    }
    if(result) {
        for(uint i=0,i_end=ICompUniverse<T>::exclusions->size();i<i_end;i++) {
            if(ICompUniverse<T>::exclusions->at(i)->divides(monomial) && ICompUniverse<T>::exclusions->at(i)->compare(monomial)!=0)
                return false;
        }
    }
    return result;
}

template<typename T>
void SpecificCompUniverseNoOrderPos<T>::add(uint64_t pos)
{
    NOT_IMPLEMENTED;
}

template<typename T>
void SpecificCompUniverseNoOrderPos<T>::addBorder()
{
    OwningVector<IPolynomial<T>*>* newElements = new OwningVector<IPolynomial<T>*>();
    for(uint i=0,end_i=U->size();i<end_i;i++) {
        Polynomial<T>* p = new Polynomial<T>(ICompUniverse<T>::indet);
        for(uint k=0,end_k=U->getIndet();k<end_k;k++) {
            IMonomial* m = U->at(i)->getMonomial()->copy();
            m = m->extend(k,1);

            bool excluded = false;
            for(uint i=0;i<ICompUniverse<T>::exclusions->size() && !excluded;i++) {
                excluded |= ICompUniverse<T>::exclusions->at(i)->divides(m);
                excluded |= (ICompUniverse<T>::exclusions->at(i)->compare(m)==0);
            }

            if(!excluded)
                p->push(new Term<T>(1,m));
            else
                m->del();
        }
        newElements->push_back(p);
    }
    ICompUniverse<T>::add(newElements);
    delete newElements;
}

template<typename T>
void SpecificCompUniverseNoOrderPos<T>::add(IMonomial* monomial)
{
    for(uint i=0;i<U->size();i++) {
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

template<typename T>
bool SpecificCompUniverseNoOrderPos<T>::beyondLastElement(IMonomial* monomial) const
{
    for(uint i=0,end=U->size();i<end;i++) {
        if(monomial->divides(U->at(i)->getMonomial()))
            return false;
    }
    return true;
}

template class SpecificCompUniverseNoOrderPos<uint64_t>;
template class SpecificCompUniverseNoOrderPos<int64_t>;

} // namespace borderbasis
