#include "include/i_polynomial.h"

#include <stack>
#include "hash/xxhash.h"

namespace polynomial {

template<typename T>
IPolynomial<T>::IPolynomial()
{

}

template<typename T>
IPolynomial<T>::~IPolynomial()
{

}

template<typename T>
int IPolynomial<T>::compare(const IPolynomial* other) const
{
    int result = 0;

    for(uint i=0,k=0;i<size() && k<other->size();i++,k++)
    {
        result = at(i)->compare(other->at(k));
        if(result!=0)
            return result;
    }

    if(size()==other->size())
        return 0;
    return size()-other->size();
}

template<typename T>
void IPolynomial<T>::hash(uint64_t* out) const
{
    void* state[2];
    state[0] = XXH64_init(0xdeadbeef);
    state[1] = XXH64_init(0xbaadf00d);
    uint repSize = size();
    T factor = 0;
    uint64_t monomial = 0;
    char tmpBuf[18];

    for(uint i=0;i<repSize;i++) {
        factor = at(i)->getCoef();
        const IMonomial* cm = at(i)->getMonomial();
        if(cm->supportsGetPos()) {
            monomial = cm->getPos();
        } else {
            uint indet = getIndet();
            char* tmpBuf2 = new char[indet];
            for(uint k=0;k<indet;k++) {
                tmpBuf2[k] = (char)cm->at(k);
            }
            void* stateTmp = XXH64_init(0xdeadbeef);
            XXH64_update(stateTmp,tmpBuf2,indet);
            monomial = XXH64_digest(stateTmp);
            delete tmpBuf2;
        }
        *((T*)tmpBuf) = factor;
        *((uint64_t*)(tmpBuf+sizeof(T)+1)) = monomial;
        tmpBuf[sizeof(T)] = '*';
        tmpBuf[sizeof(T)+9] = '*';
        XXH64_update (state[i%2],tmpBuf,sizeof(T)+10);
    }
    out[0] = (uint64_t)XXH64_digest(state[0]);
    out[1] = (uint64_t)XXH64_digest(state[1]);
}

template<typename T>
void IPolynomial<T>::subtract(IPolynomial<T>* other,IField<T>* f)
{
    uint curPos2 = 0;
    uint limit1 = other->size();
    Term<T>* term2 = at(0);

    for(uint curPos1=0;curPos1<limit1;curPos1++) {
        Term<T>* term1 = other->at(curPos1);
        int cmp = 1;
        for(;cmp>0 && curPos2<size();curPos2++) {
            term2 = at(curPos2);
            cmp = term2->compare(term1);
        }
        if(cmp==0) {
            // we found an equal monomial => subtract
            T newCoef = f->subtract(term2->getCoef(),term1->getCoef());
            if(newCoef==0) {
                curPos2--;
                remove(curPos2);
            } else {
                term2->setCoef(newCoef);
            }
        } else {
            // no equal monomial, just insert -1
            Term<T>* newTerm = new Term<T>();
            newTerm->setCoef(f->subtract(0,1));
            newTerm->setMonomial(term1->getMonomial()->copy());
            push(newTerm);
        }
    }
}

template<typename T>
TAKE_OWN IPolynomial<T>* IPolynomial<T>::getLinearReducible(int* index,const bool* indexMap,const IField<T>* f) const
{
    return NULL; // aka don't try to reduce
}

template<typename T>
void IPolynomial<T>::multiply(uint64_t constant, const IField<T>* f)
{
    for(uint i=0;i<size();i++) {
        at(i)->setCoef(f->multiply(at(i)->getCoef(),constant));
    }
}

template<typename T>
void IPolynomial<T>::substitute(int indet,const IPolynomial<T>* replacement,const IField<T>* f)
{
    NOT_IMPLEMENTED;
}

template class IPolynomial<uint64_t>;
template class IPolynomial<int64_t>;

} // namespace polynomial
