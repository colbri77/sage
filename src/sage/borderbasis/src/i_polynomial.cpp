#include "include/i_polynomial.h"

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
uint64_t IPolynomial<T>::hash() const
{
    void* state = XXH32_init(0xdeadbeef);
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
            char* tmpBuf2 = new char[getIndet()];
            uint indet = getIndet();
            for(uint k=0;k<indet;k++) {
                tmpBuf2[k] = (char)cm->at(k);
            }
            void* state2 = XXH32_init(0xdeadbeef);
            XXH32_update(state2,tmpBuf2,getIndet());
            monomial = XXH32_digest(state2);
            delete tmpBuf2;
        }
        *((T*)tmpBuf) = factor;
        *((uint64_t*)(tmpBuf+sizeof(T)+1)) = monomial;
        tmpBuf[sizeof(T)] = '*';
        tmpBuf[sizeof(T)+9] = '*';
        XXH32_update (state,tmpBuf,sizeof(T)+10);
    }
    return (uint64_t)XXH32_digest (state);
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

template class IPolynomial<uint64_t>;
template class IPolynomial<int64_t>;

} // namespace polynomial
