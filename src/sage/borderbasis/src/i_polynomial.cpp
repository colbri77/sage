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
        monomial = at(i)->getMonomial()->getPos();
        *((T*)tmpBuf) = factor;
        *((uint64_t*)(tmpBuf+sizeof(T)+1)) = monomial;
        tmpBuf[sizeof(T)] = '*';
        tmpBuf[sizeof(T)+9] = '*';
        XXH32_update (state,tmpBuf,sizeof(T)+10);
    }
    return (uint64_t)XXH32_digest (state);
}

template class IPolynomial<uint64_t>;
template class IPolynomial<int64_t>;

} // namespace polynomial
