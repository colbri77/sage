#include "include/i_monomial.h"

#include <cstring>

namespace polynomial {

IMonomial::IMonomial(IMonomial::termOrder termOrdering)
: termOrdering(termOrdering)
{

}

IMonomial::~IMonomial()
{

}

int IMonomial::compare(const IMonomial* other) const
{
    ENSURE(termOrdering == other->termOrdering, "IMonomial::compare(other): monomials have different term ordering");
    ENSURE(getIndet() == other->getIndet(), "IMonomial::compare(other): monomials have a different amount of variables");
    // pos is uint64_t => simple subtraction might cause an overflow
    return getPos()>other->getPos() ? 1 :
            (getPos()<other->getPos() ? -1 : 0);
}

void IMonomial::del()
{
    delete this;
}

uint IMonomial::getLV() const
{
    for(uint i=0,i_end=getIndet();i<i_end;i++) {
        if(at(i)>0)
            return i;
    }
    ASSERT_NOT_REACHED;
    return 0;
}

int IMonomial::getSingleVarIndex() const
{
    int result = 0;
    for(uint i=0;i<getIndet();i++) {
        if(at(i)!=0) {
            if(result!=-1)
                return -1;
            else
                result = i;
        }
    }
    return result;
}

} // namespace polynomial
