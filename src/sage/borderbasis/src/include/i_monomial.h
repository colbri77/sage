#ifndef __I_MONOMIAL_H__
#define __I_MONOMIAL_H__

#include "definitions.h"
#include "fastFlexibleArray.h"

namespace polynomial {

class IMonomial : public FFArrayElement
{
public:
    enum termOrder {
        DEGLEX=0,
        DEGREVLEX=1
    };

    termOrder termOrdering;

    IMonomial(termOrder termOrdering);

    virtual const uint& at(uint const& index) const = 0;
    // set() creates a (probably) new monomial, the old one should be considered deleted.
    virtual TAKE_OWN IMonomial* set(uint index, uint value) = 0;

    virtual uint getIndet() const = 0;
    virtual uint getDegree() const = 0;
    // extend() creates a (probably) new monomial, the old one should be considered deleted.
    virtual TAKE_OWN IMonomial* extend(uint index, int value) = 0;
    virtual TAKE_OWN IMonomial* copy() const = 0;
    virtual TAKE_OWN IMonomial* next() const = 0;
    virtual bool divides(const IMonomial* numerator) const = 0;
    virtual bool isBorderOf(const IMonomial* monomial) const = 0;

    virtual bool supportsGetPos() const = 0;
    virtual uint64_t getPos() const = 0;

    virtual void del();
    virtual int compare(const IMonomial* other) const;
    virtual uint getLV() const;
    virtual int getSingleVarIndex() const;

protected:
    virtual ~IMonomial();
};

} // namespace polynomial

#endif // __I_MONOMIAL_H__
