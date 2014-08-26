#ifndef __I_MONOMIAL_H__
#define __I_MONOMIAL_H__

#include "definitions.h"

namespace polynomial {

class IMonomial
{
public:
    enum termOrder {
        DEGLEX
    };

    const termOrder termOrdering;

    IMonomial(termOrder termOrdering);
    virtual ~IMonomial();

    virtual const uint& at(uint const& index) const = 0;
    virtual void set(uint index, uint value) = 0;

    virtual uint getIndet() const = 0;
    virtual uint getDegree() const = 0;
    virtual void extend(uint index, int value) = 0;
    virtual TAKE_OWN IMonomial* copy() const = 0;
    virtual TAKE_OWN IMonomial* next() const = 0;
    virtual bool divides(const IMonomial* numerator) const = 0;
    virtual bool isBorderOf(const IMonomial* monomial) const = 0;

    virtual bool supportsGetPos() const = 0;
    virtual uint64_t getPos() const = 0;

    int compare(const IMonomial* other) const;
};

} // namespace polynomial

#endif // __I_MONOMIAL_H__
