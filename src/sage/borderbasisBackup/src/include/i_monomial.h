#ifndef __I_MONOMIAL_H__
#define __I_MONOMIAL_H__

#include "definitions.h"

namespace polynomial {

template<typename T>
class IMonomial
{
public:
    enum termOrder {
        DEGLEX
    };

    const termOrder termOrdering;

    IMonomial<T>(termOrder termOrdering);
    virtual ~IMonomial<T>();

    virtual const uint& at(uint const& index) const = 0;
    virtual void set(uint index, uint value) = 0;

    virtual uint getIndet() const = 0;
    virtual uint getDegree() const = 0;
    virtual void extend(uint index, int value) = 0;
    virtual TAKE_OWN IMonomial<T>* copy() const = 0;
    virtual TAKE_OWN IMonomial<T>* next() const = 0;
    virtual bool divides(const IMonomial<T>* numerator) const = 0;
    virtual bool isBorderOf(const IMonomial<T>* monomial) const = 0;

    virtual uint64_t getPos() const = 0;

    int compare(const IMonomial<T>* other) const;
};

} // namespace polynomial

#endif // __I_MONOMIAL_H__
