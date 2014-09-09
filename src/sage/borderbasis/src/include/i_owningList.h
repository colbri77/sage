#ifndef __I_OWNINGLIST_H__
#define __I_OWNINGLIST_H__

#include "definitions.h"

namespace base {

template <class T>
class IOwningList
{
public:
    IOwningList(){}
    virtual ~IOwningList(){}

    virtual void pop_back() = 0;
    virtual void pop() = 0;
    virtual T pop_lift() = 0;
    virtual void clear() = 0;
    virtual void clear_keep() = 0;
    virtual void remove(uint pos) = 0;
    virtual T lift(uint pos) = 0;
    virtual void push_front(T val) = 0;
    virtual void push_back(T val) = 0;
    virtual const T& at(uint pos) const = 0;
    virtual T& at(uint pos) = 0;
    virtual uint size() const = 0;
};

} // namespace base

#endif // __I_OWNINGLIST_H__
