#ifndef __OWNINGVECTOR_H__
#define __OWNINGVECTOR_H__

#include <vector>
#include "i_owningList.h"

namespace base {

template <class T>
class OwningVector : public vector<T>, public IOwningList<T>
{
    public:
        OwningVector() : vector<T>() {}
        OwningVector(uint n) : vector<T>(n) {}
        OwningVector(const OwningVector& x) : vector<T>(x) {}
        virtual ~OwningVector() {clear();}

        virtual void pop_back()
        {
            delete vector<T>::back();
            vector<T>::pop_back();
        }

        virtual void pop()
        {
            delete vector<T>::back();
            vector<T>::pop_back();
        }

        virtual T pop_lift()
        {
            T result = vector<T>::back();
            vector<T>::pop_back();
            return result;
        }

        virtual void clear()
        {
            uint c = vector<T>::size();
            for(uint i=0;i<c;i++)
                delete vector<T>::at(i);
            vector<T>::clear();
        }

        virtual void push_front(T val)
        {
            vector<T>::insert(vector<T>::begin(),val);
        }

        virtual void push_back(T val)
        {
            vector<T>::push_back(val);
        }

        virtual const T& at(uint pos) const
        {
            return vector<T>::at(pos);
        }

        virtual T& at(uint pos)
        {
            return vector<T>::at(pos);
        }

        virtual uint size() const
        {
            return vector<T>::size();
        }

        virtual void remove(uint pos)
        {
            delete vector<T>::at(pos);
            typename OwningVector<T>::iterator it = vector<T>::begin();
            it += pos;
            vector<T>::erase(it);
        }

        virtual T lift(uint pos)
        {
            T result = vector<T>::at(pos);
            typename OwningVector<T>::iterator it = vector<T>::begin();
            it += pos;
            vector<T>::erase(it);
            return result;
        }
};

}

#endif // __OWNINGVECTOR_H__
