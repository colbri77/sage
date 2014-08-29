#ifndef __FASTFLEXIBLEARRAY_H__
#define __FASTFLEXIBLEARRAY_H__

#include "definitions.h"

#define FFA_BLOCKEXP 10    // measured in 2^x
#define FFA_BLOCKSIZE (1<<FFA_BLOCKEXP)
#define FFA_BLOCKOUT (FFA_BLOCKSIZE-1)

namespace base {

class FFArrayElement
{
public:
    FFArrayElement(){}
    virtual ~FFArrayElement(){}
};

class FastFlexibleArray
{
public:
    FastFlexibleArray();
    virtual ~FastFlexibleArray();

    void add(uint64_t index,FFArrayElement* elem);
    FFArrayElement* get(uint64_t index) const;

private:
    void** nullPtrList;
    void** list;

    void clear(void** list,uint depth);
    void** newList();
};

} // namespace base

#endif // __FASTFLEXIBLEARRAY_H__
