#include "include/fastFlexibleArray.h"

namespace base {

FastFlexibleArray::FastFlexibleArray()
: nullPtrList(new void*[FFA_BLOCKSIZE]),
list(NULL)
{
    for(uint i=0;i<FFA_BLOCKSIZE;i++)
        nullPtrList[i] = nullPtrList;
    list = newList();
}

FastFlexibleArray::~FastFlexibleArray()
{
    clear(list,4);
    delete list;
    delete nullPtrList;
}

void FastFlexibleArray::add(uint64_t index,FFArrayElement* elem)
{
    void** l = (void**)(list[index>>(3*FFA_BLOCKEXP)]);
    l = (void**)(l[(index>>(FFA_BLOCKEXP<<1))&FFA_BLOCKOUT]);
    l = (void**)(l[(index>>FFA_BLOCKEXP)&FFA_BLOCKOUT]);

    if(l==nullPtrList) {
        uint pos = index>>(3*FFA_BLOCKEXP);
        if(list[pos]==nullPtrList)
            list[pos] = newList();
        l = (void**)list[pos];
        pos = (index>>(FFA_BLOCKEXP<<1))&FFA_BLOCKOUT;
        if(l[pos]==nullPtrList)
            l[pos] = newList();
        l = (void**)l[pos];
        pos = (index>>FFA_BLOCKEXP)&FFA_BLOCKOUT;
        if(l[pos]==nullPtrList)
            l[pos] = newList();
        l = (void**)l[pos];
        l[index&FFA_BLOCKOUT] = (void*)elem;
    } else {
        l[index&FFA_BLOCKOUT] = (void*)elem;
    }
}

FFArrayElement* FastFlexibleArray::get(uint64_t index) const
{
    void* result = list[index>>(3*FFA_BLOCKEXP)];
    result = ((void**)result)[(index>>(FFA_BLOCKEXP<<1))&FFA_BLOCKOUT];
    result = ((void**)result)[(index>>FFA_BLOCKEXP)&FFA_BLOCKOUT];
    result = ((void**)result)[index&FFA_BLOCKOUT];

    return (FFArrayElement*)(result==nullPtrList ? NULL : result);
}

void FastFlexibleArray::clear(void** list,uint depth)
{
    depth--;
    for(uint i=0;i<FFA_BLOCKSIZE;i++) {
        void* ptr = list[i];
        if(ptr!=nullPtrList) {
            if(depth>0) {
                clear((void**)ptr,depth);
            } else {
                delete ((FFArrayElement*)ptr);
            }
        }
    }
}

void** FastFlexibleArray::newList()
{
    void** result = new void*[FFA_BLOCKSIZE];
    for(uint i=0;i<FFA_BLOCKSIZE;i++)
        result[i] = nullPtrList;
    return result;
}

} // namespace base
