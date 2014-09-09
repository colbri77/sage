#ifndef __BMAP128_H__
#define __BMAP128_H__

#include <map>
#include "definitions.h"

namespace base {

class BMap128
{
public:
    BMap128();
    ~BMap128();

    bool contains(const uint64_t* index) const;
    bool get(const uint64_t* index) const;
    void set(const uint64_t* index,bool value);
    void clear();

private:
    map<uint64_t,map<uint64_t,bool> > _map;
};

} // namespace base

#endif // __MAP128_H__
