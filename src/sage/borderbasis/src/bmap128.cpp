#include "include/bmap128.h"

namespace base {

BMap128::BMap128()
: _map(map<uint64_t,map<uint64_t,bool> >())
{
}

BMap128::~BMap128()
{
}

bool BMap128::contains(const uint64_t* index) const
{
    map<uint64_t,map<uint64_t,bool> >::const_iterator it1 = _map.find(index[0]);
    if(it1==_map.end())
        return false;

    map<uint64_t,bool>::const_iterator it2 = it1->second.find(index[1]);
    if(it2==it1->second.end())
        return false;

    return true;
}

bool BMap128::get(const uint64_t* index) const
{
    const map<uint64_t,bool>& tmpMap =  _map.at(index[0]);
    return tmpMap.at(index[1]);
}

void BMap128::set(const uint64_t* index,bool value)
{
    map<uint64_t,map<uint64_t,bool> >::iterator it1 = _map.find(index[0]);
    if(it1==_map.end())
        _map[index[0]] = map<uint64_t,bool>();

    _map[index[0]][index[1]] = value;
}

void BMap128::clear()
{
    _map.clear();
}

} // namespace base
