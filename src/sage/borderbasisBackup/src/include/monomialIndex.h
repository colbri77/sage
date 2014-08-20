#ifndef __MONOMIALINDEX_H__
#define __MONOMIALINDEX_H__

#include <map>
#include "i_owningList.h"
#include "i_polynomial.h"

namespace borderbasis {

template<typename T>
class MonomialIndex
{
public:
    MonomialIndex(const IOwningList<KEEP_REF IPolynomial<T>*>* generators);
    virtual ~MonomialIndex();

    uint getColumns() const;
    uint toIndex(const IMonomial<T>* monomial) const;
    TAKE_OWN IMonomial<T>* toMonomial(uint index) const;

private:
    union cTableData {
        IMonomial<T>* monomial;
        uint integer;
    };

    static uint processors;

    cTableData* cTableFast;
    map<uint,cTableData>* cTableMemory;
    cTableData* cTableRev;
    uint tableColumns;
};

} // namespace borderbasis

#endif // __MONOMIALINDEX_H__
