#include "include/monomialIndex.h"

#include <omp.h>

#define MEMORY_BORDER 0xa0000000    /*Nr of integers, amount is 10GB*/

namespace borderbasis {

template<typename T>
uint MonomialIndex<T>::processors = omp_get_num_procs();

template<typename T>
MonomialIndex<T>::MonomialIndex(const IOwningList<KEEP_REF IPolynomial<T>*>* generators)
: cTableFast(NULL),
cTableMemory(NULL),
cTableRev(NULL),
tableColumns(0)
{
    uint64_t maxPos = 0;
    uint i = 0;
    uint limit = generators->size();
    uint bs = max((uint)1,limit/MonomialIndex<T>::processors);

    // 1. Find the biggest monomial, it defines how much memory we have to use during the index generation
    #pragma omp parallel private(i) shared(generators,limit,maxPos,bs)
    {
        #pragma omp for schedule(dynamic,bs) nowait
        for(i=0;i<limit;i++) {
            if(generators->at(i)->size()==0) continue;
            uint64_t uTmp = generators->at(i)->at(0)->getMonomial()->getPos();
            #pragma omp critical
            {
                if(uTmp>maxPos)
                    maxPos = uTmp;
            }
        }
    }

    // 2. Decide if we can create a fast mapping or if that would need to much memory
    bool useMemoryExpensive = (maxPos<(uint64_t)MEMORY_BORDER);

    if(useMemoryExpensive) {
        cTableFast = new cTableData[maxPos+1]();
    } else {
        cTableMemory = new map<uint,cTableData>();
    }

    limit = generators->size();
    bs = max((uint)1,limit/MonomialIndex<T>::processors);

    // 3. Fill the table - each unique monomial is placed at its exact position
    #pragma omp parallel private(i) shared(limit,bs,generators)
    {
        #pragma omp for schedule(dynamic,bs) nowait
        for(i=0;i<limit;i++) {
            for(uint k=0,ts=generators->at(i)->size();k<ts;k++) {
                IMonomial* t = generators->at(i)->at(k)->getMonomial();
                uint64_t iTemp = t->getPos();
                #pragma omp critical
                {
                    if(useMemoryExpensive) {
                        if(cTableFast[iTemp].integer==0) {
                            cTableFast[iTemp].monomial = t;
                            tableColumns++;
                        }
                    } else {
                        if(cTableMemory->find(iTemp)==cTableMemory->end()) {
                            cTableMemory->at(iTemp).monomial = t;
                            tableColumns++;
                        }
                    }
                }
            }
        }
    }

    // 4. Since we've countet the unique monomials, we know how big the reverse table must be
    cTableRev =  new cTableData[tableColumns];

    // 5. Replace the monomials in the index list by their positions, store them in the reverse list instead
    uint runner=0;
    if(useMemoryExpensive) {
        for(int k=maxPos;k>=0;k--) {
            if(cTableFast[k].monomial!=0) {
                cTableRev[runner].monomial = cTableFast[k].monomial;
                cTableFast[k].integer = runner;
                runner++;
            }
        }
    } else {
        for(typename map<uint,cTableData>::reverse_iterator it = cTableMemory->rbegin();it!=cTableMemory->rend();it++) {
            cTableRev[runner].monomial = it->second.monomial;
            it->second.integer = runner;
            runner++;
        }
    }
}

template<typename T>
MonomialIndex<T>::~MonomialIndex()
{
    DEL_SAFE(cTableFast);
    DEL_SAFE(cTableMemory);
    delete cTableRev;
}

template<typename T>
uint MonomialIndex<T>::toIndex(const IMonomial* monomial) const
{
    if(cTableFast!=NULL)
        return cTableFast[monomial->getPos()].integer;
    else
        return cTableMemory->at(monomial->getPos()).integer;
}

template<typename T>
TAKE_OWN IMonomial* MonomialIndex<T>::toMonomial(uint index) const
{
    return cTableRev[index].monomial->copy();
}

template<typename T>
uint MonomialIndex<T>::getColumns() const
{
    return tableColumns;
}

template class MonomialIndex<uint64_t>;
template class MonomialIndex<int64_t>;

} // namespace borderbasis
