#include "include/borderBasisTools.h"

#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <stack>
#include <stdlib.h>
#include <map>
#include <omp.h>
#include "include/owningVector.h"
#include "include/monomialIndex.h"
#include "include/degLexMonomial.h"

namespace borderbasis {

template<typename T>
uint BorderBasisTools<T>::processors = omp_get_num_procs();

template<typename T>
BorderBasisTools<T>::BorderBasisTools(IMatrixFactory<T>* matrixFactory,
                     PolynomialFactory<T>* polFactory,
                     MonomialFactory* monFactory,
                     uint indeterminates,
                     OptLevel optimization)
: matrixFactory(matrixFactory),
polFactory(polFactory),
monFactory(monFactory),
indet(indeterminates),
statistics(new Statistics()),
optimization(optimization),
universe(NULL)
{
    if(optimization==NONE)
        universe = new LinearCompUniverse<T>(indet);
    else
        universe = new SpecificCompUniverse<T>(indet);
}

template<typename T>
BorderBasisTools<T>::~BorderBasisTools()
{
    delete universe;
    //delete matrixFactory;
    //delete polFactory;
    //delete monFactory;
    delete statistics;
}

template<typename T>
void BorderBasisTools<T>::getStatistics(Statistics* out) const
{
    out->maxMatrix.columns = statistics->maxMatrix.columns;
    out->maxMatrix.rows = statistics->maxMatrix.rows;
}

template<typename T>
void BorderBasisTools<T>::calculateBasis(const IOwningList<IPolynomial<T>*>* in,
                            IOwningList<IPolynomial<T>*>* out,
                            IPolynomial<T>* orderIdeal)
{
    statistics->start();
    out->clear();
    universe->clear();

    // initialize working vector
    OwningVector<IPolynomial<T>*> tmpVec = OwningVector<IPolynomial<T>*>();
    for(uint i=0,ts=in->size();i<ts;i++) {
        if(!in->at(i)->isZero()) {
            tmpVec.push_back(in->at(i)->copy());
        }
    }

    // start the algorithm
    // 1. Initialize the computational universe
    universe->add(&tmpVec);

    for(bool firstRun=true; true; firstRun=false) {
        // 2. Extend the polynomial set according to the computational universe
        extend(&tmpVec,!firstRun);

        // 3. The result is already in row echelon form
        // 4. Read the candidate order ideal from the set
        getOrderIdeal(&tmpVec,orderIdeal);

        // 5. Check and possibly extend the computational universe and (maybe) repeat
        if(checkOrderIdeal(orderIdeal))
            break;
    }

    // remove polynomials that are not border bases
    while(tmpVec.size()>0) {
        IPolynomial<T>* p = tmpVec.pop_lift();
        bool isBorderBasis = false;
        for(uint i=0,ts=orderIdeal->size();i<ts;i++) {
            IMonomial* t = orderIdeal->at(i)->getMonomial();
            if(p->at(0)->getMonomial()->isBorderOf(t)) {
                isBorderBasis = true;
                break;
            }
        }
        if(isBorderBasis)
            out->push_back(p);
        else
            delete p;
    }

    statistics->stop();
}

template<typename T>
bool BorderBasisTools<T>::checkOrderIdeal(const IPolynomial<T>* orderIdeal)
{
    bool result = true;

    if(orderIdeal->size()==(uint)0)
            return true;

    if(optimization==NONE) {
        if(orderIdeal->at(0)->getMonomial()->getDegree()>=universe->getMaxDegree()) {
            universe->addBorder();
            result = false;
        }
    }
    else if(optimization==ENHANCED) {
        for(uint i=0;i<indet && result;i++) {
            IPolynomial<T>* pNew = orderIdeal->copy();
            pNew->incrementAtIndet(i);
            if(!universe->contains(pNew)) {
                universe->addBorder();
                result = false;
            }
            delete pNew;
        }
    }
    else if(optimization==OPTIMISTIC || optimization==EXPERIMENTAL) {
        OwningVector<IPolynomial<T>*>* ovTemp = new OwningVector<IPolynomial<T>*>();
        for(uint i=0;i<indet;i++) {
            IPolynomial<T>* pNew = orderIdeal->copy();
            pNew->incrementAtIndet(i);
            if(!universe->contains(pNew)) {
                ovTemp->push_back(pNew);
            } else {
                delete pNew;
            }
        }
        if(ovTemp->size()>0) {
            universe->add(ovTemp);
            result = false;
        }
        delete ovTemp;
    }

    return result;
}

template<typename T>
void BorderBasisTools<T>::toSimpleBasis(IOwningList<IPolynomial<T>*>* in,bool full)
{
    MonomialIndex<T>* index = new MonomialIndex<T>(in);

    uint columns = index->getColumns();
    uint rows = in->size();

    IMatrix<T>* matrix = matrixFactory->create(rows,columns);

    statistics->logMatrix(rows,columns);

    uint i = 0, r = 0;
    uint bs = max((uint)1,rows/BorderBasisTools<T>::processors);

    // fill the matrix
    #pragma omp parallel private(i,r) shared(rows,in,matrix,bs)
    {
        #pragma omp for schedule(dynamic,bs) nowait
        for(r=0;r<rows;r++) {
            IPolynomial<T>* p = in->at(r);
            uint ts=p->size();
            for(i=0;i<ts;i++) {
                IMonomial* t = p->at(i)->getMonomial();
                uint cTmp = index->toIndex(t);
                uint64_t vTmp = p->at(i)->getCoef();
                matrix->set(r,cTmp,vTmp);
            }
        }
    }

    // transformation in row echelon form
    matrix->toRowEchelon(true);

    // store polynomials in "in" temporarily, the MonomialIndex relies on their existence.
    OwningVector<IPolynomial<T>*> tmpStorage = OwningVector<IPolynomial<T>*>();
    for(i=in->size();i>0;i--)
        tmpStorage.push_back(in->pop_lift());

    // extract the new polynomials from the matrix
    bs = max((uint)1,rows/BorderBasisTools<T>::processors);
    #pragma omp parallel private(r) shared(rows,columns,bs,matrix,in)
    {
        #pragma omp for ordered schedule(dynamic,bs) nowait
        for(r=0;r<rows;r++) {
            uint curPos = 0;
            IPolynomial<T>* p = polFactory->create(indet);
            for(uint c=r;c<columns;c++) {
                if(matrix->get(r,c)!=0) {
                    IMonomial* t = index->toMonomial(c);
                    p->push_back(new Term<T>(matrix->get(r,c),t));
                }
            }
            if(!p->isZero()) {
                #pragma omp ordered
                {
                    in->push_back(p);
                }
            }
            else
                delete p;
        }
    }

    delete matrix;
    delete index;
}

template<typename T>
void BorderBasisTools<T>::extend(IOwningList<IPolynomial<T>*>* in,bool isBasis)
{
    // 1. If the list does not already describe a basis, convert it to one
    if(!isBasis)
        toSimpleBasis(in,true);

    // 2. We create a map to store hashes of all processed polynomials so that we can ignore
    //    them when they come back up later. If the value in the list is false, then they are
    //    in the "to do"-list, but have not been processed yet - they don't need to be added.
    //    I the value is true, they have been handled completely and can be ignored.
    map<uint64_t,bool> polMap = map<uint64_t,bool>();
    uint64_t hash = 0;
    map<uint64_t,bool>::iterator it;

    while(true) {
        uint end = in->size();
        uint lastOriginalPolynomial = end-1;

        // 3. For each polynomial in the original list:
        for(uint i=0;i<end;i++) {
            IPolynomial<T>* currentPol = in->at(i);

            // 3.1 We check if the "offspring" would still be in the computational universe
            // ---this is a shortcut for the NONE-Optimization
            if(optimization==NONE) {
                if(currentPol->at(0)->getMonomial()->getDegree()>=universe->getMaxDegree())
                    continue;   // 3.1.1 If not, continue with the next polynomial
            }

            // 3.2 Check if this polynomial has already been handled completely
            hash = currentPol->hash();
            it=polMap.find(hash);
            if(it!=polMap.end() && it->second) {
                continue; // 3.2.1 If this has already been handled, continue with the next polynomial
            } else {
                polMap[hash] = true;
            }

            // 3.3 Generate "offsprings"
            for(uint k=0;k<indet;k++) {
                IPolynomial<T>* p = currentPol->copy();
                p->incrementAtIndet(k);
                // 3.4 Check if the leading term is in the universe. On OptLevel NONE this is already sure.
                if(optimization==NONE || universe->contains(p->at(0)->getMonomial())) {
                    hash = p->hash();
                    it = polMap.find(hash);
                    // 3.5 If offspring has never been seen before, add it to the list
                    if(it==polMap.end()) {
                        polMap[hash] = false;
                        in->push_back(p);
                    } else {
                        delete p;
                    }
                } else {
                    delete p;
                }
            }

            // 3.6 On OptLevel EXPERIMENTAL, we try to extend the universe as much as possible immediately
            if(optimization==EXPERIMENTAL) {
                end = in->size();
            }
        }

        // 4. If the OptLevel is high enough, we have to extend the comp. universe
        if(optimization!=NONE)
            universe->add(in,lastOriginalPolynomial);

        // 5. We now have the extended list, calculate a new basis of it.
        toSimpleBasis(in,true);
        // 6. If the size didn't change, its still the same basis and we're done.
        if(in->size()==lastOriginalPolynomial+1)
            break;
    }
}

template<typename T>
void BorderBasisTools<T>::getOrderIdeal(IOwningList<IPolynomial<T>*>* in,IPolynomial<T>* out)
{
    out->clear();

    IMonomial* t = monFactory->create(indet);
    IMonomial* tTemp = NULL;

    uint lastElemPos = universe->getMaxPos();

    for(int i=in->size()-1;i>=0;i--) {
        IMonomial* tLead = in->at(i)->at(0)->getMonomial();
        while(tLead->compare(t)>0) {
            bool inUniverse = universe->contains(t);
            if(inUniverse)
                out->push(new Term<T>(1,t));
            tTemp = t;
            t = t->next();
            if(!inUniverse)
                delete tTemp;
        }
        tTemp = t;
        t = t->next();
        delete tTemp;
    }
    while(t->getPos()<=lastElemPos) {
        bool inUniverse = universe->contains(t);
        if(inUniverse)
            out->push(new Term<T>(1,t));
        tTemp = t;
        t = t->next();
        if(!inUniverse)
            delete tTemp;
    }
    delete t;
}

template class BorderBasisTools<uint64_t>;
template class BorderBasisTools<int64_t>;

} // namespace borderbasis
