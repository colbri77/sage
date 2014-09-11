#include "include/borderBasisTools.h"

#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <stack>
#include <stdlib.h>
#include <omp.h>
#include <cstring>
#include "include/bmap128.h"
#include "include/owningVector.h"
#include "include/degLexMonomial.h"
#include "include/field.h"
#include "include/monomialIndex.h"

namespace borderbasis {

template<typename T>
uint BorderBasisTools<T>::processors = omp_get_num_procs();

template<typename T>
BorderBasisTools<T>::BorderBasisTools(IField<T>* field,
                     PolynomialFactory<T>* polFactory,
                     IMonomialFactory* monFactory,
                     uint indeterminates,
                     OptLevel optimization)
: field(field),
matrixFactory(NULL),
polFactory(polFactory),
monFactory(monFactory),
indet(indeterminates),
statistics(new Statistics()),
optimization(optimization),
universe(NULL),
getPosSupport(monFactory->supportsGetPos())
{
    if(getPosSupport) {
        switch(optimization) {
        case NONE: universe = new LinearCompUniverse<T>(indet); break;
        case ENHANCED: universe = new SpecificCompUniverse<T>(indet); break;
        case OPTIMISTIC: universe = new SpecificCompUniverseNoBorderLog<T>(indet); break;
        case EXPERIMENTAL: universe = new SpecificCompUniverseNoBorderLog<T>(indet); break;
        case MUTANT: universe = new SpecificCompUniverse<T>(indet); break;
        case IMPROVED_MUTANT: universe = new SpecificCompUniverse<T>(indet); break;
        case IMPROVED_MUTANT_LINEAR: universe = new LinearCompUniverse<T>(indet); break;
        case IMPROVED_MUTANT_OPTIMISTIC: universe = new SpecificCompUniverseNoBorderLog<T>(indet); break;
        case IMPROVED_MUTANT_EXPERIMENTAL: universe = new SpecificCompUniverseNoBorderLog<T>(indet); break;
        default: ASSERT_NOT_REACHED;
        }
    } else {
        if(optimization==NONE)
            universe = new LinearCompUniverse<T>(indet);
        else
            universe = new SpecificCompUniverseNoOrderPos<T>(indet);
    }
}

template<typename T>
BorderBasisTools<T>::BorderBasisTools(int dummyNeccessaryForCython,
                     IMatrixFactory<T>* matrixFactory,
                     PolynomialFactory<T>* polFactory,
                     IMonomialFactory* monFactory,
                     uint indeterminates,
                     OptLevel optimization)
: field(NULL),
matrixFactory(matrixFactory),
polFactory(polFactory),
monFactory(monFactory),
indet(indeterminates),
statistics(new Statistics()),
optimization(optimization),
universe(NULL),
getPosSupport(monFactory->supportsGetPos())
{
    if(getPosSupport) {
        switch(optimization) {
        case NONE: universe = new LinearCompUniverse<T>(indet); break;
        case ENHANCED: universe = new SpecificCompUniverse<T>(indet); break;
        case OPTIMISTIC: universe = new SpecificCompUniverseNoBorderLog<T>(indet); break;
        case EXPERIMENTAL: universe = new SpecificCompUniverseNoBorderLog<T>(indet); break;
        case MUTANT: universe = new SpecificCompUniverse<T>(indet); break; 
        case IMPROVED_MUTANT: universe = new SpecificCompUniverse<T>(indet); break;
        case IMPROVED_MUTANT_LINEAR: universe = new LinearCompUniverse<T>(indet); break;
        case IMPROVED_MUTANT_OPTIMISTIC: universe = new SpecificCompUniverseNoBorderLog<T>(indet); break;
        case IMPROVED_MUTANT_EXPERIMENTAL: universe = new SpecificCompUniverseNoBorderLog<T>(indet); break;
        default: ASSERT_NOT_REACHED;
        }
    } else {
        if(optimization==NONE)
            universe = new LinearCompUniverse<T>(indet);
        else
            universe = new SpecificCompUniverseNoOrderPos<T>(indet);
    }
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
    out->max_comparisons_in_reduction = statistics->max_comparisons_in_reduction;
    out->maxMatrix.rows = statistics->maxMatrix.rows;
    out->maxMatrix.columns = statistics->maxMatrix.columns;
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

    // 2. Initialize state for Mutant algorithms
    MutantState* mstate = new MutantState(indet);

    for(bool firstRun=true; true; firstRun=false) {
        // 3. Extend the polynomial set according to the computational universe
        if(optimization==MUTANT || optimization==IMPROVED_MUTANT ||
           optimization==IMPROVED_MUTANT_LINEAR || optimization==IMPROVED_MUTANT_OPTIMISTIC) {
            extendMutant(&tmpVec,!firstRun,mstate);
        } else {
            extend(&tmpVec,!firstRun);
        }

        // 4. The result is already in row echelon form
        // 5. Read the candidate order ideal from the set
        getOrderIdeal(&tmpVec, orderIdeal);

        // 6. Check and possibly extend the computational universe and (maybe) repeat
        if(checkOrderIdeal(orderIdeal,mstate))
            break;
    }

    delete mstate;

    reduceFinal(&tmpVec);

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
bool BorderBasisTools<T>::checkOrderIdeal(const IPolynomial<T>* orderIdeal,MutantState* mstate)
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
    else if(optimization==MUTANT) {
        for(uint i=0;i<indet && result;i++) {
            IPolynomial<T>* pNew = orderIdeal->copy();
            pNew->incrementAtIndet(i);
            if(!universe->contains(pNew)) {
                universe->addBorder();
                mstate->d_min = universe->getMaxDegree();
                mstate->d_max = mstate->d_min;
                result = false;
            }
            delete pNew;
        }
    }
    else if(optimization==IMPROVED_MUTANT || optimization==IMPROVED_MUTANT_LINEAR) {
        bool x_empty = true;
        for(uint i=0;i<indet;i++) {
            if(mstate->X_[i]) {
                x_empty = false;
                break;
            }
        }
        for(uint i=0;i<indet && result;i++) {
            IPolynomial<T>* pNew = orderIdeal->copy();
            pNew->incrementAtIndet(i);
            if(!universe->contains(pNew)) {
                if(x_empty) {
                    universe->addBorder();
                    mstate->d_min = universe->getMaxDegree();
                    mstate->d_max = mstate->d_min;
                }
                result = false;
            }
            delete pNew;
        }
    }
    else if(optimization==IMPROVED_MUTANT_OPTIMISTIC || optimization==IMPROVED_MUTANT_EXPERIMENTAL) {
        bool x_empty = true;
        for(uint i=0;i<indet;i++) {
            if(mstate->X_[i]) {
                x_empty = false;
                break;
            }
        }
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
            if(x_empty) {
                universe->add(ovTemp);
                mstate->d_min = universe->getMaxDegree();
                mstate->d_max = mstate->d_min;
            }
            result = false;
        }
        delete ovTemp;
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
void BorderBasisTools<T>::addAndReduce(IOwningList<IPolynomial<T>*>* in,int pos)
{
    uint64_t cmpCounter = 0;
    uint i = 0;
    IPolynomial<T>* f = NULL;

    //reduction1:
        uint H = pos;
        uint q = 0;

    reduction2:
        if(H == in->size()) {
            if(cmpCounter>statistics->max_comparisons_in_reduction)
                statistics->max_comparisons_in_reduction = cmpCounter;
            return;
        }

    //reduction3:
        f = in->at(H);  // not really removed here, just do that if necessary later
        H++;
        i = 1;
        uint64_t coef = f->at(0)->getCoef();
        // reduce the polynomials leading terms coefficient to one
        if(coef!=1) {
            for(uint i=0;i<f->size();i++) {
                Term<T>* term = f->at(i);
                uint64_t newCoef = field->divide(term->getCoef(),coef);
                if(newCoef==0) {
                    f->remove(i);
                    i--;
                } else {
                    term->setCoef(newCoef);
                }
            }
        }

    reduction4:
        if(f->isZero() || i>pos+q)
            goto reduction7;

    //reduction5:
        if(f->at(0)->getMonomial()->compare(in->at(i-1)->at(0)->getMonomial())==0) {
            cmpCounter++;
            f->subtract(in->at(i-1),field);
            i = 1;
            goto reduction4;
        }

    //reduction6:
        i++;
        goto reduction4;

    reduction7:
        if(f->isZero()) {
            H--;
            in->remove(H);
        } else {
            q++;
        }
        goto reduction2;
}

template<typename T>
void BorderBasisTools<T>::reduceFinal(IOwningList<IPolynomial<T>*>* in)
{
    // if we don't have a field we use matrices, and they already transform into row echelon form.
    if(!field)
        return;

    bool* processed = new bool[in->size()]();

    for(uint procCtr=0;procCtr<in->size();procCtr++) {
        uint posMin = 0;
        for(uint i=0;i<in->size();i++) {
            if(!processed[i] &&
               in->at(i)->at(0)->getMonomial()->compare(in->at(posMin)->at(0)->getMonomial())<0)
                posMin = i;
        }
        processed[posMin] = true;
        IMonomial* mini = in->at(posMin)->at(0)->getMonomial();

        for(uint i=0;i<in->size();i++) {
            if(!processed[i]) {
                for(uint k=0;k<in->at(i)->size();k++) {
                    if(in->at(i)->at(k)->getMonomial()->compare(mini)==0) {
                        in->at(i)->subtract(in->at(posMin),field);
                    }
                }
            }
        }
    }

    delete processed;
}

template<typename T>
void BorderBasisTools<T>::extend(IOwningList<IPolynomial<T>*>* in,bool isBasis)
{
    // 1. If the list does not already describe a basis, convert it to one
    if(!isBasis) {
        if(field) addAndReduce(in,0);
        else toSimpleBasis(in,true);
    }

    // 2. We create a map to store hashes of all processed polynomials so that we can ignore
    //    them when they come back up later. If the value in the list is false, then they are
    //    in the "to do"-list, but have not been processed yet - they don't need to be added.
    //    I the value is true, they have been handled completely and can be ignored.
    BMap128 polMap = BMap128();
    uint64_t hash[2] = {0,0};

    while(true) {
        uint end = in->size();
        uint lastOriginalPolynomial = end-1;

        // 3. For each polynomial in the original list:
        for(uint i=0;i<end;i++) {
            IPolynomial<T>* currentPol = in->at(i);
            if(currentPol->isZero())
                continue;

            // 3.1 on "experimental", we have to make sure this one is still in the computational universe
            if(optimization==EXPERIMENTAL) {
                if(!universe->contains(currentPol->at(0)->getMonomial()))
                    continue;
            }

            // 3.2 Check if this polynomial has already been handled completely
            currentPol->hash(hash);
            if(polMap.contains(hash) && polMap.get(hash)) {
                continue; // 3.2.1 If this has already been handled, continue with the next polynomial
            } else {
                polMap.set(hash,true);
            }

            // 3.3 Generate "offsprings"
            for(uint k=0;k<indet;k++) {
                IPolynomial<T>* p = currentPol->copy();
                p->incrementAtIndet(k);
                p->hash(hash);
                // 3.5 If offspring has never been seen before, add it to the list
                if(!polMap.contains(hash)) {
                    polMap.set(hash,false);
                    in->push_back(p);
                } else {
                    delete p;
                }
            }

            // 3.6 On OptLevel EXPERIMENTAL, we try to extend the universe as much as possible immediately
            if(optimization==EXPERIMENTAL) {
                end = in->size();
            }
        }
        
        // 4. We now have the extended list, calculate a new basis of it.
        if(field) addAndReduce(in,lastOriginalPolynomial+1);
        else toSimpleBasis(in,false);
        
        OwningVector<IPolynomial<T>*> ovTemp = OwningVector<IPolynomial<T>*>();

        // 5. remove elements that have a leading term outside the universe
        int limit = (field ? lastOriginalPolynomial : -1);
        for(int i=(int)in->size()-1;i>limit;i--) {
            if(!universe->contains(in->at(i)->at(0)->getMonomial()))
               //in->remove(i);
                ovTemp.push_back(in->lift(i));
        }

        // 6. If the OptLevel is high enough, we have to extend the comp. universe
        if(optimization!=NONE)
            universe->add(in,field ? lastOriginalPolynomial : 0);

        // 7. check if polynomials in ovTemp are valid now
        uint newPolCtr = 0;
        do {
            newPolCtr = 0;
            for(uint i=0;i<ovTemp.size();i++) {
                if(universe->contains(ovTemp.at(i)->at(0)->getMonomial())) {
                    in->push_back(ovTemp.lift(i));
                    newPolCtr++;
                    i--;
                }
            }
            universe->add(in,in->size()-newPolCtr-1);
        } while(newPolCtr>0);


        // 8. If the size didn't change, its still the same basis and we're done.
        if(in->size()==lastOriginalPolynomial+1)
            break;
    }
}

template<typename T>
void BorderBasisTools<T>::extendMutant(IOwningList<IPolynomial<T>*>* in,bool isBasis,MutantState* mstate)
{
    uint H = 0;
    uint64_t hash[2] = {0,0};
    uint d_elim = 0;
    stack<IPolynomial<T>*> M = stack<IPolynomial<T>*>();
    IPolynomial<T>* currentPol = NULL;
    OwningVector<IPolynomial<T>*> W_ = OwningVector<IPolynomial<T>*>();
    OwningVector<IPolynomial<T>*> W = OwningVector<IPolynomial<T>*>();

    // never called before, execute step 1
    if(!isBasis) {
        if(field) addAndReduce(in,0);
        else toSimpleBasis(in,true);

        mstate->G->clear();
        mstate->P_mutant->clear();
        mstate->d_max = 0;
        mstate->d_min = 0x7fffffff;
        for(uint i=0;i<in->size();i++) {
            in->at(i)->hash(hash);
            mstate->VHash->set(hash,true);
            mstate->G->push_back(in->at(i)->copy());
            uint degree = in->at(i)->at(0)->getMonomial()->getDegree();
            if(mstate->d_max < degree)
                mstate->d_max = degree;
            if(mstate->d_min > degree)
                mstate->d_min = degree;
        }
        for(uint i=0;i<indet;i++)
            mstate->X_[i] = false;
    }

    mutantS2:
        if(optimization==MUTANT) {
            H = mstate->G->size();
            for(uint i=0,i_end=mstate->G->size();i<i_end;i++) {
                currentPol = mstate->G->at(i);
                currentPol->hash(hash);
                if(currentPol->at(0)->getMonomial()->getDegree()==mstate->d_min &&
                   !mstate->P_mutant->contains(hash)) {

                    mstate->P_mutant->set(hash,true);
                    d_elim = mstate->d_min+1;

                    for(uint k=0;k<indet;k++) {
                        IPolynomial<T>* p = currentPol->copy();
                        p->incrementAtIndet(k);
                        mstate->G->push_back(p);
                    }
                }
            }
        }
        else if(optimization==IMPROVED_MUTANT || optimization==IMPROVED_MUTANT_LINEAR || 
                optimization==IMPROVED_MUTANT_OPTIMISTIC || optimization==IMPROVED_MUTANT_EXPERIMENTAL) {
            //mutantS2a:
            int xl = -1;
            for(int i=indet-1;i>=0;i--) {
                if(mstate->X_[i]) {
                    mstate->X_[i] = false;
                    xl = i;
                    break;
                }
            }
            if(xl==-1) {
                for(uint i=0,i_end=mstate->G->size();i<i_end;i++) {
                    currentPol = mstate->G->at(i);
                    if(currentPol->at(0)->getMonomial()->getDegree()==mstate->d_min) {
                        mstate->X_[currentPol->at(0)->getMonomial()->getLV()] = true;
                    }
                }
                d_elim = mstate->d_min + 1;
                for(int i=indet-1;i>=0;i--) {
                    if(mstate->X_[i]) {
                        mstate->X_[i] = false;
                        xl = i;
                        break;
                    }
                }
            }

            //mutantS2b:
            H = mstate->G->size();
            for(uint i=0,i_end=mstate->G->size();i<i_end;i++) {
                currentPol = mstate->G->at(i);
                currentPol->hash(hash);
                IMonomial* monomial = currentPol->at(0)->getMonomial();
                if(monomial->getDegree()==mstate->d_min &&
                   monomial->getLV()==xl &&
                   (optimization != IMPROVED_MUTANT_EXPERIMENTAL || universe->contains(monomial)) &&
                   !mstate->P_mutant->contains(hash)) {

                    mstate->P_mutant->set(hash,true);
                    d_elim = mstate->d_min+1;

                    for(uint k=0;k<indet;k++) {
                        IPolynomial<T>* p = currentPol->copy();
                        p->incrementAtIndet(k);
                        mstate->G->push_back(p);
                    }
                }
            }
        }

    mutantS3:
        if(field) addAndReduce(mstate->G,H);
        else toSimpleBasis(mstate->G,false);
    
    //mutantS4:
        for(uint i=0,i_end=mstate->G->size();i<i_end;i++) {
            if(mstate->G->at(i)->at(0)->getMonomial()->getDegree()<d_elim)
                M.push(mstate->G->at(i));
        }

    //mutantS5:
    //mutantS5_:
        int necessary = (int)M.size();
        if(optimization==IMPROVED_MUTANT || optimization==IMPROVED_MUTANT_LINEAR || 
           optimization==IMPROVED_MUTANT_OPTIMISTIC || optimization==IMPROVED_MUTANT_EXPERIMENTAL) {
            // caluclate the "necessary" amount of polynomials
            uint k = 0x7fffffff;
            uint Q = 0;
            stack<IPolynomial<T>*> sTmp = stack<IPolynomial<T>*>();
            
            while(M.size()>0) {
                currentPol = M.top();
                M.pop();
                sTmp.push(currentPol);
            }
            k = d_elim-1;
            while(!sTmp.empty()) {
                currentPol = sTmp.top();
                M.push(currentPol);
                sTmp.pop();
            }

            for(uint i=0,i_end=mstate->G->size();i<i_end;i++) {
                uint degree = mstate->G->at(i)->at(0)->getMonomial()->getDegree();
                if(degree <= k+1)
                    Q++;
            }

            int64_t Sk = 0;
            int64_t lastValue = 1;
            for(uint l=1;l<=k+1;l++) {
                lastValue *= (indet-l+1);
                lastValue /= l;
                Sk += lastValue;
            }

            int nc = (Sk-Q)/indet+1;
            if(nc<=0) necessary = 1;
            else if(nc<necessary) necessary = nc;
        }

        H = mstate->G->size();
        bool isEmpty = true;
        uint d_elim_new = 0x7fffffff;
        for(;!M.empty();) {
            currentPol = M.top();
            M.pop();
            currentPol->hash(hash);
            if(!mstate->P_mutant->contains(hash)) {
                necessary--;
                isEmpty = false;
                if(currentPol->at(0)->getMonomial()->getDegree()<d_elim_new)
                    d_elim_new = currentPol->at(0)->getMonomial()->getDegree();
                if(necessary>=0) {
                    mstate->P_mutant->set(hash,true);
                    for(uint k=0;k<indet;k++) {
                        IPolynomial<T>* p = currentPol->copy();
                        p->incrementAtIndet(k);
                        mstate->G->push_back(p);
                    }
                }
            }
        }
        while(!M.empty())
            M.pop();
        if(!isEmpty) {
            d_elim = d_elim_new + 1;
            goto mutantS3;
        }

    //mutantS6:
        if(d_elim <= mstate->d_min) {
            d_elim++;
            goto mutantS3;
        }

    mutantS7:
        W_.clear_keep();
        for(uint i=0,i_end=mstate->G->size();i<i_end;i++) {
            mstate->G->at(i)->hash(hash);
            if(!mstate->VHash->contains(hash)) {
                W_.push_back(mstate->G->at(i));
            }
        }

    //mutantS8:
        W.clear_keep();
        for(uint i=0;i<W_.size();i++) {
            if(universe->contains(W_[i]->at(0)->getMonomial())) {
                W.push_back(W_[i]);
                W_.lift(i);
                i--;
            }
        }

    //mutantS9:
        if(W.size()>0) {
            universe->add(&W);
            //mutantS10:
            while(W.size()>0) {
                IPolynomial<T>* p = W.pop_lift();
                in->push_back(p->copy());
                p->hash(hash);
                mstate->VHash->set(hash,true);
            }
            goto mutantS7; 
        }
    //mutantS11:
        if(mstate->d_min<mstate->d_max) {
            mstate->d_min++;
            goto mutantS2;
        }
    //mutantS12:
        W.clear_keep();
        W_.clear_keep();

        if(!field) toSimpleBasis(in,false);
}

template<typename T>
void BorderBasisTools<T>::getOrderIdeal(IOwningList<IPolynomial<T>*>* in,IPolynomial<T>* out)
{
    out->clear();

    IMonomial* t = monFactory->create();
    IMonomial* tTemp = NULL;

    IPolynomial<T>* p = polFactory->create(indet);
    // collect leading monomials
    for(uint i=0,end_i=in->size();i<end_i;i++) {
        p->push(new Term<T>(1,in->at(i)->at(0)->getMonomial()->copy()));
    }

    for(int i=p->size()-1;i>=0;i--) {
        IMonomial* tLead = p->at(i)->getMonomial();
        while(tLead->compare(t)>0) {
            bool inUniverse = universe->contains(t);
            if(inUniverse)
                out->push(new Term<T>(1,t));
            tTemp = t;
            t = t->next();
            if(tTemp == t)
                break;
            if(!inUniverse)
                tTemp->del();
        }
        tTemp = t;
        t = t->next();
        if(tTemp == t)
            break;
        tTemp->del();
    }
    while(!universe->beyondLastElement(t) && tTemp != t) {
        bool inUniverse = universe->contains(t);
        if(inUniverse)
            out->push(new Term<T>(1,t));
        tTemp = t;
        t = t->next();
        if(tTemp == t)
            break;
        if(!inUniverse)
            tTemp->del();
    }
    delete p;
    t->del();
}

template class BorderBasisTools<uint64_t>;
template class BorderBasisTools<int64_t>;

} // namespace borderbasis
