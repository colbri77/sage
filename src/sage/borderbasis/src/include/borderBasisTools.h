#ifndef __BORDERBASISTOOLS_H__
#define __BORDERBASISTOOLS_H__

#include "i_polynomial.h"
#include "i_matrixFactory.h"
#include "polynomialFactory.h"
#include "monomialFactory.h"
#include "i_owningList.h"
#include "compUniverse.h"
#include "statistics.h"
#include "field.h"
#include "bmap128.h"

namespace borderbasis {

enum OptLevel {
	NONE,
	ENHANCED,
	OPTIMISTIC,
	EXPERIMENTAL,
	MUTANT,
	IMPROVED_MUTANT,
	IMPROVED_MUTANT_LINEAR,
	IMPROVED_MUTANT_OPTIMISTIC,
};

class BBConfig {
    public:
        // Cython is extremely picky, struct is not working, only class with explicit setters
	BBConfig(uint indeterminates,OptLevel optimization,void* field,void* matrixFactory,void* polFactory, void* monFactory, bool use_pol_exclusion,bool use_variable_exclusion,bool* variable_exclusions)
	: indeterminates(indeterminates), optimization(optimization), field(field), matrixFactory(matrixFactory), polFactory(polFactory),
	  monFactory(monFactory),use_pol_exclusion(use_pol_exclusion), use_variable_exclusion(use_variable_exclusion),
	  variable_exclusions(variable_exclusions){}
    uint indeterminates;
	OptLevel optimization;
	void* field;
	void* matrixFactory;
    void* polFactory;
	void* monFactory;
	bool use_pol_exclusion;
	bool use_variable_exclusion;
	bool* variable_exclusions;
};

template<typename T>
class BorderBasisTools
{
public:
    BorderBasisTools(BBConfig* config);

    virtual ~BorderBasisTools();

    void calculateBasis(const IOwningList<IPolynomial<T>*>* in,
                            IOwningList<IPolynomial<T>*>* out,
                            IPolynomial<T>* orderIdeal);

    void calculateBasis(const void* in,void* out,void* orderIdeal) {
        // cython has problems with double templates like <IPolynomial<T>*>, therefore this wrapper
        calculateBasis((const IOwningList<IPolynomial<T>*>*)in,
                       (IOwningList<IPolynomial<T>*>*)out,
                       (IPolynomial<T>*)orderIdeal);
    }

    void getStatistics(Statistics* out) const;

private:
    IField<T>* field;
    IMatrixFactory<T>* matrixFactory;
    PolynomialFactory<T>* polFactory;
    IMonomialFactory* monFactory;
    uint indet;
    Statistics* statistics;
    OptLevel optimization;
    ICompUniverse<T>* universe;
    bool getPosSupport;
    BBConfig* config;
    bool* excludedList;
    uint excludedListLen;

    class MutantState {
    public:
        MutantState(uint indet)
            : P_mutant(new BMap128()),
            G(new OwningVector<IPolynomial<T>*>()),
            VHash(new BMap128()),
            d_min(0),
            d_max(0),
            X_(new bool[indet]()) {}
        ~MutantState(){
            delete P_mutant;
            delete G;
            delete VHash;
            delete X_;
        }
        BMap128* P_mutant;
        OwningVector<IPolynomial<T>*>* G; // The position where the new Polynomials have been added
        BMap128* VHash;
        uint d_min;
        uint d_max;
        bool* X_;
    };

    static uint processors;

    void toSimpleBasis(IOwningList<IPolynomial<T>*>* in,bool full);
    void extend(IOwningList<IPolynomial<T>*>* in,bool isBasis);
    void extendMutant(IOwningList<IPolynomial<T>*>* in,bool isBasis,MutantState* mstate);
    void getOrderIdeal(IOwningList<IPolynomial<T>*>* in,IPolynomial<T>* out);
    bool checkOrderIdeal(const IPolynomial<T>* orderIdeal,MutantState* mstate);
    void addAndReduce(IOwningList<IPolynomial<T>*>* in,int pos);
    void reduceFinal(IOwningList<IPolynomial<T>*>* in);
};

}

#endif // __BORDERBASISTOOLS_H__
