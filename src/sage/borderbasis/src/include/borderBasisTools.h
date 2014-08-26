#ifndef __BORDERBASISTOOLS_H__
#define __BORDERBASISTOOLS_H__

#include "i_polynomial.h"
#include "i_matrixFactory.h"
#include "polynomialFactory.h"
#include "monomialFactory.h"
#include "i_owningList.h"
#include "compUniverse.h"
#include "statistics.h"

namespace borderbasis {

enum OptLevel {
	NONE,
	ENHANCED,
	OPTIMISTIC,
	EXPERIMENTAL
};

template<typename T>
class BorderBasisTools
{
public:
    BorderBasisTools(IMatrixFactory<T>* matrixFactory,
                     PolynomialFactory<T>* polFactory,
                     MonomialFactory* monFactory,
                     uint indeterminates,
                     OptLevel optimization);
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
    IMatrixFactory<T>* matrixFactory;
    PolynomialFactory<T>* polFactory;
    MonomialFactory* monFactory;
    uint indet;
    Statistics* statistics;
    OptLevel optimization;
    ICompUniverse<T>* universe;
    bool getPosSupport;

    static uint processors;

    void toSimpleBasis(IOwningList<IPolynomial<T>*>* in,bool full);
    void extend(IOwningList<IPolynomial<T>*>* in,bool isBasis);
    void getOrderIdeal(IOwningList<IPolynomial<T>*>* in,IPolynomial<T>* out);
    bool checkOrderIdeal(const IPolynomial<T>* orderIdeal);
};

}

#endif // __BORDERBASISTOOLS_H__
