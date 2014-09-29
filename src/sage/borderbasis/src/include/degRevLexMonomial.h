#ifndef __DEGREVLEXMONOMIAL_H__
#define __DEGREVLEXMONOMIAL_H__

#include "degLexMonomial.h"

namespace polynomial {

class DegRevLexMonomial : public DegLexMonomial
{
protected:
    friend class MonomialFactoryDegRevLex;
    friend class FastFlexibleArray;

    virtual DegLexMonomial* create(uint indet, FastFlexibleArray* monomBox) const OVERRIDE;

    DegRevLexMonomial(uint64_t pos, uint indet, FastFlexibleArray* monomBox);
    DegRevLexMonomial(uint indet, FastFlexibleArray* monomBox);
    virtual ~DegRevLexMonomial();

    virtual void recalcPos() OVERRIDE;
    virtual void initFromPos(uint64_t pos) OVERRIDE;
};

class DegRevLexMonomialGF2 : public DegLexMonomialGF2
{
protected:
    friend class MonomialFactoryDegRevLexGF2;
    friend class FastFlexibleArray;

    virtual DegLexMonomial* create(uint indet, FastFlexibleArray* monomBox) const OVERRIDE;

    DegRevLexMonomialGF2(uint64_t pos, uint indet, FastFlexibleArray* monomBox,bool* excludedIndets);
    DegRevLexMonomialGF2(uint indet, FastFlexibleArray* monomBox,bool* excludedIndets);
    virtual ~DegRevLexMonomialGF2();

    virtual void recalcPos() OVERRIDE;
    virtual void initFromPos(uint64_t pos) OVERRIDE;
};

} // namespace polynomial

#endif // __DEGREVLEXMONOMIAL_H__

