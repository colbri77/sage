#ifndef __DEGLEXMONOMIAL_H__
#define __DEGLEXMONOMIAL_H__

#include "i_monomial.h"
#include "fastFlexibleArray.h"

namespace polynomial {

class MonomialFactoryDegLex;
class MonomialFactoryNoOrderPos;

class DegLexMonomial : public IMonomial
{
public:
    virtual const uint& at(uint const& index) const OVERRIDE;
    virtual TAKE_OWN IMonomial* set(uint index, uint value) OVERRIDE;

    virtual uint getIndet() const OVERRIDE;
    virtual uint getDegree() const OVERRIDE;
    virtual TAKE_OWN IMonomial* extend(uint index, int value) OVERRIDE;
    virtual TAKE_OWN IMonomial* copy() const OVERRIDE;
    virtual TAKE_OWN IMonomial* next() const OVERRIDE;
    virtual bool divides(const IMonomial* numerator) const OVERRIDE;
    virtual bool isBorderOf(const IMonomial* monomial) const OVERRIDE;

    virtual bool supportsGetPos() const OVERRIDE;
    virtual uint64_t getPos() const OVERRIDE;

    virtual void del() OVERRIDE;

private:
    uint* rep;
    uint indet;
    uint64_t pos;
    uint degree;
    FastFlexibleArray* monomBox;

    friend class MonomialFactoryDegLex;
    friend class FastFlexibleArray;

    DegLexMonomial(uint64_t pos, uint indet, FastFlexibleArray* monomBox);
    DegLexMonomial(uint indet, FastFlexibleArray* monomBox);
    virtual ~DegLexMonomial();

    void recalcPos();
    void recalcDegree();
    void initFromPos(uint64_t pos);
};

class DegLexMonomialNoOrderPos : public IMonomial
{
public:
    virtual const uint& at(uint const& index) const OVERRIDE;
    virtual TAKE_OWN IMonomial* set(uint index, uint value) OVERRIDE;

    virtual uint getIndet() const OVERRIDE;
    virtual uint getDegree() const OVERRIDE;
    virtual TAKE_OWN IMonomial* extend(uint index, int value) OVERRIDE;
    virtual TAKE_OWN IMonomial* copy() const OVERRIDE;
    virtual TAKE_OWN IMonomial* next() const OVERRIDE;
    virtual bool divides(const IMonomial* numerator) const OVERRIDE;
    virtual bool isBorderOf(const IMonomial* monomial) const OVERRIDE;

    virtual bool supportsGetPos() const OVERRIDE;
    virtual uint64_t getPos() const OVERRIDE;

    virtual int compare(const IMonomial* other) const OVERRIDE;

private:
    uint* rep;
    uint indet;
    uint degree;

    friend class MonomialFactoryNoOrderPos;

    DegLexMonomialNoOrderPos(uint indet);
    virtual ~DegLexMonomialNoOrderPos();

    void recalcDegree();
};


} // namespace polynomial

#endif // __DEGLEXMONOMIAL_H__
