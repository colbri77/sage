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

protected:
    friend class MonomialFactoryDegLex;
    friend class FastFlexibleArray;

    DegLexMonomial(uint64_t pos, uint indet, FastFlexibleArray* monomBox);
    DegLexMonomial(uint indet, FastFlexibleArray* monomBox);
    virtual ~DegLexMonomial();

    virtual DegLexMonomial* create(uint indet, FastFlexibleArray* monomBox) const;

    uint* rep;
    uint indet;
    uint64_t pos;
    uint degree;
    FastFlexibleArray* monomBox;

    void recalcPos();
    void recalcDegree();
    void initFromPos(uint64_t pos);
};

class DegLexMonomialGF2 : public DegLexMonomial
{
public:
    virtual TAKE_OWN IMonomial* set(uint index, uint value) OVERRIDE;
    virtual TAKE_OWN IMonomial* extend(uint index, int value) OVERRIDE;
    virtual TAKE_OWN IMonomial* next() const OVERRIDE;

protected:
    friend class MonomialFactoryDegLexGF2;
    friend class FastFlexibleArray;

    virtual DegLexMonomial* create(uint indet, FastFlexibleArray* monomBox) const OVERRIDE;

    DegLexMonomialGF2(uint64_t pos, uint indet, FastFlexibleArray* monomBox);
    DegLexMonomialGF2(uint indet, FastFlexibleArray* monomBox);
    virtual ~DegLexMonomialGF2();
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

    uint* rep;
    uint indet;
    uint degree;

    friend class MonomialFactoryNoOrderPos;

    virtual DegLexMonomialNoOrderPos* create(uint indet) const;

    DegLexMonomialNoOrderPos(uint indet);
    virtual ~DegLexMonomialNoOrderPos();

    void recalcDegree();
};

class DegLexMonomialNoOrderPosGF2 : public DegLexMonomialNoOrderPos
{
public:
    virtual TAKE_OWN IMonomial* set(uint index, uint value) OVERRIDE;
    virtual TAKE_OWN IMonomial* extend(uint index, int value) OVERRIDE;
    virtual TAKE_OWN IMonomial* next() const OVERRIDE;

protected:
    friend class MonomialFactoryNoOrderPosGF2;

    virtual DegLexMonomialNoOrderPos* create(uint indet) const OVERRIDE;

    DegLexMonomialNoOrderPosGF2(uint indet);
    virtual ~DegLexMonomialNoOrderPosGF2();
};


} // namespace polynomial

#endif // __DEGLEXMONOMIAL_H__
