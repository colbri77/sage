#ifndef __MONOMIALFACTORY_H__
#define __MONOMIALFACTORY_H__

#include "i_monomial.h"
#include "fastFlexibleArray.h"

namespace polynomial {

class IMonomialFactory
{
public:
    IMonomialFactory(){}
    virtual ~IMonomialFactory(){}

    virtual bool supportsGetPos() const = 0;
    virtual TAKE_OWN IMonomial* create() const = 0;
    virtual TAKE_OWN IMonomial* create(uint64_t pos) const = 0;
    virtual void excludeIndet(uint index) {NOT_IMPLEMENTED;}
};

class MonomialFactoryNoOrderPos : public IMonomialFactory
{
public:
    MonomialFactoryNoOrderPos(uint indet);
    virtual ~MonomialFactoryNoOrderPos();

    virtual bool supportsGetPos() const OVERRIDE;

    virtual TAKE_OWN IMonomial* create() const OVERRIDE;
    virtual TAKE_OWN IMonomial* create(uint64_t pos) const OVERRIDE;

protected:
    uint indet;
};

class MonomialFactoryDegLex : public IMonomialFactory
{
public:
    MonomialFactoryDegLex(uint indet);
    virtual ~MonomialFactoryDegLex();

    virtual bool supportsGetPos() const OVERRIDE;

    virtual TAKE_OWN IMonomial* create() const OVERRIDE;
    virtual TAKE_OWN IMonomial* create(uint64_t pos) const OVERRIDE;

protected:
    MonomialFactoryDegLex();

    FastFlexibleArray* monomials;
    uint indet;
};

class MonomialFactoryNoOrderPosGF2 : public MonomialFactoryNoOrderPos
{
public:
    MonomialFactoryNoOrderPosGF2(uint indet);
    virtual ~MonomialFactoryNoOrderPosGF2();

    virtual TAKE_OWN IMonomial* create() const OVERRIDE;
    virtual TAKE_OWN IMonomial* create(uint64_t pos) const OVERRIDE;

    virtual void excludeIndet(uint index) OVERRIDE;

private:
    bool* excludedIndets;
};

class MonomialFactoryDegLexGF2 : public MonomialFactoryDegLex
{
public:
    MonomialFactoryDegLexGF2(uint indet);
    virtual ~MonomialFactoryDegLexGF2();

    virtual TAKE_OWN IMonomial* create(uint64_t pos) const OVERRIDE;

    virtual void excludeIndet(uint index) OVERRIDE;

protected:
    MonomialFactoryDegLexGF2();

    bool* excludedIndets;
};

class MonomialFactoryDegRevLex : public MonomialFactoryDegLex
{
public:
    MonomialFactoryDegRevLex(uint indet);
    virtual ~MonomialFactoryDegRevLex();

    virtual TAKE_OWN IMonomial* create(uint64_t pos) const OVERRIDE;
};

class MonomialFactoryDegRevLexGF2 : public MonomialFactoryDegLexGF2
{
public:
    MonomialFactoryDegRevLexGF2(uint indet);
    virtual ~MonomialFactoryDegRevLexGF2();

    virtual TAKE_OWN IMonomial* create(uint64_t pos) const OVERRIDE;
};


} // namespace polynomial

#endif // __MONOMIALFACTORY_H__
