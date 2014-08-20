#ifndef __DEGLEXMONOMIAL_H__
#define __DEGLEXMONOMIAL_H__

#include "i_monomial.h"

namespace polynomial {

class DegLexMonomial : public IMonomial<uint64_t>
{
public:
    DegLexMonomial(uint64_t pos, uint indet);
    DegLexMonomial(uint indet);
    DegLexMonomial(uint values[], uint indet);
    virtual ~DegLexMonomial();

    virtual const uint& at(uint const& index) const OVERRIDE;
    virtual void set(uint index, uint value) OVERRIDE;

    virtual uint getIndet() const OVERRIDE;
    virtual uint getDegree() const OVERRIDE;
    virtual void extend(uint index, int value) OVERRIDE;
    virtual TAKE_OWN IMonomial<uint64_t>* copy() const OVERRIDE;
    virtual TAKE_OWN IMonomial<uint64_t>* next() const OVERRIDE;
    virtual bool divides(const IMonomial<uint64_t>* numerator) const OVERRIDE;
    virtual bool isBorderOf(const IMonomial<uint64_t>* monomial) const OVERRIDE;

    virtual uint64_t getPos() const OVERRIDE;

private:
    uint* rep;
    uint indet;
    uint64_t pos;
    uint degree;

    void recalcPos();
    void recalcDegree();
    void initFromPos(uint64_t pos);
};

} // namespace polynomial

#endif // __DEGLEXMONOMIAL_H__
