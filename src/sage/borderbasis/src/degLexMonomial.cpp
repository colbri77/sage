#include "include/degLexMonomial.h"

#include <cstring>
#include <stdio.h>
#include <iostream>
namespace polynomial {
void printM(IMonomial* m) {
if(m->getDegree()==0) {cout << "1"; return;}
for(uint i=0;i<m->getIndet();i++) {
    if(m->at(i)>1)
        printf("%c^%d",'x'+i,m->at(i));
    else if(m->at(i)>0)
        printf("%c",'x'+i);
}
}
//----- DegLexMonomial ------------------------------------

DegLexMonomial::DegLexMonomial(uint64_t pos, uint indet,FastFlexibleArray* monomBox)
: IMonomial(IMonomial::DEGLEX),
rep(new uint[indet]()),
indet(indet),
pos(pos),
degree(0),
monomBox(monomBox)
{
    initFromPos(pos);
}

DegLexMonomial::DegLexMonomial(uint indet,FastFlexibleArray* monomBox)
: IMonomial(IMonomial::DEGLEX),
rep(new uint[indet]()),
indet(indet),
pos(0),
degree(0),
monomBox(monomBox)
{

}

DegLexMonomial::~DegLexMonomial()
{
    delete rep;
}

const uint& DegLexMonomial::at(uint const& index) const
{
    ENSURE(index < indet, "DegLexMonomial::at(index): index is out of bounds");
    return rep[index];
}

TAKE_OWN IMonomial* DegLexMonomial::set(uint index, uint value)
{
    ENSURE(index < indet, "DegLexMonomial::set(index): index is out of bounds");
    uint repBackup = rep[index];
    uint degBackup = degree;
    uint64_t posBackup = pos;

    degree += (value-rep[index]);
    rep[index] = value;
    recalcPos();

    uint64_t newPos = pos;

    rep[index] = repBackup;
    degree = degBackup;
    pos = posBackup;

    DegLexMonomial* result = (DegLexMonomial*)(monomBox->get(newPos));

    if(result == NULL) {
        result = create(indet,monomBox);

        for(uint i=0;i<indet;i++)
            result->rep[i] = rep[i];
        result->degree = degree + (value-rep[index]);
        result->rep[index] = value;
        result->pos = newPos;

        monomBox->add(newPos,result);
    }

    return (IMonomial*)result;
}

uint DegLexMonomial::getIndet() const
{
    return indet;
}

uint DegLexMonomial::getDegree() const
{
    return degree;
}

DegLexMonomial* DegLexMonomial::create(uint indet, FastFlexibleArray* monomBox) const
{
    return new DegLexMonomial(indet,monomBox);
}

TAKE_OWN IMonomial* DegLexMonomial::extend(uint index, int value)
{
    ENSURE(index < indet, "DegLexMonomial::extend(index,value): index is out of bounds");
    ENSURE(value > 0 || rep[index] >= (uint)(-value), "DegLexMonomial::extend(index,value): value would make index negative");

    uint64_t posBackup = pos;
    uint repBackup = rep[index];
    uint degBackup = degree;

    rep[index] += value;
    degree += value;
    recalcPos();

    uint64_t newPos = pos;
    pos = posBackup;
    degree = degBackup;
    rep[index] = repBackup;

    DegLexMonomial* result = (DegLexMonomial*)(monomBox->get(newPos));

    if(result == NULL) {
        result = create(indet,monomBox);

        for(uint i=0;i<indet;i++)
            result->rep[i] = rep[i];
        result->rep[index] += value;
        result->degree = degree + value;
        result->pos = newPos;

        monomBox->add(newPos,result);
    }

    return (IMonomial*)result;
}

TAKE_OWN IMonomial* DegLexMonomial::copy() const
{
    return (IMonomial*)this;
}

TAKE_OWN IMonomial* DegLexMonomial::next() const
{
    DegLexMonomial* result = (DegLexMonomial*)(monomBox->get(pos+1));

    if(result == NULL) {
        result = create(indet,monomBox);

        for(uint i=0;i<indet;i++)
            result->rep[i] = rep[i];

        if(degree==0) {
            result->rep[indet-1] = 1;
        } else if(indet==1) {
            result->rep[0]++;
        } else if(rep[indet-1]>0) {
            result->rep[indet-1]--;
            result->rep[indet-2]++;
        } else {
            uint nextNonZero = indet-2;
            for(;rep[nextNonZero]==0;nextNonZero--);
            if(nextNonZero==0) {
                result->rep[0] = 0;
                result->rep[indet-1] = degree+1;
            } else {
                result->rep[indet-1] = rep[nextNonZero]-1;
                result->rep[nextNonZero] = 0;
                result->rep[nextNonZero-1] = rep[nextNonZero-1]+1;
            }
        }

        result->pos = pos+1;
        result->recalcDegree();

        monomBox->add(pos+1,result);
    }

    return (IMonomial*)result;
}

bool DegLexMonomial::divides(const IMonomial* numerator) const
{
    ENSURE(indet == numerator->getIndet(), "DegLexMonomial::divides(numerator): monomials have different amounts of variables");
    for(uint i=0;i<indet;i++) {
        if(rep[i]>numerator->at(i))
            return false;
    }
    return true;
}

bool DegLexMonomial::isBorderOf(const IMonomial* monomial) const
{
    ENSURE(indet == monomial->getIndet(), "DegLexMonomial::isBorderOf(monomial): monomials have different amounts of variables");
    if(degree!=monomial->getDegree()+1)
        return false;
    bool foundPos = false;
    for(uint i=0;i<indet;i++) {
        if(rep[i]!=monomial->at(i)) {
            if(foundPos)
                return false;
            else
                foundPos = true;
        }
    }
    return true;
}

bool DegLexMonomial::supportsGetPos() const
{
    return true;
}

uint64_t DegLexMonomial::getPos() const
{
    return pos;
}

void DegLexMonomial::recalcPos()
{
    if(degree==0) {
        pos = 0;
        return;
    }

    uint64_t alphaLast = rep[0];

    uint64_t sigma = indet+degree-1;
    uint64_t gamma = degree;
    uint64_t base = indet;
    uint64_t tmp = indet;

    if(degree==1) {
            base = 0;
    } else {
        for(uint i=1;i<degree-1;i++) {
            tmp *= (indet+i);
            tmp /= (i+1);
            base += tmp;
        }
    }

    uint64_t I = 0;
    uint64_t fLast = 0;

    for(uint i=1;i<indet;i++) {
        for(uint j=0;j<rep[i-1]+2;j++) {
            if(i==1 && j==0) {
                uint64_t a = 1,b = 1;
                for(uint k=1;k<=degree;k++) {
                    a *= (sigma-k);
                    b *= k;
                }
                fLast = a/b;
            }
            else if(j==0) {
                alphaLast = rep[i-2];
                fLast *= (sigma-gamma-1);
                gamma -= alphaLast;
                if(gamma==0) break; // we're done
                fLast /= gamma;
                sigma -= (alphaLast+1);
            }
            else {
                // for verification: this is not the calculation f(j+1)=xyz*f(j), but rather f(j)=xyz*f(j-1)
                fLast *= (gamma-j+1);
                if(sigma-j==0) break;   // we're done - this happens after the last round of i
                fLast /= (sigma-j);
            }
            if(j<rep[i-1])
                I += fLast;
        }
        if(gamma==0) break; // we're done
    }

    pos = base + I + 1;
}

void DegLexMonomial::recalcDegree()
{
    degree = 0;
    for(uint i=0;i<indet;i++)
        degree += rep[i];
}

void DegLexMonomial::initFromPos(uint64_t pos)
{
    if(pos==0)
        return;

    degree = 1;

    uint64_t alphaLast = 0;

    uint64_t base = indet;
    uint64_t tmp = indet;

    for(uint i=1;base<pos;i++,degree++) {
        tmp *= (indet+i);
        tmp /= (i+1);
        base += tmp;
    }
    base -= tmp;

    uint64_t sigma = indet+degree-1;
    uint64_t gamma = degree;
    uint64_t I = pos - base - 1;
    uint64_t ITest = 0;
    uint64_t fLast = 0;
    uint runningDegreeCount = 0;

    for(uint i=1;i<indet;i++) {
        int keepRunning = 2;
        for(uint j=0;keepRunning>0;j++) {
            if(i==1 && j==0) {
                uint64_t a = 1,b = 1;
                for(uint k=1;k<=degree;k++) {
                    a *= (sigma-k);
                    b *= k;
                }
                fLast = a/b;
            }
            else if(j==0) {
                alphaLast = rep[i-2];
                fLast *= (sigma-gamma-1);
                gamma -= alphaLast;
                if(gamma==0) break; // we're done
                fLast /= gamma;
                sigma -= (alphaLast+1);
            }
            else {
                fLast *= (gamma-j+1);
                if(sigma-j==0) break;   // we're done - this happens after the last round of i
                fLast /= (sigma-j);
            }
            if(keepRunning==2) {
                ITest += fLast;
                if(ITest > I) {
                    keepRunning--;
                    ITest -= fLast;
                    rep[i-1] = j;
                    runningDegreeCount += j;
                }
            } else {
                keepRunning--;
            }
        }
        if(gamma==0) break; // we're done
    }

    rep[indet-1] = degree-runningDegreeCount;
}

void DegLexMonomial::del()
{
    // left blank on purpose
}

//----- DegLexMonomialGF2 -------------------------------------------
DegLexMonomialGF2::DegLexMonomialGF2(uint64_t pos, uint indet, FastFlexibleArray* monomBox)
: DegLexMonomial(pos,indet,monomBox)
{

}

DegLexMonomialGF2::DegLexMonomialGF2(uint indet, FastFlexibleArray* monomBox)
: DegLexMonomial(indet,monomBox)
{

}

DegLexMonomialGF2::~DegLexMonomialGF2()
{

}

DegLexMonomial* DegLexMonomialGF2::create(uint indet, FastFlexibleArray* monomBox) const
{
    return new DegLexMonomialGF2(indet,monomBox);
}

TAKE_OWN IMonomial* DegLexMonomialGF2::set(uint index, uint value)
{
    return DegLexMonomial::set(index,value>0 ? 1 : 0);
}

TAKE_OWN IMonomial* DegLexMonomialGF2::extend(uint index, int value)
{
    uint oldValue = DegLexMonomial::at(index);
    if(oldValue==0 && value>0)
        return DegLexMonomial::extend(index,1);
    else if(oldValue>0 && value<0)
        return DegLexMonomial::extend(index,-1);
    return this;
}

TAKE_OWN IMonomial* DegLexMonomialGF2::next() const
{
    DegLexMonomialGF2* result = NULL;
    bool fitting = false;

    for(uint64_t curPos=pos+1;!fitting;curPos++) {
        result = (DegLexMonomialGF2*)(monomBox->get(curPos));

        if(result == NULL) {
            result = new DegLexMonomialGF2(curPos,indet,monomBox);
            monomBox->add(curPos,result);
        }

        fitting = true;
        for(uint i=0;i<indet && fitting;i++) {
            if(result->at(i)>1)
                fitting = false;
        }
        if(result->getDegree()>indet)
            ASSERT_NOT_REACHED;
    }
    return (IMonomial*)result;
}

//----- DegLexMonomialNoOrderPos ------------------------------------

DegLexMonomialNoOrderPos::DegLexMonomialNoOrderPos(uint indet)
: IMonomial(IMonomial::DEGLEX),
rep(new uint[indet]()),
indet(indet),
degree(0)
{

}

DegLexMonomialNoOrderPos::~DegLexMonomialNoOrderPos()
{
    delete rep;
}

DegLexMonomialNoOrderPos* DegLexMonomialNoOrderPos::create(uint indet) const
{
    return new DegLexMonomialNoOrderPos(indet);
}

const uint& DegLexMonomialNoOrderPos::at(uint const& index) const
{
    ENSURE(index < indet, "DegLexMonomial::at(index): index is out of bounds");
    return rep[index];
}

TAKE_OWN IMonomial* DegLexMonomialNoOrderPos::set(uint index, uint value)
{
    ENSURE(index < indet, "DegLexMonomial::set(index): index is out of bounds");
    degree += (value-rep[index]);
    rep[index] = value;
    return this;
}

uint DegLexMonomialNoOrderPos::getIndet() const
{
    return indet;
}

uint DegLexMonomialNoOrderPos::getDegree() const
{
    return degree;
}

TAKE_OWN IMonomial* DegLexMonomialNoOrderPos::extend(uint index, int value)
{
    ENSURE(index < indet, "DegLexMonomial::extend(index,value): index is out of bounds");
    ENSURE(value > 0 || rep[index] >= (uint)(-value), "DegLexMonomial::extend(index,value): value would make index negative");
    rep[index] += value;
    degree += value;
    return this;
}

TAKE_OWN IMonomial* DegLexMonomialNoOrderPos::copy() const
{
    DegLexMonomialNoOrderPos* result = create(indet);

    for(uint i=0;i<indet;i++)
        result->rep[i] = rep[i];

    result->degree = degree;
    return result;
}

TAKE_OWN IMonomial* DegLexMonomialNoOrderPos::next() const
{
    DegLexMonomialNoOrderPos* result = (DegLexMonomialNoOrderPos*)copy();

    if(degree==0) {
        result->rep[indet-1] = 1;
    } else if(indet==1) {
        result->rep[0]++;
    } else if(rep[indet-1]>0) {
        result->rep[indet-1]--;
        result->rep[indet-2]++;
    } else {
        uint nextNonZero = indet-2;
        for(;rep[nextNonZero]==0;nextNonZero--);
        if(nextNonZero==0) {
            result->rep[0] = 0;
            result->rep[indet-1] = degree+1;
        } else {
            result->rep[indet-1] = rep[nextNonZero]-1;
            result->rep[nextNonZero] = 0;
            result->rep[nextNonZero-1] = rep[nextNonZero-1]+1;
        }
    }

    result->recalcDegree();

    return result;
}

bool DegLexMonomialNoOrderPos::divides(const IMonomial* numerator) const
{
    ENSURE(indet == numerator->getIndet(), "DegLexMonomial::divides(numerator): monomials have different amounts of variables");
    for(uint i=0;i<indet;i++) {
        if(rep[i]>numerator->at(i))
            return false;
    }
    return true;
}

bool DegLexMonomialNoOrderPos::isBorderOf(const IMonomial* monomial) const
{
    ENSURE(indet == monomial->getIndet(), "DegLexMonomial::isBorderOf(monomial): monomials have different amounts of variables");
    if(degree!=monomial->getDegree()+1)
        return false;
    bool foundPos = false;
    for(uint i=0;i<indet;i++) {
        if(rep[i]!=monomial->at(i)) {
            if(foundPos)
                return false;
            else
                foundPos = true;
        }
    }
    return true;
}

bool DegLexMonomialNoOrderPos::supportsGetPos() const
{
    return false;
}

uint64_t DegLexMonomialNoOrderPos::getPos() const
{
    NOT_IMPLEMENTED;
    return 0;
}

void DegLexMonomialNoOrderPos::recalcDegree()
{
    degree = 0;
    for(uint i=0;i<indet;i++)
        degree += rep[i];
}

int DegLexMonomialNoOrderPos::compare(const IMonomial* other) const
{
    ENSURE(termOrdering == other->termOrdering, "IMonomial::compare(other): monomials have different term ordering");
    ENSURE(getIndet() == other->getIndet(), "IMonomial::compare(other): monomials have a different amount of variables");

    if(degree!=other->getDegree())
        return ((int)degree)-((int)other->getDegree());

    for(uint i=0;i<indet;i++) {
        int a = (int)rep[i];
        int b = (int)other->at(i);
        if(a!=b)
            return a-b;
    }

    return 0;
}

//----- DegLexMonomialNoOrderPosGF2 -------------------------------------------

DegLexMonomialNoOrderPosGF2::DegLexMonomialNoOrderPosGF2(uint indet)
: DegLexMonomialNoOrderPos(indet)
{

}

DegLexMonomialNoOrderPosGF2::~DegLexMonomialNoOrderPosGF2()
{

}

DegLexMonomialNoOrderPos* DegLexMonomialNoOrderPosGF2::create(uint indet) const
{
    return new DegLexMonomialNoOrderPosGF2(indet);
}


TAKE_OWN IMonomial* DegLexMonomialNoOrderPosGF2::set(uint index, uint value)
{
    return DegLexMonomialNoOrderPos::set(index,value>0 ? 1 : 0);
}

TAKE_OWN IMonomial* DegLexMonomialNoOrderPosGF2::extend(uint index, int value)
{
    uint oldValue = DegLexMonomialNoOrderPos::at(index);
    if(oldValue==0 && value>0)
        return DegLexMonomialNoOrderPos::extend(index,1);
    else if(oldValue>0 && value<0)
        return DegLexMonomialNoOrderPos::extend(index,-1);
    return this;
}

TAKE_OWN IMonomial* DegLexMonomialNoOrderPosGF2::next() const
{
    DegLexMonomialNoOrderPosGF2* result = (DegLexMonomialNoOrderPosGF2*)copy();

    if(degree==0) {
        result->rep[indet-1] = 1;
    } else {
        int oneCtr = 0;
        int i = (int)indet-1;
        for(;i>=0;i--) {
            if(result->rep[i]==1) {
                oneCtr++;
            } else if(oneCtr>0) {
                oneCtr--;
                break;
            }
        }
        if(i<0)
            ASSERT_NOT_REACHED;
        result->rep[i]++;
        for(int k=(int)indet-1;k>i;k++,oneCtr--) {
            if(oneCtr>0)
                result->rep[k] = 1;
            else
                result->rep[k] = 0;
        }
    }

    result->recalcDegree();

    return result;
}

} // namespace polynomial
