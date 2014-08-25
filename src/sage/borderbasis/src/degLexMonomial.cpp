#include "include/degLexMonomial.h"

#include <cstring>

namespace polynomial {

DegLexMonomial::DegLexMonomial(uint64_t pos, uint indet)
: IMonomial(IMonomial::DEGLEX),
rep(new uint[indet]()),
indet(indet),
pos(0),
degree(0)
{
    initFromPos(pos);
}

DegLexMonomial::DegLexMonomial(uint indet)
: IMonomial(IMonomial::DEGLEX),
rep(new uint[indet]()),
indet(indet),
pos(0),
degree(0)
{

}

DegLexMonomial::DegLexMonomial(uint values[], uint indet)
: IMonomial(IMonomial::DEGLEX),
rep(new uint[indet]),
indet(indet),
pos(0),
degree(0)
{
    memcpy(rep,values,sizeof(uint)*indet);
    recalcDegree();
    recalcPos();
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

void DegLexMonomial::set(uint index, uint value)
{
    ENSURE(index < indet, "DegLexMonomial::set(index): index is out of bounds");
    degree += (value-rep[index]);
    rep[index] = value;
    recalcPos();
}

uint DegLexMonomial::getIndet() const
{
    return indet;
}

uint DegLexMonomial::getDegree() const
{
    return degree;
}

void DegLexMonomial::extend(uint index, int value)
{
    ENSURE(index < indet, "DegLexMonomial::extend(index,value): index is out of bounds");
    ENSURE(value > 0 || rep[index] >= (uint)(-value), "DegLexMonomial::extend(index,value): value would make index negative");
    rep[index] += value;
    degree += value;
    recalcPos();
}

TAKE_OWN IMonomial* DegLexMonomial::copy() const
{
    DegLexMonomial* result = new DegLexMonomial(indet);

    for(uint i=0;i<indet;i++)
        result->rep[i] = rep[i];

    result->pos = pos;
    result->degree = degree;
    return result;
}

TAKE_OWN IMonomial* DegLexMonomial::next() const
{
    DegLexMonomial* result = (DegLexMonomial*)copy();

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

    return result;
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

} // namespace polynomial
