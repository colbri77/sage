#include "include/degRevLexMonomial.h"

namespace polynomial {

//-----DegRevLexMonomial----------------------------------
DegRevLexMonomial::DegRevLexMonomial(uint64_t pos, uint indet, FastFlexibleArray* monomBox)
: DegLexMonomial(indet,monomBox)
{
    IMonomial::termOrdering = IMonomial::DEGREVLEX;
    DegLexMonomial::pos = pos;
    initFromPos(pos);
    DegLexMonomial::recalcSingleVarIndex();
}

DegRevLexMonomial::DegRevLexMonomial(uint indet, FastFlexibleArray* monomBox)
: DegLexMonomial(indet,monomBox)
{
    IMonomial::termOrdering = IMonomial::DEGREVLEX;
}

DegRevLexMonomial::~DegRevLexMonomial()
{

}

DegLexMonomial* DegRevLexMonomial::create(uint indet, FastFlexibleArray* monomBox) const
{
    return new DegRevLexMonomial(indet,monomBox);
}

void DegRevLexMonomial::recalcPos()
{
    if(DegLexMonomial::degree==0) {
        pos = 0;
        return;
    }

    uint64_t alphaLast = DegLexMonomial::rep[indet-1];

    uint64_t sigma = DegLexMonomial::indet+DegLexMonomial::degree-1;
    uint64_t gamma = DegLexMonomial::degree;
    uint64_t base = DegLexMonomial::indet;
    uint64_t tmp = DegLexMonomial::indet;

    for(uint i=1;i<DegLexMonomial::degree;i++) {
        tmp *= (DegLexMonomial::indet+i);
        tmp /= (i+1);
        base += tmp;
    }

    uint64_t I = 0;
    uint64_t fLast = 0;

    for(uint i=1;i<DegLexMonomial::indet;i++) {
        for(uint j=0;j<DegLexMonomial::rep[DegLexMonomial::indet-i]+2;j++) {
            if(i==1 && j==0) {
                uint64_t a = 1,b = 1;
                for(uint k=1;k<=DegLexMonomial::degree;k++) {
                    a *= (sigma-k);
                    b *= k;
                }
                fLast = a/b;
            }
            else if(j==0) {
                alphaLast = DegLexMonomial::rep[DegLexMonomial::indet-i+1];
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
            if(j<DegLexMonomial::rep[DegLexMonomial::indet-i])
                I += fLast;
        }
        if(gamma==0) break; // we're done
    }

    DegLexMonomial::pos = base - I;
}

void DegRevLexMonomial::initFromPos(uint64_t pos)
{
    if(pos==0)
        return;

    DegLexMonomial::degree = 1;

    uint64_t alphaLast = 0;

    uint64_t base = DegLexMonomial::indet;
    uint64_t tmp = DegLexMonomial::indet;

    for(uint i=1;base<pos;i++,DegLexMonomial::degree++) {
        tmp *= (DegLexMonomial::indet+i);
        tmp /= (i+1);
        base += tmp;
    }

    uint64_t sigma = DegLexMonomial::indet+DegLexMonomial::degree-1;
    uint64_t gamma = DegLexMonomial::degree;
    uint64_t I = base - pos;
    uint64_t ITest = 0;
    uint64_t fLast = 0;
    uint runningDegreeCount = 0;

    for(uint i=1;i<DegLexMonomial::indet;i++) {
        int keepRunning = 2;
        for(uint j=0;keepRunning>0;j++) {
            if(i==1 && j==0) {
                uint64_t a = 1,b = 1;
                for(uint k=1;k<=DegLexMonomial::degree;k++) {
                    a *= (sigma-k);
                    b *= k;
                }
                fLast = a/b;
            }
            else if(j==0) {
                alphaLast = DegLexMonomial::rep[DegLexMonomial::indet-i+1];
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
                    DegLexMonomial::rep[DegLexMonomial::indet-i] = j;
                    runningDegreeCount += j;
                }
            } else {
                keepRunning--;
            }
        }
        if(gamma==0) break; // we're done
    }

    DegLexMonomial::rep[0] = DegLexMonomial::degree-runningDegreeCount;
}

//-----DegRevLexMonomialGF2-------------------------------

DegRevLexMonomialGF2::DegRevLexMonomialGF2(uint64_t pos, uint indet, FastFlexibleArray* monomBox,bool* excludedIndets)
: DegLexMonomialGF2(indet,monomBox,excludedIndets)
{
    IMonomial::termOrdering = IMonomial::DEGREVLEX;
    DegLexMonomialGF2::pos = pos;
    initFromPos(pos);
    DegLexMonomialGF2::recalcSingleVarIndex();
}

DegRevLexMonomialGF2::DegRevLexMonomialGF2(uint indet, FastFlexibleArray* monomBox,bool* excludedIndets)
: DegLexMonomialGF2(indet,monomBox,excludedIndets)
{
    IMonomial::termOrdering = IMonomial::DEGREVLEX;
}

DegRevLexMonomialGF2::~DegRevLexMonomialGF2()
{

}

DegLexMonomial* DegRevLexMonomialGF2::create(uint indet, FastFlexibleArray* monomBox) const
{
    return new DegRevLexMonomialGF2(indet,monomBox,DegLexMonomialGF2::excludedIndets);
}

void DegRevLexMonomialGF2::recalcPos()
{
    if(DegLexMonomial::degree==0) {
        pos = 0;
        return;
    }

    uint64_t alphaLast = DegLexMonomial::rep[indet-1];

    uint64_t sigma = DegLexMonomial::indet+DegLexMonomial::degree-1;
    uint64_t gamma = DegLexMonomial::degree;
    uint64_t base = DegLexMonomial::indet;
    uint64_t tmp = DegLexMonomial::indet;

    for(uint i=1;i<DegLexMonomial::degree;i++) {
        tmp *= (DegLexMonomial::indet+i);
        tmp /= (i+1);
        base += tmp;
    }

    uint64_t I = 0;
    uint64_t fLast = 0;

    for(uint i=1;i<DegLexMonomial::indet;i++) {
        for(uint j=0;j<DegLexMonomial::rep[DegLexMonomial::indet-i]+2;j++) {
            if(i==1 && j==0) {
                uint64_t a = 1,b = 1;
                for(uint k=1;k<=DegLexMonomial::degree;k++) {
                    a *= (sigma-k);
                    b *= k;
                }
                fLast = a/b;
            }
            else if(j==0) {
                alphaLast = DegLexMonomial::rep[DegLexMonomial::indet-i+1];
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
            if(j<DegLexMonomial::rep[DegLexMonomial::indet-i])
                I += fLast;
        }
        if(gamma==0) break; // we're done
    }

    DegLexMonomial::pos = base - I;
}

void DegRevLexMonomialGF2::initFromPos(uint64_t pos)
{
    if(pos==0)
        return;

    DegLexMonomial::degree = 1;

    uint64_t alphaLast = 0;

    uint64_t base = DegLexMonomial::indet;
    uint64_t tmp = DegLexMonomial::indet;

    for(uint i=1;base<pos;i++,DegLexMonomial::degree++) {
        tmp *= (DegLexMonomial::indet+i);
        tmp /= (i+1);
        base += tmp;
    }

    uint64_t sigma = DegLexMonomial::indet+DegLexMonomial::degree-1;
    uint64_t gamma = DegLexMonomial::degree;
    uint64_t I = base - pos;
    uint64_t ITest = 0;
    uint64_t fLast = 0;
    uint runningDegreeCount = 0;

    for(uint i=1;i<DegLexMonomial::indet;i++) {
        int keepRunning = 2;
        for(uint j=0;keepRunning>0;j++) {
            if(i==1 && j==0) {
                uint64_t a = 1,b = 1;
                for(uint k=1;k<=DegLexMonomial::degree;k++) {
                    a *= (sigma-k);
                    b *= k;
                }
                fLast = a/b;
            }
            else if(j==0) {
                alphaLast = DegLexMonomial::rep[DegLexMonomial::indet-i+1];
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
                    DegLexMonomial::rep[DegLexMonomial::indet-i] = j;
                    runningDegreeCount += j;
                }
            } else {
                keepRunning--;
            }
        }
        if(gamma==0) break; // we're done
    }

    DegLexMonomial::rep[0] = DegLexMonomial::degree-runningDegreeCount;
}

} // namespace polynomial

