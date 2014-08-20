#include "include/statistics.h"

namespace borderbasis {

Statistics::Statistics()
{
    maxMatrix.rows = 0;
    maxMatrix.columns = 0;
}

Statistics::~Statistics()
{

}

void Statistics::start()
{
    maxMatrix.rows = 0;
    maxMatrix.columns = 0;
}

void Statistics::stop()
{

}

void Statistics::logMatrix(uint r,uint c)
{
    // convert to uint64_t to prevent integer overflow on multiplication
    uint64_t mr = maxMatrix.rows;
    uint64_t mc = maxMatrix.columns;
    uint64_t cr = r;
    uint64_t cc = c;
    if(mr*mc < cr*cc) {
        maxMatrix.rows = r;
        maxMatrix.columns = c;
    }
}

} // namespace borderbasis
