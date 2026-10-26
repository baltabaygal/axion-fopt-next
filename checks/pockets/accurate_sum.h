#ifndef POCKETS_ACCURATE_SUM_H
#define POCKETS_ACCURATE_SUM_H
#include <cmath>

// Neumaier compensated host-side reduction. Field evolution is untouched.
// Needed when an O(N^3) diagnostic sum is compared with pairwise NumPy sums
// at a fixed 1e-12 relative accounting tolerance on large grids.
class AccurateSum {
    double value_=0.,correction_=0.;
public:
    AccurateSum& operator+=(double increment) {
        const double next=value_+increment;
        correction_+=std::abs(value_)>=std::abs(increment)
            ? (value_-next)+increment : (increment-next)+value_;
        value_=next;
        return *this;
    }
    operator double() const {return value_+correction_;}
};
#endif
