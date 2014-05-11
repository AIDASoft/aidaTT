#ifndef IFITTINGALGORITHM_HH
#define IFITTINGALGORITHM_HH

#include "trajectory.hh"

namespace aidaTT
{
    class trajectory; // needs forward declaration to build

    class IFittingAlgorithm
    {
        public:
            virtual bool fit(const trajectory&) const       = 0;
            virtual unsigned int getNDF() const             = 0;
            virtual double getChiSquare() const             = 0;
            virtual double lostWeight() const               = 0;
    };
}

#endif // IFITTINGALGORITHM_HH
