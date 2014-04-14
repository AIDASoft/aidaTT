#ifndef IFITTINGALGORITHM_HH
#define IFITTINGALGORITHM_HH

#include "trajectory.hh"

namespace aidaTT
{
    class trajectory; // needs forward declaration to build

    class IFittingAlgorithm
    {
        public:
            virtual bool const fit(const trajectory&)       = 0;
            virtual unsigned int const getNDF()             = 0;
            virtual double const getChiSquare()              = 0;
            virtual double const lostWeight()                = 0;
    };
}

#endif // IFITTINGALGORITHM_HH
