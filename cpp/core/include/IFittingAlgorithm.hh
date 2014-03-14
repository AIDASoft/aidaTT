#ifndef IFITTINGALGORITHM_HH
#define IFITTINGALGORITHM_HH

#include "trajectory.hh"

namespace aidaTT
{
    class trajectory; // needs forward declaration to build

    class IFittingAlgorithm
    {
        public:
            virtual bool initializeFitter(const trajectory&) = 0; 
            virtual bool const fit() = 0;
            virtual unsigned int const getNDF() = 0;
            virtual double const getChiSquare() = 0;
    };
}

#endif // IFITTINGALGORITHM_HH
