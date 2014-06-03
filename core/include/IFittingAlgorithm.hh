#ifndef IFITTINGALGORITHM_HH
#define IFITTINGALGORITHM_HH

#include "trajectory.hh"
#include "fitResults.hh"

namespace aidaTT
{
    class trajectory; // needs forward declaration to build

    class IFittingAlgorithm
    {
        public:
            virtual bool fit(const trajectory&)   = 0;
            virtual const fitResults& getResults() const  = 0;
    };
}

#endif // IFITTINGALGORITHM_HH
