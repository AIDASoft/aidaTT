#ifndef IFITTINGALGORITHM_HH
#define IFITTINGALGORITHM_HH

#include "fitResults.hh"

namespace aidaTT
{
    class trajectory; // needs forward declaration to build

    class IFittingAlgorithm
    {
        public:
            virtual bool fit(const trajectory&)   = 0;
            virtual const fitResults& getResults() const  = 0;
            virtual ~IFittingAlgorithm(){}
    };
}

#endif // IFITTINGALGORITHM_HH
