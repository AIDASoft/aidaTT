#ifndef IFITTINGALGORITHM_HH
#define IFITTINGALGORITHM_HH

#include "trajectory.hh"

namespace aidaTT
{
    class trajectory; // needs forward declaration to build

    class IFittingAlgorithm
    {
        public:
            bool const fit(const trajectory&);
            unsigned int const getNDF();
            double const getChiSquare();
    };

}

#endif // IFITTINGALGORITHM_HH
