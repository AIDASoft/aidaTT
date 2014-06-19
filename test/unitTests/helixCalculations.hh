#ifndef HELIXCALCULATIONS_HH
#define HELIXCALCULATIONS_HH

/// simple track parameter initialization and read back test
#include "trackParameters.hh"

#include "UnitTest.hh"
#include <vector>

class helixCalculations : public UnitTesting::UnitTest
{
    public:
        helixCalculations();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

        aidaTT::trackParameters* _one;
        aidaTT::trackParameters* _two;
        aidaTT::trackParameters* _three;

        aidaTT::fullCovariance* _covmatrix;
        aidaTT::Vector3D* _refpoint ;
        aidaTT::Vector5* _helix ;

};
#endif // HELIXCALCULATIONS_HH
