#ifndef MATRIXOPS_HH
#define MATRIXOPS_HH

/// initialization, access, write tests for matrix abstraction class
#include "fiveByFiveMatrix.hh"

#include "UnitTest.hh"
#include <vector>

class matrixOps : public UnitTesting::UnitTest
{
    public:
        matrixOps();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

        aidaTT::fiveByFiveMatrix* _matrix1;
        aidaTT::fiveByFiveMatrix* _matrix2;
        aidaTT::fiveByFiveMatrix* _matrix3;
        std::vector<double> _values;

};
#endif // MATRIXOPS_HH
