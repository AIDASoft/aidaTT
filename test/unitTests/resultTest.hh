#ifndef RESULTTEST_HH
#define RESULTTEST_HH

/// initialization, access, write tests for matrix and vector abstraction class
#include "fitResults.hh"

#include "UnitTest.hh"

class resultTest : public UnitTesting::UnitTest
{
    public:
        resultTest();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

        aidaTT::fitResults* _fR;
};
#endif // RESULTTEST_HH
