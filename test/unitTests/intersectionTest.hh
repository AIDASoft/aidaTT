#ifndef INTERSECTIONTEST_HH
#define INTERSECTIONTEST_HH

/// simple track parameter initialization and read back test
#include "intersections.hh"

#include "UnitTest.hh"
#include <vector>

class intersectionTest : public UnitTesting::UnitTest
{
    public:
        intersectionTest();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

        aidaTT::circle* c1;
        aidaTT::circle* c2;
        aidaTT::circle* c3;
        aidaTT::straightLine* sL1;
        aidaTT::straightLine* sL2;
        aidaTT::straightLine* sL3;
        aidaTT::straightLine* sL4;

};
#endif // INTERSECTIONTEST_HH
