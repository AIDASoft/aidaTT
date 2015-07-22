#ifndef FINALTRACKTEST_HH
#define FINALTRACKTEST_HH

/// initialization, access, write tests for matrix and vector abstraction class

#include "AidaTT.hh"
#include "DD4hep/DD4hepUnits.h"
#include "lcio.h"
#include <IMPL/TrackImpl.h>
#include "EVENT/TrackState.h"
#include "UnitTest.hh"
#include "LCIOPersistency.hh"

class finalTrackTest : public UnitTesting::UnitTest
{
    public:
        finalTrackTest();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

};
#endif // FINALTRACKTEST_HH
