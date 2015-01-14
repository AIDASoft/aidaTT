#ifndef INITIALTRACKTEST_HH
#define INITIALTRACKTEST_HH

/// initialization, access, write tests for matrix and vector abstraction class

#include "AidaTT.hh"
#include "DD4hep/DD4hepUnits.h"
#include "lcio.h"
#include <IMPL/TrackImpl.h>
#include "EVENT/TrackState.h"
#include "UnitTest.hh"
#include "LCIOPersistency.hh"

class initialTrackTest : public UnitTesting::UnitTest
{
    public:
        initialTrackTest();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

};
#endif // INITIALTRACKTEST_HH
