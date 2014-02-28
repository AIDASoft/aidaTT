#include "UnitTestSuite.hh"
#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;
using namespace UnitTesting;

void UnitTestSuite::addTest(UnitTest* t) throw(UnitTestSuiteError)
{
    // Verify test is valid and has a stream:
    if(t == 0)
        throw UnitTestSuiteError("Null test in UnitTestSuite::addTest");
    else if(osptr && !t->getStream())
        t->setStream(osptr);
    tests.push_back(t);
    t->reset();
}



void UnitTestSuite::addSuite(const UnitTestSuite& s)
{
    for(size_t i = 0; i < s.tests.size(); ++i)
        {
            assert(tests[i]);
            addTest(s.tests[i]);
        }
}



void UnitTestSuite::free()
{
    for(size_t i = 0; i < tests.size(); ++i)
        {
            delete tests[i];
            tests[i] = 0;
        }
}



void UnitTestSuite::run()
{
    reset();
    for(size_t i = 0; i < tests.size(); ++i)
        {
            assert(tests[i]);
            tests[i]->intro();
            tests[i]->run();
            tests[i]->report();
        }
}



void UnitTestSuite::intro() const
{
    if(osptr)
        {
            *osptr << "******************************************************************" << endl;
            *osptr << "******************* AIDATT U N I T test suite ********************" << endl;
            *osptr << "******************************************************************" << endl;
        }
}



long UnitTestSuite::report() const
{
    if(osptr)
        {
            long totFail = 0, totSuccess = 0;
            *osptr << "************************************************************************************************" << endl;
            size_t i;
            for(i = 0; i < tests.size(); ++i)
                {
                    assert(tests[i]);
                    *osptr << " Test " << tests[i]->name() << ": ran with " << tests[i]->getNumPassed() << " success(es) and " << tests[i]->getNumFailed() << " fail(s)." << endl;

                    totFail    += tests[i]->getNumFailed();
                    totSuccess += tests[i]->getNumPassed();
                }
            *osptr << "************************************************************************************************" << endl;
            *osptr << "************************************************************************************************" << endl;

            *osptr << " SUMMARY: " << totFail << " FAILS and " << totSuccess << " successes in total." << endl;
            *osptr << "************************************************************************************************" << endl;
            *osptr << "************************************************************************************************" << endl;

            return totFail;
        }
    else
        return getNumFailed();
}



long UnitTestSuite::getNumPassed() const
{
    long totPass = 0;
    for(size_t i = 0; i < tests.size(); ++i)
        {
            assert(tests[i]);
            totPass += tests[i]->getNumPassed();
        }
    return totPass;
}



long UnitTestSuite::getNumFailed() const
{
    long totFail = 0;
    for(size_t i = 0; i < tests.size(); ++i)
        {
            assert(tests[i]);
            totFail += tests[i]->getNumFailed();
        }
    return totFail;
}



void UnitTestSuite::reset()
{
    for(size_t i = 0; i < tests.size(); ++i)
        {
            assert(tests[i]);
            tests[i]->reset();
        }
}
