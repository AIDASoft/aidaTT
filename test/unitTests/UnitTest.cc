#include "UnitTest.hh"
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <cmath>

using namespace std;
using namespace UnitTesting;

void UnitTest::do_test(bool cond, const std::string& lbl, const char* fname, long lineno)
{
    if(!cond)
        do_fail(lbl, fname, lineno);
    else
        succeed_();
}



void UnitTest::do_fail(const std::string& lbl, const char* /* fname */ , long lineno)
{
    ++nFail;
    if(osptr)
        {
            *osptr << "failure: " << setw(76) << std::left << lbl << std::right << " [line " << setw(3) << lineno << "]\n";
        }
}



bool UnitTest::floatCompare(float x1, float x2)
{
    bool ret = (fabs(x1 - x2) < 1.e-9); // fix me, empirical number
    if(!ret) cout << "warning floatCompare: x1=" << x1 << " x2=" << x2 << ", diff=" << x1 - x2 << endl;
    return ret;
}



bool UnitTest::roughFloatCompare(float x1, float x2)
{
    bool ret = (fabs(x1 - x2) < 1.e-5); // fix me, empirical number
    if(!ret) cout << "warning roughFloatCompare: x1=" << x1 << " x2=" << x2 << ", diff=" << x1 - x2 << endl;
    return ret;
}



const std::string & UnitTest::name() const
{
    if(fName == "")
        fName = std::string(typeid(*this).name());
    return fName;
}



long UnitTest::report() const
{
    return nFail;
}



void UnitTest::intro() const
{
    if(osptr)
        {
            *osptr << "************************************************************************************************" << endl;
            *osptr << "* Starting U N I T test : " << name() << " (file " << fFileName << ")" << endl;
            *osptr << "************************************************************************************************" << endl;
        }
}
