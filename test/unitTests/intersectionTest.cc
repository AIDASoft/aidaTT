#include "intersectionTest.hh"

using namespace std;
using namespace aidaTT;

intersectionTest::intersectionTest() : UnitTest("IntersectionCalculationTest", __FILE__)
{
    c1 = new circle(0., 0., 1.);
    c2 = new circle(1., 1., 1.);
    sL1 = new straightLine(0., 1., .5);
    sL2 = new straightLine(1., 1., 0.);
}



void intersectionTest::_test()
{
    /// check circle and circle
    intersections sectCC = intersectCircleCircle(*c1, *c2);
    test_(sectCC.number() == 2);
    test_(floatCompare(sectCC[0].first, 1.));
    test_(floatCompare(sectCC[0].second, 0.));
    test_(floatCompare(sectCC[1].first, 0.));
    test_(floatCompare(sectCC[1].second, 1.));

    /// check circle and line
    intersections sectCL = intersectCircleStraightLine(*c1, *sL1);
    test_(sectCL.number() == 2);
    test_(floatCompare(sectCL[0].first, sqrt(.75)));
    test_(floatCompare(sectCL[0].second, 0.5));
    test_(floatCompare(sectCL[1].first, -sqrt(.75)));
    test_(floatCompare(sectCL[1].second, 0.5));

    /// check line and line
    intersections sectLL = intersectStraightLineStraightLine(*sL2, *sL1);
    test_(sectLL.number() == 1);
    test_(floatCompare(sectLL[0].first, -0.5));
    test_(floatCompare(sectLL[0].second, 0.5));
    //    test_(floatCompare(_one->getYPerp(),     _one->getTrackParameters()(4)));

}



void intersectionTest::run()
{
    _test();
}
