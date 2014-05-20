#include "intersectionTest.hh"

using namespace std;
using namespace aidaTT;

intersectionTest::intersectionTest() : UnitTest("IntersectionCalculationTest", __FILE__)
{
    c1 = new circle(0., 0., 1.);
    sL1 = new straightLine(0., 1., .5);
}



void intersectionTest::_test()
{

    intersections sect = intersectCircleStraightLine(*c1, *sL1);
    test_(sect.number() == 2);
    test_(floatCompare(sect[0].first, sqrt(.75)));
    test_(floatCompare(sect[0].second, 0.5));
    test_(floatCompare(sect[1].first, -sqrt(.75)));
    test_(floatCompare(sect[1].second, 0.5));
//    test_(floatCompare(_one->getYPerp(),     _one->getTrackParameters()(4)));

}



void intersectionTest::run()
{
    _test();
}
