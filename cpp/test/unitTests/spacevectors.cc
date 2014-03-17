#include "spacevectors.hh"
#include <vector>

using namespace std;
using namespace aidaTT;

spacevectors::spacevectors() : UnitTest("SpaceVectors", __FILE__), _a(3.1), _b(2.7), _c(1.414)
{
    _f[0] =  .1;
    _f[1] =  .2;
    _f[2] =  .3;
    _d[0] =  1.1;
    _d[1] =  1.2;
    _d[2] =  1.3;

    // the test values for initialisation
    _vecA = new Vector3D();
    _vecB = new Vector3D(_f);
    _vecC = new Vector3D(_d);
    _vecD = new Vector3D(_a, _b, _c);
}



void spacevectors::_test()
{
    test_(floatCompare(_vecA->x() , 0.));
    test_(floatCompare(_vecA->y() , 0.));
    test_(floatCompare(_vecA->z() , 0.));

    test_(floatCompare(_vecB->x(), _f[0]));
    test_(floatCompare(_vecB->y(), _f[1]));
    test_(floatCompare(_vecB->z(), _f[2]));
    test_(floatCompare((*_vecB)[0], _f[0]));
    test_(floatCompare((*_vecB)[1], _f[1]));
    test_(floatCompare((*_vecB)[2], _f[2]));

    test_(floatCompare(_vecC->x(), _d[0]));
    test_(floatCompare(_vecC->y(), _d[1]));
    test_(floatCompare(_vecC->z(), _d[2]));

    test_(floatCompare(_vecD->x(), _a));
    test_(floatCompare(_vecD->y(), _b));
    test_(floatCompare(_vecD->z(), _c));

    test_(floatCompare(_vecA->rho() , 0.));
    test_(floatCompare(_vecB->rho() , .22360679774997896964));
    test_(floatCompare(_vecC->trans() , 1.62788205960997063873));


    (*_vecB)[0] = _a;
    (*_vecB)[1] = _b;
    test_(floatCompare(_vecB->x(), _a));
    test_(floatCompare(_vecB->y(), _b));
    test_(floatCompare(_vecB->rho() , sqrt(_a * _a + _b * _b)));
    //~ test_( floatCompare(a , b));
}



void spacevectors::run()
{
    _test();
}
