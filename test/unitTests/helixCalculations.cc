#include "helixCalculations.hh"
#include "helixUtils.hh"

#include <vector>

using namespace std;
using namespace aidaTT;

helixCalculations::helixCalculations() : UnitTest("HelixCalculations", __FILE__)
{

    // first: fix sizes of parameters
    typedef const unsigned int cui;
    cui ms = 25, ps = 3, hs = 5;

    // define values
    const double cm[ms] = { 1., 0. , 0., 0., 0.,
                            0., 1. , 0., 0., 0.,
                            0., 0. , 1., 0., 0.,
                            0., 0. , 0., 1., 0.,
                            0., 0. , 0., 0., 1.
                          };

    const double rp[ps] = {2., 3., 5.};
    // helix parameters: omega = -.1 , tanLambda = 1., phi0 = M_PI4, d0 = 0., z0 = 0.
    const double hp[hs] = {.1, 1., M_PI_4, 0., 0.};

    vector<double> covmatrix;
    vector<double> refpoint;
    vector<double> helix;
    covmatrix.resize(ms);
    refpoint.resize(ps);
    helix.resize(hs);

    covmatrix.assign(cm, cm + ms); // vector fill
    helix.assign(hp, hp + hs); // vector fill
    _covmatrix = new fullCovariance(covmatrix);
    _refpoint = new Vector3D(rp);
    _helix = new Vector5(helix);
}



void helixCalculations::_test()
{
    _one = new trackParameters();
    _two = new trackParameters();
    _two->setTrackParameters(*_helix);

    _three = new trackParameters();


// zero testing
    test_(floatCompare(calculateRadius(*_one), 0.));
    test_(floatCompare(calculateXCenter(*_one), 0.));
    test_(floatCompare(calculateYCenter(*_one), 0.));
    test_(floatCompare(calculatePhifromXY(0., 0., *_one), 0.));
    test_(floatCompare(calculateSfromXY(0., 0., *_one), 0.));
    test_(floatCompare(calculateXfromS(0., *_one), 0.));
    test_(floatCompare(calculateYfromS(0., *_one), 0.));
    test_(floatCompare(calculateZfromS(0., *_one), 0.));

    // test a special track
    test_(floatCompare(calculateRadius(*_two), 10.));
    test_(floatCompare(calculateXCenter(*_two), 7.07106781186547524401));
    test_(floatCompare(calculateYCenter(*_two), -7.07106781186547524401));
    test_(floatCompare(calculatePhifromXY(0., 0., *_two), M_PI_4));

    test_(floatCompare(calculateSfromXY(7.07106781186547524401, 2.92893218813452475599, *_two), 10. * M_PI_4));
    test_(floatCompare(calculateXfromS(10. * M_PI_4, *_two), 7.07106781186547524401));
    test_(floatCompare(calculateYfromS(10. * M_PI_4, *_two),  2.92893218813452475599));
    test_(floatCompare(calculateZfromS(1.234567 / 10., *_two), 1.234567 / 10.));


}



void helixCalculations::run()
{
    _test();
}
