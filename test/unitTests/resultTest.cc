#include "resultTest.hh"

using namespace std;
using namespace aidaTT;

resultTest::resultTest() : UnitTest("FitResultTest", __FILE__)
{
    _fR = new fitResults;
}



void resultTest::_test()
{
    // test_( floatCompare ( (*_vec4)(j), (*_vec5)(j) ) );
    test_(_fR->areValid() ==  false) ;
    test_(_fR->ndf() ==  0) ;
    test_(floatCompare(_fR->chiSquare(), 0.));
    test_(floatCompare(_fR->weightLost(), 0.));

    const trackParameters& tp = _fR->estimatedParameters();
    test_(floatCompare(tp.parameters()(0), 0.));
    test_(floatCompare(tp.parameters()(1), 0.));
    test_(floatCompare(tp.parameters()(2), 0.));
    test_(floatCompare(tp.parameters()(3), 0.));
    test_(floatCompare(tp.parameters()(4), 0.));



    // create some test track parameters
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
    const double hp[hs] = {1., 0.02, 1.2, .2, .3};

    vector<double> covmatrix;
    vector<double> refpoint;
    vector<double> helix;
    covmatrix.resize(ms);
    refpoint.resize(ps);
    helix.resize(hs);

    covmatrix.assign(cm, cm + ms); // vector fill
    helix.assign(hp, hp + hs); // vector fill

    fullCovariance* _covmatrix = new fullCovariance(covmatrix);
    Vector3D* _refpoint = new Vector3D(rp);
    Vector5* _helix = new Vector5(helix);

    trackParameters* TP = new trackParameters;
    TP->setTrackParameters(*_helix, *_covmatrix,  *_refpoint);

    const unsigned int n = 64;
    const double wl = 12.34;
    const double chs = 225.256;

    _fR->setResults(true, chs, n, wl, *TP);

    test_(_fR->areValid() ==  true) ;
    test_(_fR->ndf() ==  n) ;
    test_(floatCompare(_fR->chiSquare(), chs));
    test_(floatCompare(_fR->weightLost(), wl));

    const trackParameters& tp2 = _fR->estimatedParameters();
    test_(floatCompare(tp2.parameters()(0), hp[0]));
    test_(floatCompare(tp2.parameters()(1), hp[1]));
    test_(floatCompare(tp2.parameters()(2), hp[2]));
    test_(floatCompare(tp2.parameters()(3), hp[3]));
    test_(floatCompare(tp2.parameters()(4), hp[4]));

}



void resultTest::run()
{
    _test();
}
