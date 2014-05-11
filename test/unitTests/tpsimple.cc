#include "tpsimple.hh"

using namespace std;
using namespace aidaTT;

tpsimple::tpsimple() : UnitTest("SimpleTrackParameters", __FILE__), a(1.2), b(2.3), c(3.4), d(4.5), e(5.6)
{
    // the test values for initialisation

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
    const double hp[hs] = {1., 0.02, 1.2, .2, .3};

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



void tpsimple::_test()
{
    // create new trackparameters
    _one = new trackParameters();

    // check whether all functions return the right values
    test_(floatCompare(_one->getQoverP() , 0.));
    test_(floatCompare(_one->getLambda() , 0.));
    test_(floatCompare(_one->getPhi() , 0.));
    test_(floatCompare(_one->getXPerp() , 0.));
    test_(floatCompare(_one->getYPerp() , 0.));

    for(size_t i = 0; i < 5; ++i)
        test_(floatCompare(_one->getTrackParameters()(i) , 0.));

    for(size_t i = 0; i < 3; ++i)
        test_(floatCompare(_one->getReferencePoint()[i], 0.));

    for(size_t i = 0; i < 5; ++i)
        for(size_t j = 0; j < 5; ++j)
            {
                test_(floatCompare(_one->getCovarianceMatrix()(i, j), 0.));
            }
    // now use the setters
    _one->setTrackParameters(*_helix, *_covmatrix,  *_refpoint);

    // check the new values
    for(size_t i = 0; i < 5; ++i)
        test_(floatCompare(_one->getTrackParameters()(i) , (*_helix)(i)));

    for(size_t i = 0; i < 3; ++i)
        test_(floatCompare(_one->getReferencePoint()[i], (*_refpoint)[i]));

    for(size_t i = 0; i < 5; ++i)
        for(size_t j = 0; j < 5; ++j)
            {
                if(i != j)
                    test_(floatCompare(_one->getCovarianceMatrix()(i, j), 0.));
                else
                    test_(floatCompare(_one->getCovarianceMatrix()(i, j), 1.));

            }

    // check whether all functions return the right values
    test_(floatCompare(_one->getQoverP(),     _one->getTrackParameters()(0)));
    test_(floatCompare(_one->getLambda(), _one->getTrackParameters()(1)));
    test_(floatCompare(_one->getPhi(),   _one->getTrackParameters()(2)));
    test_(floatCompare(_one->getXPerp(),     _one->getTrackParameters()(3)));
    test_(floatCompare(_one->getYPerp(),     _one->getTrackParameters()(4)));

}



void tpsimple::run()
{
    _test();
}
