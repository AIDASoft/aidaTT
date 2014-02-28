#include "tpsimple.hh"

using namespace std;
using namespace aidaTT;

tpsimple::tpsimple() : UnitTest("SimpleTrackParameters", __FILE__), a(1.2), b(2.3), c(3.4), d(4.5), e(5.6)
{
    // the test values for initialisation

    // first: fix sizes of parameters
    typedef const unsigned int cui;
    cui ms = 15, ps = 3, hs = 5;

    // define values
    const double cm[ms] = {1., 0., 1., 0., 0., 1., 0. , 0. , 0., 1., 0., 0., 0., 0., 1. };
    const double rp[ps] = {2., 3., 5.};
    const double hp[hs] = {1., 0.02, 1.2, .2, .3};

    _covmatrix.resize(15);
    _refpoint.resize(3) ;
    _helix.resize(5);
    _covmatrix.assign(cm, cm + ms); // vector fill
    _refpoint.assign(rp, rp + ps); // vector fill
    _helix.assign(hp, hp + hs); // vector fill
}



void tpsimple::_test()
{
    // create new trackparameters
    _one = new trackParameters();

    // check whether all functions return the right values
    test_(floatCompare(_one->getOmega() , 0.));
    test_(floatCompare(_one->getTanLambda() , 0.));
    test_(floatCompare(_one->getPhiZero() , 0.));
    test_(floatCompare(_one->getDZero() , 0.));
    test_(floatCompare(_one->getZZero() , 0.));

    for(size_t i = 0; i < 5; ++i)
        test_(floatCompare(_one->getTrackParameters().at(i) , 0.));

    for(size_t i = 0; i < 3; ++i)
        test_(floatCompare(_one->getReferencePoint().at(i), 0.));

    for(size_t i = 0; i < 15; ++i)
        test_(floatCompare(_one->getCovarianceMatrix().at(i), 0.));

    // now use the setters
    _one->setTrackParameters(_helix, _refpoint, _covmatrix);

    // check the new values
    for(size_t i = 0; i < 5; ++i)
        test_(floatCompare(_one->getTrackParameters().at(i) , _helix.at(i)));

    for(size_t i = 0; i < 3; ++i)
        test_(floatCompare(_one->getReferencePoint().at(i), _refpoint.at(i)));

    for(size_t i = 0; i < 15; ++i)
        test_(floatCompare(_one->getCovarianceMatrix().at(i), _covmatrix.at(i)));

    // check whether all functions return the right values
    test_(floatCompare(_one->getOmega(),     _one->getTrackParameters().at(0)));
    test_(floatCompare(_one->getTanLambda(), _one->getTrackParameters().at(1)));
    test_(floatCompare(_one->getPhiZero(),   _one->getTrackParameters().at(2)));
    test_(floatCompare(_one->getDZero(),     _one->getTrackParameters().at(3)));
    test_(floatCompare(_one->getZZero(),     _one->getTrackParameters().at(4)));

}



void tpsimple::run()
{
    _test();
}
