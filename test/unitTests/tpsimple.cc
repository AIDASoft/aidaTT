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
    _two = new trackParameters();
    _three = new trackParameters();
    _four = new trackParameters();

    // check whether all functions return the right values
    test_(floatCompare((*_one)(0) , 0.));
    test_(floatCompare((*_one)(1) , 0.));
    test_(floatCompare((*_one)(2) , 0.));
    test_(floatCompare((*_one)(3) , 0.));
    test_(floatCompare((*_one)(4) , 0.));


    _two->setReferencePoint(*_refpoint);
    _three->setTrackParameters(*_helix);
    _four->setCovarianceMatrix(*_covmatrix);

    for(size_t i = 0; i < 5; ++i)
        {
            test_(floatCompare(_one->parameters()(i) , 0.));
            test_(floatCompare(_two->parameters()(i) , 0.));
            test_(floatCompare(_four->parameters()(i) , 0.));
        }
    for(size_t i = 0; i < 3; ++i)
        {
            test_(floatCompare(_one->referencePoint()[i], 0.));
            test_(floatCompare(_three->referencePoint()[i], 0.));
            test_(floatCompare(_four->referencePoint()[i], 0.));
        }

    for(size_t i = 0; i < 5; ++i)
        for(size_t j = 0; j < 5; ++j)
            {
                test_(floatCompare(_one->covarianceMatrix()(i, j), 0.));
                test_(floatCompare(_two->covarianceMatrix()(i, j), 0.));
                test_(floatCompare(_three->covarianceMatrix()(i, j), 0.));
            }
    // now use the setters
    _one->setTrackParameters(*_helix, *_covmatrix,  *_refpoint);

    // check the new values
    for(size_t i = 0; i < 5; ++i)
        {
            test_(floatCompare(_one->parameters()(i) , (*_helix)(i)));
            test_(floatCompare(_three->parameters()(i) , (*_helix)(i)));
        }
    for(size_t i = 0; i < 3; ++i)
        {
            test_(floatCompare(_one->referencePoint()[i], (*_refpoint)[i]));
            test_(floatCompare(_two->referencePoint()[i], (*_refpoint)[i]));
        }

    for(size_t i = 0; i < 5; ++i)
        for(size_t j = 0; j < 5; ++j)
            {
                if(i != j)
                    {
                        test_(floatCompare(_one->covarianceMatrix()(i, j), 0.));
                        test_(floatCompare(_four->covarianceMatrix()(i, j), 0.));
                    }
                else
                    {
                        test_(floatCompare(_one->covarianceMatrix()(i, j), 1.));
                        test_(floatCompare(_four->covarianceMatrix()(i, j), 1.));
                    }
            }

    // check whether all functions return the right values
    test_(floatCompare((*_one)(0),     _one->parameters()(0)));
    test_(floatCompare((*_one)(1), _one->parameters()(1)));
    test_(floatCompare((*_one)(2),   _one->parameters()(2)));
    test_(floatCompare((*_one)(3),     _one->parameters()(3)));
    test_(floatCompare((*_one)(4),     _one->parameters()(4)));

}



void tpsimple::run()
{
    _test();
}
