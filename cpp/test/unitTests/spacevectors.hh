#ifndef SPACEVECTORS_HH
#define SPACEVECTORS_HH

/// initialization, access, write tests for Vector3D
#include "Vector3D.hh"
#include "UnitTest.hh"

class spacevectors : public UnitTesting::UnitTest
{
    public:
        spacevectors();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

        aidaTT::Vector3D* _vecA;
        aidaTT::Vector3D* _vecB;
        aidaTT::Vector3D* _vecC;
        aidaTT::Vector3D* _vecD;
        float _f[3];
        double _d[3];
        double _a, _b, _c;

};
#endif // SPACEVECTORS_HH
