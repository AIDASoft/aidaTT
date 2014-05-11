#ifndef CONSTANTBFIELD_HH
#define CONSTANTBFIELD_HH

#include "ConstantSolenoidBField.hh"

#include "UnitTest.hh"
#include "Vector3D.hh"

class constantbfield : public UnitTesting::UnitTest
{
    public:
        constantbfield();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

        aidaTT::ConstantSolenoidBField* _bf1;
        aidaTT::ConstantSolenoidBField* _bf2;
        typedef const double cdbl;
        cdbl a, b;
        aidaTT::Vector3D _point ;
        aidaTT::Vector3D _val1 ;
        aidaTT::Vector3D _val2 ;
};
#endif // CONSTANTBFIELD_HH
