#include "constantbfield.hh"

using namespace aidaTT;

constantbfield::constantbfield() : UnitTest("ConstantBField", __FILE__), a(-1.2), b(2.3)
{
    _point[0] = 0.;
    _point[1] = 3.;
    _point[2] = 2.;

    _val1[0] = 0.;
    _val1[1] = 0.;
    _val1[2] = a;

    _val2[0] = 0.;
    _val2[1] = 0.;
    _val2[2] = b;
}



void constantbfield::_test()
{
    // create new trackparameters
    _bf1 = new ConstantSolenoidBField(a);
    _bf2 = new ConstantSolenoidBField(b);

    // check whether all functions return the right values
    test_(floatCompare(_bf1->Bx(_point) , 0.));
    test_(floatCompare(_bf1->By(_point) , 0.));
    test_(floatCompare(_bf1->Bz(_point) , a));

    test_(floatCompare(_bf2->Bx(_point) , 0.));
    test_(floatCompare(_bf2->By(_point) , 0.));
    test_(floatCompare(_bf2->Bz(_point) , b));

    test_(floatCompare(_bf1->BField(_point).x() , 0.));
    test_(floatCompare(_bf1->BField(_point).y() , 0.));
    test_(floatCompare(_bf1->BField(_point).z() , a));

    test_(floatCompare(_bf2->BField(_point).x() , 0.));
    test_(floatCompare(_bf2->BField(_point).y() , 0.));
    test_(floatCompare(_bf2->BField(_point).z() , b));

    test_(floatCompare(_bf1->BField(2., 3.14, -23.4).x() , 0.));
    test_(floatCompare(_bf1->BField(2., 3.14, -23.4).y() , 0.));
    test_(floatCompare(_bf1->BField(2., 3.14, -23.4).z() , a));

}



void constantbfield::run()
{
    _test();
}
