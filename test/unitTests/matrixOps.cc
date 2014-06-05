#include "matrixOps.hh"
#include <vector>
#include <iostream>
using namespace std;
using namespace aidaTT;

matrixOps::matrixOps() : UnitTest("MatrixOperations", __FILE__)
{
    // the test values for initialisation
    double vals[] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25.};
    _values.resize(25);
    _values.assign(vals, vals + 25);
    _matrix1 = new fiveByFiveMatrix();
    _matrix2 = new fiveByFiveMatrix(_values);

    double vv[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
    _vvv.resize(5);
    _vvv.assign(vv, vv + 5);

    _vec1 = new Vector5(1.1, 2.2, 3.3, 4.4, 5.5);
    _vec2 = new Vector5();
    _vec3 = new Vector5(*_vec1);
    //_vec4 = new Vector5(vv);
    _vec5 = new Vector5(_vvv);

}



void matrixOps::_test()
{
    for(unsigned int i = 0; i < 5; ++i)
        for(unsigned int j = 0; j < 5; ++j)
            test_(floatCompare((*_matrix1)(i, j) , 0.));

    for(unsigned int row = 0; row < 5; ++row)
        {
            for(unsigned int column = 0; column < 5; ++column)
                {
                    unsigned int index = row * 5 + column;
                    (*_matrix1)(row, column) = _values.at(index);
                    test_(floatCompare((*_matrix2)(row, column) , _values.at(index)));
                }
        }

    _matrix1->Unit();
    for(unsigned int i = 0; i < 5; ++i)
        for(unsigned int j = 0; j < 5; ++j)
            {
                if(i == j)
                    test_(floatCompare((*_matrix1)(i, j) , 1.));
                else
                    test_(floatCompare((*_matrix1)(i, j) , 0.));
            }

    for(unsigned int row = 0; row < 5; ++row)
        {
            for(unsigned int column = 0; column < 5; ++column)
                {
                    unsigned int index = row * 5 + column;
                    test_(floatCompare((*_matrix2)(row, column) , _values.at(index)));
                }
        }

    fiveByFiveMatrix productMM = (*_matrix2)  * (*_matrix2);

    test_(floatCompare(productMM(0, 0) , 215));
    test_(floatCompare(productMM(0, 1) , 230));
    test_(floatCompare(productMM(0, 2) , 245));
    test_(floatCompare(productMM(0, 3) , 260));
    test_(floatCompare(productMM(0, 4) , 275));

    test_(floatCompare(productMM(1, 0) , 490));
    test_(floatCompare(productMM(1, 1) , 530));
    test_(floatCompare(productMM(1, 2) , 570));
    test_(floatCompare(productMM(1, 3) , 610));
    test_(floatCompare(productMM(1, 4) , 650));

    test_(floatCompare(productMM(2, 0) , 765));
    test_(floatCompare(productMM(2, 1) , 830));
    test_(floatCompare(productMM(2, 2) , 895));
    test_(floatCompare(productMM(2, 3) , 960));
    test_(floatCompare(productMM(2, 4) , 1025));

    test_(floatCompare(productMM(3, 0) , 1040));
    test_(floatCompare(productMM(3, 1) , 1130));
    test_(floatCompare(productMM(3, 2) , 1220));
    test_(floatCompare(productMM(3, 3) , 1310));
    test_(floatCompare(productMM(3, 4) , 1400));

    test_(floatCompare(productMM(4, 0) , 1315));
    test_(floatCompare(productMM(4, 1) , 1430));
    test_(floatCompare(productMM(4, 2) , 1545));
    test_(floatCompare(productMM(4, 3) , 1660));
    test_(floatCompare(productMM(4, 4) , 1775));

    Vector5 productMV = (* _matrix2)  * (*_vec1);

    test_(floatCompare(productMV(0) , 60.5));
    test_(floatCompare(productMV(1) , 143.));
    test_(floatCompare(productMV(2) , 225.5));
    test_(floatCompare(productMV(3) , 308.));
    test_(floatCompare(productMV(4) , 390.5));

    for(unsigned int j = 0; j < 5; ++j)
        {
            test_(floatCompare((*_vec2)(j), 0.));
            test_(floatCompare((*_vec1)(j), (*_vec3)(j)));
            test_(floatCompare((*_vec5)(j), _vvv[j]));
            // test_( floatCompare ( (*_vec4)(j), (*_vec5)(j) ) );
        }
}



void matrixOps::run()
{
    _test();
}
