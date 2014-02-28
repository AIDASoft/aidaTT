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
    _values.assign(vals, vals + 25 );
    _matrix1 = new fiveByFiveMatrix();
    _matrix2 = new fiveByFiveMatrix(_values);

}



void matrixOps::_test()
{
    for (unsigned int i = 0; i < 5; ++i)
        for(unsigned int j = 0; j < 5; ++j)
            test_( floatCompare( (*_matrix1)(i,j) ,0. ));

    for(unsigned int row = 0; row < 5; ++row)
    {
        for(unsigned int column = 0; column < 5; ++column)
            {
                unsigned int index = row * 5 + column;
                (*_matrix1)(row,column) = _values.at(index);
                test_( floatCompare( (*_matrix2)(row,column) ,_values.at(index) ));
            }
    }
    
    _matrix1->Unit();
    for (unsigned int i = 0; i < 5; ++i)
        for(unsigned int j = 0; j < 5; ++j)
            {
                if( i == j)
                    test_( floatCompare( (*_matrix1)(i,j) ,1. ));
                else
                    test_( floatCompare( (*_matrix1)(i,j) ,0. ));
            }
    
    for(unsigned int row = 0; row < 5; ++row)
    {
        for(unsigned int column = 0; column < 5; ++column)
            {
                unsigned int index = row * 5 + column;
                test_( floatCompare( (*_matrix2)(row,column) ,_values.at(index) ));
            }
    }
    
    
}



void matrixOps::run()
{
    _test();
}
