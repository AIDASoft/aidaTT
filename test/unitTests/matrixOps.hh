#ifndef MATRIXOPS_HH
#define MATRIXOPS_HH

/// initialization, access, write tests for matrix and vector abstraction class
#include "fiveByFiveMatrix.hh"
#include "Vector5.hh"

#include "UnitTest.hh"
#include <vector>

class matrixOps : public UnitTesting::UnitTest
{
    public:
        matrixOps();
        void run();

    private:
        // the test calls in different blocks
        // the distinctions are arbitrary:
        void _test();

        aidaTT::fiveByFiveMatrix* _matrix1;
        aidaTT::fiveByFiveMatrix* _matrix2;
        aidaTT::fiveByFiveMatrix* _matrix3;
        aidaTT::Vector5* _vec1;
        aidaTT::Vector5* _vec2;
        aidaTT::Vector5* _vec3;
        aidaTT::Vector5* _vec4;
        aidaTT::Vector5* _vec5;
        std::vector<double> _values;
        std::vector<double> _vvv;

};
#endif // MATRIXOPS_HH
