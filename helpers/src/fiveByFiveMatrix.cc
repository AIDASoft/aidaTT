#include "fiveByFiveMatrix.hh"
#include <stdexcept>
#include <cmath>

#include <gsl/gsl_blas.h>

namespace aidaTT
{
    fiveByFiveMatrix::fiveByFiveMatrix()
    {
        _matrix = gsl_matrix_alloc(5, 5);
        gsl_matrix_set_zero(_matrix);
    }



    fiveByFiveMatrix::fiveByFiveMatrix(const fiveByFiveMatrix& fbfm)
    {
        _matrix = gsl_matrix_alloc(5, 5);
        gsl_matrix_memcpy(_matrix, fbfm._matrix);
    }



    fiveByFiveMatrix::fiveByFiveMatrix(const std::vector<double>& initValues)
    {
        _matrix = gsl_matrix_alloc(5, 5);
        gsl_matrix_set_zero(_matrix);

        if(initValues.size() != 25)
            throw std::length_error("Wrong length of vector for matrix initialization.");
        else
            {
                for(std::vector<double>::const_iterator elem = initValues.begin(), start = initValues.begin(), last = initValues.end(); elem < last; ++elem)
                    {
                        unsigned int index = elem - start ; // the index of the current vector element
                        unsigned int column = index % 5; // integer modulo
                        unsigned int row = index / 5; // integer division
                        gsl_matrix_set(_matrix, row, column, *elem);
                    }
            }
    }



    fiveByFiveMatrix::~fiveByFiveMatrix()
    {
        gsl_matrix_free(_matrix);
    }



    fiveByFiveMatrix& fiveByFiveMatrix::operator=(const fiveByFiveMatrix& fbfm)
    {
        if(this == &fbfm)
            return *this;

       	gsl_matrix_free(_matrix);
        _matrix = gsl_matrix_alloc(5, 5);
        gsl_matrix_memcpy(_matrix, fbfm._matrix);
        return *this;
    }



    double fiveByFiveMatrix::operator()(unsigned int row, unsigned int column) const
    {
        // check range !
        if((row > 4) || (column > 4))
            throw std::invalid_argument("Wrong index when accessing matrix elements for reading.");
        else
            return gsl_matrix_get(_matrix, row, column);
    }



    double& fiveByFiveMatrix::operator()(unsigned int row, unsigned int column)
    {
        // check range !
        if((row > 4) || (column > 4))
            throw std::invalid_argument("Wrong index when accessing matrix elements for writing.");
        else
            return *gsl_matrix_ptr(_matrix, row, column);
    }



    fiveByFiveMatrix fiveByFiveMatrix::operator*(const fiveByFiveMatrix& rm)
    {
        fiveByFiveMatrix C;
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, this->_matrix, rm._matrix, 0.0, C._matrix);

        return C;
    }



    Vector5 fiveByFiveMatrix::operator*(const Vector5& rv)
    {
        Vector5 C;
        gsl_blas_dgemv(CblasNoTrans, 1., this->_matrix, rv._vector, 1., C._vector);
        return C;
    }




    double* fiveByFiveMatrix::array() const
    {
        return _matrix->data;
    }
}
