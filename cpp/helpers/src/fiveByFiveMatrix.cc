#include "fiveByFiveMatrix.hh"
#include <stdexcept>
#include <iostream>
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



    fiveByFiveMatrix fiveByFiveMatrix::operator=(const fiveByFiveMatrix& fbfm)
    {
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


    double* fiveByFiveMatrix::array() const
    {
        double retVal[25] =
        {
            gsl_matrix_get(_matrix, 0, 0), gsl_matrix_get(_matrix, 0, 1), gsl_matrix_get(_matrix, 0, 2), gsl_matrix_get(_matrix, 0, 3), gsl_matrix_get(_matrix, 0, 4),
            gsl_matrix_get(_matrix, 1, 0), gsl_matrix_get(_matrix, 1, 1), gsl_matrix_get(_matrix, 1, 2), gsl_matrix_get(_matrix, 1, 3), gsl_matrix_get(_matrix, 1, 4),
            gsl_matrix_get(_matrix, 2, 0), gsl_matrix_get(_matrix, 2, 1), gsl_matrix_get(_matrix, 2, 2), gsl_matrix_get(_matrix, 2, 3), gsl_matrix_get(_matrix, 2, 4),
            gsl_matrix_get(_matrix, 3, 0), gsl_matrix_get(_matrix, 3, 1), gsl_matrix_get(_matrix, 3, 2), gsl_matrix_get(_matrix, 3, 3), gsl_matrix_get(_matrix, 3, 4),
            gsl_matrix_get(_matrix, 4, 0), gsl_matrix_get(_matrix, 4, 1), gsl_matrix_get(_matrix, 4, 2), gsl_matrix_get(_matrix, 4, 3), gsl_matrix_get(_matrix, 4, 4)
        };
        return retVal;
    }

}
