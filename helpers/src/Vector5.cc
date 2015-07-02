#include "Vector5.hh"
#include <stdexcept>
#include <iostream>
namespace aidaTT
{
    Vector5::Vector5()
    {
        _vector = gsl_vector_alloc(5);
        gsl_vector_set_zero(_vector);
    }



    Vector5::Vector5(const Vector5& vec)
    {
        _vector = gsl_vector_alloc(5);
        gsl_vector_memcpy(_vector, vec._vector);
    }



    Vector5::Vector5(const std::vector<double>& initValues)
    {
        _vector = gsl_vector_alloc(5);
        gsl_vector_set_zero(_vector);

        if(initValues.size() != 5)
            throw std::length_error("Wrong length of vector for Vector5 initialization.");
        else
            {
                for(std::vector<double>::const_iterator elem = initValues.begin(), start = initValues.begin(), last = initValues.end(); elem < last; ++elem)
                    {
                        unsigned int index = elem - start ; // the index of the current vector element
                        gsl_vector_set(_vector, index, *elem);
                    }
            }
    }



    Vector5::Vector5(const double* array)
    {
        std::vector<double> tmpVector;
        try
            {
                tmpVector.assign(array, array + 5);
            }
        catch(...)
            {
                throw std::length_error("Wrong length of array for Vector5 initialization.");
            }

    }


    Vector5::Vector5(double i, double j, double k, double m, double n)
    {
        _vector = gsl_vector_alloc(5);
        gsl_vector_set(_vector, 0, i);
        gsl_vector_set(_vector, 1, j);
        gsl_vector_set(_vector, 2, k);
        gsl_vector_set(_vector, 3, m);
        gsl_vector_set(_vector, 4, n);
    }



    Vector5::~Vector5()
    {
        gsl_vector_free(_vector);
    }



    Vector5& Vector5::operator=(const Vector5& vec)
    {
        if(this == &vec)
            return *this;

	// gsl_vector_free(_vector);
        // _vector = gsl_vector_alloc(5);
        gsl_vector_memcpy(_vector, vec._vector);
        return (*this);
    }



    Vector5& Vector5::operator+(const Vector5& vec)
    {
        gsl_vector_add(_vector, vec._vector);

        return (*this);
    }



    double Vector5::operator()(unsigned int index) const
    {
        // check range !
        if(index > 4)
            throw std::invalid_argument("Wrong index when accessing vector elements for reading.");
        else
            return gsl_vector_get(_vector, index);
    }



    double& Vector5::operator()(unsigned int index)
    {
        // check range !
        if(index > 4)
            throw std::invalid_argument("Wrong index when accessing vector elements for writing.");
        else
            return *gsl_vector_ptr(_vector, index);
    }
}
