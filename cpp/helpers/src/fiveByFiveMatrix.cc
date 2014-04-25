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



    fiveByFiveMatrix fiveByFiveMatrix::operator=(const fiveByFiveMatrix& fbfm)
    {
        if(this == &fbfm)
            return *this;

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




    const double* const fiveByFiveMatrix::array() const
    {
        return _matrix->data;
    }



    /// calculate and return transformation matrix from curivilinear to perigee track parameters (at reference point)
    fiveByFiveMatrix curvilinearToPerigeeJacobian(const Vector5& clParams, const Vector3D& bfield)
    {
        const double qop    = clParams(0);
        const double lambda = clParams(1);
        const double phi0   = clParams(2);

        // define local curvilinear coordinate system: U = Z x T / |Z x T|, V = T x U
        const double cosPhi = cos(phi0);
        const double sinPhi = sin(phi0);
        const double tanLambda = tan(lambda);
        const double cosLambda = 1. / sqrt(1. + tanLambda * tanLambda);
        const double sinLambda = sin(lambda);

        Vector3D T(cosPhi * cosLambda, sinPhi * cosLambda, sinLambda);
        Vector3D U(-sinPhi, cosPhi, 0.);
        Vector3D V(-cosPhi * sinLambda, -sinPhi * sinLambda, cosLambda);

        // helper vectors for local (perigee) system
        // J = -U , K = Z, I = J x K
        Vector3D J(sinPhi, -cosPhi, 0.);
        Vector3D K(0., 0., 1.);
        Vector3D I = J.cross(K);
        // set the right variables for helix parametrisation by wittek/strandlie
        // H is direction of magnetic field, N is H x T/a, a = |H x T|
        Vector3D H(bfield.unit()); // magnetic field direction
        Vector3D aN = H.cross(T); // a*N

        // formal helpers for later evaluation of jacobian transformation from curvilinear track parameters to perigee frame
        const double ui = U.dot(I), anv = aN.dot(V), ti = T.dot(I);
        const double vi = V.dot(I), anu = aN.dot(U), vk = V.dot(K);

        fiveByFiveMatrix jacobian;
        jacobian.Unit();

        // differences to unit matrix:
        const double Q = - qop * bfield.r(); // -Bz*q/p
        const double qbar = qop * bfield.z();
        jacobian(0, 0) = - bfield.z() / cosLambda;
        jacobian(0, 1) = -qbar * tanLambda / cosLambda;
        jacobian(0, 3) = qbar * Q * tanLambda * ui * anv / cosLambda / ti;
        jacobian(0, 4) = qbar * Q * tanLambda * vi * anv / cosLambda / ti;

        jacobian(1, 1) = -1.;
        jacobian(1, 3) = Q * ui * anv / ti;
        jacobian(1, 4) = Q * vi * anv / ti;

        jacobian(2, 3) = -Q * ui * anu / cosLambda / ti;
        jacobian(2, 4) = -Q * vi * anu / cosLambda / ti;

        jacobian(3, 3) = vk / ti;

        jacobian(4, 4) = -1. / ti;

        return jacobian;
    }



    fiveByFiveMatrix perigeeToCurvilinearJacobian(const Vector5& perParams, const Vector3D& bfield)
    {
        const double kappa   = perParams(0);
        const double theta   = perParams(1);
        const double phi     = perParams(2);

        // define local curvilinear coordinate system: U = Z x T / |Z x T|, V = T x U
        const double qop = - kappa * sin(theta) / bfield.z();
        const double lambda = M_PI_2 - theta;
        const double cosPhi = cos(phi);
        const double sinPhi = sin(phi);
        const double tanLambda = tan(lambda);
        const double cosLambda = 1. / sqrt(1. + tanLambda * tanLambda);
        const double sinLambda = sin(lambda);

        Vector3D T(cosPhi * cosLambda, sinPhi * cosLambda, sinLambda);
        Vector3D U(-sinPhi, cosPhi, 0.);
        Vector3D V(-cosPhi * sinLambda, -sinPhi * sinLambda, cosLambda);

        // helper vectors for local (perigee) system
        // J = -U , K = Z, I = J x K
        Vector3D J(sinPhi, -cosPhi, 0.);
        Vector3D K(0., 0., 1.);
        // set the right variables for helix parametrisation by wittek/strandlie
        // H is direction of magnetic field, N is H x T/a, a = |H x T|
        Vector3D H(bfield.unit()); // magnetic field direction
        Vector3D aN = H.cross(T); // a*N

        // formal helpers for later evaluation of jacobian transformation from curvilinear track parameters to perigee frame
        const double anv = aN.dot(V), tj = T.dot(J), tk = T.dot(K);
        const double anu = aN.dot(U), vk = V.dot(K);

        fiveByFiveMatrix jacobian;
        jacobian.Unit();

        // differences to unit matrix:
        const double Q = - qop * bfield.r(); // -Bz*q/p
        jacobian(0, 0) = - sin(theta) / bfield.z() ;
        jacobian(0, 1) = qop / tanLambda ;

        jacobian(1, 1) = -1.;
        jacobian(1, 3) = - Q * tj * anv;
        jacobian(1, 4) = - Q * tk * anv;

        jacobian(2, 3) = -Q * tj * anu / cosLambda;
        jacobian(2, 4) = -Q * tk * anu / cosLambda;

        jacobian(3, 3) = -1 ;

        jacobian(4, 4) = vk;

        return jacobian;
    }



    fiveByFiveMatrix perigeeToILDJacobian(const Vector5& perParams)
    {
        const double tanLambda = tan(M_PI_2 - perParams(1));

        fiveByFiveMatrix per2ILDjacobian;
        per2ILDjacobian.Unit();

        // differences to unit matrix:
        per2ILDjacobian(0, 0) = -1.;
        per2ILDjacobian(1, 1) = -(1. + tanLambda * tanLambda);
        per2ILDjacobian(3, 3) = -1.;

        return per2ILDjacobian;
    }



    fiveByFiveMatrix ildToPerigeeJacobian(const Vector5& ildParams)
    {
        const double tanLambda = ildParams(1);

        fiveByFiveMatrix ild2PERjacobian;
        ild2PERjacobian.Unit();

        // differences to unit matrix:
        ild2PERjacobian(0, 0) = -1.;
        ild2PERjacobian(1, 1) = -1 / (1. + tanLambda * tanLambda);
        ild2PERjacobian(3, 3) = -1.;

        return ild2PERjacobian;
    }



    fiveByFiveMatrix curvilinearToILDJacobian(const Vector5& curvilinearJacobian)
    {
        return fiveByFiveMatrix();
    }



    fiveByFiveMatrix ildToCurvilinearJacobian(const Vector5& ILDJacobian)
    {
        return fiveByFiveMatrix();
    }
}
