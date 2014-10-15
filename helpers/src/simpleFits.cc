#include "simpleFits.hh"
#include <cmath>
#include <gsl/gsl_blas.h>

using namespace std;

namespace aidaTT
{


    simpleFitXY::simpleFitXY(bool flag, double xr, double yr) :
        _curved(flag), _npar((flag) ? 3 : 2), _xRef(xr), _yRef(yr), _numPoints(0),
        _sx(0.), _sy(0.), _sxx(0.), _sxy(0.), _syy(0.), _sw(0.), _sr(0.), _sxr(0.), _syr(0.), _srr(0.)
    {
        if(_curved)
            {
                _parameters = gsl_vector_alloc(3);
                _covariance = gsl_matrix_alloc(3, 3);
                gsl_matrix_set_zero(_covariance);
            }
        else
            {
                _parameters = gsl_vector_alloc(2);
                _covariance = gsl_matrix_alloc(2, 2);
                gsl_matrix_set_zero(_covariance);
            }

    }



/// Add point.
    /**
     * \param [in]  x       X of point
     * \param [in]  y       Y of point
     * \param [in]  w       weight of point
     */
    void simpleFitXY::addPoint(double x, double y, double w)
    {
        _numPoints++;
        double xl = x - _xRef;
        double yl = y - _yRef;
        _sw += w;
        _sx += w * xl;
        _sy += w * yl;
        _sxx += w * xl * xl;
        _sxy += w * xl * yl;
        _syy += w * yl * yl;

        if(_curved)
            {
                double r2 = xl * xl + yl * yl;
                _sr += w * r2;
                _sxr += w * r2 * xl;
                _syr += w * r2 * yl;
                _srr += w * r2 * r2;
            }
    }



/// Perform fit.
    /**
     * \param [out]  Chi2     chi2 of fit
     * \param [out]  nPoints  number of points
     * \return number of fit parameters (2 or 3)
     */
    int simpleFitXY::fit(double& Chi2, int& nPoints)
    {
        // averages
        double ax = _sx / _sw;
        double ay = _sy / _sw;
        double ar = _sr / _sw;
        double axx = _sxx / _sw;
        double ayy = _syy / _sw;
        double axy = _sxy / _sw;
        double axr = _sxr / _sw;
        double ayr = _syr / _sw;
        double arr = _srr / _sw;
        // variances
        double cxx = axx - ax * ax;
        double cyy = ayy - ay * ay;
        double cxy = axy - ax * ay;
        double cxr = axr - ax * ar;
        double cyr = ayr - ay * ar;
        double crr = arr - ar * ar;

        double q1, q2;
        if(_curved)
            {
                q1 = crr * cxy - cxr * cyr;
                q2 = crr * (cxx - cyy) - cxr * cxr + cyr * cyr;
            }
        else
            {
                q1 = cxy;
                q2 = cxx - cyy;
            }
        double phi = 0.5 * atan2(2. * q1, q2);
        double sinphi = sin(phi);
        double cosphi = cos(phi);

        // compare phi with initial track direction
        if(cosphi * (ax + _xRef) + sinphi * (ay + _yRef) < 0.)
            {
                // reverse direction
                phi -= (phi > 0.) ? M_PI : -M_PI;
                cosphi = -cosphi;
                sinphi = -sinphi;
            }

        if(_curved)
            {
                double kappa = (sinphi * cxr - cosphi * cyr) / crr;
                double delta = -kappa * ar + sinphi * ax - cosphi * ay;
                // track parameters
                double rho = -2. * kappa / sqrt(1. - 4. * delta * kappa);
                double d = 2. * delta / (1. + sqrt(1. - 4. * delta * kappa));
                gsl_vector_set(_parameters, 0, rho);
                gsl_vector_set(_parameters, 1, phi);
                gsl_vector_set(_parameters, 2, d);

                // chi2
                double u = 1. - rho * d;

                Chi2 = _sw * u * u * (sinphi * sinphi * cxx - 2. * sinphi * cosphi * cxy + cosphi * cosphi * cyy - kappa * kappa * crr);

                // calculate covariance matrix
                double sa = sinphi * _sx - cosphi * _sy;
                double sb = cosphi * _sx + sinphi * _sy;
                double sc = (sinphi * sinphi - cosphi * cosphi) * _sxy + sinphi * cosphi * (_sxx - _syy);
                double sd = sinphi * _sxr - cosphi * _syr;
                double saa = sinphi * sinphi * _sxx - 2. * sinphi * cosphi * _sxy + cosphi * cosphi * _syy;

                gsl_matrix_set(_covariance, 0, 0,  0.25 * _srr - d * (sd - d * (saa + 0.5 * _sr - d * (sa - 0.25 * d * _sw))));
                gsl_matrix_set(_covariance, 0, 1,  u * (0.5 * (cosphi * _sxr + sinphi * _syr) - d * (sc - 0.5 * d * sb)));
                gsl_matrix_set(_covariance, 1, 0,  gsl_matrix_get(_covariance, 0, 1));
                gsl_matrix_set(_covariance, 1, 1,  u * u * (cosphi * cosphi * _sxx + 2. * cosphi * sinphi * _sxy + sinphi * sinphi * _syy));
                gsl_matrix_set(_covariance, 0, 2,  rho * (-0.5 * sd + d * saa) - 0.5 * u * _sr + 0.5 * d * ((3 * u - 1.) * sa - u * d * _sw));
                gsl_matrix_set(_covariance, 2, 0,  gsl_matrix_get(_covariance, 0, 2));
                gsl_matrix_set(_covariance, 1, 2,  -u * (rho * sc + u * sb));
                gsl_matrix_set(_covariance, 2, 1,  gsl_matrix_get(_covariance, 1, 2));
                gsl_matrix_set(_covariance, 2, 2,  rho * (rho * saa + 2 * u * sa) + u * u * _sw);
            }
        else
            {
                // track parameters
                double d = sinphi * ax - cosphi * ay;
                gsl_vector_set(_parameters, 0, phi);
                gsl_vector_set(_parameters, 1, d);

                // chi2
                Chi2 = _sw * (sinphi * sinphi * cxx - 2. * sinphi * cosphi * cxy + cosphi * cosphi * cyy);

                // calculate covariance matrix
                gsl_matrix_set(_covariance, 0, 0, cosphi * cosphi * _sxx + 2. * cosphi * sinphi * _sxy + sinphi * sinphi * _syy);
                gsl_matrix_set(_covariance, 0, 1, -(cosphi * _sx + sinphi * _sy));
                gsl_matrix_set(_covariance, 1, 0, gsl_matrix_get(_covariance, 0, 1));
                gsl_matrix_set(_covariance, 1, 1, _sw);
            }

/// TODO: invert covariance    _covariance.Invert();

        nPoints = _numPoints;
        return _npar;
    }

/// Get parameters vector.
    /**
     * \return parameter vector
     */
    gsl_vector* simpleFitXY::getPar() const
    {
        return _parameters;
    }

/// Get covariance matrix.
    /**
     * \return covariance matrix
     */
    gsl_matrix* simpleFitXY::getCov() const
    {
        return _covariance;
    }



/// Constructor for simple fit in ZS.
    simpleFitZS::simpleFitZS() :
        _npar(2), _numPoints(0), _sx(0.), _sy(0.), _sxx(0.), _sxy(0.), _syy(0.), _sw(0.)
    {
        _parameters = gsl_vector_alloc(2);
        _covariance = gsl_matrix_alloc(2, 2);
    }

/// Add point.
    /**
     * \param [in]  x       arc-length S of point
     * \param [in]  y       Z of point
     * \param [in]  w       weight of point
     */
    void simpleFitZS::addPoint(double x, double y, double w)
    {
        _numPoints++;

        // x is S, y is Z
        _sw += w;
        _sx += w * x;
        _sy += w * y;
        _sxx += w * x * x;
        _sxy += w * x * y;
        _syy += w * y * y;
    }

/// Perform fit.
    /**
     * \param [out]  Chi2     chi2 of fit
     * \param [out]  nPoints  number of points
     * \return number of fit parameters (2)
     */
    int simpleFitZS::fit(double& Chi2, int& nPoints)
    {
        // linear equation system A*x = b
        gsl_vector* bVec = gsl_vector_alloc(2);
        gsl_vector_set(bVec, 0, _sxy);
        gsl_vector_set(bVec, 1, _sy);

        gsl_matrix_set(_covariance, 0, 0, _sxx);
        gsl_matrix_set(_covariance, 0, 1, _sx);
        gsl_matrix_set(_covariance, 1, 0, _sx);
        gsl_matrix_set(_covariance, 1, 1, _sw);

        // solve
//!!!   aMat.Invert();

        gsl_blas_dgemv(CblasNoTrans, 1., _covariance, bVec, 1., _parameters);

        // chi2
        Chi2 =  pow(gsl_vector_get(_parameters, 0), 2) * _sxx +
                pow(gsl_vector_get(_parameters, 1), 2) * _sw  + _syy
                + 2. * gsl_vector_get(_parameters, 0) * gsl_vector_get(_parameters, 1) * _sx
                - 2. * gsl_vector_get(_parameters, 0) * _sxy
                - 2. * gsl_vector_get(_parameters, 1) * _sy;

        nPoints = _numPoints;
        return _npar;
    }

/// Get parameters vector.
    /**
     * \return parameter vector
     */
    gsl_vector* simpleFitZS::getPar() const
    {
        return _parameters;
    }

/// Get covariance matrix.
    /**
     * \return covariance matrix
     */
    gsl_matrix* simpleFitZS::getCov() const
    {
        return _covariance;
    }


} // close namespace aidaTT
