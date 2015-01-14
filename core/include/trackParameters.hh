#ifndef TRACKPARAMETERS_H
#define TRACKPARAMETERS_H

#include <vector>
#include <bitset>
#include <stdexcept>
#include <ostream>

#include "fiveByFiveMatrix.hh"
#include "Vector3D.hh"
#include "Vector5.hh"

namespace aidaTT
{
    /// trackParameters:: the main class for track parameter definition and interaction
    /***
     *  @version $Id
     *  @author Ch. Rosemann, DESY
     *
     * A track is fully characterized by 23 floating point numbers.
     * Five values parametrize a helix, their covariance matrix has 15 entries; w.r.t. a reference point (three values).
     *
     * The choice for interfacing aidaTT to the outer world is the L3 parametrization, a particular perigee representation:
     *  [ Omega, tan(lambda), phi_0, d_0, z_0, ], e.g. explained in Alcaraz, L3 Internal Note 1666 (1995)
     * A similar choice, but in a different order is given in a publicly accessible note:
     *  Kraemer, Track parameters in LCIO [http://www-flc.desy.de/lcnotes/notes/LC-DET-2006-004.pdf]
     *
     ***/

    /* TODO:
     *  implement copy construction & assignment
     */

    class trackParametrization;

    class trackParameters
    {
        public:
            /// default constructor, sets everything to zero
            trackParameters();

            // define copy ctor and assignment
            trackParameters(const trackParameters&);
            trackParameters& operator=(const trackParameters&);

            ///~ getter functions
            const Vector5& parameters() const
            {
                return _helixparams;
            }

            /// access parameters to set values
            Vector5& parameters()
            {
                return _helixparams;
            }

            ///~ direct read access to the individual elements by index
            double operator()(unsigned int) const;

            ///~ get the reference point in global, cartesian coordinates (x,y,z)
            const Vector3D& referencePoint() const
            {
                return _refpoint;
            }

            /// access the reference point in global, cartesian coordinates (x,y,z)
            Vector3D& referencePoint()
            {
                return _refpoint;
            }

            /// read access covariance matrix
            const fullCovariance& covarianceMatrix() const
            {
                return _covmatrix;
            };

            /// write access covariance matrix
            fullCovariance& covarianceMatrix()
            {
                return _covmatrix;
            };

            ///~ setter functions
            ///~ set everything at once
            void setTrackParameters(const Vector5&, const fullCovariance&, const Vector3D&);

            ///~ only set the helix parameters
            void setTrackParameters(const Vector5&);

            ///~ direct write access to the individual track parameters by index
            double& operator()(unsigned int);

            /// set the reference point, defaults to the nominal center
            void setReferencePoint(const Vector3D&);

            ///~ set the covariance matrix
            void setCovarianceMatrix(const fiveByFiveMatrix&);

        private:
            fiveByFiveMatrix _covmatrix;
            Vector5        _helixparams;
            Vector3D       _refpoint;
    };



    inline std::ostream & operator << (std::ostream & os, const trackParameters& tp)
    {
        os << " [parameters: " << tp.parameters() << "]" << std::endl << " { covariance: " << tp.covarianceMatrix() << " } , ( reference point: " << tp.referencePoint() << " )";
        return os ;
    }
}
#endif // TRACKPARAMETERS_H
