#ifndef TRACKPARAMETERS_H
#define TRACKPARAMETERS_H

#include <vector>
#include <bitset>
#include <stdexcept>

#include "fiveByFiveMatrix.hh"
#include "Vector3D.hh"
#include "Vector5.hh"

namespace aidaTT
{

    /*** trackParameters:: the main class for track parameter definition and interaction.
     *
     *  @version $Id
     *  @author Ch. Rosemann, DESY
     *
     * A track is fully characterized by 23 floating point numbers.
     * Five values parametrize a helix, their covariance matrix has 15 entries; w.r.t. a reference point (three values).
     *
     * There are a lot of different choices for parametrization, the choice depending on preference.
     *
     * Any other parametrization can be deduced from this by using the appropriate transformation matrix.
     * There are helper functions defined and implemented (in fiveByFiveMatrix.hh) to switch to/from/in between curvilinear and
     *      - perigee track parametrization [ kappa, theta, Phi, epsilon, z_p ], given by Billoir & Qian, NIM A449 (2000) 344
     *      - L3 perigee track parametrization [ Omega, tan(lambda), phi_0, d_0, z_0, ], given by Alcaraz, L3 Internal Note 1666 (1995)
     *
     * The default helix parametrization for propagation is the curvilinear track parametrization, drawn from Strandlie & Wittek, NIM A 566 (2006), 678ff.
     * This is: [ q/p, lambda, Phi, x_perp, y_perp ]
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

            ///~ encode the type of parametrization
            const trackParametrization& parametrization() const;

            ///~ getter functions
            Vector5 getTrackParameters() const;

            ///~ direct read access to the individual elements by index
            double operator()(unsigned int) const;

            ///~ get the reference point in global, cartesian coordinates (x,y,z)
            Vector3D getReferencePoint() const
            {
                return _refpoint;
            };

            fullCovariance getCovarianceMatrix() const
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
            void setCovarianceMatrix(const fullCovariance&);

            void print() const;

        private:
            fullCovariance _covmatrix;
            Vector5        _helixparams;
            Vector3D       _refpoint;
    };



    class trackParametrization
    {
        public:
            /// enum for the type
            enum
            {
                curvilinear = 0, perigee, ildl3
            };

            trackParametrization() : _idbits(0) {}

            /// set the trackparametrization to the given choice
            trackParametrization(unsigned int parametrization) : _idbits(0)
            {
                _idbits.set(parametrization, true) ;
            }

            bool isCurvilinear() const
            {
                return _idbits[trackParametrization::curvilinear];
            };
            bool isPerigee() const
            {
                return _idbits[trackParametrization::perigee];
            };
            bool isILD() const
            {
                return _idbits[trackParametrization::ildl3];
            };
            bool isL3() const
            {
                return _idbits[trackParametrization::ildl3];
            };

            void setParametrization(unsigned int parametrization)
            {
                _idbits.set(parametrization , true) ;
            }

        protected:

            mutable std::bitset<32> _idbits;

    };
}
#endif // TRACKPARAMETERS_H
