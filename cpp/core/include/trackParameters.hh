#ifndef TRACKPARAMETERS_H
#define TRACKPARAMETERS_H

#include <vector>

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
     * The default interface uses (for now) the L3/LCIO perigee parametrization: [\Omega, \tan\lambda, \Phi_0, d_0, z_0]
     *
     * Any other parametrization can be deduced from this by using the appropriate transformation matrix.
     *
     * The default internal helix parametrization for propagation is the curvilinear track parametrization, drawn from Strandlie & Wittek, NIM A 566 (2006), 678ff.
     * This is: [ q/p, lambda, Phi, x_perp, y_perp ]
     *
     * Implemented are:
     *      - perigee track parametrization [ kappa, theta, Phi, epsilon, z_p ], given by Billoir & Qian, NIM A449 (2000) 344
     *      - L3 perigee track parametrization [ Omega, tan(lambda), phi_0, d_0m z_0, ], given by Alcaraz, L3 Internal Note 1666 (1995)
     ***/

    /* TODO: decide the nicenesss of using. internal switching?
     *  use only strandlie&wittek parametrization !?
     *  implement the automatic conversion with matrices!
     *  implement copy construction & assignment
     *
     */

    class trackParameters
    {
        public:
            /// default constructor, sets everything to zero
            trackParameters();

            // define copy ctor and assignment
            trackParameters(const trackParameters&);
            trackParameters operator=(const trackParameters&);

            ///~ getter functions
            std::vector<double> getTrackParameters() const;

            ///~ the individual values
            double getOmega()     const
            {
                return _omega;
            };
            double getTanLambda() const
            {
                return _tanlambda;
            };
            double getPhiZero()   const
            {
                return _phizero;
            };
            double getDZero()     const
            {
                return _dzero;
            };
            double getZZero()     const
            {
                return _zzero;
            };

            ///~ get the reference point in global, cartesian coordinates (x,y,z)
            std::vector<double> getReferencePoint() const
            {
                return _refpoint;
            };

            ///~ get the covariance matrix in lower triangle form
            ///~ index-like structure:
            ///~ [  0   1   3   6  10 ]
            ///~ [  1   2   4   7  11 ]
            ///~ [  3   4   5   8  12 ]
            ///~ [  6   7   8   9  13 ]
            ///~ [ 10  11  12  13  14 ]

            std::vector<double> getCovarianceMatrix() const
            {
                return _covmatrix;
            } ;

            ///~ setter functions
            ///~ set everything at once
            void setTrackParameters(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);

            ///~ only set the helix parameters
            void setTrackParameters(const std::vector<double>&);
            void setOmega(const double& x)
            {
                _omega = x;
                _helixparams.at(0) = _omega;
            };
            void setTanLambda(const double& x)
            {
                _tanlambda = x;
                _helixparams.at(1) = _tanlambda;
            };
            void setPhiZero(const double& x)
            {
                _phizero = x;
                _helixparams.at(2) = _phizero;
            };
            void setDZero(const double& x)
            {
                _dzero = x;
                _helixparams.at(3) = _dzero;
            };
            void setZZero(const double& x)
            {
                _zzero = x;
                _helixparams.at(4) = _zzero;
            };

            /// set the reference point, defaults to the nominal center
            void setReferencePoint(const std::vector<double>&);

            ///~ set the covariance matrix
            void setCovarianceMatrix(const std::vector<double>& cm)
            {
                _covmatrix = cm;
            };

            void print() const;

        private:
            double _omega;
            double _tanlambda;
            double _phizero;
            double _dzero;
            double _zzero;

            double _xref;
            double _yref;
            double _zref;

            std::vector<double> _covmatrix;
            std::vector<double> _helixparams;
            std::vector<double> _refpoint;
    };

}
#endif // TRACKPARAMETERS_H
