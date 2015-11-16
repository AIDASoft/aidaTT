#ifndef SIMPLEFITS_HH
#define SIMPLEFITS_HH

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/// simple fits in xy and sz projections

namespace aidaTT
{

/// Simple fit in XY.
    /**
     * Fit circle ( V.Karimaeki: "Effective circle fitting for particle trajectories" ) or straight line
     */
    class simpleFitXY
    {
        public:
            simpleFitXY(bool, double, double);
            void addPoint(double, double, double);
            int fit(double&, int&);
            gsl_vector* getPar() const;
            gsl_matrix* getCov() const;

        private:
            /// flag for curved (circle) track
            const bool _curved;
            /// number of track parameters (3: circle, 2: line)
            const int _npar;
            /// X of reference point
            const double _xRef;
            /// Y of reference point
            const double _yRef;
            /// number of points (hits) used
            int _numPoints;
            /// weighted sum(x)
            double _sx;
            /// weighted sum(y)
            double _sy;
            /// weighted sum(x*x)
            double _sxx;
            /// weighted sum(x*y)
            double _sxy;
            /// weighted sum(y*y)
            double _syy;
            /// sum of weights
            double _sw;
            /// weighted sum(r*r)
            double _sr;
            /// weighted sum(x*r*r)
            double _sxr;
            /// weighted sum(y*r*r)
            double _syr;
            /// weighted sum(r*r*r*r)
            double _srr;
            /// parameter vector
            gsl_vector* _parameters;
            /// covariance matrix
            gsl_matrix* _covariance;
    };

/// Simple fit in ZS.
    /**
     * Fit straight line.
     */
    class simpleFitZS
    {
        public:
            simpleFitZS();
            void addPoint(double, double, double);
            int fit(double&, int&);
            gsl_vector* getPar() const;
            gsl_matrix* getCov() const;

        private:
            /// number of track parameters (2)
            const int _npar;
            /// number of points (hits) used
            int _numPoints;
            /// weighted sum(s)
            double _sx;
            /// weighted sum(z)
            double _sy;
            /// weighted sum(s*s)
            double _sxx;
            /// weighted sum(s*z)
            double _sxy;
            /// weighted sum(z*z)
            double _syy;
            /// sum of weights
            double _sw;
            /// parameter vector
            gsl_vector* _parameters;
            /// covariance matrix
            gsl_matrix* _covariance;
    };
}
#endif //  SIMPLEFITS_HH
