//~ #include <cmath>
//~ 
//~ #include <gsl/gsl_matrix.h>
//~ 
//~ #include "IBField.hh"
//~ 
//~ namespace aidaTT
//~ {
//~ 
//~ 
//~ /// Simple Helix.
    //~ /**
     //~ * Utilities for simple helix track model build from perigee parameters, assuming constant magnetic field in Z-direction.
     //~ *
     //~ * Propagation and transformation jacobians according to Strandlie & Wittek (NIM A 566(2006) 687-698)
     //~ *
     //~ *  @author C. Kleinwort, DESY (130801)
     //~ */
    //~ class simpleHelix
    //~ {
        //~ public:
            //~ simpleHelix(const double, const double, const double, const double, const double, const double*, const double);
            //~ simpleHelix(const double*, const double*, const double);
            //~ void dump() const;
            //~ double getArcLengthXY(const double*) const;
            //~ void getZSDirection(double*) const;
            //~ void getPosAtArcLength(const double, double*) const;
            //~ bool getExpectedPlanePos(const double*, const double, double&, double&, double&, double&, double&) const;
//~ 
            //~ fiveByFiveMatrix analyticalHelixJacobian(const double, const double) const;
            //~ fiveByFiveMatrix simplifiedHelixJacobian(const double, const double) const;
            //~ fiveByFiveMatrix curvilinearToPerigeeJacobian() const;
//~ 
            //~ fiveByFiveMatrix perigeeToILDJacobian() const;
            //~ fiveByFiveMatrix perigeeToLCIOJacobian() const;
            //~ fiveByFiveMatrix helixToLCIOJacobian() const;
            //~ void moveTo(const double, const double, const double, double*);
            //~ void getStateAt(const double, const double, const double, fiveByFiveMatrixSym&, Vector3D&, fiveByFiveMatrixSym&);
//~ 
        //~ protected:
            //~ /// X of reference point
            //~ const double _xref;
            //~ /// Y of reference point
            //~ const double _yref;
            //~ /// Z of reference point
            //~ const double _zref;
            //~ /// Z component of magnetic field (*c)
            //~ const double _bzc;
            //~ /// 1/R
            //~ const double _rinv;
            //~ /// flight direction at point of closest approach (in XY)
            //~ const double _phi0;
            //~ /// distance of closest approach in (XY)
            //~ const double _dca;
            //~ /// dZ/ds
            //~ const double _dzds;
            //~ /// Z position at distance of closest approach
            //~ const double _z0;
            //~ /// XY circle parameter: X position of center / R
            //~ const double _xRelCenter;
            //~ /// XY circle parameter: Y position of center / R
            //~ const double _yRelCenter;
    //~ };
//~ 
//~ /// Local Helix.
    //~ /**
     //~ * Utilities for local helix track model build from perigee parameters and magnetic field at the PCA
     //~ * (point of closest approach) to the reference point, assuming constant magnetic field in arbitrary direction.
     //~ *
     //~ * In a local system a simple helix is used. The offset of the local system is the PCA of the track (to the reference point)
     //~ * and the directions are (p,n,h), h is the direction of the magnetic field, n = h x t / |h x t| with track direction t and
     //~ * p = n x h. The bending plane is defined by p and n and p is the local track direction at the PCA in this plane.
     //~ *
     //~ * Implemented is a stepwise propagation in a changing (only locally constant) magnetic field.
     //~ *
     //~ * Propagation and transformation jacobians according to Strandlie & Wittek (NIM A 566(2006) 687-698)
     //~ *
     //~ *  @author C. Kleinwort, DESY (131112)
     //~ */
    //~ class localHelix
    //~ {
        //~ public:
            //~ localHelix(const double*, const double*, const double);
            //~ void dump() const;
            //~ fiveByFiveMatrix propagateTo(const double*, const double);
            //~ simpleHelix getSimpleHelix() const;
//~ 
        //~ private:
            //~ /// c * scaling factor (0/1) for magnetic field
            //~ double _bfac;
            //~ /// B-field.
            //~ Vector3D _bFieldc;
            //~ /// offset of local system
            //~ Vector3D _offset;
            //~ /// rotation to local system
            //~ fiveByFiveMatrix _rotation;
            //~ /// rotation from local system
            //~ fiveByFiveMatrix _rotInv;
            //~ /// (global) track direction
            //~ Vector3D _trackDir;
            //~ /// (local) cos(lambda)
            //~ double _cosLambda;
            //~ /// Q/P = const !
            //~ double _qbyp;
            //~ /// local helix parameter
            //~ double _localPar[5];
            //~ /// local reference point
            //~ double _localRef[3];
//~ 
            //~ Vector3D getBFieldc(const double*) const;
            //~ void calcLocal();
            //~ void calcGlobal();
    //~ };
//~ 
//~ /// Simple fit in XY.
    //~ /**
     //~ * Fit circle (Karimaki) or straight line.
     //~ */
    //~ class simpleFitXY
    //~ {
        //~ public:
            //~ simpleFitXY(bool, double, double);
            //~ void addPoint(double, double, double);
            //~ int fit(double&, int&);
            //~ Vector3D getPar() const;
            //~ fiveByFiveMatrixSym getCov() const;
//~ 
        //~ private:
            //~ /// flag for curved (circle) track
            //~ const bool _curved;
            //~ /// number of track parameters (3: circle, 2: line)
            //~ const int _npar;
            //~ /// X of reference point
            //~ const double _xRef;
            //~ /// Y of reference point
            //~ const double _yRef;
            //~ /// number of points (hits) used
            //~ int _numPoints;
            //~ /// weighted sum(x)
            //~ double _sx;
            //~ /// weighted sum(y)
            //~ double _sy;
            //~ /// weighted sum(x*x)
            //~ double _sxx;
            //~ /// weighted sum(x*y)
            //~ double _sxy;
            //~ /// weighted sum(y*y)
            //~ double _syy;
            //~ /// sum of weights
            //~ double _sw;
            //~ /// weighted sum(r*r)
            //~ double _sr;
            //~ /// weighted sum(x*r*r)
            //~ double _sxr;
            //~ /// weighted sum(y*r*r)
            //~ double _syr;
            //~ /// weighted sum(r*r*r*r)
            //~ double _srr;
            //~ /// parameter vector
            //~ Vector3D _parameters;
            //~ /// covariance matrix
            //~ fiveByFiveMatrixSym _covariance;
    //~ };
//~ 
//~ /// Simple fit in ZS.
    //~ /**
     //~ * Fit straight line.
     //~ */
    //~ class simpleFitZS
    //~ {
        //~ public:
            //~ simpleFitZS();
            //~ void addPoint(double, double, double);
            //~ int fit(double&, int&);
            //~ Vector3D getPar() const;
            //~ fiveByFiveMatrixSym getCov() const;
//~ 
        //~ private:
            //~ /// number of track parameters (2)
            //~ const int _npar;
            //~ /// number of points (hits) used
            //~ int _numPoints;
            //~ /// weighted sum(s)
            //~ double _sx;
            //~ /// weighted sum(z)
            //~ double _sy;
            //~ /// weighted sum(s*s)
            //~ double _sxx;
            //~ /// weighted sum(s*z)
            //~ double _sxy;
            //~ /// weighted sum(z*z)
            //~ double _syy;
            //~ /// sum of weights
            //~ double _sw;
            //~ /// parameter vector
            //~ Vector3D _parameters;
            //~ /// covariance matrix
            //~ fiveByFiveMatrixSym _covariance;
    //~ };
//~ 
//~ 
//~ }
