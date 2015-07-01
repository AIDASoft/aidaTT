#include "utilities.hh"
#include "helixGymnastics.hh"

#include "aidaTT-Units.hh"

namespace aidaTT
{

    double calculateQoverP(const trackParameters& tp, double bfield)
    {
        if(bfield != 0.)
            return (- cos(calculateLambda(tp)) * calculateCurvature(tp) / (bfield * convertBr2P_cm));
        else
            return 0.;
    }

    Vector3D pointOnTrajectory(const trackParameters& tp,  double s) {

      Vector3D p ;
      const Vector3D& rp = tp.referencePoint() ;

      double omega = tp.parameters()( OMEGA ) ;  
      double phi0  = tp.parameters()( PHI0  ) ;
      double tanl  = tp.parameters()( TANL  ) ;
      double d0    = tp.parameters()( D0    ) ;
      double z0    = tp.parameters()( Z0    ) ;

      double sinphi = sin( phi0 ) ;
      double cosphi = cos( phi0 ) ;

      p.x() = rp.x() - d0 * sinphi + (1./omega) * ( sinphi - sin( phi0 - s * omega ) ) ;
      
      p.y() = rp.y() + d0 * cosphi - (1./omega) * ( cosphi - cos( phi0 - s * omega ) ) ;
	  
      p.z() = rp.z() + z0 + s * tanl ;
	  
      return  p ;
    }



    void calculateStartHelix(const Vector3D& x1, const Vector3D& x2,   const Vector3D& x3 , trackParameters& tp , bool backward)
    {

        //-------------------------------------------------------------------------------------------------
        // modified version of original code copied from KalTest::THelicalTrack (2003/10/03  K.Fujii  )
        //------------------------------------------------------------------------------------------------

        static const Vector3D ez(0., 0., 1.);

        Vector3D x12 = x2 - x1 ;
        Vector3D x13 = x3 - x1 ;
        Vector3D x23 = x3 - x2 ;

        x12.z() = 0.;
        x13.z() = 0.;
        x23.z() = 0.;

        double x12mag = x12.r();
        x12 = x12.unit();
        double x13mag = x13.r();
        x13 = x13.unit();
        double x23mag = x23.r();
        x23 = x23.unit();

        double sinHalfPhi23 = x12.cross(x13).z();

        double cosHalfPhi23 = 0.5 * (x13mag / x12mag + (1. - x23mag / x12mag) * (x12mag + x23mag) / x13mag);
        double halfPhi23 = std::atan2(sinHalfPhi23, cosHalfPhi23);

        double r  = -0.5 * x23mag / sinHalfPhi23;

        Vector3D xc = 0.5 * (x2 + x3) + r * cosHalfPhi23 * x23.cross(ez);

        if(backward) r = -r;

        tp.parameters()(OMEGA)  =  1. / r ;
        tp.parameters()(TANL)  = (x2.z() - x3.z()) / (r * 2 * halfPhi23) ;
        tp.parameters()(PHI0)  =  std::atan2(r * (xc.y() - x1.y()),  r * (xc.x() - x1.x())) + M_PI / 2. ;
        tp.parameters()(D0)  =  0. ;
        tp.parameters()(Z0)  =  0. ;

        tp.referencePoint() = x1 ;

    }

    double moveHelixTo(trackParameters& tp,  const Vector3D& refNew)
    {

        //-------------------------------------------------------------------------------------------------
        // modified version of original code copied from KalTest::THelicalTrack (2003/10/03  K.Fujii  )
        //------------------------------------------------------------------------------------------------

        // ---------------------------------------------------
        // (0) Preparation
        // ---------------------------------------------------
        //   Define some numerical constants.
        //

        static const double kPi    = M_PI ;
        static const double kTwoPi = 2.0 * kPi;

        //   Copy helix parmeters to local variables

        double dr    = - tp(D0) ;
        double fi0   =   tp(PHI0) - M_PI / 2. ;

        while(fi0 < 0.)      fi0 += kTwoPi;
        while(fi0 > kTwoPi)  fi0 -= kTwoPi;

        double cpa   =  tp(OMEGA)  ;    //fixme: ???
        double dz    =  tp(Z0) ;
        double tnl   =  tp(TANL) ;

        double x0    = tp.referencePoint().x();
        double y0    = tp.referencePoint().y();
        double z0    = tp.referencePoint().z();

        double xv    = refNew.x();
        double yv    = refNew.y();
        double zv    = refNew.z();

        // ---------------------------------------------------
        // (1) Calculate a' = f_k-1(a_k-1)
        // ---------------------------------------------------
        //        a' = (dr', fi0', cpa', dz', tnl')
        //        a  = (dr , fi0 , cpa , dz , tnl )
        //
        double r     = 1. / cpa;
        double rdr   = r + dr;
        double csf0  = std::cos(fi0);
        double snf0  = std::sqrt(std::max(0.0, (1.0 - csf0) * (1.0 + csf0)));
        if(fi0 > kPi) snf0 = -snf0;

        double xc    = x0 + rdr * csf0;
        double yc    = y0 + rdr * snf0;

        // std::cout << " ***++ xc = " << xc << ", yc = " << yc << " XCenter " << calculateXCenter( tp ) << " YCenter " << calculateYCenter( tp ) << std::endl ;

        double fi0p  = 0.;

        if(cpa > 0.) fi0p = std::atan2((yc - yv), (xc - xv));
        if(cpa < 0.) fi0p = std::atan2((yv - yc), (xv - xc));

        while(fi0p < 0.)      fi0p += kTwoPi;
        while(fi0p > kTwoPi)  fi0p -= kTwoPi;

        double csf   = std::cos(fi0p);
        double snf   = std::sqrt(std::max(0.0, (1.0 - csf) * (1.0 + csf)));
        if(fi0p > kPi) snf = -snf;

        double anrm  = 1.0 / std::sqrt(csf * csf + snf * snf);
        csf  *= anrm;
        snf  *= anrm;
        //~ double csfd  = csf * csf0 + snf * snf0;
        //~ double snfd  = snf * csf0 - csf * snf0;

        double fid   = fi0p - fi0;
        while(fid < 0)      fid += kTwoPi;
        while(fid > kTwoPi) fid -= kTwoPi;
        if(fid > kPi)    fid -= kTwoPi;

        double drp   = (xc - xv) * csf + (yc - yv) * snf - r;
        double dzp   = z0 - zv + dz - r * tnl * fid;

        //fg:
        fi0p += M_PI / 2. ;
        while(fi0p < -kPi)  fi0p += kTwoPi;
        while(fi0p >  kPi)  fi0p -= kTwoPi;

        tp(D0)    = - drp ;
        tp(PHI0)  =   fi0p ;
        tp(OMEGA) =   cpa ;
        tp(Z0)    =   dzp ;
        tp(TANL)  =   tnl;

        tp.setReferencePoint(refNew) ;

        return fid ;
    }



    /// Calculate transformation matrix from curvilinear track parameter corrections to perigee parametrization
    /// \param [return] 5x5 matrix -- the transformation
    /// \param [input] trackParameters -- the reference parameters, at/towards which the corrections are to be applied
    /// \param [input] Vector5 curvilinear corrections -- as calculated e.g. from a fit
    /// \param [input] Vector3D BField -- needed to evaluate expressions
    fiveByFiveMatrix curvilinearToPerigeeJacobian(const trackParameters& tP, const Vector3D& bfield)
    {
        const double qop    = calculateQoverP(tP, bfield.z());
        const double lambda = calculateLambda(tP);
        const double phi0   = calculatePhi0(tP);

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

	//std::cout << " z coordinate of magnetic field " << bfield.z() << " magnitude of magnetic field " << bfield.r() << std::endl ;

        /// TODO: more comments on the conversion of the actual "B" field values
        /// actually (!) they are defined WITH the conversion to 1/R already /intended/
        const double Q    = - qop * bfield.r() * convertBr2P_cm ; // -B*c*q/p
	//const double Q    = - qop * bfield.z() * convertBr2P_cm ; 
        const double qbar =   qop * bfield.z() * convertBr2P_cm;
        jacobian(0, 0) = - bfield.z() * convertBr2P_cm / cosLambda;
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


    /* comment out for now -- ever needed??
        /// Calculate transformation matrix from perigee to curvilinear track parametrization
        /// \param [return] 5x5 matrix -- the wanted transformation
        /// \param [input] trackParameters -- the reference parameters, at/towards which the corrections are to be applied
        /// \param [input] Vector5 perigee track parameter corrections
        /// \param [input] Vector3D BField -- needed to evaluate expressions
        fiveByFiveMatrix perigeeToCurvilinearJacobian(const trackParameters& tP, const Vector5& perParams, const Vector3D& bfield)
        {
            /// FIXME !!// this is probably broken
            const double kappa   = perParams(0);
            const double theta   = perParams(1);
            const double phi     = perParams(2);

            // **************************************************************************
            //fg: make sure this is not called before it is fixed:
            std::cerr  << " **** incomplete method fiveByFiveMatrix perigeeToCurvilinearJacobian() (utilities.cc) called - fix it ! " << std::endl ;
            exit(1) ;
            // **************************************************************************


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
            const double Q = - qop * bfield.r() * convertBr2P_cm; // -Bz*q/p
            jacobian(0, 0) = - sin(theta) / (bfield.z() * convertBr2P_cm) ;
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
    */


    /// Calculate the simple transformation matrix from perigee to L3 parametrization
    /// (kappa, theta, phi, epsilon, z_p) -> ( Omega, tanLambda, phi0, d0, z0)
    /// with d0 = -epsilon, phi0 = phi, omega = -kappa, , z0 = z_p, lambda = pi/2 - theta
    /// \param [return] 5x5 matrix -- the transformation
    /// \param [input] trackParameters -- the reference parameters, at/towards which the corrections are to be applied
    fiveByFiveMatrix perigeeToL3Jacobian(const trackParameters& tP)
    {
        //// FIXME !
        const double tanLambda = calculateTanLambda(tP);

        fiveByFiveMatrix p2l3;

        // differences to zero matrix:
        p2l3(OMEGA, perigeeKAPPA) = -1.; // omega = -kappa
        p2l3(TANL,  perigeeTHETA) = -1. * (1. + tanLambda * tanLambda);
        //p2l3(TANL,  perigeeTHETA) = -1. / (1. + tanLambda * tanLambda); // lambda = pi/2 - theta -> theta( tanLambda = pi/2 - arctan (tanLambda)
        p2l3(PHI0, perigeePHI)    = +1.; // phi0 = phi
        p2l3(D0, perigeeEPSILON)  = -1.; // d0 = -epsilon
        p2l3(Z0, perigeeZP)       = +1.; //  z0 = z_p

        return p2l3;
    }



    /// Calculate the simple transformation matrix from L3 to perigee parametrization
    /// ( Omega, phi0, d0, z0, tanLambda) -> (kappa, theta, phi, epsilon, z_p)
    /// with d0 = -epsilon, phi0 = phi, omega = -kappa, , z0 = z_p, lambda = pi/2 - theta
    /// \param [return] 5x5 matrix -- the transformation
    /// \param [input] trackParameters -- the reference parameters, at/towards which the corrections are to be applied
    fiveByFiveMatrix L3ToPerigeeJacobian(const trackParameters& tP)
    {
        const double tanLambda = calculateTanLambda(tP) ;

        fiveByFiveMatrix l32p;

        // differences to zero matrix: sign differences and the usage of tan Lambda insted of the polar angle
        l32p(0, 0) = -1.; // kappa = - Omega
        l32p(4, 1) = -1. / (1. + tanLambda * tanLambda); // theta = pi/2 - lambda
        l32p(1, 2) = +1.; // phi = phi0
        l32p(2, 3) = -1.; // epsilon = - d0
        l32p(3, 4) = +1.; // z_p = z0

        return l32p;
    }



    /// Calculate transformation matrix from curvilinear track parametrization to L3 parametrization
    /// -- this joins the two transformations from curvilinear to perigee and then from perigee to L3
    /// \param [return] 5x5 matrix -- the transformation
    /// \param [input] trackParameters -- the reference parameters, at/towards which the corrections are to be applied
    /// \param [input] Vector5 curvilinear track parameter corrections
    /// \param [input] Vector3D BField -- needed to evaluate expressions
    fiveByFiveMatrix curvilinearToL3Jacobian(const trackParameters& tP, const Vector3D& bfield)
    {
        // first calculate the perigee transform and parameter values
        fiveByFiveMatrix cl2p = curvilinearToPerigeeJacobian(tP, bfield) ;

        // now the second transformation
        fiveByFiveMatrix p2L3 = perigeeToL3Jacobian(tP);

        // and return the product
        return p2L3 * cl2p;
    }


    /* not needed ???
        /// Calculate transformation matrix from L3 track parametrization to curvilinear parametrization
        /// -- this joins the two transformations from L3 to perigee and then perigee to curvilinear
        /// \param [return] 5x5 matrix -- the transformation
        /// \param [input] trackParameters -- the reference parameters, at/towards which the corrections are to be applied
        /// \param [input] Vector5 L3 parameter corrections
        /// \param [input] Vector3D BField -- needed to evaluate expressions
        fiveByFiveMatrix L3ToCurvilinearJacobian(const trackParameters& tP, const Vector3D& bfield)
        {
            // first move from L3 to perigee
            fiveByFiveMatrix l32p = L3ToPerigeeJacobian(tP);

            // now calculate the second transformation matrix
            fiveByFiveMatrix p2cl = perigeeToCurvilinearJacobian(tP, bfield);

            // return the full product
            return p2cl * l32p;
        }
        */
}
