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
        double csfd  = csf * csf0 + snf * snf0;
        double snfd  = snf * csf0 - csf * snf0;

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

    /// calculate and return transformation matrix from curivilinear to perigee track parameters (at reference point)
    fiveByFiveMatrix curvilinearToPerigeeJacobian(const trackParameters& tP, const Vector5& clParams, const Vector3D& bfield)
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

        /// TODO: more comments on the conversion of the actual "B" field values
        /// actually (!) they are defined WITH the conversion to 1/R already /intended/
        const double Q    = - qop * bfield.r() * convertBr2P_cm ; // -B*c*q/p
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



    fiveByFiveMatrix perigeeToCurvilinearJacobian(const trackParameters& tP, const Vector5& perParams, const Vector3D& bfield)
    {
        /// !!!! WRONG - TODO!
        const double kappa   = perParams(0);
        const double theta   = perParams(1);
        const double phi     = perParams(2);

        //**************************************************************************
        //fg: make sure this is not called before it is fixed:
        std::cerr  << " **** incomplete method fiveByFiveMatrix perigeeToCurvilinearJacobian() (utilities.cc) called - fix it ! " << std::endl ;
        exit(1) ;
        //**************************************************************************


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



    fiveByFiveMatrix perigeeToILDJacobian(const trackParameters& tP, const Vector5& perParams)
    {
      const double tanLambda = calculateTanLambda(tP) + perParams( TANL ) ;

        fiveByFiveMatrix per2ILDjacobian;
        per2ILDjacobian.Unit();

        // differences to unit matrix:
        per2ILDjacobian(0, 0) = -1.;
        per2ILDjacobian(1, 1) = -(1. + tanLambda * tanLambda);
        per2ILDjacobian(3, 3) = -1.;

        return per2ILDjacobian;
    }



    fiveByFiveMatrix ildToPerigeeJacobian(const trackParameters& tP, const Vector5& ildParams)
    {
      const double tanLambda = calculateTanLambda(tP) ;

        fiveByFiveMatrix ild2PERjacobian;
        ild2PERjacobian.Unit();

        // differences to unit matrix:
        ild2PERjacobian(0, 0) = -1.;
        ild2PERjacobian(1, 1) = -1 / (1. + tanLambda * tanLambda);
        ild2PERjacobian(3, 3) = -1.;

        return ild2PERjacobian;
    }



    fiveByFiveMatrix curvilinearToILDJacobian(const trackParameters& tP, const Vector5&  clParams, const Vector3D& bfield)
    {
      
      fiveByFiveMatrix cur2per = curvilinearToPerigeeJacobian(tP, clParams, bfield) ;
	
        Vector5 perParams =  cur2per * clParams;

        return perigeeToILDJacobian(tP, perParams) * cur2per ;
    }



    fiveByFiveMatrix ildToCurvilinearJacobian(const trackParameters& tP, const Vector5& ildParams, const Vector3D& bfield)
    {
        Vector5 perParams = ildToPerigeeJacobian(tP, ildParams) * ildParams;

        return perigeeToCurvilinearJacobian(tP, perParams, bfield);
    }
}
