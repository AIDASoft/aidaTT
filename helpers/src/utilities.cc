#include "utilities.hh"

namespace aidaTT
{

    /// SUPER DUMMY !!! only L3 type: [ Omega, tan(lambda), phi_0, d_0, z_0 ]

  //fg: this is what we have in CEDViewer :
  // double xs = ts->getReferencePoint()[0] -  ts->getD0() * sin( ts->getPhi() ) ;
  // double ys = ts->getReferencePoint()[1] +  ts->getD0() * cos( ts->getPhi() ) ;
  // double zs = ts->getReferencePoint()[2] +  ts->getZ0() ;

  // x0 and y0: p.c.a. coordinates w.r.t reference point

    double calculateX0(const trackParameters& tp)
    {
        /// L3 type only
      //      return cos(tp.parameters()(2)) * tp.parameters()(3)  / 10.;
      return sin( - tp.parameters()(2)) * tp.parameters()(3)  / 10.;
    }



    double calculateY0(const trackParameters& tp)
    {
        /// L3 type only
      //return sin(tp.parameters()(2)) * tp.parameters()(3)  / 10.;
      return cos(tp.parameters()(2)) * tp.parameters()(3)  / 10.;
    }



    double calculatePhi0(const trackParameters& tp)
    {
        /// L3 type only
        return tp.parameters()(2);
    }



    double  calculateTanLambda(const trackParameters& tp)
    {
        /// L3 type only
        return tp.parameters()(1);
    }



    double  calculateLambda(const trackParameters& tp)
    {
        /// L3 type only
        return atan(tp.parameters()(1));
    }



    double calculateZ0(const trackParameters& tp)
    {
        /// L3 type only
        return tp.parameters()(4) / 10.;
    }



    double calculateDistanceFromPCA(const trackParameters& tp)
    {
        /// L3 type only
        return tp.parameters()(3) / 10.;
    }



    double calculateD0(const trackParameters& tp)
    {
        /// L3 type only
        return tp.parameters()(3) / 10.;
    }



    double calculateCurvature(const trackParameters& tp)
    {
        /// L3 type only
        return tp.parameters()(0) * 10.;
    }



    double calculateOmega(const trackParameters& tp)
    {
        /// L3 type only
        return tp.parameters()(0) * 10.;
    }



    double calculateQoverP(const trackParameters& tp, double bfield)
    {
        if(bfield != 0.)
            return (- cos(calculateLambda(tp)) * calculateCurvature(tp) / bfield);
        else
            return 0.;
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



    fiveByFiveMatrix perigeeToCurvilinearJacobian(const trackParameters& tP, const Vector5& perParams, const Vector3D& bfield)
    {
        /// !!!! WRONG - TODO!
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



    fiveByFiveMatrix perigeeToILDJacobian(const trackParameters& tP, const Vector5& perParams)
    {
        const double tanLambda = calculateTanLambda(tP);

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
        const double tanLambda = calculateTanLambda(tP);

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
        Vector5 perParams = curvilinearToPerigeeJacobian(tP, clParams, bfield) * clParams;

        return perigeeToILDJacobian(tP, perParams);
    }



    fiveByFiveMatrix ildToCurvilinearJacobian(const trackParameters& tP, const Vector5& ildParams, const Vector3D& bfield)
    {
        Vector5 perParams = ildToPerigeeJacobian(tP, ildParams) * ildParams;

        return perigeeToCurvilinearJacobian(tP, perParams, bfield);
    }
}
