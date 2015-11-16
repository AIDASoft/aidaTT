#include "utilities.hh"
#include "helixUtils.hh"

#include "aidaTT-Units.hh"

namespace aidaTT
{

    double calculateQoverP(const trackParameters& tp, double bfield)
    {
      //std::cout << " MAGNETIC FIELD CONVERSION FACTOR " <<  convertBr2P_cm << std::endl ;
      if(bfield != 0.)
            return (cos(calculateLambda(tp)) * calculateCurvature(tp) / (bfield * convertBr2P_cm));
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
	const double omega_test = calculateCurvature(tP);

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
        //const double Q    = - qop * bfield.r() * convertBr2P_cm ; // -B*c*q/p
	//const double Q    = - qop * bfield.z() * convertBr2P_cm ; 
        //const double qbar =   qop * bfield.z() * convertBr2P_cm;

	const double Q = cosLambda * omega_test ;
        const double qbar = -1.0*Q ;

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
	
	/*	
	p2l3(2, 0) = -1.;
	p2l3(4, 1) = -1.* (1. + tanLambda * tanLambda);
	p2l3(1, 2) = 1.;
	p2l3(0, 3) = -1.;
	p2l3(3, 4) = 1. ;
	*/
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
	//return cl2p;
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

 std::pair<Vector3D, Vector3D>* calculateLocalCurvilinearSystem(double s, const trackParameters& tP)
  {
    const double omega = calculateCurvature(tP);
    const double phi0 = calculatePhi0(tP);
    const double lambda = calculateLambda(tP);

    const double u0 = - sin(phi0 - omega * s);
    const double u1 =   cos(phi0 - omega * s);
    const double u2 = 0.;

    const double v0 = - cos(phi0 - omega * s) * sin(lambda);
    const double v1 = - sin(phi0 - omega * s) * sin(lambda);
    const double v2 = cos(lambda);

    return new std::pair<Vector3D, Vector3D> (Vector3D(u0, u1, u2), Vector3D(v0, v1, v2));
  }



  std::vector<double>* calculateLocalToMeasurementProjectionMatrix(const Vector3D& clU,  const Vector3D& clV, const std::vector<Vector3D>& measDirs)
  {
    // calculate the projection matrix from the local curvilinear system to the measurement system
    // done in two steps: first compute the easier measurement to local projection, then invert the result
    // note: ALL measurements are taken as two dimensional
    //  dim(M) = 1, then an arbitrary direction is taken, which will allow inversion of the matrix
    //             , also the direction does not contribute, since the associated precision is zero
    // NOTE: higher dimensions (up to five) are possible, but have to be implemented

    /// the return matrix, set as four element double vector
    std::vector<double>* retVec = new std::vector<double>(4);
    double a, b, c, d;

    // two cases to cover:
    // 1D measurements -- create arbitrary orthogonal vector as fake measurement direction
    //      ! not taken into account, since precision is zero
    // 2D measurements -- straighforward

    if(measDirs.size() == 1)
      {
	// construct second orthogonal vector from meas. direction
	Vector3D mdir = measDirs[0].unit();

	Vector3D ortho(mdir.z(), mdir.z(), -mdir.x() - mdir.y());
	// check if the vector is not close to zero
	if(ortho.r() < 1e-6)
	  ortho.fill(- mdir.y() - mdir.z(), mdir.x(), mdir.x());

	a = measDirs[0] * clU;
	b = measDirs[0] * clV;
	c = ortho * clU;
	d = ortho * clV;
      }
    else if(measDirs.size() == 2)
      {
	a = measDirs[0] * clU;
	b = measDirs[0] * clV;
	c = measDirs[1] * clU;
	d = measDirs[1] * clV;
      }
    else
      throw std::invalid_argument("Measurement dimensions > 2 are not yet implemented.");

    // now invert the matrix!
    double determinant = a * d - b * c;

    if(determinant != 0.)
      {
	(*retVec)[0] =   d  / determinant;
	(*retVec)[1] = (-b) / determinant;
	(*retVec)[2] = (-c) / determinant;
	(*retVec)[3] =   a  / determinant;
      }
    else
      {
	throw std::invalid_argument("Projection matrix can't be inverted, bailing out.");
      }

    return retVec;
  }

}
