#include "analyticalPropagation.hh"
#include <cmath>

namespace aidaTT
{
bool analyticalPropagation::getJacobian(fiveByFiveMatrix& jac, double dw, double qop, const Vector3D& tStart, const Vector3D& tEnd, const Vector3D& bfield)
{
  // Get analytical helix propagator (in constant magnetic field)
/**
 * Adapted from TRPRFN.F (GEANT3) by Claus Kleinwort (DESY)
 * for curvilinear track parameters (q/p,lambda,phi,x_t,y_t).
 * \param [in] reference to (5*5) propagation matrix -- the computed object
 * \param [in] dw   (3D) arc length to end point
 * \param [in] qop q/p (signed inverse momentum)
 * \param [in] tStart  track direction at start point
 * \param [in] tEnd   track direction at end point
 * \param [in] bfield   B*c (magnetic field *c)
 * \return bool for success
 */

	const double qp = -bfield.r(); // -|B*c|
	const double q = qp * qop; // Q
	if (q == 0.) {
		// line
		jac(3,2) = dw * sqrt(tStart[0] * tStart[0] + tStart[1] * tStart[1]);
		;
		jac(4,1) = dw;
	} else {
		// helix
		// at start
		const double coslambdaStart = sqrt(tStart[0] * tStart[0] + tStart[1] * tStart[1]);
		// at end
		const double coslambdaEnd = sqrt(tEnd[0] * tEnd[0] + tEnd[1] * tEnd[1]);
		const double coslambdaEndInv = 1. / coslambdaEnd;
		// magnetic field direction
		Vector3D hn(bfield);
        hn.unit();
		// (signed) momentum
		const double pav = 1.0 / qop;
		//
		const double theta = q * dw;
		const double sint = sin(theta);
		const double cost = cos(theta);
		const double gamma = hn.dot(tEnd); // H*T
		const Vector3D an1 = hn.cross(tStart); // HxT0
		const Vector3D an2 = hn.cross(tEnd); // HxT
		// U0, V0
		const double au1 = 1. / tStart.trans();
		const Vector3D u1(-au1 * tStart[1], au1 * tStart[0], 0.);
		const Vector3D v1(-tStart[2] * u1[1], tStart[2] * u1[0], tStart[0] * u1[1] - tStart[1] * u1[0]);
		// U, V
		const double au2 = 1. / tEnd.trans();
		const Vector3D u2(-au2 * tEnd[1], au2 * tEnd[0], 0.);
		const Vector3D v2(-tEnd[2] * u2[1], tEnd[2] * u2[0], tEnd[0] * u2[1] - tEnd[1] * u2[0]);
		//
		const double anv = -hn.dot(u2); // N*V=-H*U
		const double anu = hn.dot(v2);  // N*U= H*V
		const double omcost = 1. - cost;
		const double tmsint = theta - sint;
		// M0-M
		const Vector3D dx(-(gamma * tmsint * hn[0] + sint * tStart[0] + omcost * an1[0]) / q,
				-(gamma * tmsint * hn[1] + sint * tStart[1] + omcost * an1[1]) / q,
				-(gamma * tmsint * hn[2] + sint * tStart[2] + omcost * an1[2]) / q);
		// HxU0
		const Vector3D hu1 = hn.cross(u1);
		// HxV0
		const Vector3D hv1 = hn.cross(v1);
		// some dot products
		const double u1u2 = u1.dot(u2), u1v2 = u1.dot(v2), v1u2 = v1.dot(u2), v1v2 = v1.dot(v2);
		const double hu1u2 = hu1.dot(u2), hu1v2 = hu1.dot(v2), hv1u2 = hv1.dot(u2), hv1v2 = hv1.dot(v2);
		const double hnu1 = hn.dot(u1), hnv1 = hn.dot(v1), hnu2 = hn.dot(u2), hnv2 = hn.dot(v2);
		const double tEndu1 = tEnd.dot(u1), tEndv1 = tEnd.dot(v1);
		const double tEnddx = tEnd.dot(dx), u2dx = u2.dot(dx), v2dx = v2.dot(dx);
		const double an2u1 = an2.dot(u1), an2v1 = an2.dot(v1);
		// jacobian
		// 1/P
		jac(0,0) = 1.;
		// Lambda
		jac(1,0) = -qp * anv * tEnddx;
		jac(1,1) = cost * v1v2 + sint * hv1v2 + omcost * hnv1 * hnv2 + anv * (-sint * tEndv1 + omcost * an2v1 - gamma * tmsint * hnv1);
		jac(1,2) = coslambdaStart
				* (cost * u1v2 + sint * hu1v2 + omcost * hnu1 * hnv2 + anv * (-sint * tEndu1 + omcost * an2u1 - gamma * tmsint * hnu1));
		jac(1,3) = -q * anv * tEndu1;
		jac(1,4) = -q * anv * tEndv1;
		// Phi
		jac(2,0) = -qp * anu * tEnddx * coslambdaEndInv;
		jac(2,1) = coslambdaEndInv
				* (cost * v1u2 + sint * hv1u2 + omcost * hnv1 * hnu2 + anu * (-sint * tEndv1 + omcost * an2v1 - gamma * tmsint * hnv1));
		jac(2,2) = coslambdaEndInv * coslambdaStart
				* (cost * u1u2 + sint * hu1u2 + omcost * hnu1 * hnu2 + anu * (-sint * tEndu1 + omcost * an2u1 - gamma * tmsint * hnu1));
		jac(2,3) = -q * anu * tEndu1 * coslambdaEndInv;
		jac(2,4) = -q * anu * tEndv1 * coslambdaEndInv;
		// Xt
		jac(3,0) = pav * u2dx;
		jac(3,1) = (sint * v1u2 + omcost * hv1u2 + tmsint * hnu2 * hnv1) / q;
		jac(3,2) = (sint * u1u2 + omcost * hu1u2 + tmsint * hnu2 * hnu1) * coslambdaStart / q;
		jac(3,3) = u1u2;
		jac(3,4) = v1u2;
		// Yt
		jac(4,0) = pav * v2dx;
		jac(4,1) = (sint * v1v2 + omcost * hv1v2 + tmsint * hnv2 * hnv1) / q;
		jac(4,2) = (sint * u1v2 + omcost * hu1v2 + tmsint * hnv2 * hnu1) * coslambdaStart / q;
		jac(4,3) = u1v2;
		jac(4,4) = v1v2;
	}
    return true;
}
}
