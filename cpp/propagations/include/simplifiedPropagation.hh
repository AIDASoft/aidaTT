#ifndef SIMPLIFIEDPROPAGATION_HH
#define SIMPLIFIEDPROPAGATION_HH

#include "IPropagation.hh"

/** simplified helix propagator in constant solenoidal magnetic field
 * calculated in curvilinear track parameters (q/p,lambda,phi,x_t,y_t)
 *
 * The propagator is quadratic in the arc length difference.
 * The only other used information is the direction of the track at the start point.
**/

namespace aidaTT
{

    class simplifiedPropagation : public IPropagation
    {
            bool  getJacobian(fiveByFiveMatrix& jacobian, double dw, double qop, const Vector3D& tstart, const Vector3D& tend, const Vector3D& bfield);
    };

}

#endif //SIMPLIFIEDPROPAGATION_HH
