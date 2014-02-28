#ifndef IPROPAGATION_HH
#define IPROPAGATION_HH

#include "fiveByFiveMatrix.hh"
#include "Vector3D.hh"

/* the propagation method abstraction is one of the key concepts in the tracking toolkit
 *
 * in principle they are simply 5x5 matrices that describe the change of the 5 track parameters
 * between two points
 * 
 * only a single call is possible: getJacobian(...)
 *  
 * the parameters are:
 * \param [in] fiveByFiveMatrix - the jacobian (5*5 matrix) that is calculated in the function call
 * \param [in] double::dw - the 3d path length between the two points 
 * \param [in] double::qop q/p - the signed inverse momentum
 * \param [in] Vector3D::tstart - the track direction at the starting point
 * \param [in] Vector3D::tend - the track direction at the end point
 * \param [in] Vector3D::bfield - the magnetic field vector
 * \returns 
 *      bool value if the calculation was successful
 * input:
 *  the current set of track parameters (even if only a subset is needed)
 *  the position of the two points, expressed in the distance between them -- in arc length!
 * 
 * NOTE: currently a single bfield value is used -- so it is assumed to be constant for the time!
 *
 * 
 */

namespace aidaTT
{
    class IPropagation
    {
            bool  getJacobian(fiveByFiveMatrix& jac, double dw, double qop, const Vector3D& tstart, const Vector3D& tend, const Vector3D& bfield);
    };
}
#endif // IPROPAGATION_HH
