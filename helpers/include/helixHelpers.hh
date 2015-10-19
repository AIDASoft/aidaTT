#ifndef HELIXHELPERS_HH
#define HELIXHELPERS_HH

#include "trackParameters.hh"
#include "Vector3D.hh"
#include "fiveByFiveMatrix.hh"
#include "utilities.hh"

/*** 

 */


namespace aidaTT
{

    /// compute start parameters for a helix from three points, e.g. first, last and middle point )
    void calculateStartHelix(const Vector3D& x1, const Vector3D& x2,   const Vector3D& x3 , trackParameters& tp , bool backward = false) ;

    /// move the helix parameters to a new reference point
    /// FIXME: covariance matrix is not yet updated !
    double moveHelixTo(trackParameters& tp,  const Vector3D& ref) ;


}

#endif // HELIXHELPERS_HH
