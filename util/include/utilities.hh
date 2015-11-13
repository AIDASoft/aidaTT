#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "trackParameters.hh"
#include "fiveByFiveMatrix.hh"
#include "Vector3D.hh"
#include "Vector5.hh"

/*** this is a collection of helper functionality
 * the current name "utilities" is a misnomer
 * everything revolves around the transformation or conversion of track parameters
 * -> !!! it should be explicit that the track parameters ARE FIXED in l3 convention
 */


namespace aidaTT
{

    enum
    {
        perigeeKAPPA = 0,
        perigeeTHETA,
        perigeePHI,
        perigeeEPSILON,
        perigeeZP
    };


    /// helper function to calculate the needed/wanted values from the current parametrization
    inline double calculateX0(const trackParameters& tp)
    {
        return (sin(-tp.parameters()(PHI0)) * tp.parameters()(D0)) + tp.referencePoint().x();
    }
    inline double calculateY0(const trackParameters& tp)
    {
        return (cos(tp.parameters()(PHI0)) * tp.parameters()(D0)) + tp.referencePoint().y();
    }

    /// compute the point on the trajectory for a given value of s
    Vector3D pointOnTrajectory(const trackParameters& tp,  double s) ;

    double calculateQoverP(const trackParameters& , double BField);


/// a number of track parametrization conversion matrices
/// externally/interface to outside is "L3" parametrization
/// internally used is curvilinear track parametrization
/// intermediate to them is a different kind of "perigee" parametrization

    fiveByFiveMatrix curvilinearToPerigeeJacobian(const trackParameters&, const Vector3D&);
    //fiveByFiveMatrix perigeeToCurvilinearJacobian(const trackParameters&, const Vector5&, const Vector3D&);

    fiveByFiveMatrix perigeeToL3Jacobian(const trackParameters&);
    fiveByFiveMatrix L3ToPerigeeJacobian(const trackParameters&);

    fiveByFiveMatrix curvilinearToL3Jacobian(const trackParameters&, const Vector3D&);
    //fiveByFiveMatrix L3ToCurvilinearJacobian(const trackParameters&, const Vector3D&);
}

#endif // UTITILITIES_HH
