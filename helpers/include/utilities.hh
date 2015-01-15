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
 *
 * the two helix helpers should be moved out into a separate class (and verified btw.)
 */


namespace aidaTT
{

    /** Enum for accessing the track parameters in L3/LCIO perigee convention */
    enum
    {
        OMEGA = 0,
        TANL,
        PHI0,
        D0,
        Z0
    } ;

    enum
    {
        perigeeKAPPA = 0,
        perigeeTHETA,
        perigeePHI,
        perigeeEPSILON,
        perigeeZP
    };


    /// backward compatible accessors
    inline double calculateOmega(const trackParameters& tp)
    {
        return tp.parameters()(OMEGA);
    }
    inline double calculateTanLambda(const trackParameters& tp)
    {
        return tp.parameters()(TANL);
    }
    inline double calculatePhi0(const trackParameters& tp)
    {
        return tp.parameters()(PHI0);
    }
    inline double calculateD0(const trackParameters& tp)
    {
        return tp.parameters()(D0);
    }
    inline double calculateZ0(const trackParameters& tp)
    {
        return tp.parameters()(Z0);
    }

    inline double calculateCurvature(const trackParameters& tp)
    {
        return tp.parameters()(OMEGA);
    }
    inline double calculateLambda(const trackParameters& tp)
    {
        return atan(tp.parameters()(TANL));
    }

    /// helper function to calculate the needed/wanted values from the current parametrization
    inline double calculateX0(const trackParameters& tp)
    {
        return (sin(-tp.parameters()(PHI0)) * tp.parameters()(D0)) + tp.referencePoint().x();
    }
    inline double calculateY0(const trackParameters& tp)
    {
        return (cos(tp.parameters()(PHI0)) * tp.parameters()(D0)) + tp.referencePoint().y();
    }

    double calculateQoverP(const trackParameters& , double BField);

    /// compute start parameters for a helix from three points, e.g. first, last and middle point )
    void calculateStartHelix(const Vector3D& x1, const Vector3D& x2,   const Vector3D& x3 , trackParameters& tp , bool backward = false) ;

    /// move the helix parameters to a new reference point
    /// FIXME: covariance matrix is not yet updated !
    double moveHelixTo(trackParameters& tp,  const Vector3D& ref) ;


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
