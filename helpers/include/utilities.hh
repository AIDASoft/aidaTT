#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "trackParameters.hh"
#include "fiveByFiveMatrix.hh"
#include "Vector3D.hh"
#include "Vector5.hh"

namespace aidaTT
{
    enum {
      OMEGA = 0, 
      TANL, 
      PHI0,
      D0,
      Z0 
    } ;

  
  
	/// helper function to calculate the needed/wanted values from the current parametrization
    double calculateX0(const trackParameters&);
    double calculateY0(const trackParameters&);
    double calculatePhi0(const trackParameters&);
    double calculateTanLambda(const trackParameters&);
    double calculateLambda(const trackParameters&);
    double calculateZ0(const trackParameters&);
    double calculateDistanceFromPCA(const trackParameters&);
    double calculateD0(const trackParameters&);
    double calculateCurvature(const trackParameters&);
    double calculateOmega(const trackParameters&);
    double calculateQoverP(const trackParameters& , double BField);

  /// compute start parameters for a helix from three points, e.g. first, last and middle point ) 
  void calculateStartHelix( const Vector3D& x1, const Vector3D& x2,   const Vector3D& x3 , trackParameters& tp , bool backward=false) ;

  /// move the helix parameters to a new reference point - fixme: covariance matrix is not yet updated !
  double moveHelixTo( trackParameters& tp,  const Vector3D& ref ) ;

    fiveByFiveMatrix curvilinearToPerigeeJacobian(const trackParameters&, const Vector5&, const Vector3D&);
    fiveByFiveMatrix perigeeToCurvilinearJacobian(const trackParameters&, const Vector5&, const Vector3D&);

    fiveByFiveMatrix perigeeToILDJacobian(const trackParameters&, const Vector5&);
    fiveByFiveMatrix ildToPerigeeJacobian(const trackParameters&, const Vector5&);

    fiveByFiveMatrix curvilinearToILDJacobian(const trackParameters&, const Vector5&, const Vector3D&);
    fiveByFiveMatrix ildToCurvilinearJacobian(const trackParameters&, const Vector5&, const Vector3D&);
}

#endif // UTITILITIES_HH
