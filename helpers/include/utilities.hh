#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "trackParameters.hh"
#include "fiveByFiveMatrix.hh"
#include "Vector3D.hh"
#include "Vector5.hh"

namespace aidaTT
{
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

    fiveByFiveMatrix curvilinearToPerigeeJacobian(const trackParameters&, const Vector5&, const Vector3D&);
    fiveByFiveMatrix perigeeToCurvilinearJacobian(const trackParameters&, const Vector5&, const Vector3D&);

    fiveByFiveMatrix perigeeToILDJacobian(const trackParameters&, const Vector5&);
    fiveByFiveMatrix ildToPerigeeJacobian(const trackParameters&, const Vector5&);

    fiveByFiveMatrix curvilinearToILDJacobian(const trackParameters&, const Vector5&, const Vector3D&);
    fiveByFiveMatrix ildToCurvilinearJacobian(const trackParameters&, const Vector5&, const Vector3D&);
}

#endif // UTITILITIES_HH
