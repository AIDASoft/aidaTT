#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "trackParameters.hh"
namespace aidaTT
{
    double calculateX0(const trackParameters&);
    double calculateY0(const trackParameters&);
    double calculatePhi0(const trackParameters&);
    double calculateTanLambda(const trackParameters&);
    double calculateLambda(const trackParameters&);
    double calculateZ0(const trackParameters&);
    double calculateDistanceFromPCA(const trackParameters&);
    double calculateCurvature(const trackParameters&);
    double calculateQoverP(const trackParameters& , double BField);
}

#endif // UTITILITIES_HH
