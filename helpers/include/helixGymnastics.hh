#ifndef HELIXGYMNASTICS_HH
#define HELIXGYMNASTICS_HH

#include "trackParameters.hh"

namespace aidaTT
{
    double calculateRadius(const trackParameters&) ;
    double calculateXCenter(const trackParameters&) ;
    double calculateYCenter(const trackParameters&) ;

    double calculatePhifromXY(double, double, const trackParameters&) ;

    double calculateSfromXY(double, double, const trackParameters&) ;
    double calculateSfromXY(std::pair<double, double>, const trackParameters&);

    double calculateXfromS(double, const trackParameters&) ;
    double calculateYfromS(double, const trackParameters&) ;
    double calculateZfromS(double, const trackParameters&) ;

    Vector3D calculateTangent(double, const trackParameters&);

    std::pair<Vector3D, Vector3D>* calculateLocalCurvilinearSystem(double, const trackParameters&);

    std::vector<double>* calculateLocalToMeasurementProjectionMatrix(const Vector3D&, const Vector3D&, const std::vector<Vector3D>&);

}
#endif // HELIXGYMNASTICS_HH
