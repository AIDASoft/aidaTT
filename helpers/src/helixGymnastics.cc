#include "helixGymnastics.hh"
#include "utilities.hh"


namespace aidaTT
{

    double calculateRadius(const trackParameters& tP)
    {
        const double curvature = calculateCurvature(tP);

        if(curvature != 0.)
            return 1. / curvature;
        return 0.;
    }



    double calculateXCenter(const trackParameters& tP)
    {
        const double curvature = calculateCurvature(tP);
        const double dzero = calculateDistanceFromPCA(tP);
        const double phi0 = calculatePhi0(tP);

        if(curvature != 0.)
            return tP.referencePoint().x() + (1. / curvature - dzero) * sin(phi0);
        return 0.;
    }



    double calculateYCenter(const trackParameters& tP)
    {
        const double curvature = calculateCurvature(tP);
        const double dzero = calculateDistanceFromPCA(tP);
        const double phi0 = calculatePhi0(tP);

        if(curvature != 0.)
            return tP.referencePoint().y() - (1. / curvature - dzero) * cos(phi0);
        return 0.;
    }



    double  calculatePhifromXY(double x, double y, const trackParameters& tP)
    {
        // x0 and y0: p.c.a. coordinates w.r.t reference point
        const double x0   = calculateX0(tP);
        const double y0   = calculateY0(tP);
        const double phi0 = calculatePhi0(tP);
        const double curvature = calculateCurvature(tP);

        return atan2(sin(phi0) - curvature * (x - x0), cos(phi0) + curvature * (y - y0));
    }



    double  calculateSfromXY(double x, double y, const trackParameters& tP)
    {
        // x0 and y0: p.c.a. coordinates w.r.t reference point
        const double x0   = calculateX0(tP);
        const double y0   = calculateY0(tP);
        const double phi0 = calculatePhi0(tP);

        double phi = calculatePhifromXY(x, y, tP);
        double dphi = phi - phi0;

        if(phi != 0.)
            return ((x - x0) * cos(phi0) + (y - y0) * sin(phi0)) / (sin(dphi) / dphi) ;
        else
            return 0.;
    }



    double calculateSfromXY(std::pair<double, double> k, const trackParameters& tp)
    {
        return calculateSfromXY(k.first, k.second, tp);
    }



    double  calculateXfromS(double s, const trackParameters& tP)
    {
        const double x0   = calculateX0(tP);
        const double phi0 = calculatePhi0(tP);
        const double curvature = calculateCurvature(tP);
        if(curvature != 0.)
            return (x0 + 2. / curvature * sin(curvature * s / 2.) * cos(phi0 - curvature * s / 2.));
        else
            return (x0 + s * cos(phi0));
    }



    double  calculateYfromS(double s, const trackParameters& tP)
    {
        const double y0   = calculateY0(tP);
        const double phi0 = calculatePhi0(tP);
        const double curvature = calculateCurvature(tP);
        if(curvature != 0.)
            return (y0 + 2. / curvature * sin(curvature * s / 2.) * sin(phi0 - curvature * s / 2.));
        else
            return (y0 + s * cos(phi0));
    }



    double  calculateZfromS(double s, const trackParameters& tP)
    {
        const double tanLambda = calculateTanLambda(tP);
        const double z0        = calculateZ0(tP);
        return (z0 + s * tanLambda);
    }



    Vector3D calculateTangent(double s, const trackParameters& tP)
    {
        const double omega = calculateCurvature(tP);
        const double phi0 = calculatePhi0(tP);
        const double lambda = calculateLambda(tP);

        double t0 = cos(phi0 - omega * s) * cos(lambda);
        double t1 = sin(phi0 - omega * s) * cos(lambda);
        double t2 = sin(lambda);

        return Vector3D(t0, t1, t2);
    }



    std::pair<Vector3D, Vector3D>* calculateLocalCurvilinearSystem(double s, const trackParameters& tP)
    {
        const double omega = calculateCurvature(tP);
        const double phi0 = calculatePhi0(tP);
        const double lambda = calculateLambda(tP);

        const double u0 = - sin(phi0 - omega * s);
        const double u1 =   cos(phi0 - omega * s);
        const double u2 = 0.;

        const double v0 = - cos(phi0 - omega * s) * sin(lambda);
        const double v1 = - sin(phi0 - omega * s) * sin(lambda);
        const double v2 = cos(lambda);

        return new std::pair<Vector3D, Vector3D> (Vector3D(u0, u1, u2), Vector3D(v0, v1, v2));
    }



    std::vector<Vector3D>* calculateLocalToMeasurementProjectionMatrix(const Vector3D& clU,  const Vector3D& clV, const std::vector<Vector3D>& measDirs)
    {
        std::vector<Vector3D>* retVec = new std::vector<Vector3D>;

        for(std::vector<Vector3D>::const_iterator measDir = measDirs.begin(), last = measDirs.end(); measDir < last; ++measDir)
            {
                Vector3D projection = ((*measDir) * clU) * clU + ((*measDir) * clV) * clV;
                retVec->push_back(projection);
            }
        return retVec;
    }

}
