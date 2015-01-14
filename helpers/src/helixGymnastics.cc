#include "helixGymnastics.hh"
#include "utilities.hh"


namespace aidaTT
{

    double calculateRadius(const trackParameters& tP)
    {
        const double curvature = calculateCurvature(tP);

        if(curvature != 0.)
            return fabs(1. / curvature);
        return 0.;
    }



    double calculateXCenter(const trackParameters& tP)
    {
        const double omega = calculateOmega(tP);
        if(omega == 0.)
            return 0.;

        //fg: need signed radius here (1/omega)
        const double radius = 1. / omega;
        const double dzero = calculateD0(tP);
        const double phi0 = calculatePhi0(tP);

        return tP.referencePoint().x()  + (radius - dzero) * sin(phi0);
    }



    double calculateYCenter(const trackParameters& tP)
    {
        const double omega = calculateOmega(tP);
        if(omega == 0.)
            return 0.;

        const double radius = 1. / omega;
        const double dzero = calculateD0(tP);
        const double phi0 = calculatePhi0(tP);

        return tP.referencePoint().y() - (radius - dzero) * cos(phi0);
    }



    double  calculatePhifromXY(double x, double y, const trackParameters& tP)
    {
        // x0 and y0: p.c.a. coordinates - not w.r.t reference point
        const double x0   = calculateX0(tP);
        const double y0   = calculateY0(tP);
        const double phi0 = calculatePhi0(tP);
        const double curvature = calculateCurvature(tP);

        return atan2(sin(phi0) - curvature * (x - x0), cos(phi0) + curvature * (y - y0));
    }



    double  calculateSfromXY(double x, double y, const trackParameters& tP)
    {
        // x0 and y0: p.c.a. coordinates w.r.t reference point
        const double x0   = calculateX0(tP) ;
        const double y0   = calculateY0(tP) ;
        const double phi0 = calculatePhi0(tP);

        double phi = calculatePhifromXY(x, y, tP);
        double dphi = phi - phi0;

        if(dphi != 0.)
            {

                double s = ((x - x0) * cos(phi0) + (y - y0) * sin(phi0)) / (sin(dphi) / dphi)  ;

                // std::cout <<  " *** calculateSfromXY :  s = " << s
                //        << " dphi*R = " << dphi / std::fabs( calculateCurvature( tP ) )
                //        << " dphi = " << dphi
                //        << std::endl ;

                return s ;

            }
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
        return (z0 + tP.referencePoint().z() + s * tanLambda);
        //  return (z0 + s * tanLambda);
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



    std::vector<double>* calculateLocalToMeasurementProjectionMatrix(const Vector3D& clU,  const Vector3D& clV, const std::vector<Vector3D>& measDirs)
    {
        // calculate the projection matrix from the local curvilinear system to the measurement system
        // done in two steps: first compute the easier measurement to local projection, then invert the result
        // note: ALL measurements are taken as two dimensional
        //  dim(M) = 1, then an arbitrary direction is taken, which will allow inversion of the matrix
        //             , also the direction does not contribute, since the associated precision is zero
        // NOTE: higher dimensions (up to five) are possible, but have to be implemented

        /// the return matrix, set as four element double vector
        std::vector<double>* retVec = new std::vector<double>(4);
        double a, b, c, d;

        // two cases to cover:
        // 1D measurements -- create arbitrary orthogonal vector as fake measurement direction
        //      ! not taken into account, since precision is zero
        // 2D measurements -- straighforward

        if(measDirs.size() == 1)
            {
                // construct second orthogonal vector from meas. direction
                Vector3D mdir = measDirs[0].unit();

                Vector3D ortho(mdir.z(), mdir.z(), -mdir.x() - mdir.y());
                // check if the vector is not close to zero
                if(ortho.r() < 1e-6)
                    ortho.fill(- mdir.y() - mdir.z(), mdir.x(), mdir.x());

                a = measDirs[0] * clU;
                b = measDirs[0] * clV;
                c = ortho * clU;
                d = ortho * clV;
            }
        else if(measDirs.size() == 2)
            {
                a = measDirs[0] * clU;
                b = measDirs[0] * clV;
                c = measDirs[1] * clU;
                d = measDirs[1] * clV;
            }
        else
            throw std::invalid_argument("Measurement dimensions > 2 are not yet implemented.");

        // now invert the matrix!
        double determinant = a * d - b * c;

        if(determinant != 0.)
            {
                (*retVec)[0] = 1. / determinant * d;
                (*retVec)[1] = 1. / determinant * (-b);
                (*retVec)[2] = 1. / determinant * (-c);
                (*retVec)[3] = 1. / determinant * a;
            }
        else
            {
                throw std::invalid_argument("Projection matrix can't be inverted, bailing out.");
            }

        return retVec;
    }

}
