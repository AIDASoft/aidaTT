#include "trajectory.hh"

#include <stdexcept>
#include <utility>
#include <iostream>
#include <algorithm>

#include "intersections.hh"
#include "utilities.hh"
#include "fiveByFiveMatrix.hh"

namespace aidaTT
{
    // constructor (A) for a trajectory -- track found, but needs fitting
    trajectory::trajectory(const trackParameters& tp, IFittingAlgorithm* fa,
                           const IBField* bfield, IPropagation* pm, const IGeometry* geom)
        : _referenceParameters(tp) , _fittingAlgorithm(fa) , _bfield(bfield),  _propagation(pm), _geometry(geom)
    {
        _initialTrajectoryElements.clear();

        /// TODO : implement different BField!
        _bfieldZ = _bfield->Bz(_referenceParameters.referencePoint());
    }



    //~ the minimal useful constructor
    trajectory::trajectory(const trackParameters& tp, const IGeometry* geom) : _referenceParameters(tp), _fittingAlgorithm(NULL) , _bfield(NULL), _propagation(NULL), _geometry(geom)
    {
        _initialTrajectoryElements.clear();
    }



    //~ copy constructor -- doesn't copy any internals!
    trajectory::trajectory(const trajectory& traj) : _referenceParameters(traj._referenceParameters),
        _fittingAlgorithm(traj._fittingAlgorithm), _bfield(traj._bfield), _propagation(traj._propagation), _geometry(traj._geometry)
    {
        _initialTrajectoryElements.clear();
        _intersectionsList.clear();
    }



    trajectory::~trajectory()
    {
        for(std::vector<trajectoryElement*>::iterator element = _initialTrajectoryElements.begin(), last = _initialTrajectoryElements.end(); element < last; ++element)
            delete *element;
    }


    const std::vector<trajectoryElement*>& trajectory::trajectoryElements() const
    {
        return _initialTrajectoryElements;
    }



    const std::vector<std::pair<double, const ISurface*> >& trajectory::getIntersectionsWithSurfaces(const std::list<const aidaTT::ISurface*>& surfaces)
    {
        /// master method for intersection calculation, subdelegates. steps:
        /// 1. calculate all intersections
        /// 2. check for every intersection if inside the bounds
        /// 3. chose the one with the smaller s (!)
        /// 4. build the vector with pairs of s and surface

        for(std::list<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf)
            {
                double s = 0.; /// random value for init, never used
                bool intersects = _calculateIntersectionWithSurface(*surf, s);

                if(intersects)
                    _intersectionsList.push_back(std::make_pair(s, (*surf)));
            }

        return _intersectionsList;
    }



    const fitResults& trajectory::getFitResults()
    {
        return _fittingAlgorithm->getResults();
    }



    bool trajectory::_calculateIntersectionWithSurface(const ISurface* surf, double& s, Vector2D* localUV)
    {
        /// currently three different types of surfaces are available
        if(surf->type().isZCylinder())
            return _intersectsWithinZCylinderBounds(surf, s, localUV);
        else if(surf->type().isZPlane())
            return _intersectWithinZPlaneBounds(surf, s, localUV);
        else if(surf->type().isZDisk())
            return _intersectWithinZDiskBounds(surf, s, localUV);
        else
            throw std::invalid_argument("[aidaTT::trajectory::getIntersectionWithSurfaces] Unknown surface type!");
    }



    bool trajectory::_intersectsWithinZCylinderBounds(const ISurface* surf, double& s, Vector2D* localUV)
    {
        //// TODO:: MISSING IMPLEMENTATION
        return false;
    }



    bool trajectory::_intersectWithinZPlaneBounds(const ISurface* surf, double& s, Vector2D* localUV)
    {
        const Vector3D refpoint = _referenceParameters.referencePoint();

        // the straight line: normals plus distance; distance must be positive !
        const double nx = surf->normal().x();
        const double ny = surf->normal().y();
        const double dist = fabs(surf->distance(refpoint));
        straightLine line(nx, ny, dist);
        // create circle
        const double radius  = _calculateRadius();
        const double xcenter = _calculateXCenter();
        const double ycenter = _calculateYCenter();
        circle circ(xcenter, ycenter, radius);

        intersections candidates = intersectCircleStraightLine(circ, line);

        if(candidates.number() < 1)
            return false;
        else if(candidates.number() == 1)
            {
                const double S = _calculateSfromXY(candidates[0]);
                const double Z = _calculateZfromS(S);
                Vector3D thePlace(candidates[0].first, candidates[0].second, Z);
                bool inside = surf->insideBounds(thePlace);

                if(inside)
                    {
                        s = S;
                        return true;
                    }
                else
                    return false;
            }

        ///  else -- the standard case: two solutions index 0 and 1
        /// calculate all values first, then evaluate
        const double X0 = candidates[0].first;
        const double Y0 = candidates[0].second;
        const double S0 = _calculateSfromXY(X0, Y0);
        const double Z0 = _calculateZfromS(S0);

        const double X1 = candidates[1].first;
        const double Y1 = candidates[1].second;
        const double S1 = _calculateSfromXY(X1, Y1);
        const double Z1 = _calculateZfromS(S1);
        Vector3D sol0(X0, Y0, Z0);
        Vector3D sol1(X1, Y1, Z1);

        const bool insideFirst  = surf->insideBounds(sol0);
        const bool insideSecond = surf->insideBounds(sol1);

        if((!insideFirst && !insideSecond) || (S0 < 0. && S1 < 0.))      // discard negative or no solution
            return false;
        else if(insideFirst && S0 > 0. && !insideSecond)
            {
                s = S0;
                if(localUV != NULL)
                    _calculateLocalCoordinates(surf, sol0, localUV);
                return true;
            }
        else if(!insideFirst && insideSecond &&  S1 > 0.)
            {
                s = S1;
                if(localUV != NULL)
                    _calculateLocalCoordinates(surf, sol1, localUV);

                return true;
            }
        else // both are valid , choose the smaller solution
            {
                if(S0 > 0. && S1 < 0.)
                    {
                        s = S0;
                        if(localUV != NULL)
                            _calculateLocalCoordinates(surf, sol0, localUV);
                    }
                else if(S0 < 0. && S1 > 0.)
                    {
                        s = S1;
                        if(localUV != NULL)
                            _calculateLocalCoordinates(surf, sol1, localUV);
                    }
                else // both are positive
                    {
                        if(S0 < S1)
                            {
                                s = S0;
                                if(localUV != NULL)
                                    _calculateLocalCoordinates(surf, sol0, localUV);
                            }
                        else
                            {
                                s = S1;
                                if(localUV != NULL)
                                    _calculateLocalCoordinates(surf, sol1, localUV);
                            }
                    }
                return true;
            }
    }



    bool trajectory::_intersectWithinZDiskBounds(const ISurface* surf, double& s, Vector2D* localUV)
    {
        // the z position of the plane
        double planePositionZ = surf->origin().z();
        s = (planePositionZ - calculateZ0(_referenceParameters)) / calculateTanLambda(_referenceParameters);
        double x = _calculateXfromS(s);
        double y = _calculateYfromS(s);

        Vector3D thePlace = Vector3D(x, y, planePositionZ);

        if(surf->insideBounds(thePlace))
            {
                if(localUV != NULL)
                    _calculateLocalCoordinates(surf, thePlace, localUV);
                return true;
            }
        else
            return false;
    }



    void trajectory::_calculateLocalCoordinates(const ISurface* surf, const Vector3D& position, Vector2D* localUV)
    {
        Vector2D local = surf->globalToLocal(position);
        localUV->_u = local.u();
        localUV->_v = local.v();
    }



    double  trajectory::_calculatePhifromXY(double x, double y) const
    {
        // x0 and y0: p.c.a. coordinates w.r.t reference point
        static const double x0   = calculateX0(_referenceParameters);
        static const double y0   = calculateY0(_referenceParameters);
        static const double phi0 = calculatePhi0(_referenceParameters);
        static const double curvature = calculateCurvature(_referenceParameters);

        return atan2(sin(phi0) - curvature * (x - x0), cos(phi0) + curvature * (y - y0));
    }



    double  trajectory::_calculateSfromXY(double x, double y) const
    {
        // x0 and y0: p.c.a. coordinates w.r.t reference point
        static const double x0   = calculateX0(_referenceParameters);
        static const double y0   = calculateY0(_referenceParameters);
        static const double phi0 = calculatePhi0(_referenceParameters);

        double phi = _calculatePhifromXY(x, y);
        return ((x - x0) * cos(phi0) + (y - y0) * sin(phi0)) / (sin(phi) / phi) ;
    }



    double  trajectory::_calculateXfromS(double s) const
    {
        static const double x0   = calculateX0(_referenceParameters);
        static const double phi0 = calculatePhi0(_referenceParameters);
        static const double curvature = calculateCurvature(_referenceParameters);
        if(curvature != 0.)
            return (x0 + 2. / curvature * sin(curvature * s / 2.) * cos(phi0 - curvature * s / 2.));
        else
            return (x0 + s * cos(phi0));
    }



    double  trajectory::_calculateYfromS(double s) const
    {
        static const double y0   = calculateY0(_referenceParameters);
        static const double phi0 = calculatePhi0(_referenceParameters);
        static const double curvature = calculateCurvature(_referenceParameters);
        if(curvature != 0.)
            return (y0 + 2. / curvature * sin(curvature * s / 2.) * sin(phi0 - curvature * s / 2.));
        else
            return (y0 + s * cos(phi0));
    }



    double  trajectory::_calculateZfromS(double s) const
    {
        static const double tanLambda = calculateTanLambda(_referenceParameters);
        static const double z0        = calculateZ0(_referenceParameters);
        return z0 + s * tanLambda;
    }



    Vector3D trajectory::_calculateTangent(double s)
    {
        static const double omega = calculateCurvature(_referenceParameters);
        static const double phi0 = calculatePhi0(_referenceParameters);
        static const double lambda = calculateLambda(_referenceParameters);

        double t0 = cos(phi0 - omega * s) * cos(lambda);
        double t1 = sin(phi0 - omega * s) * cos(lambda);
        double t2 = sin(lambda);

        return Vector3D(t0, t1, t2);
    }


    std::pair<Vector3D, Vector3D>* trajectory::_calculateLocalCurvilinearSystem(double s)
    {
        static const double omega = calculateCurvature(_referenceParameters);
        static const double phi0 = calculatePhi0(_referenceParameters);
        static const double lambda = calculateLambda(_referenceParameters);

        const double u0 = - sin(phi0 - omega * s);
        const double u1 =   cos(phi0 - omega * s);
        const double u2 = 0.;

        const double v0 = - cos(phi0 - omega * s) * sin(lambda);
        const double v1 = - sin(phi0 - omega * s) * sin(lambda);
        const double v2 = cos(lambda);

        return new std::pair<Vector3D, Vector3D> (Vector3D(u0, u1, u2), Vector3D(v0, v1, v2));
    }




    double trajectory::_calculateRadius() const
    {
        static const double curvature = calculateCurvature(_referenceParameters);

        if(curvature != 0.)
            return 1. / curvature;
        return 0.;
    }



    double trajectory::_calculateXCenter() const
    {
        static const double curvature = calculateCurvature(_referenceParameters);
        static const double dzero = calculateDistanceFromPCA(_referenceParameters);
        static const double phi0 = calculatePhi0(_referenceParameters);

        if(curvature != 0.)
            return _referenceParameters.referencePoint().x() + (1. / curvature - dzero) * sin(phi0);
        return 0.;
    }



    double trajectory::_calculateYCenter() const
    {
        static const double curvature = calculateCurvature(_referenceParameters);
        static const double dzero = calculateDistanceFromPCA(_referenceParameters);
        static const double phi0 = calculatePhi0(_referenceParameters);

        if(curvature != 0.)
            return _referenceParameters.referencePoint().y() - (1. / curvature - dzero) * cos(phi0);
        return 0.;
    }


    void trajectory::addMeasurement(const Vector3D& position, const std::vector<double>& resolution, const ISurface& surface, void* id)
    {
        Vector2D* referenceUV = new Vector2D();
        double s =  0;
        _calculateIntersectionWithSurface(&surface, s, referenceUV);

        Vector2D* measuredUV = new Vector2D(surface.globalToLocal(position));

        std::pair<Vector2D*, Vector2D*>* theUVs = new std::pair<Vector2D*, Vector2D*> (measuredUV, referenceUV);

        std::vector<Vector3D>* measDir = new std::vector<Vector3D>;
        measDir->push_back(surface.u(position));
        measDir->push_back(surface.v(position));

        _initialTrajectoryElements.push_back(new trajectoryElement(s, surface, measDir, resolution, theUVs, _calculateLocalCurvilinearSystem(s), id));
    }



    void trajectory::addElement(const Vector3D& point, void* id)
    {
        double s =  _calculateSfromXY(point.x(), point.y());
        _initialTrajectoryElements.push_back(new trajectoryElement(s, id));
    }



    void trajectory::addElement(const Vector3D& point, const ISurface& surface, void* id)
    {
        double s =  _calculateSfromXY(point.x(), point.y());
// REDO!        _initialTrajectoryElements.push_back(new trajectoryElement(s, surface, id));
    }



    static    bool compareTrajectoryElements(trajectoryElement* one, trajectoryElement* two)
    {
        return (one->arcLength() < two->arcLength());
    }



    void trajectory::prepareForFitting()
    {
        ///~ first sort the trajectory elements by arclength
        sort(_initialTrajectoryElements.begin(), _initialTrajectoryElements.end(), compareTrajectoryElements);

        static const double cosLambda = cos(calculateLambda(_referenceParameters));

        ///~ TODO: this uses a static b field for now ...
        const double qbyp      = calculateQoverP(_referenceParameters, _bfieldZ) ;
        const Vector3D BField(0., 0., _bfieldZ);

        ///~ calculate and add the jacobians

        /// the first jacobian is useless, just use an empty 5x5 matrix
        if(_initialTrajectoryElements.size() > 0)
            {
                fiveByFiveMatrix* j = new fiveByFiveMatrix;
                (_initialTrajectoryElements.at(0))->setJacobian(j);
            }
        /// now the rest
        for(std::vector<trajectoryElement*>::iterator element = _initialTrajectoryElements.begin() + 1, last = _initialTrajectoryElements.end(); element < last; ++element)
            {
                ///~ obtain the two arc lengths
                double prevS = (*(element - 1))->arcLength();
                double currS = (*element)->arcLength();
                ///~ calculate 3D arclength -> divide 2D arc length by cos lambda
                double dw = (currS - prevS) / cosLambda;

                Vector3D tstart = _calculateTangent(prevS);
                Vector3D tend   = _calculateTangent(currS);

                fiveByFiveMatrix* jacob = new fiveByFiveMatrix;
                _propagation->getJacobian(*jacob, dw, qbyp, tstart, tend, BField);
                (*element)->setJacobian(jacob);
            }
    }



    bool trajectory::fit()
    {
        return _fittingAlgorithm->fit(*this);
    }


}
