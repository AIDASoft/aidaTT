#include "trajectory.hh"

#include <stdexcept>
#include <utility>
#include <iostream>

#include "intersections.hh"
#include "utilities.hh"

namespace aidaTT
{
// the only constructor for a trajectory: known track parameters and trajectory elements
    trajectory::trajectory(const trackParameters& tp,
                           const std::vector<trajectoryElement>& elems,
                           const IFittingAlgorithm* fa, const  IPropagation* pm, const IGeometry* geom)
        : _referenceParameters(tp) , _fittingAlgorithm(fa) , _propagation(pm), _geometry(geom)
    {
        // copy the elements
        // TODO like this? trajectoryElements should know about their position?
        //            _initialTrajectoryElements.resize(elems.size());
        //~ for(vector<trajectoryElement>::const_iterator elem = elems.begin(), last = elems.end(); elem < last; ++elem)
        //~ _initialTrajectoryElements.push_back(*elem);
    }



    //~ the minimal useful constructor
    trajectory::trajectory(const trackParameters& tp, const IGeometry* geom) : _referenceParameters(tp), _fittingAlgorithm(NULL) , _propagation(NULL), _geometry(geom)
    {}



    trajectory::~trajectory()
    {}


    const std::vector<trajectoryElement>& trajectory::getTrajectoryElements() const
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
                double __s = 0.; /// random value for init, never used
                bool __intersects = false;

                /// currently three different types of surfaces are available
                if((*surf)->type().isZCylinder())
                    __intersects = _intersectsWithinZCylinderBounds(*surf, __s);
                else if((*surf)->type().isZPlane())
                    __intersects = _intersectWithinZPlaneBounds(*surf, __s);
                else if((*surf)->type().isZDisk())
                    __intersects = _intersectWithinZDiskBounds(*surf, __s);
                else
                    throw std::invalid_argument("[aidaTT::trajectory::getIntersectionWithSurfaces] Unknown surface type!");

                if(__intersects)
                    _intersectionsList.push_back(std::make_pair(__s, (*surf)));
            }

        return _intersectionsList;
    }



    bool trajectory::_intersectsWithinZCylinderBounds(const ISurface* surf, double& s)
    {
        return false;
    }



    bool trajectory::_intersectWithinZPlaneBounds(const ISurface* surf, double& s)
    {
        const Vector3D refpoint = _referenceParameters.referencePoint();
        // the straight line: normals plus distance
        const double nx = surf->normal().x();
        const double ny = surf->normal().y();

        const double dist = surf->distance(refpoint);

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
                std::cout << " s is " << S << " [x, y , z] = [" << candidates[0].first << " , " <<  candidates[0].second  << " , " << Z << "]" << std::endl;
                Vector3D thePlace(candidates[0].first * .1, candidates[0].second * .1, Z * .1);
                std::cout << " pass vector to function: " << thePlace << std::endl;
                bool inside = surf->insideBounds(thePlace, 10000000000000000000.);
                std::cout << " the distance is supposedly : " << surf->distance(thePlace) << std::endl;

                if(inside)
                    {
                        std::cout << " YES ! " << std::endl;
                        s = S;
                        return true;
                    }
                else
                    {
                        std::cout << " out of bounds " << std::endl;
                        return false;
                    }
            }

        ///  else -- the standard case: two solutions index 0 and 1
        /// calculate all values first, then evaluate
        const double X0 = candidates[0].first;
        const double Y0 = candidates[0].second;
        const double X1 = candidates[1].first;
        const double Y1 = candidates[1].second;
        const double S0 = _calculateSfromXY(X0, Y0);
        const double S1 = _calculateSfromXY(X1, Y1);
        const double Z0 = _calculateZfromS(S0);
        const double Z1 = _calculateZfromS(S1);
        const bool insideFirst  =  surf->insideBounds(Vector3D(X0, Y0, Z0));
        const bool insideSecond = surf->insideBounds(Vector3D(X1, Y1, Z1));
        if(!insideFirst && !insideSecond)
            return false;
        else if(insideFirst && !insideSecond)
            {
                s = S0;
                return true;
            }
        else if(!insideFirst && insideSecond)
            {
                s = S1;
                return true;
            }
        else // both are valid , choose the smaller solution
            {
                s = S0;
                if(S1 < S0)
                    s = S1;
                return true;
            }
    }



    bool trajectory::_intersectWithinZDiskBounds(const ISurface* surf, double& s)
    {
        // the z position of the plane
        //double __planePositionZ = surf->origin().z();

        return false;
    }



    double  trajectory::_calculateSfromXY(double x, double y) const
    {
        // x0 and y0: p.c.a. coordinates w.r.t reference point
        const double x0   = calculateX0(_referenceParameters);
        const double y0   = calculateY0(_referenceParameters);
        const double phi0 = calculatePhi0(_referenceParameters);
        return (x - x0) * cos(phi0) + (y - y0) * sin(phi0);
    }



    double  trajectory::_calculateZfromS(double s) const
    {
        const double tanLambda = calculateTanLambda(_referenceParameters);
        const double z0        = calculateZ0(_referenceParameters);
        return z0 + s * tanLambda;
    }



    double trajectory::_calculateRadius() const
    {
        return 1000.;
    }



    double trajectory::_calculateXCenter() const
    {
        return 10.;
    }



    double trajectory::_calculateYCenter() const
    {
        return 0.;
    }

}
