#include "trajectory.hh"
#include <stdexcept>
#include <utility>


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


    const std::vector<trajectoryElement>& trajectory::getTrajectoryElements() const
    {
        return _initialTrajectoryElements;
    }



    const std::vector<std::pair<double, const ISurface*> >& trajectory::getIntersectionsWithSurfaces(const std::list<const aidaTT::ISurface*>& surfaces)
    {
        /// 1. calculate all intersections
        /// 2. chose the one with the smaller s (!)
        /// 3. check for every intersection if inside the bounds
        /// 4. build the vector with pairs of s and surface

        //~ for(std::list<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf)
        //~ {
        //~ double __s = -2.;
        //~ bool __intersects = false;
        //~ if((*surf)->type().isZCylinder())
        //~ __intersects = _intersectsWithinZCylinderBounds(*surf, __s);
        //~ else if((*surf)->type().isZPlane())
        //~ __intersects = _intersectWithinZPlaneBounds(*surf, __s);
        //~ else if((*surf)->type().isZDisk())
        //~ __intersects = _intersectWithinZDiskBounds(*surf, __s);
        //~ else
        //~ throw std::invalid_argument("[aidaTT::trajectory::getIntersectionWithSurfaces] Unknown surface type!");
//~
        //~ if( __intersects)
        //~ _intersectionsList.push_back(make_pair(__s, (*surf)));
        //~ }

        return _intersectionsList;
    }


    //~
    //~ bool trajectory::_intersectsWithinZCylinderBounds(const ISurface* surf, double& s)
    //~ {
    //~ return false;
    //~ }
    //~
    //~
    //~
    //~ bool trajectory::_intersectWithinZPlaneBounds(const ISurface* surf, double& s)
    //~ {
    //~ // the normal directions in x and y
    //~ const double __nx = surf->normal().x();
    //~ const double __ny = surf->normal().y();
//~
    //~ return false;
    //~ }
    //~
    //~
    //~
    //~ bool trajectory::_intersectWithinZDiskBounds(const ISurface* surf, double& s)
    //~ {
    //~ // the z position of the plane
    //~ double __planePositionZ = surf->origin().z();
//~
    //~ return false;
    //~ }
}
