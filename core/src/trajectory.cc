#include "trajectory.hh"

#include <stdexcept>
#include <utility>
#include <iostream>
#include <algorithm>

#include "intersections.hh"
#include "utilities.hh"
#include "fiveByFiveMatrix.hh"

#ifdef AIDATT_USE_DD4HEP
#include "DD4hep/DD4hepUnits.h"
#endif // AIDATT_USE_DD4HEP

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
        /// 3. choose the one with the smaller s (!), if needed
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



    bool trajectory::_calculateIntersectionWithSurface(const ISurface* surf, double& s, Vector2D* localUV, Vector3D* xx)
    {
        /// currently three different types of surfaces are available
        if(surf->type().isZCylinder())
            return _intersectsWithinZCylinderBounds(surf, s, localUV, xx);
        else if(surf->type().isZPlane())
            return _intersectWithinZPlaneBounds(surf, s, localUV, xx);
        else if(surf->type().isZDisk())
            return _intersectWithinZDiskBounds(surf, s, localUV, xx);
        else
            throw std::invalid_argument("[aidaTT::trajectory::getIntersectionWithSurfaces] Unknown surface type!");
    }



    bool trajectory::_intersectsWithinZCylinderBounds(const ISurface* surf, double& s, Vector2D* localUV, Vector3D* xx)
    {
      //// see: L3 internal note 1666 "Helicoidal tracks", J.Alcaraz

      if( ! surf->type().isParallelToZ() ){

	throw std::runtime_error( " **** _intersectsWithinZCylinderBounds: only works for xylinders parallel to the z-axis " ) ;
      }

      const double omega = calculateOmega( _referenceParameters);
      const double phi0  = calculatePhi0(  _referenceParameters);

      const Vector3D& rp =  _referenceParameters.referencePoint() ;
      const double x0    = rp.x()  ;
      const double y0    = rp.y()  ;
      
      const double rho = dynamic_cast<const ICylinder*>(surf)->radius() ;

      const Vector3D&  cylc = surf->origin() - rho * surf->normal() ;
      const double xrho = cylc.x()  ;
      const double yrho = cylc.y()  ;

      //--------

      const double sinph = sin( phi0 ) ;
      const double cosph = cos( phi0 ) ;

      const double dx = xrho - x0 ;
      const double dy = yrho - y0 ;

      double sox = sinph - omega * dx  ;
      double coy = cosph + omega * dy ;
      
      double gamma = ( 2*dx * sinph - 2 * dy * cosph - omega * rho * rho  - omega * ( dx*dx + dy*dy ) ) ;
      gamma /= ( 2 * rho * sqrt( sox * sox + coy * coy) ) ;

      if( std::fabs( gamma ) > 1. )  // no solution  ( could have faster check at beginning  ...  ) 
	return false ;

      const double phirho = atan2( sox , coy ) ;

      const double asing = asin( gamma ) ;

      const double phic0 = asing + phirho ;

      double phic1 = ( asing  > 0. ?  M_PI - asing :  - M_PI - asing  ) ;
      phic1 += phirho ;

      const double X0 = xrho + rho * cos( phic0 )  ;
      const double Y0 = yrho + rho * sin( phic0 )  ;
      const double X1 = xrho + rho * cos( phic1 )  ;
      const double Y1 = yrho + rho * sin( phic1 )  ;
      
      const double S0 = calculateSfromXY( std::make_pair( X0 , Y0 ) , _referenceParameters);
      const double S1 = calculateSfromXY( std::make_pair( X1 , Y1 ) , _referenceParameters);

      const double Z0 = calculateZfromS(S0, _referenceParameters);
      const double Z1 = calculateZfromS(S1, _referenceParameters);

      Vector3D sol0(X0, Y0, Z0);
      Vector3D sol1(X1, Y1, Z1);

      const bool insideFirst  = surf->insideBounds(sol0);
      const bool insideSecond = surf->insideBounds(sol1);


      if((!insideFirst && !insideSecond))   // || (S0 < 0. && S1 < 0.))      //do not  discard negative or no solution
	return false;

      else if(insideFirst && !insideSecond)
	{
	  s = S0;
	  if(localUV != NULL)
	    _calculateLocalCoordinates(surf, sol0, localUV);

	  if(xx) xx->fill(sol0) ;

	  //        return true;
	}
      else if(!insideFirst && insideSecond)
	{
	  s = S1;
	  if(localUV != NULL)
	    _calculateLocalCoordinates(surf, sol1, localUV);

	  if(xx) xx->fill(sol1) ;

	  // return true;
	}
      else // both are valid , choose the smaller absolute solution
	{
	  if(std::fabs(S0)  < std::fabs(S1))
	    {
	      s = S0;
	      if(localUV != NULL)
		_calculateLocalCoordinates(surf, sol0, localUV);

	      if(xx) xx->fill(sol0) ;
	    }
	  else
	    {
	      s = S1;
	      if(localUV != NULL)
		_calculateLocalCoordinates(surf, sol1, localUV);

	      if(xx) xx->fill(sol1) ;
	    }
	}

      return true;


    }



    bool trajectory::_intersectWithinZPlaneBounds(const ISurface* surf, double& s, Vector2D* localUV, Vector3D* xx)
    {
        // the straight line: normals plus distance; distance must be positive !
        const double nx = surf->normal().x();
        const double ny = surf->normal().y();

        // calculate distance from origin
        const double dist = fabs(surf->distance(Vector3D()));

        // define straight line from this
        straightLine line(nx, ny, dist);
        // create circle
        const double radius  = calculateRadius(_referenceParameters);
        const double xcenter = calculateXCenter(_referenceParameters);
        const double ycenter = calculateYCenter(_referenceParameters);
        circle circ(xcenter, ycenter, radius);

        //  std::cout << " ++ trajectory::_intersectWithinZPlaneBounds: circle center:  ( " << xcenter << ", " << ycenter << ")" << std::endl ;

        intersections candidates = intersectCircleStraightLine(circ, line);

        if(candidates.number() < 1)
            return false;
        else if(candidates.number() == 1)
            {
                const double S = calculateSfromXY(candidates[0], _referenceParameters);
                const double Z = calculateZfromS(S, _referenceParameters);
                Vector3D thePlace(candidates[0].first, candidates[0].second, Z);
                bool inside = surf->insideBounds(thePlace);

                if(inside)
                    {
                        s = S;
                        if(localUV != NULL)
                            _calculateLocalCoordinates(surf, thePlace, localUV);

                        if(xx) xx->fill(thePlace) ;

                        return true;
                    }
                else
                    return false;
            }

        ///  else -- the standard case: two solutions index 0 and 1
        /// calculate all values first, then evaluate
        const double X0 = candidates[0].first;
        const double Y0 = candidates[0].second;
        const double S0 = calculateSfromXY(X0, Y0, _referenceParameters);
        const double Z0 = calculateZfromS(S0, _referenceParameters);

        const double X1 = candidates[1].first;
        const double Y1 = candidates[1].second;
        const double S1 = calculateSfromXY(X1, Y1, _referenceParameters);
        const double Z1 = calculateZfromS(S1, _referenceParameters);

        Vector3D sol0(X0, Y0, Z0);
        Vector3D sol1(X1, Y1, Z1);

        const bool insideFirst  = surf->insideBounds(sol0);
        const bool insideSecond = surf->insideBounds(sol1);

        if((!insideFirst && !insideSecond))   // || (S0 < 0. && S1 < 0.))      //do not  discard negative or no solution
            return false;

        else if(insideFirst && !insideSecond)
            {
                s = S0;
                if(localUV != NULL)
                    _calculateLocalCoordinates(surf, sol0, localUV);

                if(xx) xx->fill(sol0) ;

                //        return true;
            }
        else if(!insideFirst && insideSecond)
            {
                s = S1;
                if(localUV != NULL)
                    _calculateLocalCoordinates(surf, sol1, localUV);

                if(xx) xx->fill(sol1) ;

                // return true;
            }
        else // both are valid , choose the smaller absolute solution
            {
                if(std::fabs(S0)  < std::fabs(S1))
                    {
                        s = S0;
                        if(localUV != NULL)
                            _calculateLocalCoordinates(surf, sol0, localUV);

                        if(xx) xx->fill(sol0) ;
                    }
                else
                    {
                        s = S1;
                        if(localUV != NULL)
                            _calculateLocalCoordinates(surf, sol1, localUV);

                        if(xx) xx->fill(sol1) ;
                    }
            }

        return true;

    }



    bool trajectory::_intersectWithinZDiskBounds(const ISurface* surf, double& s, Vector2D* localUV, Vector3D* xx)
    {
        // the z position of the plane
        double planePositionZ = surf->origin().z();
	double helixPositionZ = _referenceParameters.referencePoint().z() + calculateZ0( _referenceParameters ) ;

        s = ( planePositionZ - helixPositionZ ) / calculateTanLambda(_referenceParameters);

        double x = calculateXfromS(s, _referenceParameters);
        double y = calculateYfromS(s, _referenceParameters);

        Vector3D thePlace = Vector3D(x, y, planePositionZ);

        if(surf->insideBounds(thePlace))
            {
                if(localUV != NULL)
                    _calculateLocalCoordinates(surf, thePlace, localUV);

                if(xx) xx->fill(thePlace) ;

                return true;
            }
        else
            return false;
    }



    void trajectory::_calculateLocalCoordinates(const ISurface* surf, const Vector3D& position, Vector2D* localUV, Vector3D* xx)
    {
        Vector2D local = surf->globalToLocal(position);
        localUV->_u = local.u();
        localUV->_v = local.v();
    }



    void trajectory::addMeasurement(const Vector3D& position, const std::vector<double>& precision, const ISurface& surface, void* id)
    {
        /// get reference information
        Vector2D* referenceUV = new Vector2D();
        double s =  0;
        double intersects = _calculateIntersectionWithSurface(&surface, s, referenceUV);

        /// calculate measurement info
        Vector2D* measuredUV = new Vector2D(surface.globalToLocal(position));

        // combining both delivers the actual residuals: measurement MINUS reference
        const double udiff = measuredUV->u() - referenceUV->u();
        const double vdiff = measuredUV->v() - referenceUV->v();

        std::vector<double> residuals(2);
        residuals[0] = udiff;
        residuals[1] = vdiff;

        std::vector<Vector3D>* measDir = new std::vector<Vector3D>;
        measDir->push_back(surface.u(position));
        measDir->push_back(surface.v(position));

        if(intersects)
            {
                _initialTrajectoryElements.push_back(new trajectoryElement(s, surface, measDir, precision, residuals, calculateLocalCurvilinearSystem(s, _referenceParameters), id));
            }
        else
            {
                std::cout << " ERROR: hit at " << position << "  does not intersect with surface : " <<  surface
                          << "        hit will be ignored ! " << std::endl ;
            }
    }



    void trajectory::addElement(const Vector3D& point, void* id)
    {
        double s =  calculateSfromXY(point.x(), point.y(), _referenceParameters);
        _initialTrajectoryElements.push_back(new trajectoryElement(s, id));
    }



    void trajectory::addElement(const Vector3D& point, const ISurface& surface, void* id)
    {
        /// TODO :: this is not complete -- what is it needed for?
        double s =  calculateSfromXY(point.x(), point.y(), _referenceParameters);
        // suppress warning -- STILL WRONG:: TODO
        _initialTrajectoryElements.push_back(new trajectoryElement(s, id));
    }



    bool compareTrajectoryElements(trajectoryElement* one, trajectoryElement* two)
    {
        return (one->arcLength() < two->arcLength());
    }



    void trajectory::prepareForFitting()
    {
        ///~ first sort the trajectory elements by arclength
        sort(_initialTrajectoryElements.begin(), _initialTrajectoryElements.end(), compareTrajectoryElements);

        const double cosLambda = cos(calculateLambda(_referenceParameters));

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
        /// now the really interesting ones
        for(std::vector<trajectoryElement*>::iterator element = _initialTrajectoryElements.begin() + 1, last = _initialTrajectoryElements.end(); element < last; ++element)
            {
                ///~ obtain the two arc lengths
                double prevS = (*(element - 1))->arcLength();
                double currS = (*element)->arcLength();
                ///~ calculate 3D arclength -> divide 2D arc length by cos lambda
                double dw = (currS - prevS) / cosLambda;

                Vector3D tstart = calculateTangent(prevS, _referenceParameters);
                Vector3D tend   = calculateTangent(currS, _referenceParameters);

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
