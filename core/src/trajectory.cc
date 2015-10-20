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
    for(std::vector<trajectoryElement*>::iterator element = _initialTrajectoryElements.begin(), last = _initialTrajectoryElements.end(); element < last; ++element){
      //std::cout << " Now I register element with arc length " << (*element)->arcLength() << " for deletion " << std::endl ;
      delete *element;
    }
  }


  Vector3D trajectory::pointAt( double s) {
    return pointOnTrajectory( _referenceParameters, s ) ;
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
	double s = 0.;  
	
	bool intersects = _calculateIntersectionWithSurface(*surf, s );
	
	if(intersects)
	  _intersectionsList.push_back(std::make_pair(s, (*surf)));
      }
    
    return _intersectionsList;
  }
  
  
  
  const fitResults* trajectory::getFitResults(int label)
  {
    return _fittingAlgorithm->getResults(label);
  }
  
  
  
  bool trajectory::_calculateIntersectionWithSurface(const ISurface* surf, double& s, Vector2D* localUV, Vector3D* xx)
  {

    // check if we have a chached intersection
    SIMap::iterator it = _surfIntersections.find( surf ) ;
    if( it != _surfIntersections.end() ) {
      s = it->second.S ;
      if(localUV) *localUV = it->second.UV ;
      if(xx) *xx = it->second.XX ;
      return true ;
    }

    bool intersects = false ;
    
    static Vector2D _localUV ;
    static Vector3D _xx ;
    
    if( localUV == 0 ) localUV = &_localUV ;
    if( xx == 0 ) xx = &_xx ;
    
    /// currently three different types of surfaces are available
    if(surf->type().isZCylinder())
      intersects = _intersectsWithinZCylinderBounds(surf, s, localUV, xx);
    else if(surf->type().isZPlane())
      intersects = _intersectWithinZPlaneBounds(surf, s, localUV, xx);
    else if(surf->type().isZDisk())
      intersects = _intersectWithinZDiskBounds(surf, s, localUV, xx);
    

    if( intersects ){ // cache result
      _surfIntersections.insert( std::make_pair(  surf, SurfIntersection( surf, s, localUV, xx ) ) ) ;
    }

    return intersects ;
  }



  bool trajectory::_intersectsWithinZCylinderBounds(const ISurface* surf, double& s, Vector2D* localUV, Vector3D* xx)
  {
    //// see: L3 internal note 1666 "Helicoidal tracks", J.Alcaraz

    if( ! surf->type().isParallelToZ() ){

      throw std::runtime_error( " **** _intersectsWithinZCylinderBounds: only works for xylinders parallel to the z-axis " ) ;
    }

    const double omega = calculateOmega( _referenceParameters);
    const double phi0  = calculatePhi0(  _referenceParameters);
    const double d0    = calculateD0(    _referenceParameters);

    const double sinph = sin( phi0 ) ;
    const double cosph = cos( phi0 ) ;

    const Vector3D& rp =  _referenceParameters.referencePoint() ;
    const double x0    = rp.x() - d0 * sinph ;
    const double y0    = rp.y() + d0 * cosph ;
      
    const double rho = dynamic_cast<const ICylinder*>(surf)->radius() ;

    const Vector3D&  cylc = dynamic_cast<const ICylinder*>(surf)->center() ;
    const double xrho = cylc.x()  ;
    const double yrho = cylc.y()  ;

    //--------


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


    // std::cout << " trajectory::_intersectsWithinZCylinderBounds:  sol0 x = " 
    // 		<< X0 << " y = " << Y0 << " z = " << Z0 << " s = " << S0 << " isInside " << insideFirst << std::endl ;

    // std::cout << " trajectory::_intersectsWithinZCylinderBounds:  sol1 x = " 
    // 		<< X1 << " y = " << Y1 << " z = " << Z1 << " s = " << S1 << " isInside " << insideSecond << std::endl ;


    if((!insideFirst && !insideSecond))   // || (S0 < 0. && S1 < 0.))      //do not  discard negative or no solution
      return false;


#if 0 // fg: change the logic for cylinders: always only return the closest solution 
      // to have same behaviour as KalTest
      // FIXME: ideally the interscetion methods should return all solutions and have the 
      //        calling program decide which to use ...

    Vector3D& theSol = sol0 ;
    bool isInside = insideFirst ;

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

	theSol = sol1 ;
	isInside = insideSecond ;
      }
    
    return isInside ;
 
#else

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
#endif



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

    //	double epsilon = 1.e-3  ;

    const bool insideFirst  = surf->insideBounds(sol0 ) ; // , epsilon );
    const bool insideSecond = surf->insideBounds(sol1 ) ; // , epsilon );

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
    localUV->u() = local.u();
    localUV->v() = local.v();
  }


  double trajectory::computeQMS( const ISurface* surface, const trackParameters* tp){
    
    double s ; Vector2D uv ; Vector3D position ;
    double intersects = _calculateIntersectionWithSurface( surface, s , &uv, &position );
    
    if( ! intersects ) 
      return 0 ;

    aidaTT::trackParameters test_tp = ( tp ? *tp : _referenceParameters ) ;
    
    moveHelixTo( test_tp, position ) ; // move helix to the scatterer
    
    Vector5 hel = test_tp.parameters();
    
    double omega  = hel(OMEGA);
    double tnl    = hel(TANL); 
    double phi0   = hel(PHI0);
    
    
    const DDSurfaces::IMaterial& material_inn = surface->innerMaterial();
    const DDSurfaces::IMaterial& material_out = surface->outerMaterial();
    
    const double r_i = surface->innerThickness();
    const double r_o = surface->outerThickness();
    
    const double X0_o = material_out.radiationLength();
    const double X0_i = material_inn.radiationLength();
    
    double r_tot = r_i + r_o ;
    
    //calculation of effective radiation length of the surface
    //double X0_eff = 1. / ( (r_i/r_tot)/ X0_i  +  (r_o/r_tot)/ X0_o ) ;
    double X0_eff = ( r_i/X0_i + r_o/X0_o ) / ( r_i + r_o ) ; 
    
    //calculation of the path of the particle inside the material
    //compute path as projection of (straight) track to surface normal:
    DDSurfaces::Vector3D p( - std::sin( phi0 ), std::cos( phi0 ) , tnl ) ;
    DDSurfaces::Vector3D up = p.unit() ;
    
    // need to get the normal at crossing point 
    const DDSurfaces::Vector3D& n = surface->normal( position ) ;
    
    
    double cosTrk = std::fabs( up * n )  ;
    
    double path = r_i + r_o ;
    
    //note: projectedPath is already in dd4hep(TGeo) units, i.e. cm !
    path = path/cosTrk ; 
    
    double X_X0 = path * X0_eff ;
    
    double Pt = ( fabs(1./omega ) / 100.0 ) ;  // That's Pt
    double mom = Pt*TMath::Sqrt(1 + tnl*tnl);
    
    static const double mass = 0.13957018; // pion mass [GeV]
    double   beta = mom / TMath::Sqrt(mom * mom + mass * mass);
    
    double Qms = 0.0136/(mom*beta) * 1.0 * TMath::Sqrt(X_X0) * (1 + 0.0038*(TMath::Log(X_X0)));

    return Qms ;
  }


  void trajectory::addScatterer(const Vector3D& position, std::vector<double>& precision, const ISurface& surface, const trackParameters& seed_tp, void* id)
  //void trajectory::addScatterer(const Vector3D& position, TMatrixDSym& precision, const ISurface& surface,  const trackParameters& seed_tp, void* id)
  {
   
    /// get reference information
    // I am not sure what to do with these stuff
    Vector2D referenceUV ;
    double s =  0;
    double intersects = _calculateIntersectionWithSurface(&surface, s, &referenceUV);

    std::vector<Vector3D>* measDir = new std::vector<Vector3D>;
    measDir->push_back(surface.u(position));
    measDir->push_back(surface.v(position));

    std::vector<double> residuals(2);
    residuals[0] = 0;
    residuals[1] = 0;

    // I mean the stuff up here

    aidaTT::trackParameters test_tp = seed_tp ;
    /*
      Vector5 hel_IP = seed_tp.parameters();
      double omega_IP  = hel_IP(0);
      double tnl_IP    = hel_IP(1); 
      double phi0_IP   = hel_IP(2);

      std::cout << " track parameters in IP omega: " << omega_IP << " tanl " << tnl_IP << " AND phi0 " << phi0_IP << std::endl ; 
    */
    moveHelixTo( test_tp, position ) ; // move helix to the scatterer

    Vector5 hel = test_tp.parameters();

    double omega  = hel(0);
    double tnl    = hel(1); 
    double phi0   = hel(2);


    const DDSurfaces::IMaterial& material_inn = surface.innerMaterial();
    const DDSurfaces::IMaterial& material_out = surface.outerMaterial();

    const double r_i = surface.innerThickness();
    const double r_o = surface.outerThickness();

    const double X0_o = material_out.radiationLength();
    const double X0_i = material_inn.radiationLength();

    double r_tot = r_i + r_o ;

    //calculation of effective radiation length of the surface
    //double X0_eff = 1. / ( (r_i/r_tot)/ X0_i  +  (r_o/r_tot)/ X0_o ) ;
    double X0_eff = ( r_i/X0_i + r_o/X0_o ) / ( r_i + r_o ) ; 

    //calculation of the path of the particle inside the material
    //compute path as projection of (straight) track to surface normal:
    DDSurfaces::Vector3D p( - std::sin( phi0 ), std::cos( phi0 ) , tnl ) ;
    DDSurfaces::Vector3D up = p.unit() ;

    //const DDSurfaces::Vector3D& n = surface.normal() ;
    // need to get the normal at crossing point ( should be the current helix' reference point)
    const Vector3D& piv = seed_tp.referencePoint() ;
    DDSurfaces::Vector3D xx( piv.x(),piv.y(),piv.z()) ;
    const DDSurfaces::Vector3D& n = surface.normal( xx ) ;
    

    double cosTrk = std::fabs( up * n )  ;
    
    double path = r_i + r_o ;

    //note: projectedPath is already in dd4hep(TGeo) units, i.e. cm !
    path = path/cosTrk ; 

    double X_X0 = path * X0_eff ;

    double Pt = ( fabs(1./omega ) / 100.0 ) ;  // That's Pt
    double mom = Pt*TMath::Sqrt(1 + tnl*tnl);

    static const double mass = 0.13957018; // pion mass [GeV]
    double   beta = mom / TMath::Sqrt(mom * mom + mass * mass);

    double Qms = 0.0136/(mom*beta) * 1.0 * TMath::Sqrt(X_X0) * (1 + 0.0038*(TMath::Log(X_X0)));

    // std::cout << " omega par " << omega << " mom " << mom << " beta " << beta << "rinn, rout " << r_i << ", " << r_o << " X0inn, X0out " <<  X0_i << ", " << X0_o <<  " effective radiation length " << X0_eff <<  " x/X0 " << X_X0 << " path " << path << "Cosine of track angle with the surface " << cosTrk << " Qms = " << Qms << std::endl;

    precision[0] = Qms*Qms ;  precision[1] = Qms*Qms ;
    
    if(intersects)
      {
	_initialTrajectoryElements.push_back(new trajectoryElement( s , surface, measDir, precision, residuals, calculateLocalCurvilinearSystem(s, _referenceParameters), id, true ));
      } 
    else
      {
	delete measDir ;
	
	std::cout << " ERROR: hit at " << position << "  does not intersect with surface : " <<  surface
		  << "        hit will be ignored ! " << std::endl ;
      }
    
  }
  



  void trajectory::addMeasurement(const Vector3D& position, const std::vector<double>& precision, const ISurface& surface, void* id, bool isScatterer)
  {

    /// get reference information
    Vector2D referenceUV ;
    double s =  0;
    double intersects = _calculateIntersectionWithSurface(&surface, s, &referenceUV);

    if( !intersects ){
      
      std::cout << " ERROR: trajectory::addMeasurement() hit at " << position << "  does not intersect with surface : " <<  surface
		<< "        hit will be ignored ! " << std::endl ;
      
      return ;
    }
    
    /// calculate measurement info
    Vector2D measuredUV(surface.globalToLocal(position));

    // combining both delivers the actual residuals: measurement MINUS reference
    const double udiff = measuredUV.u() - referenceUV.u();
    const double vdiff = measuredUV.v() - referenceUV.v();

    std::vector<double> residuals(2);
    residuals[0] = udiff;
    residuals[1] = vdiff;

    std::vector<Vector3D>* measDir = new std::vector<Vector3D>;
    measDir->push_back(surface.u(position));
    measDir->push_back(surface.v(position));


    std::vector<double> new_prec = precision ;

    if( isScatterer ){ // also  add a scattering to the trajectory element
      
      double qms = computeQMS( &surface ) ;
      
      new_prec.push_back(  qms*qms  ) ;
    }
    
    _initialTrajectoryElements.push_back(new trajectoryElement(s, surface, measDir, new_prec, residuals, calculateLocalCurvilinearSystem(s, _referenceParameters), id , isScatterer ));
  }

  void trajectory::addScatterer( const ISurface& surface ){
    
    double s ; Vector2D uv ; Vector3D position ;
    double intersects = _calculateIntersectionWithSurface( &surface, s , &uv, &position );
    
    if( !intersects ){
      std::cout << " ERROR: trajectory::addScatterer: track does not intersect with surface : " <<  surface << std::endl ;      
      return ;
    }

    std::vector<double> residuals(2);
    residuals[0] = 0 ;
    residuals[1] = 0 ;

    std::vector<Vector3D>* measDir = new std::vector<Vector3D>;
    measDir->push_back(surface.u(position));
    measDir->push_back(surface.v(position));

    double qms = computeQMS( &surface ) ;
    
    std::vector<double> precision(1) ;
    precision[0] =  qms*qms  ;
    
    _initialTrajectoryElements.push_back(new trajectoryElement(s, surface, measDir, precision, residuals, calculateLocalCurvilinearSystem(s, _referenceParameters), 0 , true ));
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

    double prevS = 0. ;

    /// the first jacobian is useless, just use an empty 5x5 matrix
    if(_initialTrajectoryElements.size() > 0)
      {
	fiveByFiveMatrix* j = new fiveByFiveMatrix;
	j->Unit();
	(_initialTrajectoryElements.at(0))->setJacobian(j);
      }
    /// now the really interesting ones
    for(std::vector<trajectoryElement*>::iterator element = _initialTrajectoryElements.begin()+1, last = _initialTrajectoryElements.end(); element < last; ++element)
      {

	///~ obtain the two arc lengths
	//double prevS = (*(element - 1))->arcLength();
	double currS = (*element)->arcLength();
	///~ calculate 3D arclength -> divide 2D arc length by cos lambda
	double dw = (currS - prevS) / cosLambda;

	// std::cout << std::endl ;
	// std::cout << " current S " << currS << " previous S " << prevS << " currS - prevS " << dw << std::endl ;
	// std::cout << std::endl ;

	Vector3D tstart = calculateTangent(prevS, _referenceParameters);
	Vector3D tend   = calculateTangent(currS, _referenceParameters);

	fiveByFiveMatrix* jacob = new fiveByFiveMatrix;
	_propagation->getJacobian(*jacob, dw, qbyp, tstart, tend, BField);
	(*element)->setJacobian(jacob);

	prevS = currS ;

      }
  }


  bool trajectory::fit()
  {
    return _fittingAlgorithm->fit(*this);
  }
}
