#include "trajectory.hh"

#include <stdexcept>
#include <utility>
#include <iostream>
#include <algorithm>

#include "intersections.hh"
#include "utilities.hh"
#include "aidaTT-Units.hh"
#include "materialUtils.hh"


namespace aidaTT {
  
    

  trajectory::trajectory(const trackParameters& tp, IFittingAlgorithm* fa, 
			 IPropagation* pm, const IGeometry* geom) :
    _referenceParameters(tp) , _fittingAlgorithm(fa) ,  _propagation(pm), _geometry(geom), _mass( pionMass )  {
  }
  


  trajectory::trajectory(const trackParameters& tp, const IGeometry* geom) : 
    _referenceParameters(tp), _fittingAlgorithm(NULL), _propagation(NULL), _geometry(geom) , _mass( pionMass ) {
  }



  trajectory::trajectory(const trajectory& traj) : _referenceParameters(traj._referenceParameters),
						   _fittingAlgorithm(traj._fittingAlgorithm), 
						   _propagation(traj._propagation), _geometry(traj._geometry), _mass( pionMass ) {
  }



  trajectory::~trajectory() {
    
    for(std::vector<trajectoryElement*>::iterator element = _initialTrajectoryElements.begin(), 
	  last = _initialTrajectoryElements.end(); element < last; ++element){
      delete *element;
    }
  }


  const trajectoryElement* trajectory::trajectoryElementAt(double s ) {
    
    double prevS = -1.e99 ;
    const trajectoryElement* prevElement = ( _initialTrajectoryElements.size() > 0 ? 
					     _initialTrajectoryElements[0]  : 0 ) ;
    
    for(std::vector<trajectoryElement*>::iterator it = _initialTrajectoryElements.begin()+1, 
	  last = _initialTrajectoryElements.end(); it < last; ++it ){
      
      const trajectoryElement* element = *it ;
      
      double newS = element->arcLength() ;
      
      if( prevS <= s &&  s < newS ){
	break ;
      }
      prevElement = element ;
      prevS = newS ;
    }
    
    return prevElement ;
  }
  
  
  Vector3D trajectory::pointAt( double s ) {
    
    const trajectoryElement* element =  trajectoryElementAt( s ) ;
    
    if( element != 0 ) 
      return aidaTT::pointAt(  ( s - element->arcLength() ) , *element->getTrackParameters() ) ;
    
    return  aidaTT::pointAt( s, _referenceParameters ) ;
  }
  
  
  
  /// comput the tangent to the trajectory at given s
  Vector3D trajectory::tangentAt( double s ) {
    
    const trajectoryElement* element =  trajectoryElementAt( s ) ;
    
    if( element != 0 ) 
      return aidaTT::calculateTangent( ( s - element->arcLength() ),  *element->getTrackParameters() ) ;
    
    return  aidaTT::calculateTangent( s, _referenceParameters ) ;
  }




  struct SortWithS{

    bool operator()(const std::pair<double, const ISurface*>& i0,
		    const std::pair<double, const ISurface*>& i1  ){
      return i0.first < i1.first ;
    }

  };



  const IntersectionVec& trajectory::getIntersectionsWithSurfaces(const SurfaceVec& surfaces)
  {

    double maxS = M_PI * std::fabs( calculateRadius(_referenceParameters) ) ;

    for(std::vector<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; 
	surf != surfaces.end() ; ++surf) {
      
      double s = 0.;  
      
      Vector3D xx ;
      bool intersects = aidaTT::intersectWithSurface( *surf, _referenceParameters.parameters() ,
						      _referenceParameters.referencePoint(), s , xx , +1, true );
      
      // only keep intersections at positve s in the first half arc
      if( intersects && s >= 0. && s < maxS )
	_intersectionsList.push_back(std::make_pair(s, (*surf)));
    }
    

    std::sort( _intersectionsList.begin() , _intersectionsList.end() , SortWithS() ) ;

    return _intersectionsList;
  }
  
  
  const fitResults* trajectory::getFitResults(int label){

    return _fittingAlgorithm->getResults(label);
  }
  
  
  



  void trajectory::addMeasurement(const Vector3D& position, const std::vector<double>& precision, 
				  const ISurface& surface, void* id, bool isScatterer)
  {
    
    if( trajectoryElements().empty() ){ // add an initial trajectoryElement

      addElement( Vector3D(), 0 ) ;
    }


    const trackParameters* tP = trajectoryElements().back()->getTrackParameters() ;
    const double prevS        = trajectoryElements().back()->arcLength() ;

    double s = 0;
    Vector3D xx ;
    bool intersects = aidaTT::intersectWithSurface( &surface, tP->parameters() , tP->referencePoint() ,
						   s, xx, +1 , false ) ; 
    //note: if we have a measurement we do _not_ check for the bounds of the surface

    if( !intersects ){
      
      std::cout << "  ERROR: trajectory::addMeasurement() hit at " << position 
		<< "  does not intersect with surface : " <<  surface
		<< "  hit will be ignored ! " << std::endl ;
      
      return ;
    }
    
    //    std::cout << " prevS : " << prevS << " - s to next intersection: " << s << std::endl ;

    // add this trajectoryElement at the total path lengths from the IP
    s = s + prevS ;

    const Vector2D& referenceUV = surface.globalToLocal( xx ) ;

    trackParameters* trkParam = new  trackParameters( *tP ) ;
    
    // move the track paramters to the intersection point
    moveHelixTo( *trkParam, xx ) ;


    const Vector3D& mom = momentumAtPCA( *trkParam ) ;

    const Vector2D& measuredUV = surface.globalToLocal( position ) ;

    // ********************************************
    // apply the energy loss from this surface
    double energy, beta ;
    double deltaE = aidaTT::computeEnergyLoss( &surface, measuredUV , mom , energy, beta , _mass ) ; 
    trkParam->parameters()( OMEGA ) /= ( 1. - deltaE/energy ) ;
    // ********************************************


    /// calculate measurement info

    std::vector<double> residuals;
    std::vector<Vector3D>* measDir = new std::vector<Vector3D>;
    
    const double udiff = measuredUV.u() - referenceUV.u();
    residuals.push_back( udiff ) ;
    measDir->push_back( surface.u(position) );
    

    //    if( ! surface.type().isMeasurement1D()  ){
    
    const double vdiff = measuredUV.v() - referenceUV.v();
    residuals.push_back( vdiff ) ;
    measDir->push_back(surface.v(position));
    // }
    
    std::vector<double> new_prec = precision ;
    
    if( isScatterer ){ // also  add a scattering to the trajectory element
      
      
      double qms = aidaTT::computeQMS( &surface, referenceUV , mom , _mass ) ; 
      
      // projection of track onto surface vectors
      // const Vector3D& up = mom.unit() ;
      // double c1 = up * surface.u( xx );
      // double c2 = up * surface.v( xx );
      
      new_prec.push_back(  qms*qms  ) ;
      // new_prec.push_back(  c1  ) ;
      // new_prec.push_back(  c2  ) ;
      //fg: c1,c2 are scalar products of offset directions with track direction
      //    and are by construction 0. in curvilinear !
      new_prec.push_back( 0. ) ;
      new_prec.push_back( 0. ) ;
    }


    // note: need to get the curvilinear system at s==0. as this is where the local track state is defined
    _initialTrajectoryElements.push_back(new trajectoryElement(s, trkParam, surface, measDir, new_prec, residuals, calculateLocalCurvilinearSystem(0., *trkParam), id , isScatterer ));
  }

  void trajectory::addScatterer( const ISurface& surface ){
    
    if( trajectoryElements().empty() ){ // add an initial trajectoryElement

      addElement( Vector3D(), 0 ) ;
    }

    const trackParameters* tP = trajectoryElements().back()->getTrackParameters() ;
    const double prevS = trajectoryElements().back()->arcLength() ;

    double s =  0;
    Vector3D xx ;
    bool intersects = aidaTT::intersectWithSurface( &surface, tP->parameters() , tP->referencePoint() ,
						   s, xx, +1 , true ) ; 
    //note: if we have a scatterer we _do_ check for the bounds of the surface

    if( !intersects ){
      
      std::cout << "  INFO: trajectory::addScatterer() track "
		<< "  does not intersect with surface : " <<  surface
		<< std::endl ;
      return ;
    }
    
    // add this trajectoryElement at the total path lengths from the IP
    s = s + prevS ;

    trackParameters* trkParam = new  trackParameters( *tP ) ;
    
    // move the track paramters to the intersection point
    moveHelixTo( *trkParam, xx ) ;
    
    Vector2D referenceUV = surface.globalToLocal( xx ) ;

    const Vector3D& mom = momentumAtPCA( *trkParam ) ;

    // ********************************************
    // apply the energy loss from this surface
    double energy, beta ;
    double deltaE = aidaTT::computeEnergyLoss( &surface, referenceUV , mom , energy, beta , _mass ) ; 
    trkParam->parameters()( OMEGA ) /= ( 1. - deltaE/energy ) ;
    // ********************************************

    std::vector<double> residuals(2);
    residuals[0] = 0 ;
    residuals[1] = 0 ;

    std::vector<Vector3D>* measDir = new std::vector<Vector3D>;
    measDir->push_back( surface.u( xx ) );
    measDir->push_back( surface.v( xx ) );

    double qms = aidaTT::computeQMS( &surface, referenceUV , mom , _mass ) ; 
      
    // projection of track onto surface vectors
    // const Vector3D& up = mom.unit() ;
    // double c1 = up * (*measDir)[0] ; //surface.u( xx );
    // double c2 = up * (*measDir)[1] ; //surface.v( xx );

    //fg: c1,c2 are scalar products of offset directions with track direction
    //    and are by construction 0. in curvilinear !
    std::vector<double> precision(3) ;
    precision[0] =  qms*qms  ;
    precision[1] = 0. ;
    precision[2] = 0. ;
   

    // note: need to get the curvilinear system at s==0. as this is where the local track state is defined
    _initialTrajectoryElements.push_back(new trajectoryElement(s, trkParam, surface, measDir, precision, residuals, calculateLocalCurvilinearSystem(0., *trkParam), 0 , true, false ));
  }


  void trajectory::addElement(const Vector3D& point, void* id)
  {
    double s =  calculateSfromXY(point.x(), point.y(), _referenceParameters);

    //FIXME: need proper track parameters at this s ....
    trackParameters* trkParam = new  trackParameters( _referenceParameters ) ;

    _initialTrajectoryElements.push_back(new trajectoryElement(s, trkParam ,  id));
  }



  bool compareTrajectoryElements(trajectoryElement* one, trajectoryElement* two)
  {
    return (one->arcLength() < two->arcLength());
  }




  void trajectory::prepareForFitting()
  {
    ///~ first sort the trajectory elements by arclength
    sort(_initialTrajectoryElements.begin(), _initialTrajectoryElements.end(), 
	 compareTrajectoryElements);


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

	const trackParameters& trkParam = *(*element)->getTrackParameters() ;

	const double cosLambda = cos( calculateLambda( trkParam  ) );
	const Vector3D& BField = _geometry->getBField( trkParam.referencePoint() ) ;

	const double qbyp  = calculateQoverP( trkParam , BField.z() ) ;

	double currS = (*element)->arcLength();

	///~ calculate 3D arclength 
	double dw = (currS - prevS) / cosLambda;

	// Calculate the energy loss only if the trajectory element is a scatterer, meaning it has material
	double NrjLoss = 0.;
	if ((*element)->isScatterer()){
	  //calculate the energy loss
	  const aidaTT::ISurface& surf = (*element)->surface();
	  //	  NrjLoss = GetEnergyLoss(  &surf, &trkParam );

	  // check new eloss:
	  double e,b ;
	  double deltaE = aidaTT::computeEnergyLoss( &surf, trkParam , e, b, _mass ) ;
	  NrjLoss = (2.0*deltaE) / ((b*b)*e);

	  //	  std::cout << " NrjLoss : " << NrjLoss << " - NrjLossNew: " << NrjLossNew << std::endl ;
	}

	// Vector3D tstartOld = calculateTangent(prevS, trkParam);
	// Vector3D tendOld   = calculateTangent(currS, trkParam);

	Vector3D tstart = tangentAt( prevS ) ;
	Vector3D tend   = tangentAt( currS ) ;


	fiveByFiveMatrix* jacob = new fiveByFiveMatrix;
	_propagation->getJacobian(*jacob, dw, qbyp, tstart, tend, BField, NrjLoss);
	(*element)->setJacobian(jacob);

	prevS = currS ;

      }
  }


  bool trajectory::fit()
  {
    return _fittingAlgorithm->fit(*this);
  }
}
