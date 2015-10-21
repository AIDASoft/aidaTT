#ifdef USE_GBL
#include "GBLInterface.hh"
#include "utilities.hh"

namespace aidaTT
{
  GBLInterface::GBLInterface() : _trajectory(NULL), _correctionVector(NULL), _covarianceMatrix(NULL)
  {}



  GBLInterface::~GBLInterface()
  {
    if(_trajectory != NULL)
      delete _trajectory;

    _clear() ;
  }



  bool GBLInterface::fit(const trajectory& TRAJ)
  {
    /* several bits of information are needed to initialize the gbl:
     *  - a vector of GblPoints, which in turn need a p2p jacobian to be instantiated
     *  -- a measurement GblPoint needs the projection matrix, the residual vector and the precision
     *  -- a scattering GblPoint needs the precision (expected inverse standard deviation)
     */

    // clear any results from a previous fit
    _clear() ; 


    /// create vector of GBL points
    std::vector < gbl::GblPoint > theListOfPoints;

    const std::vector<trajectoryElement*>& elements = TRAJ.trajectoryElements();

    _fittedTraj = &TRAJ ;

    for(std::vector<trajectoryElement*>::const_iterator element = elements.begin(), last = elements.end(); element < last; ++element)
      {
	const fiveByFiveMatrix& jac = (*element)->jacobian();

	//std::cout << " input jacobian " << jac << std::endl ;

	///~ initialise point with jacobian from last to the current element
	gbl::GblPoint point(TMatrixD(5, 5, jac.array()));

	//std::cout << " ---  GBLInterface::fit - element : " <<  **element << std::endl ;

	//std::cout << " Now I am working on surface " << (*element)->surface();

	if((*element)->hasMeasurement())
	  {

	    std::cout << " Do I have a measurement ? " << std::endl ;

	    /// three elements are needed to add a measurement to a gblpoint:
	    /// 1) the projection matrix from the local track frame to the measurement system
	    /// 2) the residuals in the measurement directions
	    /// 3) the inverse resolution = precision of the measurements

	    const unsigned int mDim = (*element)->measurementDimension();

	    if(mDim > 2)
	      throw std::invalid_argument("Error: Currently only 1D or 2D measurements are implemented.");

	    //~ 1) projection matrix -- the basis change matrix
	    const std::vector<double>& projLocal2Meas = (*element)->localToMeasurementProjection();

	    /// convert the data from the vector into the matrix:
	    /// convention is that the first row comes first in the data
	    TMatrixD pL2M(2, 2);
	    pL2M(0, 0) = projLocal2Meas.at(0);
	    pL2M(0, 1) = projLocal2Meas.at(1);
	    pL2M(1, 0) = projLocal2Meas.at(2);
	    pL2M(1, 1) = projLocal2Meas.at(3);


	    //~ 2) the residuals in the measurement direction
	    const std::vector<double>& residuals = (*element)->measurementResiduals();

	    //~ convert the vector to an array:
	    const double* resid = residuals.data();

	    //~ 3) the precision of the measurements -- the inverse of the resolution
	    const std::vector<double>& precision = (*element)->precisions();
	    //const TMatrixDSym& precision = (*element)->precisions();


	    //~ convert the vector to an array:
	    const double* prec = precision.data();

	    // fixed size of arguments: 2D in measurements!
	    point.addMeasurement(pL2M, TVectorD(2, resid), TVectorD(2, prec));
	    //point.addMeasurement(pL2M, TVectorD(2, resid), precision);
	  }
	
	
	if((*element)->isScatterer())
	  {

	    // TODO :: implement!
	    //~ now evaluate the scattering material
	    //~ the addScatterer routine will add a thin scatterer at the given position
	    //~ a thin scatter only changes the local direction, no offset. Multiple step approach:
	    //~ -- 1.) determine whether scatterer or not
	    //~ -- 2.) if scatterer: thin or thick
	    //~ -- 3.) if thick scatterer: calculate two positions and material properties at the points!
	    //~ the scattering info enters through the inverse covariance matrix, doesn't need to be diagonalized before invocation
	    //~ point.gbl::addScatterer ( TVectorD notNeededHere, TMatrixDSym aPrecision   );

	    //~ 1) projection matrix -- the basis change matrix
	    const std::vector<double>& projLocal2Meas = (*element)->localToMeasurementProjection();


	    //std::cout << " size of projection vector " << projLocal2Meas.size() << std::endl ;

	    /// convert the data from the vector into the matrix:
	    /// convention is that the first row comes first in the data



		      
	    //~ 2) the residuals in the measurement direction
	    const std::vector<double>& residuals = (*element)->measurementResiduals();
	    //~ convert the vector to an array:
	    const double* resid = residuals.data();

	    //~ 3) precision - MPS


	    std::vector<double> precision = (*element)->precisions();

	    // fg: get the multiple scattering sigma from the precision vector - either first element or first element after measurement's precisions
	    double qms = (  (*element)->hasMeasurement() ?  precision[  (*element)->measurementDimension() ] : precision[0]   )  ;


	    double c1 = (  (*element)->hasMeasurement() ?  precision[  (*element)->measurementDimension() + 1 ] : precision[1]   )  ;

	    double c2 = (  (*element)->hasMeasurement() ?  precision[  (*element)->measurementDimension() + 2 ] : precision[2]   )  ;

	    TMatrixDSym Vk_sym(2);
	    double Scalar_value = 0 ;
	    Scalar_value = (1 - c1*c1 - c2*c2) / qms ;
	    Vk_sym(0, 0) = Scalar_value * (1 - c1*c1);
	    Vk_sym(0, 1) = -1.*Scalar_value * c1*c2;
	    Vk_sym(1, 0) = -1.*Scalar_value * c1*c2;
	    Vk_sym(1, 1) = Scalar_value * (1 - c2*c2);
	    

	    point.addScatterer( TVectorD(2, resid), Vk_sym);


	  }

	// store the point in the list that will be handed to the trajectory
	theListOfPoints.push_back(point);
      }

    std::cout << " we have added " << theListOfPoints.size() << " points to the track " << std::endl ;

    /// TODO :: check validity before continuing!

    //~ this is not elegant -- delete previous data before starting again
    if(_trajectory != NULL)
      delete _trajectory;

    _trajectory = new gbl::GblTrajectory(theListOfPoints, true); /// TODO: pass info about magnetic field

    unsigned int returnValue = _trajectory->fit(_chisquare, _ndf, _lostweight);



    //_trajectory->printTrajectory(100) ;
    //_trajectory->printPoints(100) ;


    if(returnValue == 0)
      {
	// _fillResults( TRAJ, 0 );
	return true;
      }
    else
      return false;
  }



  const fitResults* GBLInterface::_fillResults(const trajectory& TRAJ, int label) const
  {

    ResMap::const_iterator it = _theResults.find( label ) ;

    if ( it != _theResults.end() ) 
      return it->second ;


    bool v = true;
    double chs = _chisquare;
    unsigned int n = _ndf;
    double wl = _lostweight;

    TVectorD tpCorr(5);
    TMatrixDSym trackcovariance(5);

    //~ get the results at a given label in local cl track parameters
    //~ the track parameters are corrections to the curvilinear track parameters
    _trajectory->getResults( label , tpCorr, trackcovariance);

    Vector5 clCorrections(tpCorr[0], tpCorr[1], tpCorr[2], tpCorr[3], tpCorr[4]);

    fiveByFiveMatrix cl2L3Jacobian  = curvilinearToL3Jacobian(TRAJ.getInitialTrackParameters(), Vector3D(0., 0., TRAJ.Bz())) ;
    Vector5 L3corrections           =  cl2L3Jacobian * clCorrections;

 

    Vector5 fittedParameters = TRAJ.getInitialTrackParameters().parameters() + L3corrections;


    TMatrixD cl2L3_copy(5,5);
    for(int i = 0 ; i < 5 ; i++)
      {
	for(int j = 0 ; j < 5 ; j++)
	  {
	    cl2L3_copy(i, j) = cl2L3Jacobian(i, j);
	  }
      }
    /*
      for(int i = 0 ; i < 5 ; i++)
      {
      for(int j = 0 ; j < 5 ; j++)
      {
      std::cout << " original jacobian " << cl2L3Jacobian(i, j) << " copied jacobian " << cl2L3_copy(i, j) << std::endl ;
      }
      }
    */
    trackParameters tp;
    tp.setTrackParameters(fittedParameters);

    TMatrixD covarianceMatrix = trackcovariance.Similarity(cl2L3_copy);

    fiveByFiveMatrix finalCov;
    for(int i = 0 ; i < 5 ; i++)
      {
	for(int j = 0 ; j < 5 ; j++)
	  {
	    finalCov(i, j) = covarianceMatrix(i, j);
	  }
      }

    tp.setCovarianceMatrix(finalCov);

    //FIXME: the reference point is still the origin - we need to call a properly impemented moveTo 
    //       method to get it at the position of this label (hit)...
    

    fitResults* res = new fitResults( v, chs, n, wl, tp ) ;
  
    _theResults.insert( std::make_pair(  label , res )  ) ;
  
    return res ;
			
  }

  void  GBLInterface::_clear(){

    for( ResMap::iterator it = _theResults.begin() ; it != _theResults.end() ; ++it){
      delete it->second ;
    } 

    _theResults.clear() ;

    _fittedTraj = 0 ;
  }

  
}


#endif
