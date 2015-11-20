#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>
#include <list>
#include <map>

#include "trackParameters.hh"
#include "trajectoryElement.hh"
#include "helixUtils.hh"
#include "utilities.hh"
#include "IPropagation.hh"
#include "IGeometry.hh"
#include "IFittingAlgorithm.hh"

#include "fitResults.hh"

namespace aidaTT
{
  /** The main class representing a track, providing propagation and track fitting.
   *
   *  @version $Id:$
   *  @author Ch. Rosemann, F. Gaede, DESY
   */

  class IFittingAlgorithm;
  class trajectoryElement ;
  typedef std::vector<trajectoryElement*> ElementVec ;
  typedef std::vector<std::pair<double, const ISurface*> > IntersectionVec ;
  typedef std::vector<const ISurface*> SurfaceVec ;

  class trajectory {

  public:
    
    /** Create a trajectory for fitting with an initial set of track parameters, a fittingAlgorithm,
     *  a propagator and (optionally) a geometry implementation.
     */ 
    trajectory(const trackParameters& tP, IFittingAlgorithm* fitter, IPropagation* propagator, 
	       const IGeometry* geom=&IGeometry::instance() );

     /// the minimal constructor providing an intitial track state - no fitting will be possible
    /// could be used for intersection calculations
    trajectory(const trackParameters& tP, const IGeometry* geom=&IGeometry::instance() );

   ///d'tor
    ~trajectory();

    /// the intial track parameters given at construction
    const trackParameters& initialTrackParameters() const {   return _referenceParameters; };
    
    /// the geometry system used for this trajectory
    const IGeometry* geometry() const {  return _geometry;  }
    
    /// the fitting algorithm used for this trajectory
    IFittingAlgorithm* fittingAlgorithm() const {  return _fittingAlgorithm ; }
    
    /// the propagation for this trajectory
    IPropagation* propagation() const { return _propagation ; }

    /// vector of trajectory elements that have been added to the trajectory ( sorted in s )
    const ElementVec& trajectoryElements() const { return _initialTrajectoryElements; }
    

    /** add a trajectoryElement with a measurement/hit to the trajectory - identified by:
     *  a position, the precision, the surface and a user defined id/pointer
     */
    void addMeasurement(const Vector3D& position, const std::vector<double>& precision, 
			const ISurface& surface, void* id, bool isScatterer=false ) ;
    
    
    /** add a trajectoryElement that is a scatterer for the track going through the surface
     */
    void addScatterer( const ISurface& surface ) ;
    
    
    /** Helper method that computes the intersections of the track - with all surfaces 
     *  given as argument. Calculation is based on the initial track parameters.
     *  Returns a vector of pairs( arcLength, Surface*). 
     *  Can be used to get the right order for call to addScatterer and addMeasurment.
     */
    const IntersectionVec& getIntersectionsWithSurfaces(const SurfaceVec&) ;
    
    /// the point on the trajectory for a given arc length
    Vector3D pointAt( double s ) ;
    
    /// the tangent to the trajectory for a given arc length
    Vector3D tangentAt( double s ) ;
    
    /// the trajectoryElement that is valid at the given arc length, i.e.
    /// the one that is defined at the largest arc length less or equal to s
    const trajectoryElement* trajectoryElementAt( double s ) ; 
    
    ///needs to be called before the actual track fit for internal preparation
    void prepareForFitting();

    /// fit the track based on the measurements and scatterers added by the user
    bool fit();

    /// return the fit result at the given label, where label==0 corresponds to
    /// s == 0. 
    const fitResults* getFitResults(int label=0) ;


    // methods after fitting
    //~ std::vector<trajectoryElement*> getFittedTrajectoryElements() const;
    //~ std::vector<trajectoryElement*> getOutliers() const;

    //~ /// quantify the results
    //~ double getChiSquare() const;
    //~ unsigned int getNDF() const;

  private:

    /// inernal helper method for adding an initial start element to the trajectory
    void addElement(const Vector3D&, void* id);

    // disable assignment
    trajectory operator=(const trajectory&);

    // =========== members =========================================
    trackParameters     _referenceParameters;
    IFittingAlgorithm*  _fittingAlgorithm;
    IPropagation*       _propagation;
    const IGeometry* const _geometry;
    
    std::vector<trajectoryElement*>  _initialTrajectoryElements;
    std::vector<std::pair<double, const ISurface*> > _intersectionsList;
    //==============================================================
    
    
  private:
    //fg: these c'tors do not make much sense currently ( no setters for missing parameters are available )
    
    // create an empty trajectory 
    trajectory();
    
    // copy construct a trajectory -- NOT the internals
    trajectory(const trajectory&);
    
  };
}

#endif // TRAJECTORY_H
