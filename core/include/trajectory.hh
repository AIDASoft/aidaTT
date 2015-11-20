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
#include "IBField.hh"
#include "IFittingAlgorithm.hh"
#include "Vector3D.hh"

#include "fitResults.hh"

namespace aidaTT
{
  /** The main class representing a track, providing propagation and track fitting.
   *
   *  @version $Id:$
   *  @author Ch. Rosemann, F. Gaede, DESY
   */

  class IFittingAlgorithm;

  class trajectory
  {
  public:
    /// create an empty trajectory 
    trajectory();

    /// copy construct a trajectory -- NOT the internals
    trajectory(const trajectory&);

    /// constructor (A) trajectory that needs to be fitted from a reference track
    trajectory(const trackParameters&, IFittingAlgorithm*, const IBField*, IPropagation*, const IGeometry*);

    //~ the minimal useful constructor
    trajectory(const trackParameters&, const IGeometry*);

    ~trajectory();

    const trackParameters& getInitialTrackParameters() const {
      return _referenceParameters;
    };

    // void setInitialTrackParameters(const trackParameters& tp) 
    // {
    //   // clear the intersections as we have new start parameters:
    //    _surfIntersections.clear() ;
    //   _referenceParameters = tp ;
    // };

    /// methods available before fitting
    const std::vector<trajectoryElement*>& trajectoryElements() const;

    ///~ add a measurement/hit to trajectory; identified by:
    ///~     a position, the resolution, the surface and some id
    void addMeasurement(const Vector3D& position, const std::vector<double>& precision, const ISurface& surface, void* id, bool isScatterer=false ) ;

    /** Compute the multiple scattering parameter for the track at the given surface, 
     *  either using specified parameters (if given) or reference parameters.
     *  Returns 0. if no interesection is found.
     */
    double computeQMS( const ISurface* surface, double &c1, double &c2, const trackParameters* tp= 0) ;
 

    double computeDEdx( const DDSurfaces::IMaterial &mat, double mass, double mom2 );
    double GetEnergyLoss( const ISurface* surface, const trackParameters* tp );

   
    /// add a scatterer for the track going through the surface
    void addScatterer( const ISurface& surface ) ;




    ///~ TODO:: needs more thought!
    ///~ manually add an element to the trajectory; e.g. a point of interest
    ///~ optional: add a surface; e.g. which only contains material
    void addElement(const Vector3D&, void* id);
    // void addElement(const Vector3D&, const ISurface&, void* id);
    // add a scatterer
    // void addScatterer(const Vector3D& position, std::vector<double>& precision, const ISurface& surface, const trackParameters& seed_tp, void* id);

    ///~ test whether/where a surface is intersected
    bool intersectWithSurface(const ISurface* surface, Vector3D& intersect);

    const std::vector<std::pair<double, const ISurface*> >& getIntersectionsWithSurfaces(const std::vector<const ISurface*>&);

    IFittingAlgorithm* getFittingAlgorithm() const;
    IPropagation* getPropagation() const;
    IBField* getBField() const;

    /// the point on the trajectory for a given arc length
    Vector3D pointAt( double s) ;

    /// the tangent to the trajectory for a given arc length
    Vector3D tangentAt( double s ) ;

    /// the trajectoryElement that is valid at the given arc length, i.e.
    /// the one that is defined at the largest arc length less or equal to s
    const trajectoryElement* trajectoryElementAt(double s ) ; 
    

    ///~ prepare: add scattering material, sort elements, calculate and add the jacobians to all elements
    void prepareForFitting();

    ///~ the actual fit call
    bool fit();

    const fitResults* getFitResults(int label=0) ;

    /// TODO:: placeholder; only z component in constant bfield
    double Bz() const
    {
      return _bfieldZ;
    };

    // methods after fitting
    //~ std::vector<trajectoryElement*> getFittedTrajectoryElements() const;
    //~ std::vector<trajectoryElement*> getOutliers() const;

    //~ /// quantify the results
    //~ double getChiSquare() const;
    //~ unsigned int getNDF() const;

  private:

    /// disable assignment
    trajectory operator=(const trajectory&);

    // the internal parts
    trackParameters                  _referenceParameters;
    std::vector<trajectoryElement*>  _initialTrajectoryElements;

    std::vector<std::pair<double, const ISurface*> > _intersectionsList;

    IFittingAlgorithm* _fittingAlgorithm;
    IPropagation* _propagation;

    const IGeometry* const _geometry;
    const IBField* const _bfield;
    double _bfieldZ;


    // /// helper struct for caching the surface intersections
    // struct SurfIntersection{
    //   const ISurface* Surf ;
    //   double  S ;
    //   Vector2D UV ;
    //   Vector3D XX ;
    //   SurfIntersection(const ISurface* surf, double& s, Vector2D* localUV, Vector3D* xx ) : 
    // 	Surf( surf ) ,
    // 	S( s) ,
    // 	UV( *localUV )  ,
    // 	XX( *xx )   {}
    // } ;
    // /// map for caching surface intersections
    // typedef std::map<const ISurface*, SurfIntersection> SIMap ;
    
    // SIMap _surfIntersections ; 

  // public:

  //   bool _calculateIntersectionWithSurface(const ISurface*, double&, Vector2D* = NULL, Vector3D* = NULL);
  //   void _calculateLocalCoordinates(const ISurface*, const Vector3D&, Vector2D*, Vector3D* = NULL);
      
  // protected:
  //   bool _intersectsWithinZCylinderBounds(const ISurface*, double&, Vector2D* = NULL, Vector3D* = NULL);
  //   bool _intersectWithinZPlaneBounds(const ISurface*, double&, Vector2D* = NULL, Vector3D* = NULL);
  //   bool _intersectWithinZDiskBounds(const ISurface*, double&, Vector2D* = NULL, Vector3D* = NULL);

  };
}

#endif // TRAJECTORY_H
