#ifndef IGEOMETRY_HH
#define IGEOMETRY_HH

#include <vector>


#ifdef AIDATT_USE_DD4HEP
//
// Use the Surfaces and material classes (interfaces) from DD4hep 
//
#include "DDSurfaces/ISurface.h"
#include "DDSurfaces/IMaterial.h"

namespace aidaTT
{
  using DDSurfaces::ISurface ;
  using DDSurfaces::ICylinder ;
  using DDSurfaces::ICone ;
  using DDSurfaces::IMaterial ;
  using DDSurfaces::Vector2D ;
  using DDSurfaces::Vector3D ;
  using namespace DDSurfaces ;

}
#else

// include a copy of the interfaces from DD4hep
#include "AidaTTGeometry.hh"
// only meaningfull if there is an implementation
// other than DD4hep....

#endif 




namespace aidaTT {
  
  /** The geometry interface for aidaTT provides
   *  access to the tracking surfaces and the 
   *  B field.
   *  For the  implementation see DD4hepGeometry.hh.
   * 
   * @author C.Rosemann, F.gaede, DESY
   * @version $Id:$
   */
  class IGeometry{
    
  public:
    /// get a list of all surfaces in the tracking geometry 
    /// loosely sorted with radius
    virtual const std::vector<const ISurface*>& getSurfaces() const = 0;
    
    /// access the B field in Tesla at given position
    virtual Vector3D getBField( const Vector3D& xx) const = 0;
    
    /// d'tor
    virtual ~IGeometry(){}
    
    /// get the global instance of the geometry after intializing it based on the init string, e.g. a file name
    static const IGeometry& instance(const std::string& initString ) ;

    /// get the global instance of the geometry
    static const IGeometry& instance() ;

  protected:
    
    static IGeometry* _geom ;
  };
  

}
#endif // IGEOMETRY_HH
