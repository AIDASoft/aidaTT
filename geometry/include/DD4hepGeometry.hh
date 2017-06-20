#ifdef AIDATT_USE_DD4HEP

#ifndef DD4HEPGEOMETRY_HH
#define DD4HEPGEOMETRY_HH

#include "IGeometry.hh"
#include <DD4hep/Detector.h>
#include <vector>

namespace aidaTT
{
  class DD4hepGeometry : public IGeometry
  {
  public:

    DD4hepGeometry(const dd4hep::Detector& theDetector );
    
    /// get a list of all surfaces in the tracking geometry 
    /// loosely sorted with radius
    virtual const std::vector<const ISurface*>& getSurfaces() const ;
    
    /// access the B field in Tesla at given position
    virtual Vector3D getBField( const Vector3D& xx) const ;
    
    /// d'tor
    virtual ~DD4hepGeometry(){}
    

  private:
    const dd4hep::Detector& _thedetector ;

    std::vector<const ISurface* > _surfaceList;
  };
}
#endif // DD4HEPGEOMETRY_HH

#endif // AIDATT_USE_DD4HEP
