#ifdef AIDATT_USE_DD4HEP

#include "DD4hepGeometry.hh"

// DD4hep
#include "DDRec/SurfaceHelper.h"
#include "DD4hep/DD4hepUnits.h"

namespace aidaTT
{


  /// helper function for getting a sorting parameter ( basically the max. radius) for the surface
  float getSortingPolicy( const ISurface* surf ){

    //FIXME: this is a very crude sorting and needs revisiting, e.g. 
    // - what happens for disks at +z /-z ....

    if( surf->type().isZCylinder() ){

      return dynamic_cast<const dd4hep::rec::ICylinder*>(surf)->radius() ;

    } else if( surf->type().isZPlane() ){  

      return surf->origin().rho() ;

    } else if( surf->type().isZDisk() ){  

      return std::max( surf->length_along_u(), surf->length_along_v() ) / 2. ;

    } else    if( surf->type().isZCone() ){

      const ICone* c = dynamic_cast<const dd4hep::rec::ICone*>(surf) ;

      return  std::max( c->radius0() ,  c->radius1() ) ;

    } else {
      //fixme: this might be wrong for abritray surfaces 
      return surf->origin().rho() ;
    }
  }


  bool sortSurfaces( const ISurface* surf0 ,  const ISurface* surf1  ){
    
    return  getSortingPolicy(surf0) < getSortingPolicy(surf1)  ;
  }
  

  DD4hepGeometry::DD4hepGeometry(const dd4hep::Detector& thedetector ) : _thedetector( thedetector )  {
    
    const dd4hep::DetElement& det = thedetector.world() ;
    
    // create a list of all surfaces in the detector:
    dd4hep::rec::SurfaceHelper surfMan( det) ;
    
    const dd4hep::rec::SurfaceList& sL = surfMan.surfaceList() ;
    
    _surfaceList.reserve( sL.size() ) ;
    
    for(dd4hep::rec::SurfaceList::const_iterator it = sL.begin() ; it != sL.end() ; ++it){
      _surfaceList.push_back( *it );
    }
    
    std::sort( _surfaceList.begin() , _surfaceList.end() , sortSurfaces ) ; 
  }
  
  const std::vector<const ISurface*>& DD4hepGeometry::getSurfaces() const {
    return _surfaceList ;
  }
  

  Vector3D DD4hepGeometry::getBField( const Vector3D& xx) const {

    Vector3D bfield ;

    _thedetector.field().magneticField( xx.const_array() , bfield.array()  ) ;

    return (1./dd4hep::tesla) * bfield ;
  }

  


  /// =============== get the global instance of the geometry ===========================
  
  IGeometry* IGeometry::_geom = 0 ;
  
  /// return an instance of DD4hepGeometry ( initialize with initString as file name in first call)
  const IGeometry& IGeometry::instance(const std::string& initString ) {
    
    static bool first = true ;
    
    if( first ){
      
      dd4hep::Detector& thedetector = dd4hep::Detector::getInstance();
      thedetector.fromCompact( initString );
      
      IGeometry::_geom = new DD4hepGeometry( thedetector ) ;
      
      first = false;
    }

    return *_geom ;
  }
  

  /// return an instance of DD4hepGeometry assuming thedetector has been initialized elsewhere
  const IGeometry& IGeometry::instance() {
    if( _geom == 0 ) 
      _geom = new DD4hepGeometry(dd4hep::Detector::getInstance() );
    return *_geom ;
  }
}

#endif // AIDATT_USE_DD4HEP
