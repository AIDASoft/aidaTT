#include "materialUtils.hh"

#include "helixUtils.hh"
#include "IGeometry.hh"
#include "AidaTT-Units.hh"

namespace aidaTT{


  double computeQMS( const ISurface* surf, const Vector5& hp, const Vector3D& rp , double mass ) {
    
    double s =  0;
    Vector3D xx ;
    bool intersects = aidaTT::intersectWithSurface( surf, hp , rp , s, xx, 0 , true ) ; 
    
    if( !intersects ) 
      return 0. ; 
    
    Vector2D uv = surf->globalToLocal( xx ) ;
    
    double phi   = calculatePhifromXY( xx.x(), xx.y(), hp, rp) ;
    double omega = calculateOmega( hp ) ;
    double tanl  = calculateTanLambda( hp ) ;
    

    double bfieldZ  = IGeometry::instance().getBField( xx ).z() ;

    double pt = ( fabs(1./omega ) * bfieldZ * aidaTT::convertBr2P_cm  ); 

    Vector3D p( pt*std::cos( phi ), pt*std::sin( phi ) , pt*tanl ) ;
    
    return computeQMS( surf, uv, p , mass ) ;
  }
  

  double computeQMS( const ISurface* surf, const Vector2D& crossingPoint, const Vector3D& momentum , double mass ) {
    
    const DDSurfaces::IMaterial& material_inn = surf->innerMaterial();
    const DDSurfaces::IMaterial& material_out = surf->outerMaterial();
    
    const double r_i = surf->innerThickness();
    const double r_o = surf->outerThickness();
    
    const double X0_o = material_out.radiationLength();
    const double X0_i = material_inn.radiationLength();
    
    double r_tot = r_i + r_o ;
    
    //calculation of effective radiation length of the surface
    double X0_eff = ( r_i/X0_i + r_o/X0_o ) / ( r_i + r_o ) ; 
    
    //calculation of the path of the particle inside the material
    //compute path as projection of (straight) track to surface normal:
    
    DDSurfaces::Vector3D up = momentum.unit() ;
    
    const Vector3D& position = surf->localToGlobal( crossingPoint ) ;
    
    // need to get the normal at crossing point 
    const Vector3D& n = surf->normal( position ) ;
    
    double cosTrk = std::fabs( up * n )  ;
    
    double path = r_i + r_o ;
    
    path = path/cosTrk ; 
    
    double X_X0 = path * X0_eff ;
    
    double mom = momentum.r() ;
    
    double   beta = mom / std::sqrt(mom * mom + mass * mass);
    
    double Qms = 0.0136/(mom*beta) * 1.0 * std::sqrt(X_X0) * (1 + 0.038*(std::log(X_X0)));

    return Qms ;
  }


}
