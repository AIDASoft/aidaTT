#include "materialUtils.hh"

#include "helixUtils.hh"
#include "IGeometry.hh"
#include "aidaTT-Units.hh"
#include <math.h>
#include <cmath>
#include "streamlog/streamlog.h"

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
    
    const dd4hep::rec::IMaterial& material_inn = surf->innerMaterial();
    const dd4hep::rec::IMaterial& material_out = surf->outerMaterial();
    
    const double r_i = surf->innerThickness();
    const double r_o = surf->outerThickness();
    
    const double X0_o = material_out.radiationLength();
    const double X0_i = material_inn.radiationLength();
    
    double r_tot = r_i + r_o ;
    
    //calculation of effective radiation lengths per cm of the surface
    double X0_eff = ( r_i/X0_i + r_o/X0_o ) / r_tot ; 
    
    //calculation of the path of the particle inside the material
    //compute path as projection of (straight) track to surface normal:
    
    const Vector3D& up = momentum.unit() ;
    
    const Vector3D& position = surf->localToGlobal( crossingPoint ) ;
    
    // need to get the normal at crossing point 
    const Vector3D& n = surf->normal( position ) ;
    
    double cosTrk = std::fabs( up * n )  ;
    
    double path = r_i + r_o ;
    
    path = path/cosTrk ; 
    
    //    if( path > 0.1 ) path = 0.1 ; // limit the path length !!!! @@@@ ????

    double X_X0 = path * X0_eff ;
    
    double mom = momentum.r() ;
    
    double   beta = mom / std::sqrt(mom * mom + mass * mass);
    
    double Qms = 0.0136/(mom*beta) * 1.0 * std::sqrt(X_X0) * (1 + 0.038*(std::log(X_X0)));


    streamlog_out( DEBUG ) << " **QMS: surface : " << *surf << std::endl   ;
    
    streamlog_out( DEBUG2 ) << " **QMS: crossingPoint : (" <<  crossingPoint.u() << ", " 
			    << crossingPoint.v() << " ) " 
			    << " position = " << position
			    << " normal = " << n
			    << " QMS = " << Qms
			    << " path length = " << path
			    << " X0_eff : " << X0_eff
			    << std::endl ;

    if( Qms != Qms ){

      streamlog_out( ERROR ) << " **QMS: crossingPoint : (" <<  crossingPoint.u() << ", " 
			    << crossingPoint.v() << " ) " 
			    << " position = " << position
			    << " normal = " << n
			    << " QMS = " << Qms
			    << " path length = " << path
			    << " X0_eff : " << X0_eff
			    << std::endl ;

    }

    return Qms ;
  }



  double computeBetheBloch( const IMaterial &mat, double mom, double mass ){
    // -----------------------------------------
    // Bethe-Bloch eq. (Physical Review D P195.)
    // code copied from KalTest, K.Fujii, KEK
    // -----------------------------------------
    static const double kK   = 0.307075e-3;     // [GeV*cm^2]
    static const double kMe  = 0.510998902e-3;  // electron mass [GeV]
    
    double dnsty = mat.density();             // density
    double A     = mat.A();                   // atomic mass
    double Z     = mat.Z();                   // atomic number
    //double I    = Z * 1.e-8;                   // mean excitation energy [GeV]
    //double I    = (2.4 +Z) * 1.e-8;            // mean excitation energy [GeV]
    double I    = (9.76 * Z + 58.8 * std::pow(Z, -0.19)) * 1.e-9; // mean excitation energy [GeV]
    double hwp  = 28.816 * std::sqrt(dnsty * Z/A) * 1.e-9;
    double bg2  = (mom*mom) / (mass * mass);
    double gm2  = 1. + bg2;
    double meM  = kMe / mass;
    double x    = log10(std::sqrt(bg2));
    double C0   = - (2. * log(I/hwp) + 1.);
    double a    = -C0/27.;
    double del;
    if (x >= 3.)            del = 4.606 * x + C0;
    else if (0.<=x && x<3.) del = 4.606 * x + C0 + a * std::pow(3.-x, 3.);
    else                    del = 0.;
    double tmax = 2.*kMe*bg2 / (1. + meM*(2.*std::sqrt(gm2) + meM));
    double dedx = kK * Z/A * gm2/bg2 * (0.5*log(2.*kMe*bg2*tmax / (I*I))
					  - bg2/gm2 - del);

    return dedx ;
  }


  double computeEnergyLoss( const ISurface* surf, const Vector2D& crossingPoint, 
			    const Vector3D& momentum, double& energy, double& beta,
			    double mass ){
    
    const dd4hep::rec::IMaterial& material_i = surf->innerMaterial(); // crossingPoint
    const dd4hep::rec::IMaterial& material_o = surf->outerMaterial(); // crossingPoint
    
    double path_i = surf->innerThickness();
    double path_o = surf->outerThickness();
    
    const Vector3D& position = surf->localToGlobal( crossingPoint ) ;
    
    const Vector3D& up = momentum.unit() ;

    // need to get the normal at crossing point 
    const Vector3D& n = surf->normal( position ) ;

    double cosTrk = std::fabs( up * n )  ;
    
    path_i /= cosTrk ;
    path_o /= cosTrk ;
    
    double mom = momentum.r() ;
    double mom2 = mom*mom ;

    energy = std::sqrt(mom2 + mass*mass );

    beta = mom / energy ;

    double density_i = material_i.density() ;
    double density_o = material_o.density() ;
    
    double dedx_i = computeBetheBloch( material_i, mom, mass ) ;
    double dedx_o = computeBetheBloch( material_o, mom, mass ) ;
    
    double deltaE = dedx_i * density_i * path_i ;
    deltaE += dedx_o * density_o * path_o ;
 

    if( deltaE != deltaE ){

      streamlog_out( ERROR ) << " **computeEnergyLoss : (" <<  crossingPoint.u() << ", " 
			     << crossingPoint.v() << " ) " 
			     << " normal = " << n
			     << " deltaE : " << deltaE
			     << " path length = " << path_o + path_i 
			     << " momentum : " << momentum
			     << std::endl ;

    }
    

    return deltaE ;
  }


  double computeEnergyLoss( const ISurface* surf, const trackParameters& tp , 
			    double& energy, double& beta,
			    double mass ) {

    double s =  0;
    Vector3D xx ;
    bool intersects = aidaTT::intersectWithSurface( surf, tp.parameters() , tp.referencePoint() ,
						    s, xx, 0 , true ) ; 
    
    if( !intersects ) {
      // std::cout << " computeEnergyLoss - point : " << xx << " at " << s 
      // 		<< " does not intersect with surface : " << *surf << " distance : " 
      // 		<< surf->distance( xx ) << std::endl ; 
      return 0. ; 
    }

    Vector2D uv = surf->globalToLocal( xx ) ;
    
    double phi   = calculatePhifromXY( xx.x(), xx.y(), tp ) ;
    double omega = calculateOmega( tp ) ;
    double tanl  = calculateTanLambda( tp ) ;
    

    double bfieldZ  = IGeometry::instance().getBField( xx ).z() ;

    double pt = ( fabs(1./omega ) * bfieldZ * aidaTT::convertBr2P_cm  ); 

    Vector3D p( pt*std::cos( phi ), pt*std::sin( phi ) , pt*tanl ) ;
    
    return computeEnergyLoss( surf, uv, p ,energy, beta, mass ) ;
  }

}
