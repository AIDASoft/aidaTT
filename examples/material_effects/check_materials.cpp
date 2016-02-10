// aidaTT
#include "AidaTT.hh"
#include "aidaTT-Units.hh"
#include "Vector5.hh"
#include "utilities.hh"
#include "helixUtils.hh"
#include "IGeometry.hh"
#include "trackParameters.hh"

// ROOT
#include <TTree.h>
#include <TFile.h>

#include <map>
#include <sstream>

using namespace std ;
//using namespace aidaTT ;


/* Check the materials and material effects along a given direction
 * along a straight line.
 * 
 */
int main(int argc, char** argv) {
  
  if(argc < 3)
    {
      std::cout << " usage: ./check_materials ILDEx.xml theta[deg] phi[deg]" << std::endl ;
      return 1;
    }
  
  std::string inFile   =  argv[1] ;
  std::string thetaStr =  argv[2] ;
  std::string phiStr   =  argv[3] ;
  
  std::stringstream sstrTheta( thetaStr ) ;
  double theta ; 
  sstrTheta >> theta ;
  theta = theta / 180. * M_PI  ;

  std::stringstream sstrPhi( phiStr ) ;
  double phi ; 
  sstrPhi >> phi ;
  phi = phi / 180. * M_PI ;


  aidaTT::Vector3D dirVec( 1. , phi , theta,  aidaTT::Vector3D::spherical ) ;

  // create an (almost) straight track in that direction
  double omega = 1.e-12 ;
  double tanl  = tan( 0.5*M_PI - theta ) ;
  double phi0  = phi ;
  double d0    = 0. ;
  double z0    = 0. ;


  aidaTT::trackParameters tP( aidaTT::Vector5(omega, tanl, phi0 , d0 , z0 ), aidaTT::Vector3D() )  ;
  

  std::cout << " --- will check materials along : " << dirVec << std::endl ;
  
  std::cout << " --- track parameters : " << tP << std::endl ;

  const aidaTT::IGeometry& geom = aidaTT::IGeometry::instance( inFile ) ;

  const aidaTT::SurfaceVec& surfaces = geom.getSurfaces() ;

  // create a trajectory w/o fitter or propagation
  aidaTT::trajectory traj( tP, 0, 0, &geom ) ;


  aidaTT::IntersectionVec intersections = traj.getIntersectionsWithSurfaces( surfaces ) ;

  double X0_tot = 0. ;

 for(unsigned i=0,N=intersections.size() ; i<N ; ++i){
    
    
    const aidaTT::ISurface* surf = intersections[i].second ;

    double s =  0;
    aidaTT::Vector3D xx ;
    bool intersects = aidaTT::intersectWithSurface(  surf , 
						     tP.parameters() , tP.referencePoint() ,
						     s, xx, 0 , true ) ; 

    //----------------- get X0 --------------------------------------------
    const DDSurfaces::IMaterial& material_inn = surf->innerMaterial();
    const DDSurfaces::IMaterial& material_out = surf->outerMaterial();
    
    const double r_i = surf->innerThickness();
    const double r_o = surf->outerThickness();
    
    const double X0_o = material_out.radiationLength();
    const double X0_i = material_inn.radiationLength();
    
    double r_tot = r_i + r_o ;
    
    //calculation of effective radiation length of the surface
    double X0_eff = ( r_i/X0_i + r_o/X0_o ) / r_tot ; 
    
    //calculation of the path of the particle inside the material
    //compute path as projection of (straight) track to surface normal:
    
    const aidaTT::Vector3D& up = dirVec.unit() ;

        // need to get the normal at crossing point 
    const aidaTT::Vector3D& n = surf->normal( xx ) ;
    
    double cosTrk = std::fabs( up * n )  ;
    
    double path = r_i + r_o ;
    
    path = path/cosTrk ; 
    
    double X_X0 = path * X0_eff ;

    //--------------------------------------------------------------------

    X0_tot += X_X0 ;


    std::cout << " found intersection at s = " << intersections[i].first << ", " 
	      << xx << " surf: " << std::endl 
	      << *surf 
	      << " X0 : " << X_X0   << "   - total : " << X0_tot << std::endl ;;


  }


  return 0;
}
