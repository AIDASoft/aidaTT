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

  std::cout << " --- will check materials along : " << dirVec << std::endl ;


  // create an (almost) straight track in that direction
  double omega = 1.e-12 ; 
  double tanl  = M_PI - theta ;
  double phi0  = phi ;
  double d0    = 0. ;
  double z0    = 0. ;


  aidaTT::trackParameters tP( aidaTT::Vector5(omega, tanl, phi0 , d0 , z0 ), aidaTT::Vector3D() )  ;
  
  const aidaTT::IGeometry& geom = aidaTT::IGeometry::instance( inFile ) ;

  const aidaTT::SurfaceVec& surfaces = geom.getSurfaces() ;

  // create a trajectory w/o fitter or propagation
  aidaTT::trajectory traj( tP, 0, 0, &geom ) ;


  aidaTT::IntersectionVec intersections = traj.getIntersectionsWithSurfaces( surfaces ) ;


  for(unsigned i=0,N=intersections.size() ; i<N ; ++i){

    std::cout << " found intersection at s = " << intersections[i].first << " : "
	      <<  *intersections[i].second << std::endl ;

  }


  return 0;
}
