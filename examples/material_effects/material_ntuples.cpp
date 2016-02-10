// aidaTT
#include "AidaTT.hh"
#include "aidaTT-Units.hh"
#include "Vector5.hh"
#include "utilities.hh"
#include "helixUtils.hh"
#include "IGeometry.hh"
#include "trackParameters.hh"

#include "DD4hep/LCDD.h"
#include "DDRec/MaterialManager.h"

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TNtuple.h>

#include <map>
#include <sstream>

using namespace std ;
//using namespace aidaTT ;

using namespace DD4hep ;
using namespace DD4hep::DDRec ;



/* Create Ntuples with material properties along straight lines from the model
 * and from the surfaces.
 */
int main(int argc, char** argv) {
  
  if(argc < 2) {

    std::cout << " usage: ./material_ntuples ILDEx.xml " << std::endl ;
    return 1;
  }
  
  std::string inFile = argv[1] ;

  std::vector<int> thetasDeg ; // ={ 85 , 20 }  --- C++11
  thetasDeg.push_back( 90 ) ;
  thetasDeg.push_back( 85 ) ;
  thetasDeg.push_back( 60 ) ;
  thetasDeg.push_back( 45 ) ;
  thetasDeg.push_back( 30 ) ;
  thetasDeg.push_back( 20 ) ;
  thetasDeg.push_back( 15 ) ;
  thetasDeg.push_back( 10 ) ;
  thetasDeg.push_back( 7 ) ;

  std::vector<int> phisDeg ; //= { 42  } ;
  phisDeg.push_back(  0 ) ;
  phisDeg.push_back(  7 ) ;
  phisDeg.push_back( 12 ) ;
  phisDeg.push_back( 17 ) ;
  phisDeg.push_back( 25 ) ;
  phisDeg.push_back( 30 ) ;
  phisDeg.push_back( 42 ) ;
  phisDeg.push_back( 60 ) ;
  phisDeg.push_back( 71 ) ;
  phisDeg.push_back( 85 ) ;

  //----------------
  std::string varNames =  "theta:phi:epx:epy:epz:x0:lambda:eloss" ;

  TFile *ofile = new TFile( "material_ntuples.root", "RECREATE");
  TNtuple* surfTuple = new TNtuple( "surfTuple" , "material properties from surfaces", varNames.c_str() ) ;
  TNtuple* matTuple  = new TNtuple( "matTuple" , "material properties from detailed model", varNames.c_str() ) ;

  //----------------


  for( unsigned i=0,N= thetasDeg.size() ; i<N ; ++i ){
    for( unsigned j=0,M= phisDeg.size()   ; j<M ; ++j ){
      
      
      double theta = thetasDeg[i] / 180. * M_PI  ;
      double phi = phisDeg[j] / 180. * M_PI ;
      
      
      aidaTT::Vector3D dirVec( 1. , phi , theta,  aidaTT::Vector3D::spherical ) ;
      
      // create an (almost) straight track in that direction
      double omega = 1.e-12 ;
      double tanl  = tan( 0.5*M_PI - theta ) ;
      double phi0  = phi ;
      double d0    = 0. ;
      double z0    = 0. ;
      
      
      aidaTT::Vector3D p0(0.,0.,0.) ; // start point
      aidaTT::Vector3D p1  ; // end point
   
      aidaTT::trackParameters tP( aidaTT::Vector5(omega, tanl, phi0 , d0 , z0 ), p0 )  ;
      
      
      std::cout << " --- theta, phi :       " <<  thetasDeg[i] <<" , " <<  phisDeg[j]  << std::endl ;
      
      const aidaTT::IGeometry& geom = aidaTT::IGeometry::instance( inFile ) ;
      
      const aidaTT::SurfaceVec& surfaces = geom.getSurfaces() ;
      
      // create a trajectory w/o fitter or propagation
      aidaTT::trajectory traj( tP, 0, 0, &geom ) ;
      
      
      aidaTT::IntersectionVec intersections = traj.getIntersectionsWithSurfaces( surfaces ) ;
      
      double X0_tot = 0. ;
      
      // ============================ loop over surface intersections ====================
      for(unsigned ii=0,NN=intersections.size() ; ii<NN ; ++ii){
	
	
	const aidaTT::ISurface* surf = intersections[ii].second ;
	
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
	
	//"theta:phi:epx:epy:epz:x0:lambda:eloss"
	surfTuple->Fill( float(thetasDeg[i]) , float(phisDeg[j]) , xx.x() , xx.y() , xx.z() , X0_tot , 0., 0.  ) ; 

	p1 = xx ; // save last crossing point
      }


      aidaTT::Vector3D end ;
      aidaTT::Vector3D direction = ( p1 -p0 ).unit() ;


      MaterialManager matMgr;
      const MaterialVec& materials = matMgr.materialsBetween(  p0, p1 );
      double sum_x0 = 0;
      double sum_lambda = 0;
      double path_length = 0;


      for( unsigned k=0,K=materials.size();k<K;++k){
	TGeoMaterial* mat =  materials[k].first->GetMaterial();
	double length = materials[k].second;
	double nx0 = length / mat->GetRadLen();
	sum_x0 += nx0;
	double nLambda = length / mat->GetIntLen();
	sum_lambda += nLambda;
	path_length += length;
	end = path_length * direction;

	// ::printf(fmt, i+1, mat->GetName(), mat->GetZ(), mat->GetA(),
	// 	 mat->GetDensity(), mat->GetRadLen(), mat->GetIntLen(), 
	// 	 length, path_length, sum_x0, sum_lambda, end[0], end[1], end[2]);
	// //mat->Print();

	//"theta:phi:epx:epy:epz:x0:lambda:eloss"
	matTuple->Fill( float(thetasDeg[i]) , float(phisDeg[j]) , end.x() , end.y() , end.z() , sum_x0 , sum_lambda, 0.  ) ; 
      }




    } // phi
  } //theta

  ofile->Write("surfTuple");

  return 0;
}
