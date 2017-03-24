#ifdef USE_LCIO

#include "lcio.h"
#include "IO/LCReader.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/TrackerHitPlane.h"
#include "EVENT/TrackerHitZCylinder.h"
#include "EVENT/Track.h"
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

#include <IMPL/LCCollectionVec.h>
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackStateImpl.h"


// aidaTT
#include "AidaTT.hh"
#include "aidaTT-Units.hh"
#include "analyticalPropagation.hh"
#include "simplifiedPropagation.hh"
#include "GBLInterface.hh"
#include "fitResults.hh"
#include "Vector5.hh"
#include "utilities.hh"
#include "helixUtils.hh"
#include "LCIOPersistency.hh"
#include "Vector3D.hh"
#include "IGeometry.hh"

// ROOT
#include <TTree.h>
#include <TFile.h>

#include <map>

using namespace std ;
using namespace lcio;
using namespace aidaTT ;


// some global steering parameters
#define VERBOSITY streamlog::DEBUG4
//#define VERBOSITY streamlog::MESSAGE

std::string trackCollectionName("MarlinTrkTracks") ;

std::string outColName("AidaTTTracks");

const unsigned trkStateIndex = lcio::TrackState::AtIP ;

bool compute_start_helix = true ;

bool run_prefit = false ; //true ;

bool useQMS = true ;

int maxEvent = 1000 ;

//=======================================================================

aidaTT::analyticalPropagation* propagation = new aidaTT::analyticalPropagation();

aidaTT::GBLInterface* fitter = new aidaTT::GBLInterface();

typedef std::map< long, const aidaTT::ISurface* > SurfMap ;
SurfMap surfMap ;

//=======================================================================


std::string cellIDString( int detElementID) {
  lcio::BitField64 bf(  UTIL::LCTrackerCellID::encoding_string() ) ;
  bf.setValue( detElementID ) ;
  return bf.valueString() ;
}


void getHitInfo( const EVENT::TrackerHit* hit, double* hitpos, 
		 std::vector<double>& precision, const aidaTT::ISurface* surf){
  
  // get the hit position in dd4hep/aidaTT units
  for(unsigned int i = 0; i < 3; ++i) hitpos[i] = hit->getPosition()[i] * dd4hep::mm;
  
  //---- compute the precision from the hit errors
  double du,dv ;
  
  const TrackerHitPlane* planarhit = dynamic_cast<const TrackerHitPlane*>( hit );
  if( planarhit != 0 ) {
    
    du = planarhit->getdU() * dd4hep::mm  ;
    dv = planarhit->getdV() * dd4hep::mm  ;
    
  } else { // we have a TPC hit which is not yet using the CylinderTrackerHit ...
    
    const FloatVec& cov = hit->getCovMatrix();
    
    du = sqrt( cov[0] + cov[2] ) * dd4hep::mm  ;
    dv = sqrt( cov[5]          ) * dd4hep::mm  ;
  }
  
  precision.push_back( 1./ (du*du) );
  
  if( ! surf->type().isMeasurement1D()  )
    precision.push_back( 1./ (dv*dv) );
  else
    precision.push_back( 0. );
}


aidaTT::trackParameters createPreFit(aidaTT::trackParameters& tp, const TrackerHitVec& lcioHits ){
  
  // create a prefit from the hits w/o QMS and dEdx
  
  aidaTT::moveHelixTo( tp, aidaTT::Vector3D() ) ; // move to origin
  
  aidaTT::trajectory traj( tp, fitter, propagation ) ; 
  
  // try to use up to 25 or so hits....
  unsigned nHits = lcioHits.size() ;
  
  int step = ( nHits <= 25  ? 1 : int( 1.*nHits/25. )  ) ; 
  
  streamlog_out( DEBUG2 ) << " MarlinAidaTTTrack::createPreFit() :  will use every " 
			  <<  step << "-th hit for prefit ! " << std::endl ; 
  
  for(unsigned i=0 ; i < nHits ; i+= step ){
    
    EVENT::TrackerHit* hit = lcioHits[i] ;
    
    long hitid = hit->getCellID0() ;
    
    SurfMap::iterator it = surfMap.find( hitid ) ;
    
    if( it == surfMap.end() ){
      
      streamlog_out( DEBUG3 ) << " MarlinAidaTTTrack::createPreFit() : no surface found for id : " 
			      << cellIDString( hitid ) << std::endl ;
      continue;
    }
    
    const aidaTT::ISurface* surf = it->second ;
    
    double hitpos[3] ;
    std::vector<double> precision ;
    getHitInfo( hit, hitpos, precision , surf) ;
    
    streamlog_out( DEBUG2 ) << " MarlinAidaTTTrack::createPreFit() : adding hit for "
			    <<  cellIDString( hitid ) << " at " 
			    << aidaTT::Vector3D( hit->getPosition() ) << std::endl ;
    
    
    
    traj.addMeasurement( hitpos, precision, *surf, hit );
  }
  
  traj.prepareForFitting();
  
  int success = traj.fit();
  
  
  if( success ) { 
    
    const aidaTT::fitResults* result = traj.getFitResults();
    
    const aidaTT::trackParameters& newTP = result->estimatedParameters() ;
    
    streamlog_out( DEBUG4 ) << " MarlinAidaTTTrack::createPreFit() : prefit tp: " 
			    << newTP << std::endl ;
    
    
    return newTP ;
    
  } else {
    
    streamlog_out( WARNING ) << " MarlinAidaTTTrack::createPreFit() : prefit failed for tp: " 
			     << tp << std::endl ;
    
    return tp ;
  }

}





//--------------------------------------------------------------------------------------------------
/* Example program for (re) fitting LCIO tracks with aidaTT.
 * 
 */
int main(int argc, char** argv)
{

  if(argc < 3) {
    std::cout << " usage: ./lcio_tracks compact.xml input_tracks.slcio [output_file.slcio]" << std::endl ;
    return 1;
  }

#ifdef AIDATT_USE_STREAMLOG
  streamlog::out.init( std::cout , "lcio_tracks" ) ;
  streamlog::logscope scope(streamlog::out) ;
  scope.setLevel<VERBOSITY>()  ;
#endif

  std::string inFile =  argv[1] ;

  const aidaTT::IGeometry& geom = aidaTT::IGeometry::instance( inFile ) ;

  const SurfaceVec& surfaces = geom.getSurfaces() ;

  // create map of surfaces
  
  for(std::vector<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf){
    surfMap[(*surf)->id() ] = (*surf) ;
  }
  
  //*********************************************************************
  /// lcio stuff
  std::string lcioFileName = argv[2] ;

  int counter = -1 ;

  LCReader* rdr = LCFactory::getInstance()->createLCReader() ;
  rdr->open(lcioFileName) ;
  LCWriter* wrt = LCFactory::getInstance()->createLCWriter() ;

  if(argc == 4) {
    std::string outFile = argv[3];
    wrt->open(outFile) ;

  } else {
    wrt->open("aidaTT_tracks.slcio", lcio::LCIO::WRITE_NEW ) ;
  }

  LCEvent* evt = 0 ;

  UTIL::BitField64 idDecoder(LCTrackerCellID::encoding_string()) ;

  // create the propagation object
  aidaTT::analyticalPropagation* propagation = new aidaTT::analyticalPropagation();
  //aidaTT::simplifiedPropagation* propagation = new aidaTT::simplifiedPropagation();

  // create the fitter object
  aidaTT::GBLInterface* fitter = new aidaTT::GBLInterface();
  
  
  /// event loop
  while( (evt = rdr->readNextEvent()) != 0 &&  ++counter < maxEvent ) {
    
    LCCollection* trackCollection = evt->getCollection(trackCollectionName) ;
    
    // add output track collection to the event
    LCCollectionVec* outCol = new LCCollectionVec(LCIO::TRACK) ;
    
    evt->addCollection( outCol, outColName ) ;
    
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit(LCIO::TRBIT_HITS) ;
    outCol->setFlag(trkFlag.getFlag()) ;

    int nTracks = trackCollection->getNumberOfElements();
    
    // loop over all tracks in the collection
    for( unsigned i=0 ; i<nTracks ; ++i){
      
      TrackImpl* outTrk = new TrackImpl ;
      outCol->addElement(outTrk) ;

      Track* initialTrack = (Track*) trackCollection->getElementAt(i) ;
    
      aidaTT::trackParameters iTP(  aidaTT::readLCIO( initialTrack->getTrackState( trkStateIndex ) )   );  
    
      const TrackerHitVec& initialHits = initialTrack->getTrackerHits();
      unsigned nHits = initialHits.size() ;

      if( nHits < 3) {
	
	streamlog_out( DEBUG5 ) << " less than three hits - track is dropped ..." << std::endl ;
	
	continue ;
      }


      if( compute_start_helix ) { 

	//----------------------------------------------------------------------------------------------------
	aidaTT::trackParameters startHelix ;
	
	
	//--------- get the start helix from three points
	bool backwards = false ;
	
	lcio::TrackerHit* h1 = ( backwards ?  initialHits[ nHits-1 ] : initialHits[    0    ] ) ;
	lcio::TrackerHit* h2 =  initialHits[ (nHits+1) / 2 ] ;
	lcio::TrackerHit* h3 = ( backwards ?  initialHits[    0    ] : initialHits[ nHits-1 ] ) ;
	
	const double* pos1 = h1->getPosition() ;
	const double* pos2 = h2->getPosition() ;
	const double* pos3 = h3->getPosition() ;
	
	aidaTT::Vector3D x1( pos1[0] * dd4hep::mm, pos1[1] * dd4hep::mm , pos1[2] * dd4hep::mm ) ;
	aidaTT::Vector3D x2( pos2[0] * dd4hep::mm, pos2[1] * dd4hep::mm , pos2[2] * dd4hep::mm ) ;
	aidaTT::Vector3D x3( pos3[0] * dd4hep::mm, pos3[1] * dd4hep::mm , pos3[2] * dd4hep::mm ) ;
	
	calculateStartHelix( x1, x2,  x3 , startHelix , backwards ) ;
	
	moveHelixTo( startHelix, aidaTT::Vector3D(), false  ) ; // move to origin
	
	// --- set some large errors to the covariance matrix
	startHelix.covarianceMatrix().Unit() ;
	startHelix.covarianceMatrix()( aidaTT::OMEGA, aidaTT::OMEGA ) = 1.e-2 ;
	startHelix.covarianceMatrix()( aidaTT::TANL , aidaTT::TANL  ) = 1.e2 ;
	startHelix.covarianceMatrix()( aidaTT::PHI0 , aidaTT::PHI0  ) = 1.e2 ;
	startHelix.covarianceMatrix()( aidaTT::D0   , aidaTT::D0    ) = 1.e5 ;
	startHelix.covarianceMatrix()( aidaTT::Z0   , aidaTT::Z0    ) = 1.e5 ;
	
	streamlog_out( DEBUG4 ) << "  start helix from three points : " << startHelix << std::endl ;
	
	
	// use this helix as start for the fit:
	iTP = ( run_prefit ? createPreFit( startHelix , initialHits ) : startHelix  )  ;
      }
	
      TrackStateImpl* ts;
      bool success;	      
      
      aidaTT::trajectory fitTrajectory( iTP, fitter, propagation, &geom);
      
      const aidaTT::fitResults* result = 0 ; //fitTrajectory.getFitResults();
      
      streamlog_out( DEBUG1 )  << " magnetic field at origin " 
			       << fitTrajectory.geometry()->getBField( aidaTT::Vector3D() ) 
			       << std::endl ;
      

      // copy the hits in order to get strip hits in case of spacepoints
      TrackerHitVec lcioHits ;
      for(unsigned i=0 ; i < nHits ; ++i){
	TrackerHit* trkHit = initialHits[i] ;

	outTrk->addHit( trkHit  ) ;

	if( UTIL::BitSet32( trkHit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]  ){
	  
	  const EVENT::LCObjectVec rawObjects = trkHit->getRawHits();                    
	  
	  for( unsigned k=0; k< rawObjects.size(); k++ ){
	    EVENT::TrackerHit* rawHit = dynamic_cast< EVENT::TrackerHit* >( rawObjects[k] );
	    
	    lcioHits.push_back( rawHit ) ;
          }
	  
	} else { // normal non composite hit

	  lcioHits.push_back( trkHit ) ;
	}
      }


      // ==== store hits in a map and in track  ===============
      std::map< int, EVENT::TrackerHit*> hitMap ;
      for(unsigned i=0 ; i < lcioHits.size() ; ++i){
	hitMap[ lcioHits[i]->getCellID0()  ] = lcioHits[i] ;
      }
      
      //==== compute _all_ surface intersections ====================== 
      const IntersectionVec* intersections = &fitTrajectory.getIntersectionsWithSurfaces( surfaces ) ;
      
      // now we have all surface intersections of the initial track seed
      // however the hits might be on neighbouring sensors really ...

      // ... to be sorted out ...

      //========= loop over all intersections =========
      int pointLabel = 0 ; 
      for( std::vector<std::pair<double, const aidaTT::ISurface*> >::const_iterator it =  
	     intersections->begin() ; it != intersections->end() ; ++it ){
	
	const aidaTT::ISurface* surf = it->second ;
	
	EVENT::TrackerHit* hit = hitMap[ surf->id() ] ;
	
	streamlog_out(DEBUG4) << "intersection - current pointLabel : " << pointLabel  
			      << ":  at s = " << it->first <<  " surface id : " 
			      << cellIDString( surf->id()  ) << std::endl ; 

	streamlog_out(DEBUG1) << aidaTT::pointAt(it->first, iTP ) << std::endl ;

	streamlog_out(DEBUG) << *surf << std::endl ;
	
	if( hit != 0 ){ //-------- we have to add a measurement   
	  
	  double hitpos[3] ;
	  std::vector<double> precision ;
	  getHitInfo( hit, hitpos, precision , surf) ;
	  
	  fitTrajectory.addMeasurement( hitpos, precision, *surf, hit , useQMS );
	  ++pointLabel ;

	  streamlog_out(DEBUG3) << "addMeasurement called for pointLabel : " << pointLabel << std::endl ;
	  
	} else  { // we just add a scatterer
	  
	  if ( useQMS ){
	    
	    // ignore virtual surface with no material (e.g. inside the beam pipe )
	    
	    if( ! ( surf->innerMaterial().density() < 1e-6  && 
		    surf->outerMaterial().density() < 1e-6 )  ) {
	      
	      fitTrajectory.addScatterer( *surf ) ;
	      ++pointLabel ;

	      streamlog_out(DEBUG3) << " addScatterer called for pointLabel : " << pointLabel << std::endl ;
	    }
	  }
	}
      }

      
      fitTrajectory.prepareForFitting();
      
      success = fitTrajectory.fit();
      
      result = fitTrajectory.getFitResults();
      
      
      //***********************************************************************************************************

      if( ! success ) {
	streamlog_out( ERROR ) << " ********** ERROR:  Fit Failed !!!!! ****" 
			       << std::endl ;
      }
      
	      
      
      streamlog_out( DEBUG ) << " End of the loop " << std::endl ;
      streamlog_out( DEBUG ) << " initial values " << std::endl;
      streamlog_out( DEBUG ) << iTP << std::endl;
      streamlog_out( DEBUG ) << " refitted values " << std::endl;
      streamlog_out( DEBUG ) << result->estimatedParameters() << std::endl;
      
      // add Track State to track:
      ts = aidaTT::createLCIO( result->estimatedParameters() );
      //DEBUG: return seed track: ts = aidaTT::createLCIO( iTP );
      
      outTrk->setChi2( result->chiSquare() ) ;
      outTrk->setNdf( result->ndf() ) ;
      outTrk->subdetectorHitNumbers().resize(10.) ;
      
      outTrk->subdetectorHitNumbers()[0] = outTrk->getTrackerHits().size() ;
      
      float ref[3] = { 0., 0. , 0. } ;
      ts->setReferencePoint(ref);	    
      ts->setLocation(lcio::TrackState::AtIP);
      
      
      // checking the covariance matrix
      //--------------------------------------------------------------------
      std::vector<float> cm  =  ts->getCovMatrix();
      trackParameters finalAidaTP = result->estimatedParameters();
      fiveByFiveMatrix  finalAidaCovMat = finalAidaTP.covarianceMatrix();
      
      //---------------------------------------------------------------------
      
      outTrk->addTrackState(ts);
      
    }

    wrt->writeEvent(evt) ;
  }
  
  streamlog_out( DEBUG ) << " counter = " << counter << std::endl ;
  
  return 0;
}

#endif // USE_LCIO
