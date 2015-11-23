#ifdef USE_LCIO

#include "lcio.h"
#include "IO/LCReader.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackerHitImpl.h"
#include "EVENT/TrackerHitPlane.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackStateImpl.h"
#include "UTIL/ILDConf.h"
#include "UTIL/Operators.h"

#include <IMPL/LCCollectionVec.h>
#include "IMPL/TrackImpl.h"

// aidaTT
#include "AidaTT.hh"
#include "aidaTT-Units.hh"
#include "analyticalPropagation.hh"
#include "simplifiedPropagation.hh"
#include "GBLInterface.hh"
#include "fitResults.hh"
#include "Vector5.hh"

#include <map>

using namespace std ;
using namespace lcio;

/* this is an example for debugging: it reads in an lcio Track collection
 * and creates a new Track collections with hits at the comuted intersection points
 *
 */


int main(int argc, char** argv)
{
  if(argc < 3){
    std::cout << " usage: ./lcio_debug_example ILDEx.xml ILDExSimu.slcio" << std::endl ;
    return 1;
  }
  
  std::string inFile =  argv[1] ;
    
  const aidaTT::IGeometry& geom = aidaTT::IGeometry::instance() ;

  const std::vector<const aidaTT::ISurface*>& surfaces = geom.getSurfaces() ;

  // create map of surfaces
  std::map< long, const aidaTT::ISurface* > surfMap ;
  
  for(std::vector<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf){
    surfMap[(*surf)->id() ] = (*surf) ;
  }

    
  /// lcio stuff
  std::string lcioFileName = argv[2] ;
    
  LCReader* rdr = LCFactory::getInstance()->createLCReader() ;
  rdr->open(lcioFileName) ;

  LCWriter* wrt = LCFactory::getInstance()->createLCWriter() ;
  wrt->open("lcio_debug_tracks.slcio") ;
    

  LCEvent* evt = 0 ;
    
  std::string trackCollectionName = "SiTracks";
    
  UTIL::BitField64 idDecoder(ILDCellID0::encoder_string) ;
    
    
  // create the propagation object
  aidaTT::analyticalPropagation* propagation = new aidaTT::analyticalPropagation();
  //aidaTT::simplifiedPropagation* propagation = new aidaTT::simplifiedPropagation();
    
  // create the fitter object
  aidaTT::GBLInterface* fitter = new aidaTT::GBLInterface();
    
    
  /// event loop
  while((evt = rdr->readNextEvent()) != 0) {
      
    LCCollection* trackCollection = evt->getCollection(trackCollectionName) ;

    LCCollectionVec* outCol = new LCCollectionVec( LCIO::TRACK ) ;
    evt->addCollection( outCol ,  "AidaTTDebugTracks"   ) ; 
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    outCol->setFlag( trkFlag.getFlag()  ) ;

    TrackImpl* outTrk = new TrackImpl ;  
    outCol->addElement( outTrk ) ;

    LCCollectionVec* outColH = new LCCollectionVec( LCIO::TRACKERHIT ) ;
    evt->addCollection( outColH ,  "AidaTTDebugHits"   ) ; 

    int nTracks = trackCollection->getNumberOfElements();
      
    // ignore event if more or less than a single track is present
    if(nTracks != 1)
      continue;
      
    Track* initialTrack = (Track*)trackCollection->getElementAt(0);

    aidaTT::trackParameters iTP;

    // ****** use TrackState at ... *****************
    //  const TrackState* ts = initialTrack->getTrackState( TrackState::AtIP ) ;
    //  const TrackState* ts = initialTrack->getTrackState( TrackState::AtFirstHit ) ;
    const TrackState* ts = initialTrack->getTrackState( TrackState::AtLastHit ) ;
    //  const TrackState* ts = initialTrack->getTrackState( TrackState::AtCalorimeter ) ;

	    iTP.setTrackParameters(aidaTT::Vector5(ts->getOmega()/aidaTT::mm, ts->getTanLambda(), ts->getPhi(), ts->getD0()*aidaTT::mm, ts->getZ0()*aidaTT::mm ));
      


    iTP.setReferencePoint( aidaTT::Vector3D( ts->getReferencePoint()[0]*aidaTT::mm , 
					     ts->getReferencePoint()[1]*aidaTT::mm, 
					     ts->getReferencePoint()[2]*aidaTT::mm ) ) ;




    outTrk->addTrackState( new TrackStateImpl( *ts ) ) ;

    std::cout << " initialize track parameters with TrackState : " << *ts << std::endl ;

    //**********************************************

    aidaTT::trajectory fitTrajectory(iTP, fitter, propagation, &geom);
      
    std::vector<TrackerHit*> initialHits = initialTrack->getTrackerHits();
      
    for(std::vector<TrackerHit*>::iterator thit = initialHits.begin(), endIter = initialHits.end(); thit < endIter; ++thit) {
	
      long64 hitid = (*thit)->getCellID0() ;
      idDecoder.setValue(hitid) ;
	
      // if(idDecoder[ lcio::ILDCellID0::subdet] != lcio::ILDDetID::VXD)
      // 	continue;
	
      if(idDecoder[ lcio::ILDCellID0::subdet] == lcio::ILDDetID::VXD) {

	idDecoder[lcio::ILDCellID0::side] = ((*thit)->getPosition()[2]  >  0  ?   +1 : -1) ;
	
	unsigned layerID = idDecoder[lcio::ILDCellID0::layer] ;

	idDecoder[lcio::ILDCellID0::layer] = layerID + 1 ;
	
	hitid = idDecoder.lowWord() ;

      }
	
      const aidaTT::ISurface* surf = surfMap[ hitid ] ;
      
      if(surf == NULL){
	std::cerr << " lcio_debug_tracks : no surface found for id : " << idDecoder.valueString() << std::endl ;
	continue;
      }
	
      //========= compute intersection with this layer ================================

      double s = 0. ;
      aidaTT::Vector3D xxAidaTT , xx ;
      bool foundIntersect = aidaTT::intersectWithSurface( surf , fitTrajectory.initialTrackParameters(),
							  s , xxAidaTT , 0 , true );
      
      
      if( foundIntersect ){
	
	s /= aidaTT::mm  ; 
	
	xx.fill( xxAidaTT[0]/aidaTT::mm , xxAidaTT[1]/aidaTT::mm,  xxAidaTT[2]/aidaTT::mm) ;   
	

	aidaTT::Vector3D posV( (*thit)->getPosition()[0],  (*thit)->getPosition()[1], (*thit)->getPosition()[2] ) ;


	TrackerHitImpl* newHit = new TrackerHitImpl ;
	newHit->setPosition(  xx.array() ) ;
	outTrk->addHit( newHit );
	outColH->addElement( newHit ) ;


	std::cout << " ++++  intersection found for surface : " << surf << std::endl 
		  << "       at s = " << s << std::endl 
		  << "       xx  =  " << xx << std::endl 
		  << "    hitPos =  " << posV
 		  <<  std::endl ;
	
      } else {
	
	std::cout << " ++++ no intersection found for surface : " << surf << std::endl ;
      }

      
      //===============================================================================
    }

    wrt->writeEvent( evt ) ;
  }
    
    
  return 0;
}


#endif // USE_LCIO
