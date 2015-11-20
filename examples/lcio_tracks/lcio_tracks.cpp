#ifdef USE_LCIO

#include "lcio.h"
#include "IO/LCReader.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/TrackerHitPlane.h"
#include "EVENT/TrackerHitZCylinder.h"
#include "EVENT/Track.h"
#include "UTIL/ILDConf.h"

#include <IMPL/LCCollectionVec.h>
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackStateImpl.h"


// aidaTT
#include "AidaTT.hh"
#include "AidaTT-Units.hh"
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

/* this is an example how to use aidaTT with LCIO data
 *
 * in the absence of a track finding procedure, this is a refitting procedure:
 *  - read in the silicon tracks
 *  - read in initial parameters
 *  - read in hits
 *  => create trajectory object from this:
 *          dd4hep surfaces, constant magnetic field, analytical propagation, gbl fitting
 */


int main(int argc, char** argv)
{
  if(argc < 3)
    {
      std::cout << " usage: ./lcio_read_example ILDEx.xml ILDExSimu.slcio" << std::endl ;
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

  //*********************************************************************

  int counter = 0 ;

  TFile *ofile = new TFile("ofile.root","RECREATE");
  //Create tree
  TTree *t1 = new TTree("t1","t1");
  vector<double> TrackHitResidualsU ;
  t1->Branch("TrackHitResidualsU",&TrackHitResidualsU);
  vector<double> TrackHitResidualsV ;
  t1->Branch("TrackHitResidualsV",&TrackHitResidualsV);
  vector<double> pullU ;
  t1->Branch("pullU",&pullU);
  vector<double> pullV ;
  t1->Branch("pullV",&pullV);
  //int VXDlayer;
  //t1->Branch("VXDlayer",&VXDlayer,"VXDlayer/I");
  vector<int> VXDlayer ;
  t1->Branch("VXDlayer",&VXDlayer);
  vector<double> pullLCIO_U ;
  t1->Branch("pullLCIO_U",&pullLCIO_U);
  vector<double> pullLCIO_V ;
  t1->Branch("pullLCIO_V",&pullLCIO_V);
  vector<double> TrkHitRes_LCIO_U ;
  t1->Branch("TrkHitRes_LCIO_U",&TrkHitRes_LCIO_U );
  vector<double> TrkHitRes_LCIO_V ;
  t1->Branch("TrkHitRes_LCIO_V",&TrkHitRes_LCIO_V );

  //*********************************************************************

  LCReader* rdr = LCFactory::getInstance()->createLCReader() ;
  rdr->open(lcioFileName) ;
  LCWriter* wrt = LCFactory::getInstance()->createLCWriter() ;

  if(argc == 4)
    {
      std::string outFile = argv[3];
      wrt->open(outFile) ;
    }
  else
    wrt->open("innowaythisisnorway.slcio", lcio::LCIO::WRITE_NEW ) ;

  LCEvent* evt = 0 ;

  //std::string trackCollectionName = "ClupatraTracks";
  std::string trackCollectionName = "MarlinTrkTracks";


  UTIL::BitField64 idDecoder(ILDCellID0::encoder_string) ;

  // create the propagation object
  aidaTT::analyticalPropagation* propagation = new aidaTT::analyticalPropagation();
  //aidaTT::simplifiedPropagation* propagation = new aidaTT::simplifiedPropagation();

  // create the fitter object
  aidaTT::GBLInterface* fitter = new aidaTT::GBLInterface();


  /// event loop
  while((evt = rdr->readNextEvent()) != 0)
    {


      TrackHitResidualsU.clear();
      TrackHitResidualsV.clear();
      pullU.clear();
      pullV.clear();
      VXDlayer.clear();
      pullLCIO_U.clear();
      pullLCIO_V.clear();
      TrkHitRes_LCIO_U.clear(); 
      TrkHitRes_LCIO_V.clear(); 
	  
      LCCollection* trackCollection = evt->getCollection(trackCollectionName) ;
	  
      // add output track collection to the event
      LCCollectionVec* outCol = new LCCollectionVec(LCIO::TRACK) ;
      evt->addCollection(outCol ,  "AidaTTTracks") ;
      LCFlagImpl trkFlag(0) ;
      trkFlag.setBit(LCIO::TRBIT_HITS) ;
      outCol->setFlag(trkFlag.getFlag()) ;
      TrackImpl* outTrk = new TrackImpl ;
      outCol->addElement(outTrk) ;


      int nTracks = trackCollection->getNumberOfElements();

      // ignore event if more or less than a single track is present
      if(nTracks != 1)
	continue;

      Track* initialTrack = (Track*)trackCollection->getElementAt(0);

      aidaTT::trackParameters iTP(  aidaTT::readLCIO( initialTrack->getTrackState( lcio::TrackState::AtIP) )   );  

      //std::cout << "  start helix from LCIO      : " << iTP << std::endl ; 

      std::vector<TrackerHit*> initialHits = initialTrack->getTrackerHits();


      //*****************************************************************************************************
      /*
      // Simple test of MoveTo function
      //std::cout << " parameters of lcioo track at IP " << iTP << std::endl ;
      aidaTT::trackParameters iTP_AtFirstHit(  aidaTT::readLCIO( initialTrack->getTrackState( lcio::TrackState::AtFirstHit) )   );
      std::cout << " parameters of lcioo track at the first hit " << iTP_AtFirstHit << std::endl ;
      aidaTT::Vector3D firstHit = iTP_AtFirstHit.referencePoint();
      moveHelixTo( iTP, firstHit, true  ) ;
     std::cout << " parameters of lcioo track moved to the first hit" << iTP << std::endl ;

     //moveHelixTo( iTP, aidaTT::Vector3D(1.,1.,1.), true  ) ;
     //std::cout << " parameters of lcioo track at random point " << iTP << std::endl ;
     //moveHelixTo( iTP, aidaTT::Vector3D(), true  ) ;
     //std::cout << " parameters of lcioo track moved to a random point and then back to IP" << iTP << std::endl ;
     */
      //****************************************************************************************************

#define compute_start_helix 0
#if compute_start_helix //----------------------------------------------------------------------------------------------------
      aidaTT::trackParameters startHelix ;

      unsigned nHits = initialHits.size() ;
      if( nHits > 2 ) {  
	//--------- get the start helix from three points
	bool backwards = false ;

	lcio::TrackerHit* h1 = ( backwards ?  initialHits[ nHits-1 ] : initialHits[    0    ] ) ;
	lcio::TrackerHit* h2 =  initialHits[ (nHits+1) / 2 ] ;
	lcio::TrackerHit* h3 = ( backwards ?  initialHits[    0    ] : initialHits[ nHits-1 ] ) ;

	aidaTT::Vector3D x1( h1->getPosition()[0] * aidaTT::mm, h1->getPosition()[1] * aidaTT::mm , h1->getPosition()[2] * aidaTT::mm ) ;
	aidaTT::Vector3D x2( h2->getPosition()[0] * aidaTT::mm, h2->getPosition()[1] * aidaTT::mm , h2->getPosition()[2] * aidaTT::mm ) ;
	aidaTT::Vector3D x3( h3->getPosition()[0] * aidaTT::mm, h3->getPosition()[1] * aidaTT::mm , h3->getPosition()[2] * aidaTT::mm ) ;
             
	calculateStartHelix( x1, x2,  x3 , startHelix , backwards ) ;
             
	moveHelixTo( startHelix, aidaTT::Vector3D(), false  ) ; // move to origin

	// --- set some large errors to the covariance matrix
	startHelix.covarianceMatrix().Unit() ;
	startHelix.covarianceMatrix()( aidaTT::OMEGA, aidaTT::OMEGA ) = 1.e-2 ;
	startHelix.covarianceMatrix()( aidaTT::TANL , aidaTT::TANL  ) = 1.e2 ;
	startHelix.covarianceMatrix()( aidaTT::PHI0 , aidaTT::PHI0  ) = 1.e2 ;
	startHelix.covarianceMatrix()( aidaTT::D0   , aidaTT::D0    ) = 1.e5 ;
	startHelix.covarianceMatrix()( aidaTT::Z0   , aidaTT::Z0    ) = 1.e5 ;

	std::cout << "  start helix from three points : " << startHelix << std::endl ;

	// use this helix as start for the fit:
	iTP = startHelix ;

      }
      /*
	#else
	// --- set some large errors to the covariance matrix
	iTP.covarianceMatrix().Unit() ;
	iTP.covarianceMatrix()( aidaTT::OMEGA, aidaTT::OMEGA ) = 1.e-2 ;
	iTP.covarianceMatrix()( aidaTT::TANL , aidaTT::TANL  ) = 1.e2 ;
	iTP.covarianceMatrix()( aidaTT::PHI0 , aidaTT::PHI0  ) = 1.e2 ;
	iTP.covarianceMatrix()( aidaTT::D0   , aidaTT::D0    ) = 1.e5 ;
	iTP.covarianceMatrix()( aidaTT::Z0   , aidaTT::Z0    ) = 1.e5 ;
      */
#endif //----------------------------------------------------------------------------------------------------------------------






      TrackStateImpl* ts;
      
      bool success;	      
      
      aidaTT::trajectory fitTrajectory(iTP, fitter, propagation, &geom);
      const aidaTT::fitResults* result = 0 ; //fitTrajectory.getFitResults();
      
      std::cout << " magnetic field at origin " 
		<< fitTrajectory.geometry()->getBField( Vector3D() ) << std::endl ;
      
      
	    
      for(std::vector<TrackerHit*>::iterator thit = initialHits.begin(), endIter = initialHits.end(); thit < endIter; ++thit)
	{
	  long hitid = (*thit)->getCellID0() ;
	  idDecoder.setValue(hitid) ;

	  if(idDecoder[ lcio::ILDCellID0::subdet] != lcio::ILDDetID::VXD && idDecoder[ lcio::ILDCellID0::subdet] != lcio::ILDDetID::TPC)
	    continue;
	  
	  //std::cout << " hit id = " << hitid << std::endl ;
	  
	  const aidaTT::ISurface* surf = surfMap[ hitid ] ;
	  
	  if(surf == NULL)
	    {
	      std::cerr << " lcio_tracks : no surface found for id : " << idDecoder.valueString() << std::endl ;
	      continue;
	    }
	  
	  
	  double hitpos[3] = {0., 0., 0.};
	  for(unsigned int i = 0; i < 3; ++i)
	    hitpos[i] = (*thit)->getPosition()[i] * aidaTT::mm;
	  
	  //std::cout << " hit position X " << hitpos[0] << " hit position Y " << hitpos[1] << " hit position Z " << hitpos[2] << std::endl ;
	  
	  std::vector<double> precision ;
	  
	  TrackerHitPlane* planarhit = dynamic_cast<TrackerHitPlane*>(*thit);
	  
	  FloatVec TPChitCovMat = (*thit)->getCovMatrix();
	  
	  double du, dv ;
	  
	  if (planarhit != NULL) {
	    
	    du = planarhit->getdU() * aidaTT::mm  ;
	    dv = planarhit->getdV() * aidaTT::mm  ;
	    
	  } else {
	    
	    du = sqrt( TPChitCovMat[0] + TPChitCovMat[2]) * aidaTT::mm;
	    dv = sqrt( TPChitCovMat[5] ) * aidaTT::mm;
	  }
	  
	  std::cout << " u resolution " << du << " v resolution " << dv << std::endl; 

	  //we need 1./variance for the precision:
	  //what are the values for resolutiuon?
	  
	  precision.push_back( 1. /  (du*du)) ;
	  precision.push_back( 1. /  (dv*dv)) ;
	  
	  fitTrajectory.addMeasurement(hitpos, precision, *surf, (*thit), true);
	  
	  outTrk->addHit(*thit) ;
	  
	}
      
      
      for (int n=0; n < 1 ; n++){
	
	fitTrajectory.prepareForFitting();
	
	success = fitTrajectory.fit();
	
	result = fitTrajectory.getFitResults();


	//***********************************************************************************************************
	/*
	  std::cout << " loop " << n << std::endl ;
	  std::cout << " initial values " << std::endl;
	  std::cout << iTP << std::endl;
	  std::cout << " refitted values " << std::endl;
	  std::cout << result->estimatedParameters() << std::endl;
	*/
	//iTP = result->estimatedParameters();  // valid only when we make an iterative fitting 
	      

	if( ! success ) {
		
	  std::cout << " ********** ERROR:  Fit Failed !!!!! ******************************* " << std::endl ;
	}
	      
	      
      }

	         

      std::cout << " End of the loop " << std::endl ;
      std::cout << " initial values " << std::endl;
      std::cout << iTP << std::endl;
      std::cout << " refitted values " << std::endl;
      std::cout << result->estimatedParameters() << std::endl;
	    
      // add Track State to track:
      ts = aidaTT::createLCIO( result->estimatedParameters() );
	      
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
	
      wrt->writeEvent(evt) ;
    }

  ofile->Write("t1");

  std::cout << " counter = " << counter << std::endl ;

  return 0;
}

#endif // USE_LCIO
