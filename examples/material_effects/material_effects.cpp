#ifdef USE_LCIO

#include "lcio.h"
#include "IO/LCReader.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/TrackerHitPlane.h"
#include "EVENT/Track.h"
#include "UTIL/ILDConf.h"

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


/* this is an example how to use aidaTT with LCIO data when material effects are taken into account
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

    /// dd4hep stuff
    std::string inFile =  argv[1] ;

  const aidaTT::IGeometry& geom = aidaTT::IGeometry::instance(inFile) ;

  const std::vector<const aidaTT::ISurface*>& surfaces = geom.getSurfaces() ;

  // create map of surfaces
  std::map< long, const aidaTT::ISurface* > surfMap ;
  
  for(std::vector<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf){
    surfMap[(*surf)->id() ] = (*surf) ;
  }


    //create map of hits
    std::map< long, EVENT::TrackerHit* > hitMap ;

    /// lcio stuff
    std::string lcioFileName = argv[2] ;

    //*********************************************************************

    char* ofile_name = argv[4] ;

    int counter = 0 ;

    TFile *ofile = new TFile(ofile_name, "RECREATE");
    //Create tree
    TTree *t1 = new TTree("t1", "t1");
    vector<double> TrackHitResidualsU ;
    t1->Branch("TrackHitResidualsU", &TrackHitResidualsU);
    vector<double> TrackHitResidualsV ;
    t1->Branch("TrackHitResidualsV", &TrackHitResidualsV);
    vector<double> pullU ;
    t1->Branch("pullU", &pullU);
    vector<double> pullV ;
    t1->Branch("pullV", &pullV);
    //int VXDlayer;
    //t1->Branch("VXDlayer",&VXDlayer,"VXDlayer/I");
    vector<int> VXDlayer ;
    t1->Branch("VXDlayer", &VXDlayer);
    vector<double> pullLCIO_U ;
    t1->Branch("pullLCIO_U", &pullLCIO_U);
    vector<double> pullLCIO_V ;
    t1->Branch("pullLCIO_V", &pullLCIO_V);
    vector<double> TrackHitResidualsU_LCIO ;
    t1->Branch("TrackHitResidualsU_LCIO", &TrackHitResidualsU_LCIO);
    vector<double> TrackHitResidualsV_LCIO ;
    t1->Branch("TrackHitResidualsV_LCIO", &TrackHitResidualsV_LCIO);


    //*********************************************************************

    LCReader* rdr = LCFactory::getInstance()->createLCReader() ;
    rdr->open(lcioFileName) ;
    LCWriter* wrt = LCFactory::getInstance()->createLCWriter() ;

    if(argc > 3)
        {
            std::string outFile = argv[3];
            wrt->open(outFile) ;
        }
    else
        wrt->open("innowaythisisnorway.slcio", lcio::LCIO::WRITE_NEW) ;

    LCEvent* evt = 0 ;

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
	    TrackHitResidualsU_LCIO.clear();
            TrackHitResidualsV_LCIO.clear();

            LCCollection* trackCollection = evt->getCollection(trackCollectionName) ;

            // add output track collection to the event
            LCCollectionVec* outCol = new LCCollectionVec(LCIO::TRACK) ;
            evt->addCollection(outCol ,  "AidaTTTracks") ;
            LCFlagImpl trkFlag(0) ;
            trkFlag.setBit(LCIO::TRBIT_HITS) ;
            outCol->setFlag(trkFlag.getFlag()) ;
            TrackImpl* outTrk = new TrackImpl ;
            outCol->addElement(outTrk) ;
	    
	    std::vector < const aidaTT::ISurface* > measSurfs ;
	    
	    std::cout << " new event " << std::endl;

            int nTracks = trackCollection->getNumberOfElements();

            // ignore event if more or less than a single track is present
            if(nTracks != 1)
                continue;

            Track* initialTrack = (Track*)trackCollection->getElementAt(0);

            aidaTT::trackParameters iTP(aidaTT::readLCIO(initialTrack->getTrackState(lcio::TrackState::AtIP)));

            std::vector<TrackerHit*> initialHits = initialTrack->getTrackerHits();

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

	aidaTT::Vector3D x1( h1->getPosition()[0] * dd4hep::mm, h1->getPosition()[1] * dd4hep::mm , h1->getPosition()[2] * dd4hep::mm ) ;
	aidaTT::Vector3D x2( h2->getPosition()[0] * dd4hep::mm, h2->getPosition()[1] * dd4hep::mm , h2->getPosition()[2] * dd4hep::mm ) ;
	aidaTT::Vector3D x3( h3->getPosition()[0] * dd4hep::mm, h3->getPosition()[1] * dd4hep::mm , h3->getPosition()[2] * dd4hep::mm ) ;
             
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

#endif //----------------------------------------------------------------------------------------------------------------------



            TrackStateImpl* ts;

            bool success;

            aidaTT::trajectory fitTrajectory(iTP, fitter, propagation, &geom);
            //const aidaTT::fitResults* result = fitTrajectory.getFitResults();
	    const aidaTT::fitResults* result;



	    for(std::vector<TrackerHit*>::iterator thit = initialHits.begin(), endIter = initialHits.end(); thit < endIter; ++thit) {

	      //long hitid_dt = (*thit)->getCellID0() ;
	      //idDecoder.setValue(hitid_dt) ;
	      //if(idDecoder[ lcio::ILDCellID0::subdet] != lcio::ILDDetID::VXD && idDecoder[ lcio::ILDCellID0::subdet] != lcio::ILDDetID::TPC)
	      //continue;

	      if( !BitSet32( (*thit)->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ 
		long hitid = (*thit)->getCellID0() ;
		//std::cout << "  I am mapping hit " << (*thit) << " to surface " << hitid << std::endl ; 
		hitMap[hitid] = (*thit);
	      }

	      if( BitSet32( (*thit)->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
		//Split it up to 2 1-D hits
		const LCObjectVec rawObjects = (*thit)->getRawHits();			  
		for( unsigned k=0; k< rawObjects.size(); k++ ){
		  TrackerHitPlane* rawHit = dynamic_cast< TrackerHitPlane* >( rawObjects[k] );
		  if(rawHit != NULL){
		    long hitid = rawHit->getCellID0() ;
		    //std::cout << "  I am mapping hit " << rawHit << " to surface " << hitid << std::endl ; 
		    hitMap[hitid] = rawHit;
		  }
		}
	      }
	    }		
	    

	    std::vector<std::pair<double, const ISurface*> > CrossedSurfs = fitTrajectory.getIntersectionsWithSurfaces(surfaces);
	    
	    for(std::vector<std::pair<double, const aidaTT::ISurface*> >::iterator test = CrossedSurfs.begin(), fin = CrossedSurfs.end(); test < fin; ++test){
	      
	      std::cout << " *** at arc length: " << (*test).first << " found surface " << *(*test).second << std::endl;

	      //std::cout << " id of the surface " << (*test).second->id() << std::endl;
	      
	      if (hitMap[(*test).second->id()] !=NULL){

		EVENT::TrackerHit* testHit = hitMap[(*test).second->id()];
		//std::cout << " I am reading hit " << testHit << " from map " << std::endl ;
	
		double du;
		double dv;

		std::vector<double> precision;
		const aidaTT::Vector3D hitpos(testHit->getPosition()[0] * dd4hep::mm, testHit->getPosition()[1] * dd4hep::mm, testHit->getPosition()[2] * dd4hep::mm );		
		
		long ID = (*test).second->id();
		  
	  
		TrackerHitPlane* planarhit = dynamic_cast<TrackerHitPlane*>(testHit);
		if(planarhit != NULL){
		  //we need 1./variance for the precision:
		  du = planarhit->getdU() * dd4hep::mm  ;
		  dv = planarhit->getdV() * dd4hep::mm  ;
		}
		  
		  
		else {
		  const FloatVec& TPChitCovMat = testHit->getCovMatrix();
		  du = sqrt( TPChitCovMat[0] + TPChitCovMat[2]) * dd4hep::mm;
		  dv = sqrt( TPChitCovMat[5] ) * dd4hep::mm;
		}

		std::cout << " u resolution " << du << " v resolution " << dv << std::endl;

		if (dv!=0){		  
		  precision.push_back(1. / (du * du)) ;
		  precision.push_back(1. / (dv * dv)) ;
		}
		else 
		  precision.push_back(1. / (du * du)) ;

		outTrk->addHit(testHit) ;

		fitTrajectory.addMeasurement(hitpos, precision, *(*test).second, &ID, true);		
	      }
	      
		// ignore virtual surface with no material (e.g. inside the beam pipe )
	      else{
		if( ! ( (*test).second->innerMaterial().density() < 1e-6  &&
		      (*test).second->outerMaterial().density() < 1e-6 )  ) {
		  
		  fitTrajectory.addScatterer( *(*test).second ) ;
		}
	      }
	    }
	    
	    
	    fitTrajectory.prepareForFitting();
	    
	    
	    const std::vector<trajectoryElement*>& elements = fitTrajectory.trajectoryElements();
	    
	    //std::cout << " number of elements associated to the trajectory " << elements.size() << std::endl;
	    
	    success = fitTrajectory.fit();
	    
	    result = fitTrajectory.getFitResults();
	    
	    
	    //~ std::cout << " loop " << n << std::endl ;
	    //~ std::cout << " refitted values " << std::endl;
	    //~ std::cout << result->estimatedParameters() << std::endl;
	    //~
	    //iTP = result->estimatedParameters();  // valid only when we make an iterative fitting
	    
	    
	    if(! success)
	      {
		
		std::cout << " ********** ERROR:  Fit Failed !!!!! ******************************* " << std::endl ;
	      }
	    
	    if ( success ) { std::cout << " successful fit " << std::endl ; }
	    
	    
	    
	    
	    
	    
	    
            //~ std::cout << " initial values " << std::endl;
            //~ std::cout << iTP << std::endl;
            //~ std::cout << " refitted values " << std::endl;
            //~ std::cout << result->estimatedParameters() << std::endl;
            //~
            // add Track State to track:
            ts = aidaTT::createLCIO(result->estimatedParameters());
            //ts = aidaTT::createLCIO( iTP );  // only to check the initial helix

            outTrk->setChi2(result->chiSquare()) ;
            outTrk->setNdf(result->ndf()) ;
            outTrk->subdetectorHitNumbers().resize(10.) ;

            outTrk->subdetectorHitNumbers()[0] = outTrk->getTrackerHits().size() ;

            float ref[3] = { 0., 0. , 0. } ;
            ts->setReferencePoint(ref);

            ts->setLocation(lcio::TrackState::AtIP);


            outTrk->addTrackState(ts);

            wrt->writeEvent(evt) ;

	    hitMap.clear();

        }

    ofile->Write("t1");

    return 0;
}


#endif // USE_LCIO
