#ifdef AIDATT_USE_DD4HEP
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

// DD4hep
#include "DD4hepGeometry.hh"
#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/SurfaceManager.h"


// aidaTT
#include "AidaTT.hh"
#include "ConstantSolenoidBField.hh"
#include "analyticalPropagation.hh"
#include "simplifiedPropagation.hh"
#include "GBLInterface.hh"
#include "fitResults.hh"
#include "Vector5.hh"
#include "utilities.hh"
#include "LCIOPersistency.hh"

#include <map>

using namespace std ;
using namespace lcio;

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

    /// dd4hep stuff
    std::string inFile =  argv[1] ;

    /// preamble: load the geo info, get all surfaces => entry point for intersection calculation
    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
    lcdd.fromCompact(inFile);

    DD4hep::Geometry::DetElement world = lcdd.world() ;

    aidaTT::DD4hepGeometry geom(world);

    const std::list<const aidaTT::ISurface*>& surfaces = geom.getSurfaces() ;

    // create map of surfaces
    std::map< long64, const aidaTT::ISurface* > surfMap ;

    for(std::list<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf)
        {
            surfMap[(*surf)->id() ] = (*surf) ;
        }

    /// lcio stuff
    std::string lcioFileName = argv[2] ;

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

    std::string trackCollectionName = "SiTracks";

    UTIL::BitField64 idDecoder(ILDCellID0::encoder_string) ;

    // create the different objects needed for fitting
    // first a constant field parallel to z, 1T
    aidaTT::ConstantSolenoidBField*  bfield = new aidaTT::ConstantSolenoidBField(3.5);

    // create the propagation object
    aidaTT::analyticalPropagation* propagation = new aidaTT::analyticalPropagation();
    //aidaTT::simplifiedPropagation* propagation = new aidaTT::simplifiedPropagation();

    // create the fitter object
    aidaTT::GBLInterface* fitter = new aidaTT::GBLInterface();


    /// event loop
    while((evt = rdr->readNextEvent()) != 0)
        {

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

            aidaTT::trackParameters iTP(  aidaTT::readLCIO( initialTrack->getTrackState( lcio::TrackState::AtIP) )    );  

	    std::cout << "  start helix from LCIO      : " << iTP << std::endl ; 

            std::vector<TrackerHit*> initialHits = initialTrack->getTrackerHits();

#define compute_start_helix 1
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
	      
	      moveHelixTo( startHelix, aidaTT::Vector3D()  ) ; // move to origin

	      // --- set some large errors to the covariance matrix
	      startHelix.covarianceMatrix().Unit() ;
	      startHelix.covarianceMatrix()( aidaTT::OMEGA, aidaTT::OMEGA ) = 1.e-2 ;
	      startHelix.covarianceMatrix()( aidaTT::TANL , aidaTT::TANL  ) = 1.e2 ;
	      startHelix.covarianceMatrix()( aidaTT::PHI0 , aidaTT::PHI0  ) = 1.e2 ;
	      startHelix.covarianceMatrix()( aidaTT::D0   , aidaTT::D0    ) = 1.e5 ;
	      startHelix.covarianceMatrix()( aidaTT::Z0   , aidaTT::Z0    ) = 1.e5 ;

	      // std::cout << "  start helix from three points : " << startHelix << std::endl ;

	      // use this helix as start for the fit:
	      iTP = startHelix ;

	    }
#else
	    // --- set some large errors to the covariance matrix
	    iTP.covarianceMatrix().Unit() ;
	    iTP.covarianceMatrix()( aidaTT::OMEGA, aidaTT::OMEGA ) = 1.e-2 ;
	    iTP.covarianceMatrix()( aidaTT::TANL , aidaTT::TANL  ) = 1.e2 ;
	    iTP.covarianceMatrix()( aidaTT::PHI0 , aidaTT::PHI0  ) = 1.e2 ;
	    iTP.covarianceMatrix()( aidaTT::D0   , aidaTT::D0    ) = 1.e5 ;
	    iTP.covarianceMatrix()( aidaTT::Z0   , aidaTT::Z0    ) = 1.e5 ;
	    
#endif //----------------------------------------------------------------------------------------------------------------------


            aidaTT::trajectory fitTrajectory(iTP, fitter, bfield, propagation, &geom);

            for(std::vector<TrackerHit*>::iterator thit = initialHits.begin(), endIter = initialHits.end(); thit < endIter; ++thit)
                {
                    long64 hitid = (*thit)->getCellID0() ;
                    idDecoder.setValue(hitid) ;

                    if(idDecoder[ lcio::ILDCellID0::subdet] != lcio::ILDDetID::VXD)
                        continue;

                    if(idDecoder[ lcio::ILDCellID0::subdet] == lcio::ILDDetID::VXD)
                        {
                            idDecoder[lcio::ILDCellID0::side] = ((*thit)->getPosition()[2]  >  0  ?   +1 : -1) ;

                            hitid = idDecoder.lowWord() ;
                        }

                    const aidaTT::ISurface* surf = surfMap[ hitid ] ;

                    if(surf == NULL)
                        {
                            std::cerr << " lcio_tracks : no surface found for id : " << idDecoder.valueString() << std::endl ;
                            continue;
                        }

                    double hitpos[3] = {0., 0., 0.};
                    for(unsigned int i = 0; i < 3; ++i)
                        hitpos[i] = (*thit)->getPosition()[i] * dd4hep::mm;

                    std::vector<double> precision;

                    TrackerHitPlane* planarhit = dynamic_cast<TrackerHitPlane*>(*thit);
                    if(planarhit != NULL)
		      {
			//we need 1./variance for the precision: 
			double du = planarhit->getdU() * dd4hep::mm  ;
			double dv = planarhit->getdV() * dd4hep::mm  ;

			precision.push_back( 1. /  (du*du) ) ;
			precision.push_back( 1. /  (dv*dv) ) ;
			
		      }

                    fitTrajectory.addMeasurement(hitpos, precision, *surf, (*thit));

                    outTrk->addHit(*thit) ;

                }
            fitTrajectory.prepareForFitting();


            bool success = fitTrajectory.fit();
            const aidaTT::fitResults& result = fitTrajectory.getFitResults();


	    if( ! success ) {

	      std::cout << " ********** ERROR:  Fit Failed !!!!! ******************************* " << std::endl ;
	    }

            std::cout << " initial values vs. refitted: " << std::endl;
            std::cout << iTP << std::endl;
            std::cout << result.estimatedParameters() << std::endl;


            // add Track State to track:
            TrackStateImpl* ts = aidaTT::createLCIO( result.estimatedParameters()  );

	    outTrk->setChi2( result.chiSquare() ) ;
	    outTrk->setNdf( result.ndf() ) ;
	    outTrk->subdetectorHitNumbers().resize(10.) ;

	    outTrk->subdetectorHitNumbers()[0] = outTrk->getTrackerHits().size() ;


            float ref[3] = { 0., 0. , 0. } ;
            ts->setReferencePoint(ref);

            ts->setLocation(lcio::TrackState::AtIP);

            outTrk->addTrackState(ts);

            wrt->writeEvent(evt) ;
        }


    return 0;
}


#endif // USE_LCIO
#endif // USE_DD4HEP
