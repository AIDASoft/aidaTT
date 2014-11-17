#ifdef USE_DD4HEP
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
    if(argc < 2)
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

            int nTracks = trackCollection->getNumberOfElements();
            std::cout << " there are " << nTracks << " tracks in the event " << std::endl;

            // ignore event if more or less than a single track is present
            if(nTracks != 1)
                continue;

            Track* initialTrack = (Track*)trackCollection->getElementAt(0);

            aidaTT::trackParameters iTP;
            iTP.setTrackParameters(aidaTT::Vector5(initialTrack->getOmega(), initialTrack->getTanLambda(), initialTrack->getPhi(), initialTrack->getD0(), initialTrack->getZ0()));

            aidaTT::trajectory fitTrajectory(iTP, fitter, bfield, propagation, &geom);

            std::vector<TrackerHit*> initialHits = initialTrack->getTrackerHits();
            for(std::vector<TrackerHit*>::iterator thit = initialHits.begin(), endIter = initialHits.end(); thit < endIter; ++thit)
                {
                    long64 hitid = (*thit)->getCellID0() ;
                    idDecoder.setValue(hitid) ;
                    const aidaTT::ISurface* surf = surfMap[ hitid ] ;
                    std::cout << " surface is " << surf << std::endl;

                    double hitpos[3] = {0., 0., 0.};
                    for(unsigned int i = 0; i < 3; ++i)
                        hitpos[i] = (*thit)->getPosition()[i] * dd4hep::mm;
                    std::vector<double> precision;

                    TrackerHitPlane* planarhit = dynamic_cast<TrackerHitPlane*>(*thit);
                    if(planarhit != NULL)
                        {
                            precision.push_back(1. / planarhit->getdU());
                            precision.push_back(1. / planarhit->getdV());
                        }

                    fitTrajectory.addMeasurement(hitpos, precision, *surf, (*thit));
                }

            fitTrajectory.prepareForFitting();
            fitTrajectory.fit();

            const aidaTT::fitResults& result = fitTrajectory.getFitResults();

            std::cout << " estimated parameters after fitting are: " << result.estimatedParameters()(0) << "," << result.estimatedParameters()(1)  << "," <<
                      result.estimatedParameters()(2)  << "," <<  result.estimatedParameters()(3)  << "," << result.estimatedParameters()(4) << std:: endl;

        }


    return 0;
}


#endif // USE_LCIO
#endif // USE_DD4HEP
