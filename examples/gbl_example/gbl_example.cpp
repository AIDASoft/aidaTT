/*
 * this is an example of the usage
 * it's based on the example1.cpp in GBL by Claus Kleinwort
 */
#ifdef AIDATT_USE_DD4HEP
#ifdef USE_LCIO

// c++
#include <iostream>
#include <map>
#include <string>

// lcio
#include "lcio.h"
#include "IO/LCWriter.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/LCTrackerConf.h"

#include <IMPL/LCCollectionVec.h>
#include "IMPL/TrackImpl.h"

// DD4hep
#include "DD4hepGeometry.hh"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/SurfaceHelper.h"


// aidaTT
#include "AidaTT.hh"
#include "ConstantSolenoidBField.hh"
#include "analyticalPropagation.hh"
#include "simplifiedPropagation.hh"
#include "GBLInterface.hh"
#include "fitResults.hh"

using namespace std;
using namespace lcio;


int main(int argc, char** argv)
{
    if(argc < 2)
        {
            std::cout << " usage: ./gbl_example ILDEx.xml ILDExSimu.slcio" << std::endl ;
            return 1;
        }

    /// dd4hep stuff
    std::string inFile =  argv[1] ;

    /// preamble: load the geo info, get all surfaces => entry point for intersection calculation
    dd4hep::Detector& theDet = dd4hep::LCDD::getInstance();
    lcdd.fromCompact(inFile);

    dd4hep::DetElement world = theDet.world() ;

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

    std::string outFileName = "trackTest.slcio";

    LCWriter* wrt = LCFactory::getInstance()->createLCWriter() ;
    wrt->open(outFileName) ;

    LCEvent* evt = 0 ;

    std::vector< std::string > colNames ;
    colNames.push_back("VXDCollection") ;

    UTIL::BitField64 idDecoder(LCTrackerCellID::encoding_string()) ;




    // the aidaTT stuff

    /// create some bogus track parameters:
    ///~ test data is p = 3./sqrt(2) * (1 1 .1)^T
    /// => d0 = 0; tanLambda = 1/sqrt(2), phi0 = pi/4, z0 = 0
    // with B = 3.5T => omega = 3* 10^(-4) * 3.5 / sqrt(6)
    aidaTT::trackParameters bogusTP;

    //~ const double initialOmega = 3.5*1e-5;
    //~ const double initialTanLambda = 0.;
    //~ const double initialPhi = 5*M_PI_4;

    const double initialOmega = 1.5 * 1e-4;
    const double initialTanLambda = 1.;
    const double initialPhi = M_PI_2;

    bogusTP.setTrackParameters(aidaTT::Vector5(initialOmega, initialTanLambda, initialPhi, 0., 0.));

    // create the different objects needed for fitting
    // first a constant field parallel to z, 1T
    aidaTT::ConstantSolenoidBField*  bfield = new aidaTT::ConstantSolenoidBField(3.5);

    // create the propagation object
    aidaTT::analyticalPropagation* propagation = new aidaTT::analyticalPropagation();
    //aidaTT::simplifiedPropagation* propagation = new aidaTT::simplifiedPropagation();

    // create the fitter object
    aidaTT::GBLInterface* fitter = new aidaTT::GBLInterface();


    while((evt = rdr->readNextEvent()) != 0)
        {

            IMPL::LCCollectionVec* outputTrackCollection = new IMPL::LCCollectionVec(EVENT::LCIO::TRACK);
            std::string outName = "trytracks";

            std::cout << " reading an event " << std::endl;
            // now create a trajectory object to be fitted
            aidaTT::trajectory theMaster(bogusTP, fitter, bfield, propagation, &geom);

            for(unsigned icol = 0, ncol = colNames.size() ; icol < ncol ; ++icol)
                {
                    LCCollection* col = NULL;
                    try
                        {
                            col = evt->getCollection(colNames[ icol ]) ;
                        }
                    catch(DataNotAvailableException &e)
                        {
                            continue;
                        }

                    int nHit = col->getNumberOfElements() ;

                    for(int i = 0 ; i < nHit ; ++i)
                        {
                            SimTrackerHit* sHit = (SimTrackerHit*) col->getElementAt(i) ;
                            long64 id = sHit->getCellID0() ;
                            idDecoder.setValue(id) ;
                            double recalcPos[3] = {0., 0., 0.}; // could it be more ugly?

                            const aidaTT::ISurface* surf = surfMap[ id ] ;
                            if(surf->type().isSensitive())
                                {
                                    std::vector<double> precisionDummy;
                                    precisionDummy.push_back(1. / 0.001);
                                    precisionDummy.push_back(1. / 0.0012);
                                    for(unsigned int i = 0; i < 3; ++i)
                                        recalcPos[i] = sHit->getPosition()[i] * dd4hep::mm;

                                    theMaster.addMeasurement(recalcPos, precisionDummy, *surf, sHit);
                                }
                        }

                    theMaster.prepareForFitting();
                    theMaster.fit();

                    const aidaTT::fitResults& result = theMaster.getFitResults();

                    TrackImpl* theTrack = new TrackImpl();
                    theTrack->setOmega(result.estimatedParameters()(0)) ;
                    theTrack->setTanLambda(result.estimatedParameters()(1));
                    theTrack->setPhi(result.estimatedParameters()(2));
                    theTrack->setD0(result.estimatedParameters()(3));
                    theTrack->setZ0(result.estimatedParameters()(4));

                    TrackImpl* theInitialTrack = new TrackImpl();
                    theInitialTrack->setOmega(initialOmega) ;
                    theInitialTrack->setTanLambda(initialTanLambda);
                    theInitialTrack->setPhi(initialPhi);
                    theInitialTrack->setD0(0.);
                    theInitialTrack->setZ0(0.);


                    outputTrackCollection->addElement(theTrack);
                    outputTrackCollection->addElement(theInitialTrack);

                }
            evt->addCollection(outputTrackCollection , outName);

            wrt->writeEvent(evt);

        }


    return 0;
}


#endif // USE_LCIO
#endif // AIDATT_USE_DD4HEP
