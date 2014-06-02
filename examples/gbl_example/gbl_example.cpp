/*
 * this is an example of the usage
 * it's based on the example1.cpp in GBL by Claus Kleinwort
 */
#ifdef USE_DD4HEP
#ifdef USE_LCIO

// c++
#include <iostream>
#include <map>

// lcio
#include "lcio.h"
#include "IO/LCReader.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/ILDConf.h"


// DD4hep
#include "DD4hepGeometry.hh"
#include "DD4hep/LCDD.h"
#include "DDRec/SurfaceManager.h"


// aidaTT
#include "AidaTT.hh"
#include "ConstantSolenoidBField.hh"
#include "analyticalPropagation.hh"
#include "GBLInterface.hh"

using namespace std;
using namespace lcio;


int main(int argc, char** argv)
{
    if(argc < 2)
        {
            std::cout << " usage: ./example compact.xml lcio.slcio" << std::endl ;
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

    std::vector< std::string > colNames ;
    colNames.push_back("VXDCollection") ;

    UTIL::BitField64 idDecoder(ILDCellID0::encoder_string) ;




    // the aidaTT stuff

    /// create some bogus track parameters;
    aidaTT::trackParameters bogusTP;

    // create the different objects needed for fitting
    // first a constant field parallel to z, 1T
    aidaTT::ConstantSolenoidBField*  bfield = new aidaTT::ConstantSolenoidBField(1.);

    // create the propagation object
    aidaTT::analyticalPropagation* propagation = new aidaTT::analyticalPropagation();

    // create the fitter object
    aidaTT::GBLInterface* fitter = new aidaTT::GBLInterface();

    // now create a list of trajectory objects for the fitter
    aidaTT::trajectory theMaster(bogusTP, fitter, bfield, propagation, &geom);

    while((evt = rdr->readNextEvent()) != 0)
        {
            for(unsigned icol = 0, ncol = colNames.size() ; icol < ncol ; ++icol)
                {
                    LCCollection* col = evt->getCollection(colNames[ icol ]) ;
                    int nHit = col->getNumberOfElements() ;

                    for(int i = 0 ; i < nHit ; ++i)
                        {
                            SimTrackerHit* sHit = (SimTrackerHit*) col->getElementAt(i) ;
                            long64 id = sHit->getCellID0() ;
                            idDecoder.setValue(id) ;

                            const aidaTT::ISurface* surf = surfMap[ id ] ;
                            if(surf->type().isSensitive())
                                {
                                    std::cout << " thank the lord " << std::endl;
                                    std::vector<double> resolutionDummy;
                                    resolutionDummy.push_back(0.01);
                                    resolutionDummy.push_back(0.12);
                                    theMaster.addMeasurement(sHit->getPosition(), resolutionDummy, *surf, sHit);
                                }
                        }
                    theMaster.prepareForFitting();
                    theMaster.fit();
                }
        }


    return 0;
}


#endif // USE_LCIO
#endif // USE_DD4HEP
