#ifdef AIDATT_USE_DD4HEP

#include <iostream>

// aidaTT
#include "IGeometry.hh"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/SurfaceHelper.h"

int main(int argc, char** argv)
{

    if(argc < 2)
        {
            std::cout << " usage: ./geometry_test ILDEx.xml " << std::endl ;
            exit(1) ;
        }

    std::string inFile =  argv[1] ;

    dd4hep::Detector& theDet = dd4hep::Detector::getInstance();
    theDet.fromCompact(inFile);

    dd4hep::DetElement world = theDet.world() ;

    // create a list of all surfaces in the detector:
    dd4hep::rec::SurfaceHelper surfMan(world) ;

    const dd4hep::rec::SurfaceList& sL = surfMan.surfaceList() ;

    for(dd4hep::rec::SurfaceList::const_iterator it = sL.begin() ; it != sL.end() ; ++it)
        {

            aidaTT::ISurface* surf =  *it ;
            std::cout << " Id of surface: " << surf->id() << std::endl;
        }


    return 0;
}

#endif // AIDATT_USE_DD4HEP
