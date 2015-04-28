#ifdef AIDATT_USE_DD4HEP

#include <iostream>

// aidaTT
#include "IGeometry.hh"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DDRec/SurfaceHelper.h"

int main(int argc, char** argv)
{

    if(argc < 2)
        {
            std::cout << " usage: ./geometry_test ILDEx.xml " << std::endl ;
            exit(1) ;
        }

    std::string inFile =  argv[1] ;

    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
    lcdd.fromCompact(inFile);

    DD4hep::Geometry::DetElement world = lcdd.world() ;

    // create a list of all surfaces in the detector:
    DD4hep::DDRec::SurfaceHelper surfMan(world) ;

    const DD4hep::DDRec::SurfaceList& sL = surfMan.surfaceList() ;

    for(DD4hep::DDRec::SurfaceList::const_iterator it = sL.begin() ; it != sL.end() ; ++it)
        {

            aidaTT::ISurface* surf =  *it ;
            std::cout << " Id of surface: " << surf->id() << std::endl;
        }


    return 0;
}

#endif // AIDATT_USE_DD4HEP
