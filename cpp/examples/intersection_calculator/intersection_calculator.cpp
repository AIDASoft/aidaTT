#ifdef USE_DD4HEP

#include <iostream>

// aidaTT
#include "DD4hepGeometry.hh"
#include "trackParameters.hh"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DDRec/SurfaceManager.h"

int main(int argc, char** argv)
{

    if(argc < 2)
        {
            std::cout << " usage: geometry_test compact.xml" << std::endl ;
            exit(1) ;
        }

    std::string inFile =  argv[1] ;

    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
    lcdd.fromCompact(inFile);

    DD4hep::Geometry::DetElement world = lcdd.world() ;

    aidaTT::DD4hepGeometry geom(world);

    const std::list<const aidaTT::ISurface*>& surfaces = geom.getSurfaces() ;

    //~ aidaTT::trackParameters tp;
    //~
    //~ double vals[] = { -.12, 1.3, 2.4, 0.001, 0.43 };
    //~ std::vector<double> values;
    //~ values.assign(vals, vals+5);
//~
    //~ tp.setTrackParameters(aidaTT::Vector5(values));
//~
    for(std::list<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf)
        {
            std::cout << "surface with id: " << (*surf)->id() << " is of kind " << std::endl;
            std::cout << "           type: " << (*surf)->type() << std::endl;


            // getIntersectionWithSurface(
        }


    return 0;
}

#endif // USE_DD4HEP
