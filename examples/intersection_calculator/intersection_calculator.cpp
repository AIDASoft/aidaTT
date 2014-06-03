#ifdef USE_DD4HEP

#include <iostream>

// aidaTT
#include "DD4hepGeometry.hh"
#include "trackParameters.hh"
#include "trajectory.hh"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DDRec/SurfaceManager.h"

int main(int argc, char** argv)
{
    if(argc < 2)
        {
            std::cout << " usage: ./intersection_calculator ILDEx.xml " << std::endl ;
            exit(1) ;
        }

    std::string inFile =  argv[1] ;

    /// preamble: load the geo info, get all surfaces => entry point for intersection calculation
    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
    lcdd.fromCompact(inFile);

    DD4hep::Geometry::DetElement world = lcdd.world() ;

    aidaTT::DD4hepGeometry geom(world);

    const std::list<const aidaTT::ISurface*>& surfaces = geom.getSurfaces() ;


    //~ a
    double vals[] = { 0.23, 0.01, 1.57, 0.00, 0. };
    std::vector<double> values;
    values.assign(vals, vals + 5);

    aidaTT::Vector5 ildTP(values);

    std::cout << " init vector is " << ildTP << std::endl;


    aidaTT::trackParameters bla;
    bla.setTrackParameters(ildTP);

    std::cout << " parameters are " << bla << std::endl;

    aidaTT::trajectory tryout(bla, NULL);

    const std::vector<std::pair<double, const aidaTT::ISurface*> >& whatever =     tryout.getIntersectionsWithSurfaces(surfaces);

    std::cout << " length of intersections " << whatever.size() << std::endl;

    for(std::vector<std::pair<double, const aidaTT::ISurface*> >::const_iterator test = whatever.begin(), fin = whatever.end(); test < fin; ++test)
        {
            std::cout << " *** at arc length: " << (*test).first << " found surface " << *(*test).second << std::endl;
        }

    /*
     *   for(std::list<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf)
          {

              if((*surf)->type().isZCylinder())
                  {
                      std::cout << " its a cylinder parallel to z " << std::endl;
                  }
              else if((*surf)->type().isZPlane())
                  {
                      double normalX = (*surf)->normal().x();
                      double normalY = (*surf)->normal().y();
                      double normalZ = (*surf)->normal().z();
                      double xp = (*surf)->origin().x();
                      double yp = (*surf)->origin().y();

                      //~ double pv = atan2(normalX, normalY);
                      //~ double dp = (*surf)->distance(Vector3D(0., 0., 0.));
                      //~ std::cout << " the normal vector has the components nx=" << normalX << " ny=" << normalY << " nz=" << normalZ << std::endl;
                      //~ std::cout << " the atan2 of nx,ny is " <<  pv << std::endl;
                      //~ std::cout << " dp is " <<  dp << std::endl;
                      //~ std::cout << " it should be " << (*surf)->distance(Vector3D(0., 0., 0.)) << std::endl;
                      //~ std::cout << " phi0 is " << phi0 << std::endl;
                      //~ std::cout << " sin( pv - phi0) is " << sin(pv - phi0) << std::endl;
                      //~ std::cout << " omega*dp is " <<  omega*dp << std::endl;
                      //~ double omegaS = asin(sin(pv - phi0) + omega * dp) - pv + phi0;
                      //~ std::cout << " omegaS is " << omegaS << std::endl;

                  }
              else if((*surf)->type().isZDisk())
                  {
                      std::cout << " its a plane perpendicular to z " << std::endl;
                      double planePositionZ = (*surf)->origin().z();
                      std::cout << " plane is at " << planePositionZ << std::endl;
                  }
              else
                  {
                      std::cout << " No good, can't compute intersection. " << std::endl;
                  }
              //~ if( (*surf)->type().isSensitive() )

              // getIntersectionWithSurface
          }
    */

    return 0;
}

#endif // USE_DD4HEP
