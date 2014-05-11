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

    double vals[] = { -.12, 1.3, 2.4, 0.001, 0.43 };
    std::vector<double> values;
    values.assign(vals, vals + 5);

    aidaTT::Vector5 tp(values);
    //    tp.setTrackParameters(aidaTT::Vector5(values));


    std::cout << std::endl << std::endl << " there are " << surfaces.size() << " surfaces to consider " << std::endl;
    aidaTT::Vector3D bfield(0., 0., 1.);

    aidaTT::fiveByFiveMatrix cl2per = aidaTT::curvilinearToPerigeeJacobian(tp, bfield);
    aidaTT::Vector5 perTP = cl2per * tp;
    aidaTT::fiveByFiveMatrix per2ILD = aidaTT::perigeeToILDJacobian(perTP);
    aidaTT::Vector5 ildTP = per2ILD * perTP;

    std::cout << " values of the tp are:         [q/p, lambda, phi, x_t, y_t] = [ " << tp(0) << " , " << tp(1) << " , " << tp(2) << " , " << tp(3) << " , " << tp(4) << " ] " << std::endl;
    std::cout << " values of the pertp are now: [kappa, theta, phi, eps,zp ] = [ " << perTP(0) << " , " << perTP(1) << " , " << perTP(2) << " , " << perTP(3) << " , " << perTP(4) << " ] " << std::endl;
    std::cout << " values of the ildtp are now:  [Omega, tanlambda, phi0, d0, z0] = [ " << ildTP(0) << " , " << ildTP(1) << " , " << ildTP(2) << " , " << ildTP(3) << " , " << ildTP(4) << " ] " << std::endl;

    const double omega = ildTP(0);
    const double phi0  = ildTP(2);
    const double dzero = ildTP(3);
    const double xzero = dzero * cos(phi0);
    const double yzero = dzero * sin(phi0);
    aidaTT::trackParameters bla;
    bla.setTrackParameters(ildTP);
    aidaTT::trajectory tryout(bla, NULL);

    //getIntersectionsWithSurfaces(surfaces)
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
                      //~ double dp = (*surf)->distance(DDSurfaces::Vector3D(0., 0., 0.));
                      //~ std::cout << " the normal vector has the components nx=" << normalX << " ny=" << normalY << " nz=" << normalZ << std::endl;
                      //~ std::cout << " the atan2 of nx,ny is " <<  pv << std::endl;
                      //~ std::cout << " dp is " <<  dp << std::endl;
                      //~ std::cout << " it should be " << (*surf)->distance(DDSurfaces::Vector3D(0., 0., 0.)) << std::endl;
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
