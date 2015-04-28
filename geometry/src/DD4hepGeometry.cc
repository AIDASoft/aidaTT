#ifdef AIDATT_USE_DD4HEP

#include "DD4hepGeometry.hh"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DDRec/SurfaceHelper.h"

namespace aidaTT
{

    DD4hepGeometry::DD4hepGeometry(const DD4hep::Geometry::DetElement& det) : _detelem(det)
    {}



    const std::list<const ISurface*>& DD4hepGeometry::getSurfaces()
    {
        // create a list of all surfaces in the detector:
        DD4hep::DDRec::SurfaceHelper surfMan(_detelem) ;

        const DD4hep::DDRec::SurfaceList& sL = surfMan.surfaceList() ;

        for(DD4hep::DDRec::SurfaceList::const_iterator it = sL.begin() ; it != sL.end() ; ++it)
            {
                aidaTT::ISurface* surf =  *it ;
                _surfaceList.push_back(surf);
            }
        return _surfaceList;
    }
}

#endif // AIDATT_USE_DD4HEP
