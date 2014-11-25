#ifdef AIDATT_USE_DD4HEP

#ifndef DD4HEPGEOMETRY_HH
#define DD4HEPGEOMETRY_HH

#include "IGeometry.hh"
#include "DD4hep/LCDD.h"

namespace aidaTT
{
    class DD4hepGeometry : public IGeometry
    {
        public:
            DD4hepGeometry(const ::DD4hep::Geometry::DetElement&);

            const std::list<const ISurface*>& getSurfaces();

        private:
            const ::DD4hep::Geometry::DetElement& _detelem;
            std::list<const ISurface* > _surfaceList;
    };
}
#endif // DD4HEPGEOMETRY_HH

#endif // AIDATT_USE_DD4HEP
