#ifndef IGEOMETRY_HH
#define IGEOMETRY_HH

#include "Vector3D.hh"
#include <list>

/* This hides the actual implementation of geometry information.
 * Geometry is needed for the calculation of
 *  - crossing points between elements and a track
 *  - measurement directions (if a sensitive element)
 *  - material properties for multiple scattering and enery loss
 *
 * There are several elements to the geometry:
 *  - the abstract material definition
 *  - the abstract surface definition
 *
 *  @author: Christoph Rosemann, DESY and Christian Grefe, CERN
 */

#ifdef USE_DD4HEP
#include "DDSurfaces/ISurface.h"
#include "DDSurfaces/IMaterial.h"

namespace aidaTT
{
    typedef DDSurfaces::ISurface ISurface;
    typedef DDSurfaces::IMaterial IMaterial;
}
#else
namespace aidaTT
{

    class IMaterial
    {

        public:
            /// averaged proton number
            virtual double getZ() const = 0 ;

            /// averaged atomic number
            virtual double getA() const = 0 ;

            /// density - units ?
            virtual double getDensity() const = 0 ;

            /// radiation length - units ?
            virtual double getRadiationLength() const = 0 ;

            /// interaction length - units ?
            virtual double getInteractionLength() const = 0 ;

    };


    class ISurface
    {
            /// The id of a surface
            virtual long long int id() const = 0 ;

            ///~ Checks if the given point in global coordinates lies within the surface (using the given accuracy)
            virtual bool insideBounds(const Vector3D& point, double epsilon = 1.e-4) const = 0 ;

            ///~ First direction of measurement U
            virtual const Vector3D& getU(const Vector3D& point = Vector3D()) const = 0 ;

            ///~ Second direction of measurement V
            virtual const Vector3D& getV(const Vector3D& point = Vector3D()) const = 0 ;

            /// Access to the normal direction at the given point
            virtual const Vector3D& getNormal(const Vector3D& point = Vector3D()) const = 0 ;

            //~ /** Get Origin of local coordinate system on surface */
            //~ virtual const Vector3D& origin() const =0 ;

            /// Access to the material in opposite direction of the normal
            virtual const IMaterial& getInnerMaterial() const = 0 ;

            /// Access to the material in direction of the normal
            virtual const IMaterial& getOuterMaterial() const = 0 ;

            /** Thickness of inner material */
            virtual double innerThickness() const = 0 ;

            /** Thickness of outer material */
            virtual double outerThickness() const = 0 ;

            //~ /** Distance to surface */
            //~ virtual double distance(const Vector3D& point ) const =0 ;

            /// true if surface is sensitive
            virtual bool isSensitive() const = 0;

            /// true if this a cylindrical surface
            virtual bool isCylinder() const = 0;

            /// true if surface is parallel to Z
            virtual bool isParallelToZ() const = 0;

            /// true if surface is orthogonal to Z
            virtual bool isOrthogonalToZ() const = 0;

            /// true if this is a cylinder parallel to Z
            virtual bool isZCylinder() const  = 0;

            /// true if this is a plane parallel to Z
            virtual bool isZPlane() const = 0;

            /// true if this is a plane orthogonal to Z
            virtual bool isZDisk() const  = 0;


            /** True if surface is parallel to Z with accuracy epsilon - result is cached in bit SurfaceType::ParallelToZ */
            virtual bool checkParallelToZ(const ISurface& surf , double epsilon = 1.e-6) const = 0;

    };
}

#endif // USE_DD4HEP

namespace aidaTT
{
    class IGeometry
    {
        public:
            virtual const std::list<const ISurface*>& getSurfaces() = 0;
    };

}

#endif // IGEOMETRY_HH
