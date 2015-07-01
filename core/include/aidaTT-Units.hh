/* this is the declaration of units used in exchange to the outside world inside aidaTT
 *
 * the first main interaction point is the geometry, that is
 *  -- the information whether a point is on a surface
 *  -- the material information
 *  regarding the lengths, the standard unit is centimeter, in agreement with the DD4hep package and ROOT::TGeo
 *
 * the second is the declaration of track parameters which is the
 * L3 track parametrization [ Omega, tan(lambda), phi_0, d_0, z_0 ]
 *      Omega is the inverse curvature, therefore [cm^{-1}]
 *      d_0 is in [cm]
 *      z_0 is in [cm]
 *
 * This then means that the conversion from curvature in [cm^{-1}]
 * to inverse momentum in [GeV^{-1}] needs
 *      a Bfield given in [T], and
 *      a speed of light given in [cm/s]
 */

#ifdef AIDATT_USE_DD4HEP
#include "DD4hep/DD4hepUnits.h"

namespace aidaTT
{
    using  dd4hep::millimeter;
    using  dd4hep::centimeter;
    using  dd4hep::meter;

    using  dd4hep::mm;
    using  dd4hep::cm;
    using  dd4hep::m;

    using  dd4hep::megaelectronvolt;
    using  dd4hep::electronvolt;
    using  dd4hep::kiloelectronvolt;
    using  dd4hep::gigaelectronvolt;
    using  dd4hep::teraelectronvolt;
    using  dd4hep::petaelectronvolt;

    // symbols
    using  dd4hep::MeV;
    using  dd4hep::eV;
    using  dd4hep::keV;
    using  dd4hep::GeV;
    using  dd4hep::TeV;
    using  dd4hep::PeV;

  static const double convertBr2P_cm = 0.299792458 * (centimeter / meter);
  //static const double convertBr2P_cm = 1.* (centimeter / meter);
}
#else

#ifndef AIDATT_UNITS_HH
#define AIDATT_UNITS_HH
namespace aidaTT
{
    ///TODO

}
#endif // AIDATT_UNITS_HH

#endif // AIDATT_USE_DD4HEP

