#ifndef AIDATT_UNITS_HH
#define AIDATT_UNITS_HH

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

namespace aidaTT
{

#ifdef AIDATT_USE_DD4HEP
#include "DD4hep/DD4hepUnits.h"

    typedef dd4hep::millimeter millimeter;
    typedef dd4hep::centimeter centimeter;
    typedef dd4hep::meter      meter;

    typedef dd4hep::mm mm;
    typedef dd4hep::cm cm;
    typedef dd4hep::m  m;

    typedef dd4hep::megaelectronvolt = 1.e-3;
    typedef dd4hep::electronvolt = 1.e-6 * megaelectronvolt;
    typedef dd4hep::kiloelectronvolt = 1.e-3 * megaelectronvolt;
    typedef dd4hep::gigaelectronvolt = 1.e+3 * megaelectronvolt;
    typedef dd4hep::teraelectronvolt = 1.e+6 * megaelectronvolt;
    typedef dd4hep::petaelectronvolt = 1.e+9 * megaelectronvolt;

    // symbols
    typedef dd4hep::MeV = megaelectronvolt;
    typedef dd4hep::eV = electronvolt;
    typedef dd4hep::keV = kiloelectronvolt;
    typedef dd4hep::GeV = gigaelectronvolt;
    typedef dd4hep::TeV = teraelectronvolt;
    typedef dd4hep::PeV = petaelectronvolt;

    static const double convertBr2P_cm = 0.299792458 * (centimeter / meter);

#else
///TODO
#endif // AIDATT_USE_DD4HEP

}
#endif // AIDATT_UNITS_HH
