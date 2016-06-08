#ifndef materialUtils_HH
#define materialUtils_HH

#include "trackParameters.hh"
#include "IGeometry.hh"

/** Define helper functions for material effects.
 *
 *  @version $Id
 *  @author F. Gaede, DESY
 */

namespace aidaTT {

  static const double pionMass =  0.13957018 ; // [GeV]

  /** Compute the multiple scattering amplitude for a track going through
   *  the surface at the given local position and with the given momentum [GeV]
   *  Simple approximation from PDG 2012 (27.15). 
   *  If no mass is given the pion mass is assumed.
   * 
   *  @author F.Gaede, Y.Voutsinas, DESY
   *  based on original version by K.Fujii in KalTest
   */
  double computeQMS( const ISurface* surf, const Vector2D& crossingPoint, 
		     const Vector3D& mom, double mass) ; //=pionMass ) ;
  
  
  
  /** Compute the multiple scattering amplitude for a track going through
   *  the surface. The crossing point and momentum is computed from the 
   *  helix parameters and reference point.
   *  Simple approximation from PDG 2012 (27.15). 
   *  If no mass is given the pion mass is assumed.
   * 
   *  @author F.Gaede, Y.Voutsinas, DESY
   *  based on original version by K.Fujii in KalTest
   */
  double computeQMS( const ISurface* surf, const Vector5& hp, const Vector3D& rp, double mass ) ; //=pionMass ) ;
  

  /** Compute the Bethe-Bloch energy loss (per unit density and path length) 
   *  for the given material and momentum.
   *  If no mass is given the pion mass is assumed.
   *  @author F.Gaede, Y.Voutsinas, DESY
   *  based on original version by K.Fujii in KalTest
   */
  double computeBetheBloch( const IMaterial &mat, double mom, double mass ) ; //=pionMass );


  /** Compute the expected energy loss [GeV] for the particle traversing the surface
   *  at the the given (local) crossingPoint. For convenience the total energy and beta
   *  of the particle are also returned. 
   *  If no mass is given the pion mass is assumed.
   *  @author F.Gaede, Y.Voutsinas, DESY
   */
  double computeEnergyLoss( const ISurface* surf, const Vector2D& crossingPoint, 
			    const Vector3D& mom, double& energy, double& beta,
			    double mass ) ; //=pionMass  ) ;


  /** Compute the expected energy loss [GeV] for the particle traversing the surface.
   *  For convenience the total energy and beta of the particle are also returned. 
   *  If no mass is given the pion mass is assumed.
   *  @author F.Gaede, Y.Voutsinas, DESY
   */
  double computeEnergyLoss( const ISurface* surf, const trackParameters& param , 
			    double& energy, double& beta,
			    double mass ) ; //=pionMass  ) ;

  
  

}
#endif
