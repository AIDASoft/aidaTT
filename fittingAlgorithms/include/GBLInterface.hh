#ifdef USE_GBL

#ifndef GBLINTERFACE_H
#define GBLINTERFACE_H
#include "IFittingAlgorithm.hh"
#include "trajectory.hh"
#include "fiveByFiveMatrix.hh"
#include "Vector5.hh"

#include <vector>
#include <map>

// GBL:
#include "GblTrajectory.h"
#include "MilleBinary.h"


/* This is an interface class to the general broken lines package by C. Kleinwort.
 * It builds heavily on an existing implementation in the MarlinTPC package.
 *
 *  @author Ch. Rosemann, DESY (straight line), C. Kleinwort, DESY (extension to helix)
 *
 * This code interfaces the re-fitting functionality of the general broken line package by Claus Kleinwort.
 * For information on the General Broken Lines, please refer to the publication and/or the wiki pages:
 *   - C. Kleinwort, General Broken Lines as advanced track fitting method, NIM A, 673 (2012), 107-110, doi:10.1016/j.nima.2012.01.024
 *   - http://www.wiki.terascale.de/index.php/GeneralBrokenLines
 */

namespace aidaTT
{

  class GBLInterface : public IFittingAlgorithm
  {
    typedef  std::map< int, fitResults* > ResMap ;

  public:
    GBLInterface();
    ~GBLInterface();

    /// inherited methods:
    bool fit(const trajectory&);

    const fitResults* getResults(int label=0) const
    {

      const fitResults* res = 0 ;

      ResMap::const_iterator it = _theResults.find( label ) ;
      
      if( it == _theResults.end() ){
	
	res = _fillResults( *_fittedTraj , label ) ; 

      } else {

	res = it->second ;
      }
      
      return res ;
    };

  private:
    GBLInterface(const GBLInterface&);
    GBLInterface& operator=(const GBLInterface&);

    ///< GBL trajectory
    gbl::GblTrajectory* _trajectory{};
    ///< flag for curved track (helix, else straight line)
    bool _curvature{};

    ///< correction vector from GBL fit (for track parameters)
    Vector5* _correctionVector{};

    ///< covariance matrix from GBL
    fullCovariance* _covarianceMatrix{};

    ///< number of degrees of freedom in GBL fit
    int _ndf{};

    ///< chi2 from GBL fit
    double _chisquare{};

    ///< weight lost by down-weighting, not used for now!
    double _lostweight{};

    ///< arc-length of reference point
    double _refPointS{};

    ///< label of reference point
    unsigned int _refPointIndex{};

    const fitResults* _fillResults(const trajectory&, int label=0) const ;

    void _clear() ;

    const trajectory* _fittedTraj {};
    
    mutable ResMap _theResults{};

    gbl::MilleBinary* _milleBinary {};

  };

}

#endif // GBLINTERFACE_H
#endif // USE_GBL
