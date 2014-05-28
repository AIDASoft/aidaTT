#ifdef USE_GBL

#ifndef GBLINTERFACE_H
#define GBLINTERFACE_H
#include "IFittingAlgorithm.hh"

#include "fiveByFiveMatrix.hh"
#include "Vector5.hh"

#include <vector>

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
        public:
            GBLInterface();
            ~GBLInterface();

            /// inherited methods:
            bool fit(const trajectory&);

            unsigned int getNDF()       const
            {
                return _ndf;
            };
            double getChiSquare()       const
            {
                return _chisquare;
            };
            double lostWeight()         const
            {
                return _lostweight;
            };


        private:
            GBLInterface(const GBLInterface&);
            GBLInterface& operator=(const GBLInterface&);

///< GBL trajectory
            gbl::GblTrajectory* _trajectory;
            ///< flag for curved track (helix, else straight line)
            bool _curvature;

            ///< correction vector from GBL fit (for track parameters)
            Vector5* _correctionVector;

            ///< covariance matrix from GBL
            fullCovariance* _covarianceMatrix;

            ///< number of degrees of freedom in GBL fit
            int _ndf;

            ///< chi2 from GBL fit
            double _chisquare;

            ///< weight lost by down-weighting, not used for now!
            double _lostweight;

            ///< arc-length of reference point
            double _refPointS;

            ///< label of reference point
            unsigned int _refPointIndex;
    };

}

#endif // GBLINTERFACE_H
#endif // USE_GBL
