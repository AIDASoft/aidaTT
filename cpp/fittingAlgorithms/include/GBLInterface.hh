#ifdef USE_GBL

#ifndef GBLINTERFACE_H
#define GBLINTERFACE_H
#include "IFittingAlgorithm.hh"

#include "fiveByFiveMatrix.hh"
#include "Vector5.hh"

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
            bool const fit(const trajectory&)       ;
            unsigned int const getNDF()             ;
            double const getChiSquare()             ;
            double const lostWeight()               ;


        private:
            GBLInterface(const GBLInterface&);
            GBLInterface operator=(const GBLInterface&);

            gbl::GblTrajectory* _trajectory; ///< GBL trajectory
            bool _curvature; ///< flag for curved track (helix, else straight line)
            std::vector<unsigned int> _theLabels; ///< Labels of (global) MP-II parameters
            Vector5* _correctionVector; ///< correction vector from GBL fit (for track parameters)
            fullCovariance* _covarianceMatrix; ///< covariance matrix from GBL
            int _ndf; ///< number of degrees of freedom in GBL fit
            double _chisquare; ///< chi2 from GBL fit
            double _lostweight; ///< weight lost by down-weighting
            double _refPointS; ///< arc-length of reference point
            unsigned int _refPointIndex; ///< label of reference point
            double _Xcenter; ///< TPC center, X coordinate
            double _Ycenter; ///< TPC center, Y coordinate

    };

}

#endif // GBLINTERFACE_H
#endif // USE_GBL
