#ifndef API2GBL_H
#define API2GBL_H

#include "IfittingAlgorithm.h"

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

class api2gbl : public IfittingAlgorithm
{
	protected:
		std::string _inputCollectionName; ///< Name of the input collection -- track seeds
		std::string _outputTrackCollectionName; ///< Name of the output track collection
		double _bfieldScaleFactor; ///< scale factor for magnetic field (default: 1.0)
		bool _outputIsPersistent;  ///< whether the output is to be stored or not (default: true)
		std::string _fitOptions; ///< list of iterations with 'h' Huber, 't' Tukey or 'c' Cauchy down-weighting
		bool _writeMillepedeOut; ///< selects output to Millepede-II binary file
		std::string _fileNameMillepedeFile; ///< name of Millepede-II binary file
		bool _encodedModuleID; ///< Module ID is encoded in CellID0

	private:
		gbl::GblTrajectory* _trajectory; ///< GBL trajectory
// not yet:		gbl::MilleBinary* _milleBinary; ///< Millepede-II binary file
		bool _curvature; ///< flag for curved track (helix, else straight line)
		stdj::vector<unsigned int> _theLabels; ///< Labels of (global) MP-II parameters
		TVectorD* _correctionVector; ///< correction vector from GBL fit (for track parameters)
		TMatrixDSym* _covarianceMatrix; ///< covariance matrix from GBL
		int _ndf; ///< number of degrees of freedom in GBL fit
		double _chisquare; ///< chi2 from GBL fit
		double _lostweight; ///< weight lost by down-weighting
		unsigned int _refPointIndex; ///< label of reference point


		// TODO: this remains to be seen/implemented
		//void checkSorting(const std::vector<helperHit*>&) const;

};

}

#endif // API2GBL_H
