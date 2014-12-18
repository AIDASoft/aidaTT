#ifdef USE_LCIO


#include "LCIOPersistency.hh"

#include "Vector5.hh"
#include "Vector3D.hh"
#include "utilities.hh"

/* declare and define functions for LCIO interaction:
 * reading and writing LCIO::TrackState from trackParameters
 * reading and writing LCIO::Track from trajectory
 */

namespace aidaTT
{
#ifdef AIDATT_USE_DD4HEP
#include "DD4hep/DD4hepUnits.h"

    trackParameters readLCIO(const EVENT::TrackState* const ts)
    {
        trackParameters TP;
        /// track state can be read as L3 parametrization:  [ Omega, tan(lambda), phi_0, d_0, z_0, ]
	Vector5 param(ts->getOmega()/dd4hep::mm, ts->getTanLambda(), ts->getPhi(), ts->getD0()*dd4hep::mm, ts->getZ0()*dd4hep::mm);

	//debug: vary the track parameters by 1% :
	// Vector5 param( 1.01 * ts->getOmega()/dd4hep::mm, 
	// 	       1.01 * ts->getTanLambda(), 
	// 	       0.99 * ts->getPhi(), 
	// 	       0.99 * ts->getD0()*dd4hep::mm, 
	// 	       1.01  * ts->getZ0()*dd4hep::mm);

        /// get the reference point
        Vector3D refPoint(ts->getReferencePoint()[0]*dd4hep::mm , ts->getReferencePoint()[1]*dd4hep::mm , ts->getReferencePoint()[2]*dd4hep::mm);

        /// the covariance matrix is stored as lower triangle,
        ///  order of parameters is: d0, phi, omega, z0, tan(lambda)
        fullCovariance covarianceMatrix;

        // omega, omega - 0,0
        covarianceMatrix(0, 0) = ts->getCovMatrix().at(5)/(dd4hep::mm*dd4hep::mm);
        // tanLambda,tanLambda (1,1)
        covarianceMatrix(1, 1) = ts->getCovMatrix().at(14);
        // phi,phi (2,2)
        covarianceMatrix(2, 2) = ts->getCovMatrix().at(2);
        // d0,d0 (3,3)
        covarianceMatrix(3, 3) = ts->getCovMatrix().at(0)*dd4hep::mm*dd4hep::mm;
        // z0,z0 (4,4)
        covarianceMatrix(4, 4) = ts->getCovMatrix().at(9)*dd4hep::mm*dd4hep::mm;

        // omega, tanLambda (0,1) - (1,0)
        covarianceMatrix(1, 0) = ts->getCovMatrix().at(12)/dd4hep::mm;
        covarianceMatrix(0, 1) = ts->getCovMatrix().at(12)/dd4hep::mm;
        // omega, phi (0,2) - (2,0)
        covarianceMatrix(2, 0) = ts->getCovMatrix().at(4)/dd4hep::mm;
        covarianceMatrix(0, 2) = ts->getCovMatrix().at(4)/dd4hep::mm;
        //omega, d0 (0,3) - (3,0)
        covarianceMatrix(3, 0) = ts->getCovMatrix().at(3);
        covarianceMatrix(0, 3) = ts->getCovMatrix().at(3);
        // omega, z0 (0,4) - (4,0)
        covarianceMatrix(4, 0) = ts->getCovMatrix().at(8);
        covarianceMatrix(0, 4) = ts->getCovMatrix().at(8);

        // tanLambda, phi (1,2) - (2,1)
        covarianceMatrix(1, 2) = ts->getCovMatrix().at(11);
        covarianceMatrix(2, 1) = ts->getCovMatrix().at(11);
        // tanLambda, d0 (1,3) - (3,1)
        covarianceMatrix(1, 3) = ts->getCovMatrix().at(10)*dd4hep::mm;
        covarianceMatrix(3, 1) = ts->getCovMatrix().at(10)*dd4hep::mm;
        // tanLambda, z0 (1,4) - (4,1)
        covarianceMatrix(1, 4) = ts->getCovMatrix().at(13)*dd4hep::mm;
        covarianceMatrix(4, 1) = ts->getCovMatrix().at(13)*dd4hep::mm;

        // phi, d0 (2,3) - (3,2)
        covarianceMatrix(3, 2) = ts->getCovMatrix().at(1)*dd4hep::mm;
        covarianceMatrix(2, 3) = ts->getCovMatrix().at(1)*dd4hep::mm;
        // phi, z0 (2,4) - (4,2)
        covarianceMatrix(2, 4) = ts->getCovMatrix().at(7)*dd4hep::mm;
        covarianceMatrix(4, 2) = ts->getCovMatrix().at(7)*dd4hep::mm;

        // d0, z0 (3,4) - (4,3)
        covarianceMatrix(3, 4) = ts->getCovMatrix().at(6)*dd4hep::mm*dd4hep::mm;
        covarianceMatrix(4, 3) = ts->getCovMatrix().at(6)*dd4hep::mm*dd4hep::mm;

        TP.setTrackParameters(param, covarianceMatrix, refPoint);
        return  TP;
    }



    IMPL::TrackStateImpl* createLCIO(const trackParameters& tp)
    {
        Vector5 param = tp.parameters();
        Vector3D refPoint = tp.referencePoint();
        fullCovariance covarianceMatrix = tp.covarianceMatrix();

        IMPL::TrackStateImpl* tsi = new IMPL::TrackStateImpl();
        /// sets location -- yet undefined
        tsi->setLocation(0);

        /// set parameter values
        tsi->setD0( calculateD0(tp) / dd4hep::mm);
        tsi->setPhi( calculatePhi0(tp) );
        tsi->setOmega( calculateOmega(tp) * dd4hep::mm);
        tsi->setZ0( calculateZ0(tp)  / dd4hep::mm);
        tsi->setTanLambda( calculateTanLambda(tp) );

        std::vector<float> covm;

        // d0,d0 (3,3)
        covm.push_back(covarianceMatrix(3, 3) / (dd4hep::mm * dd4hep::mm) ) ;  // ts->getCovMatrix().at(0);
        // phi, d0 (2,3) - (3,2)
        covm.push_back(covarianceMatrix(3, 2) / dd4hep::mm);  // ts->getCovMatrix().at(1);
        // phi,phi (2,2)
        covm.push_back(covarianceMatrix(2, 2));  // ts->getCovMatrix().at(2);
        //omega, d0 (0,3) - (3,0)
        covm.push_back(covarianceMatrix(3, 0));  // ts->getCovMatrix().at(3);
        // omega, phi (0,2) - (2,0)
        covm.push_back(covarianceMatrix(0, 2) * dd4hep::mm);  // ts->getCovMatrix().at(4);
        // omega, omega - 0,0
        covm.push_back(covarianceMatrix(0, 0) * (dd4hep::mm * dd4hep::mm));  // ts->getCovMatrix().at(5);
        // d0, z0 (3,4) - (4,3)
        covm.push_back(covarianceMatrix(3, 4) / (dd4hep::mm * dd4hep::mm));  // ts->getCovMatrix().at(6);
        // phi, z0 (2,4) - (4,2)
        covm.push_back(covarianceMatrix(2, 4) / dd4hep::mm);  // ts->getCovMatrix().at(7);
        // omega, z0 (0,4) - (4,0)
        covm.push_back(covarianceMatrix(4, 0));  // ts->getCovMatrix().at(8);
        // z0,z0 (4,4)
        covm.push_back(covarianceMatrix(4, 4) / (dd4hep::mm * dd4hep::mm));  // ts->getCovMatrix().at(9);
        // tanLambda, d0 (1,3) - (3,1)
        covm.push_back(covarianceMatrix(1, 3) / dd4hep::mm);  // ts->getCovMatrix().at(10);
        // tanLambda, phi (1,2) - (2,1)
        covm.push_back(covarianceMatrix(1, 2));  // ts->getCovMatrix().at(11);
        // omega, tanLambda (0,1) - (1,0)
        covm.push_back(covarianceMatrix(0, 1) * dd4hep::mm);  // ts->getCovMatrix().at(12);
        // tanLambda, z0 (1,4) - (4,1)
        covm.push_back(covarianceMatrix(1, 4) / dd4hep::mm);  // ts->getCovMatrix().at(13);
        // tanLambda,tanLambda (1,1)
        covm.push_back(covarianceMatrix(1, 1));  // ts->getCovMatrix().at(14);

        tsi->setCovMatrix(covm);

        /// set reference point
        float rp[3] = { refPoint.x() / dd4hep::mm, refPoint.y() / dd4hep::mm, refPoint.z() / dd4hep::mm};
        tsi->setReferencePoint(rp);

        return tsi;
    }
#endif // AIDATT_USE_DD4HEP


//~ trajectory readLCIO(const EVENT::Track* const t)
//~ {
//~ }
//~
//~
//~
//~ IMPL::TrackImpl* createLCIO(const trajectory& traj)
//~ {
//~ }
}
#endif // USE_LCIO
