#ifndef AIDATT_HH
#define AIDATT_HH

//~ standard headers
#include <iostream>
#include <vector>

//~ aidaTT headers
#include "trajectory.hh"
#include "trajectoryElement.hh"
#include "trackParameters.hh"

/** \mainpage aidaTT -- a tracking toolkit
 *
 *  \section intro_sec Introduction
 *  This is a tracking toolkit.
 *
 *  \section Dependencies and Packaging
 *  The persistency format is LCIO, a lightweight IO format.
 *  The matrix and vector mathematics is done with GSL.
 *
 *  \section Track Fitting
 *  Implemented is the General Broken Lines package (as external dependency) by Claus Kleinwort.
 * 
 * 
 *  \section ref_sec References
 *    - V. Blobel, C. Kleinwort, F. Meier,
 *      Fast alignment of a complex tracking detector using advanced track models,
 *      Computer Phys. Communications (2011), doi:10.1016/j.cpc.2011.03.017
 *    - C. Kleinwort, General Broken Lines as advanced track fitting method,
 *      NIM A, 673 (2012), 107-110, doi:10.1016/j.nima.2012.01.024
 */


/*
 * ultimately this will be the main interface class to configure the
 * construction of a trajectory
 * - decide which methods will be used for fitting, propagation and navigation
 *
 * there are three subtasks to be completed:
 *  - the propagation is the most basic requirement for the fitting, navigation and
 *  and the evaluation of questions like: what are the errors/position at some given point
 *
 *  - the fitting is pure mathematics, the information needs to be accessed
 *  *
 *  - the navigation is a purely geometrical problem, tightly connected to the implementation
 */


namespace aidaTT
{
    class AidaTT
    {
        public:
            AidaTT();
            ~AidaTT();

//            bool  createTrajectory(std::vector<const trajectoryElement&>, const trackParameters&);
            bool  createTrajectory();

            // yet empty
            bool initializeFitter();
            bool initializePropagation();
            bool initializeGeometry();
            bool initializeBField();

        private:
            AidaTT(const AidaTT&);
            AidaTT  operator=(const AidaTT&);

    };
}
#endif // AIDATT_HH
