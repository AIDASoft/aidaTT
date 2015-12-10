#ifndef AIDATT_HH
#define AIDATT_HH

//~ standard headers
#include <iostream>
#include <vector>

//~ aidaTT headers
#include "trajectory.hh"
#include "trajectoryElement.hh"
#include "trackParameters.hh"


#ifdef AIDATT_USE_STREAMLOG
#include "streamlog/streamlog.h"
#endif

/** \mainpage aidaTT -- a tracking toolkit
 *
 *  \section intro_sec Introduction
 *  aidaTT is a tracking toolkit, developed in the context of the EU framework programme of AIDA.
 *  It is a part of work package 2, dedicated to common software tools.
 *
 * It features a completely modular design, abstracting all individual components as abstract classes:
 * Field information (IBfield), geometry implementations (IGeometry), propagation methods (IPropagation) and ftting algorithms (IFittingAlgorithm).
 *
 *  \section usage_sec Short usage desciption
 * The central user interface is the aidaTT::trajectory class.
 * For a track that is to be (re-) fitted, the construction call uses the initial track parameters:
 *  trajectory(const trackParameters&, IFittingAlgorithm*, const IBField*, IPropagation*, const IGeometry*)
 * where the chosen implementations of the abstract classes are also passed as arguments.
 *
 * Measurements are then added by
 *  void trajectory::addMeasurement(const Vector3D& position, const std::vector<double>& resolution, const ISurface& surface, void* id),
 * other points or elements like scattering material or points of interest (nominal starting point, calorimeter face,...) are added with either
 *    void trajectory::addElement(const Vector3D& position, void* id);
 *    void trajectory::addElement(const Vector3D& position, const ISurface& surface, void* id);
 *
 * In order to prepare the trajectory for fitting, that is to calculate all needed internal parts the method
 *      void trajectory::prepareForFitting()
 * is invoked.
 * After this the actual fitting call can be made with
 *      bool trajectory::fit();
 * The result from the fit is accessible in a aidaTT::fitResults object.
 *
 *
 *  \section dep_sec Dependencies and Packaging
 *  There is only one hard dependency, which is the GSL; which is used for the internal representation of vectors and matrices.
 *  [ Note: The actual internal interface is adapted, thus even this is not completely hard. ]
 *
 *  The persistency format is LCIO, a lightweight IO format.
 *  The matrix and vector mathematics is done with GSL.
 *
 *  \section fit_sec Track Fitting
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
