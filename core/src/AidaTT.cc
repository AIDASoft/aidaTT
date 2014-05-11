#include "AidaTT.hh"

#include "trackParameters.hh"
#include "trajectoryElement.hh"
#include "trajectory.hh"

//~ abstracts
#include "IGeometry.hh"
#include "IPropagation.hh"
#include "IFittingAlgorithm.hh"

//~ implementations by type
#include "analyticalPropagation.hh"
#include "simplifiedPropagation.hh"
#include "ConstantSolenoidBField.hh"
#include "GBLInterface.hh"

//~ and helper classes
#include "fiveByFiveMatrix.hh"
#include "Vector3D.hh"
#include "Vector5.hh"


namespace aidaTT
{
    AidaTT::AidaTT()
    {
        std::cout << " it works " << std::endl;
    }

    AidaTT::~AidaTT()
    {
        std::cout << " it still works "  << std::endl;
    }

    // mistake, mistake -> go back to bool to find it
    //bool AidaTT::createTrajectory(std::vector<const trajectoryElement&> vecTE, const trackParameters& tP)
    bool AidaTT::createTrajectory()
    {
        std::cout << "init traj " << std::endl;
        return false;
    }

//~ bool initializeTrajectory(std::vector<const trajectoryElements&>, const trackParameters);

    bool AidaTT::initializeFitter()
    {
        std::cout << "init fit " << std::endl;
        return false;
    }

    bool AidaTT::initializePropagation()
    {
        std::cout << "init prop " << std::endl;
        return false;
    }


}
