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




using namespace std;

namespace aidaTT
{
    AidaTT::AidaTT()
    {
        cout << " it works " << endl;
    }

    AidaTT::~AidaTT()
    {
        cout << " it still works "  << endl;
    }

    // mistake, mistake -> go back to bool to find it
    //bool AidaTT::createTrajectory(std::vector<const trajectoryElement&> vecTE, const trackParameters& tP)
    bool AidaTT::createTrajectory()
    {
        cout << "init traj " << endl;
        return false;
    }

//~ bool initializeTrajectory(std::vector<const trajectoryElements&>, const trackParameters);

    bool AidaTT::initializeFitter()
    {
        cout << "init fit " << endl;
        return false;
    }

    bool AidaTT::initializePropagation()
    {
        cout << "init prop " << endl;
        return false;
    }


}
