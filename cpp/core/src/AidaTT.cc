#include "AidaTT.hh"

#include "trackParameters.hh"
#include "trajectoryElement.hh"
#include "trajectory.hh"

#include "IGeometryAPI.hh"
#include "IMaterial.hh"
#include "IMeasurementSurface.hh"
#include "IPropagation.hh"
#include "IFittingAlgorithm.hh"

//~ #include "api2gbl.h"
//~ #include "simpleKalman.h"
//~
//~ #include "propagations/analyticalPropagation.h"
//~ #include "propagations/RungeKuttaPropagation.h"



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
