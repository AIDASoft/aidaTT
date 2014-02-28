#include "aidaTT.hh"

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
    aidaTT::aidaTT()
    {
        cout << " it works " << endl;
    }

    aidaTT::~aidaTT()
    {
        cout << " it still works "  << endl;
    }

    bool initializeTrajectory()
    {
        cout << "init traj " << endl;
        return false;
    }
//~ bool initializeTrajectory(std::vector<const trajectoryElements&>, const trackParameters);

    bool initializeFitter()
    {
        cout << "init fit " << endl;
        return false;
    }

    bool initializePropagation()
    {
        cout << "init prop " << endl;
        return false;
    }


}
