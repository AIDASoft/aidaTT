#include "trajectoryElement.hh"

namespace aidaTT
{
    ///~ standard constructor A for measurements: arc length is given, the surface it belongs to and some identification
    trajectoryElement::trajectoryElement(double arclength, const ISurface& surface, const std::vector<double>& resolutions, std::pair<Vector3D, Vector3D>* lCLS, void* id)
        : _arclength(arclength), _surface(&surface), _measurement(_surface->type().isSensitive()), _resolutions(resolutions), _localCurvilinearSystem(lCLS), _id(id)
    {
        _calculateLocalToMeasurementProjectionMatrix();
    }



    ///~ constructor B: only the arc length is given and some identification
    trajectoryElement::trajectoryElement(double arclength, void* id) : _arclength(arclength), _surface(NULL), _measurement(false), _id(id)
    {}


    ///~ constructor C: arc length and surface is given plus some identification
    trajectoryElement::trajectoryElement(double arclength, const ISurface& surface, void* id)
        : _arclength(arclength), _surface(&surface), _measurement(false), _id(id)
    {}



    trajectoryElement::~trajectoryElement()
    {
        if(_localCurvilinearSystem)
            delete _localCurvilinearSystem;
        if(_measDirections)
            delete _measDirections;
        if(_jacobianFromPrevious)
            delete _jacobianFromPrevious;
    }

// TO WRITE
//~ std::pair<trackParameters, fullCovariance> const trajectoryElement::getFullState() { return ;}



// TO WRITE
//~ trackParameters const trajectoryElement::getStateVector();



    void trajectoryElement::_calculateLocalToMeasurementProjectionMatrix()
    {
        if(!_measurement)
            return;

        const Vector3D clU = _localCurvilinearSystem->first;
        const Vector3D clV = _localCurvilinearSystem->second;

        for(std::vector<Vector3D>::iterator measDir = _measDirections->begin(), last = _measDirections->end(); measDir < last; ++measDir)
            {
                Vector3D projection = ((*measDir) * clU) * clU + ((*measDir) * clV) * clV;
                _localToMeasurementProjection.push_back(projection);
            }

    }
}
