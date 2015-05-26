#include "trajectoryElement.hh"
#include "helixGymnastics.hh"

namespace aidaTT
{
    /// standard constructor for measurements: arc length is given, the surface it belongs to;
    /// the measurement directions, resolution and residuals plus the local curvilinear system and some identification
    trajectoryElement::trajectoryElement(double arclength, const ISurface& surface, std::vector<Vector3D>* measDir, const std::vector<double>& precisions,
                                         const std::vector<double>& residuals, std::pair<Vector3D, Vector3D>* lCLS, void* id, bool isScatterer )
        : _arclength(arclength), _jacobianFromPrevious(NULL), _surface(&surface), _measurement(_surface->type().isSensitive()),
          _measDirections(measDir), _precisions(precisions), _residuals(residuals), _localCurvilinearSystem(lCLS),  _id(id), _scatterer( isScatterer )
    {
        _calculateLocalToMeasurementProjectionMatrix();
        if(_precisions.size() == 1)
            _precisions.push_back(0.);
    }



    ///~ constructor B: only the arc length is given and some identification
    trajectoryElement::trajectoryElement(double arclength, void* id) : _arclength(arclength), _jacobianFromPrevious(NULL), _surface(NULL), _measurement(false), _id(id)
    {}


    ///~ constructor C: arc length and surface is given plus some identification
    //~ trajectoryElement::trajectoryElement(double arclength, const ISurface* surface, void* id)
    //~ : _arclength(arclength), _surface(surface), _measurement(false), _id(id)
    //~ {}



    trajectoryElement::~trajectoryElement()
    {
        if(_localCurvilinearSystem)
            delete _localCurvilinearSystem;
        if(_measDirections)
            delete _measDirections;
        if(_jacobianFromPrevious)
            delete _jacobianFromPrevious;
        if(_localToMeasurementProjection)
            delete _localToMeasurementProjection;
    }



    void trajectoryElement::_calculateLocalToMeasurementProjectionMatrix()
    {
        if(!_measurement)
            return;

        _localToMeasurementProjection = calculateLocalToMeasurementProjectionMatrix(_localCurvilinearSystem->first, _localCurvilinearSystem->second, *_measDirections);
    }
}
