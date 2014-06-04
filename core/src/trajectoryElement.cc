#include "trajectoryElement.hh"

namespace aidaTT
{
    ///~ standard constructor A for measurements: arc length is given, the surface it belongs to and some identification
    trajectoryElement::trajectoryElement(double arclength, const ISurface& surface, std::vector<Vector3D>* measDir, const std::vector<double>& resolutions,
                                         std::pair<Vector2D*, Vector2D*>* uvs, std::pair<Vector3D, Vector3D>* lCLS, void* id)
        : _arclength(arclength), _jacobianFromPrevious(NULL), _surface(&surface), _measurement(_surface->type().isSensitive()),
          _measDirections(measDir), _resolutions(resolutions), _UVvalues(uvs), _localCurvilinearSystem(lCLS),  _id(id)
    {
        _calculateLocalToMeasurementProjectionMatrix();
        _calculateResiduals();
    }



    ///~ constructor B: only the arc length is given and some identification
    trajectoryElement::trajectoryElement(double arclength, void* id) : _arclength(arclength), _jacobianFromPrevious(NULL), _surface(NULL), _measurement(false), _id(id)
    {}


    ///~ constructor C: arc length and surface is given plus some identification
    //~ trajectoryElement::trajectoryElement(double arclength, const ISurface* surface, void* id)
    //~ : _arclength(arclength), _surface(surface), _measurement(false), _id(id)
    //~ {}


    ///~ calculate the residuals from the two UV vectors: measurement MINUS reference value
    void trajectoryElement::_calculateResiduals()
    {
        _residuals.push_back((_UVvalues->first)->u() - (_UVvalues->second)->u());
        _residuals.push_back((_UVvalues->first)->v() - (_UVvalues->second)->v());
    }




    trajectoryElement::~trajectoryElement()
    {
        if(_localCurvilinearSystem)
            delete _localCurvilinearSystem;
        if(_measDirections)
            delete _measDirections;
        if(_jacobianFromPrevious)
            delete _jacobianFromPrevious;
        if(_UVvalues)
            {
                delete _UVvalues->first;
                delete _UVvalues->second;
                delete _UVvalues;
            }
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
