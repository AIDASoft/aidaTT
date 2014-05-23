#include "trajectoryElement.hh"

namespace aidaTT
{
//~ class trajectoryElement::
//~ {
//~ private:
//~
//~ bool _measurement;
//~ unsigned int _measDim;
//~ double _residualU, _residualV;
//~ double _dU, _dV;
//~
//~ const Surface& _surface;
//~ double _arcLength;
//~ fiveByFiveMatrix _jacobianToNext;
//~ };

    ///~ standard constructor A for measurements: arc length is given, the surface it belongs to and some identification
    trajectoryElement::trajectoryElement(double arclength, const ISurface& surface, const std::vector<double>& resolutions, void* id)
        : _arclength(arclength), _surface(&surface), _measurement(_surface->type().isSensitive()), _resolutions(resolutions), _id(id)
    {
 
    }



    ///~ constructor B: only the arc length is given and some identification
    trajectoryElement::trajectoryElement(double arclength, void* id) : _arclength(arclength), _surface(NULL), _measurement(false), _id(id)
    {}



// TO WRITE
//~ std::pair<trackParameters, fullCovariance> const trajectoryElement::getFullState() { return ;}



// TO WRITE
//~ trackParameters const trajectoryElement::getStateVector();



///~ set the jacobian to the next element
    void trajectoryElement::setJacobian(const fiveByFiveMatrix& jacob)
    {
        /// !!! Need to rethink this:
        _jacobianFromPrevious = jacob;
    }
}
