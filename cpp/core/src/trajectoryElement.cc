#include "trajectoryElement.hh"

using namespace aidaTT;
using namespace std;
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
trajectoryElement::trajectoryElement(double arclength, const ISurface& surface, const vector<double>& resolutions, void* id)
    : _arclength(arclength), _surface(&surface), _measDim(resolutions.size()), _resolutions(resolutions), _id(id)
{ }



///~ constructor B: only the arc length is given and some identification
trajectoryElement::trajectoryElement(double arclength, void* id) : _arclength(arclength), _surface(NULL), _measDim(0), _id(id)
{}



///~ constructor C: only arc length and the jacobian to the next element plus some identification
trajectoryElement::trajectoryElement(double arclength, const fiveByFiveMatrix& jacob, void* id)
    : _arclength(arclength), _surface(NULL), _measDim(0), _jacobianToNext(jacob), _id(id)
{}



///~ constructor D: everything is already known: arc length, surface, the jacobian to the next element and some identification
trajectoryElement::trajectoryElement(double arclength,  const ISurface& surface, const vector<double>& resolutions, const fiveByFiveMatrix& jacob, void* id)
    : _arclength(arclength), _surface(&surface), _measDim(resolutions.size()), _resolutions(resolutions), _jacobianToNext(jacob), _id(id)
{}



// TO WRITE
//~ std::pair<trackParameters, fullCovariance> const trajectoryElement::getFullState() { return ;}



// TO WRITE
//~ trackParameters const trajectoryElement::getStateVector();



///~ set the jacobian to the next element
void trajectoryElement::setJacobianToNextElement(const fiveByFiveMatrix& jacob)
{
    /// !!! Need to rethink this:
    _jacobianToNext = jacob;
}
