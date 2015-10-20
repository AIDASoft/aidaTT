#ifndef TRAJECTORYELEMENT_H
#define TRAJECTORYELEMENT_H

#include "IGeometry.hh"
#include "trackParameters.hh"
#include "fiveByFiveMatrix.hh"

#include <utility>
#include <vector>
#include <ostream>

#include <TMath.h>
#include <TMatrixDSym.h>

namespace aidaTT
{


  /*
   * this is the main class that provides info on the elements that make up a track
   * - the measurements
   * - the material
   * - the information how to get to the next/previous element
   *
   * the definitions of:
   *  what a measurementsurface is
   *  what material it is
   * are external!
   *
   * an element is uniquely identified and identifiable by its arc length
   * for a fixed point of reference (!)
   *  ?? how can this be made explicit
   *
   * An element is anything that is of interest in the track fitting or extrapolation:
   *  - the key elements are the measurements, without them there is no track
   *  - also the material between the measurements or extrapolation points is crucial
   *  - also points of interest, e.g. the calorimeter surface, the beam pipe, the nominal IP, etc.
   *
   * If the trajectoryElement is nothing but a point on the trajectory, then only the jacobian to the next element will be known.
   *
   * If a trajectoryElement is a measurement, then several things must be set:
   *  - the connection to the measurement surface
   *  - the residual(s) associated with the measurement
   *  - the associated uncertainty (or -ties) matching the residual(s)
   *  - the dimensionality of the measurement (although this is clear from the surface)
   *
   * If a trajectoryElement contains material, then the following must be initialized:
   *  - the connection to the material information that allows updating
   *  - the current estimate of the radiation length and atomic charge Z, possibly averaged.
   * [ note: this follows the idea, that any thick scatterer can be approximated by two thin scatters. ]
   *
   * All information that is put into the trajectoryElement can also be retrieved.
   *
   */


  class trajectoryElement 
  {
  public:
    /// Measurement constructor: arc length, surface, measurement direction(s), precision(s) and residual(s) plus the local curvilinear system and some identification
    trajectoryElement(double, const ISurface&, std::vector<Vector3D>*, const std::vector<double>&,
		      const std::vector<double>&, std::pair<Vector3D, Vector3D>*,  void* = NULL, bool isScatterer=false);
    //alternative constructor using TMatrixDSym for the precision
    //trajectoryElement(double, const ISurface&, std::vector<Vector3D>*, const TMatrixDSym& precisions, const std::vector<double>&, std::pair<Vector3D, Vector3D>*,  void* = NULL, bool isScatterer=false);

    ///~ constructor B: only the arc length is given and some identification
    trajectoryElement(double, void* = NULL);

    ~trajectoryElement();

    /// the getting routines
    const ISurface& surface() const
    {
      return *_surface;
    };

    double arcLength() const
    {
      return _arclength;
    };

    const fiveByFiveMatrix& jacobian() const
    {
      return *_jacobianFromPrevious;
    };

    const fiveByFiveMatrix& jacobianFromPrevious() const
    {
      return *_jacobianFromPrevious;
    };

    trackParameters  fullState() const;

    bool hasMeasurement() const
    {
      return  _measurement;
    };

    bool isScatterer() const
    {
      return _scatterer;
    };

    bool isThickScatterer() const
    {
      return _thick;
    };

    // the following depend on the type of element:
    unsigned int measurementDimension() const
    {
      return _measDirections->size();
    };

    ///~ get the residual: measurement - expected position (!)
    const std::vector<double>& measurementResiduals() const
    {
      return _residuals;
    };

    const std::vector<double>& precisions() const
    //      const TMatrixDSym& precisions() const
    {
      return _precisions;
    };

    ///~ access the local curvilinear system at the given point
    const std::pair<Vector3D, Vector3D>& localCurvilinearSystem() const
    {
      return *_localCurvilinearSystem;
    };

    ///~ access the measurement direction (if available) for the corresponding surface
    const std::vector<Vector3D>& measurementDirections() const
    {
      return *_measDirections;
    };

    ///~ and finally: the projection matrix from the local track frame to the measurement system
    const std::vector<double>& localToMeasurementProjection() const
    {
      return *_localToMeasurementProjection;
    };

    ///~ set the jacobian from the previous element -- ownership of the memory is transferred
    void setJacobian(fiveByFiveMatrix* jacob)
    {
      _jacobianFromPrevious = jacob;
    }
  private:
    ///~ no construction without the arc length!
    trajectoryElement();

    ///~ no assigment, no copying
    trajectoryElement(const trajectoryElement&);
    trajectoryElement operator=(const trajectoryElement&);

    void _calculateMaterial();
    void _calculateLocalToMeasurementProjectionMatrix();

    ///~
    double _arclength;

    fiveByFiveMatrix* _jacobianFromPrevious;
    const ISurface* _surface;

    ///~ measurement variables:
    bool _measurement;
    std::vector<Vector3D>* _measDirections;
    std::vector<double> _precisions;
    //TMatrixDSym _precisions;
    std::vector<double> _residuals;

    ///~ local curvilinear system
    std::pair<Vector3D, Vector3D>* _localCurvilinearSystem;

    /// 2x2 matrix for 2D measurements, projection from local cl to measurement system
    std::vector<double>* _localToMeasurementProjection;

    ///~ scattering info
    bool _scatterer;
    bool _thick;

    const void* const _id; // just store
  };


  inline std::ostream& operator<<(  std::ostream& os, const trajectoryElement& e ) {
    
    os << " trajectory element at s= " <<  e.arcLength() << " has measurement: " <<  e.hasMeasurement() << " isScatterer : " << e.isScatterer()
       << "\n surface: " << e.surface()
      //       << "\n state: " << e.fullState()
       << "\n jacobian:                " <<  e.jacobian()
       << "\n jacobian from previous : " <<  e.jacobianFromPrevious()
       << std::endl ;

    return os ;
  }

}
#endif // TRAJECTORYELEMENT_H
