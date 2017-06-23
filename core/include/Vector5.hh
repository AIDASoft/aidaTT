#ifndef VECTOR5_HH
#define VECTOR5_HH

#include <ostream>
#include <vector>

#include <Eigen/Core>

namespace aidaTT
{

  typedef Eigen::Matrix<double, 5, 1 > Vector5d;


  /** This is a simple abstraction of a five-element vector that will work with the matrix class
   *  implemented using Eigen
   *
   */
  class Vector5 {

    friend class fiveByFiveMatrix;

  public:
    /** the default construction, it initializes all entries to zero **/
    Vector5(){  _v.Zero() ; }

    /** copy constructor **/
    Vector5(const Vector5& o) : _v( o._v ) {}

    /** the construction with a vector
     **/
    Vector5(const std::vector<double>& o) : _v( &o[0]) {}

    /** the construction with an array
     **/
    Vector5(const double* o) : _v(o) {}

    /** the construction from five doubles
     **/
    Vector5(double v0, double v1, double v2, double v3, double v4){
      _v << v0,v1,v2,v3,v4 ;
    }

    /** the destructor **/
    ~Vector5(){}

    /** assignment operator **/
    Vector5& operator=(const Vector5& o){
      if( this != &o ) {
	_v = o._v ;
      }
      return *this ;
    }

    // addition operator  - in place !?
    Vector5& operator+=(const Vector5& o){ 
      _v += o._v ;
      return *this ;
    }

    // addition w/ copy
    Vector5 operator+(const Vector5& o) const{
      return Vector5( _v + o._v ) ;
    } 


    /** Direct read access to the individual elements by index**/
    double operator()(unsigned int i) const{
      return _v(i) ;
    }

    /** direct write access to the individual elements by index**/
    double& operator()(unsigned int i){
      return _v(i) ;
    }

  protected:

    // conversion c'tor from Eigen's vector 
    Vector5(const Vector5d & o) : _v( o ) {}

  private:
    Vector5d _v ;
  };


  inline std::ostream & operator << (std::ostream & os, const Vector5& v)
  {
    os << " {" << v(0) << " , " << v(1) << " , " << v(2) << " , " << v(3) << " , " << v(4)  << "} "  ;
    return os ;
  }
}
#endif // VECTOR5_HH
