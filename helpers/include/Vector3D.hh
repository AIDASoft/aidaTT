#ifdef AIDATT_USE_DD4HEP // only use when not the DD4hep version is available
#include "DDSurfaces/Vector3D.h"

namespace aidaTT
{
    typedef DDSurfaces::Vector3D Vector3D;
}
#else
#ifndef VECTOR3D_HH
#define VECTOR3D_HH

#include <cmath>
#include <iostream>

namespace aidaTT
{
    /** Simple three dimensional vector providing the components for cartesian,
     *  cylindrical and spherical coordinate  systems - internal reperesentation is
     *  cartesian.
     *  @author F. Gaede, DESY
     */



    /* helper class by F. Gaede, DESY */

    class Vector3D
    {

        public:
            /** Default c'tor - zero vector */
            Vector3D() : _x(0.0), _y(0.0), _z(0.0) {}

            /** Copy constructor*/
            Vector3D(const Vector3D& v) : _x(v[0]), _y(v[1]), _z(v[2]) {}

            /** Constructor for float array.*/
            Vector3D(const float* v) : _x(v[0]), _y(v[1]), _z(v[2]) {}

            /** Constructor for double array.*/
            Vector3D(const double* v) : _x(v[0]), _y(v[1]), _z(v[2]) {}


            /** Templated c'tor - allows to have overloaded c'tors for different coordinates */
            template <class T>
            Vector3D(double x, double y, double z , T(&)()) ;


            /** Default coordinate system for initialization is cartesian */
            Vector3D(double x, double y, double z) :
                _x(x), _y(y), _z(z)
            {}

            //assignment operator
            Vector3D& operator=(const Vector3D& v)
            {
                _x = v[0] ;
                _y = v[1] ;
                _z = v[2] ;
                return *this ;
            }

            /// fill vector from arbitrary class that defines operator[]
            template <class T>
            inline const Vector3D& fill(const T& v)
            {
                _x = v[0] ;
                _y = v[1] ;
                _z = v[2] ;
                return *this ;
            }

            /// fill vector from double array
            inline const Vector3D& fill(const double* v)
            {
                _x = v[0] ;
                _y = v[1] ;
                _z = v[2] ;
                return *this ;
            }

            /// fill from double values
            inline const Vector3D& fill(double x, double y, double z)
            {
                _x = x ;
                _y = y ;
                _z = z ;
                return *this ;
            }


            /** Cartesian x coordinate */
            inline double x() const
            {
                return  _x ;
            }

            /** Cartesian y coordinate */
            inline double y() const
            {
                return  _y ;
            }

            /** Cartesian cartesian  z coordinate */
            inline double z() const
            {
                return  _z ;
            }

            /** Assign to  cartesian x coordinate */
            inline double& x()
            {
                return  _x ;
            }

            /**  Assign to  cartesian y coordinate */
            inline double& y()
            {
                return  _y ;
            }

            /**  Assign to cartesian z coordinate */
            inline double& z()
            {
                return  _z ;
            }


            /** Accessing x,y,z with bracket operator */
            inline double operator[](int i) const
            {
                switch(i)
                    {
                        case 0:
                            return _x ;
                            break ;
                        case 1:
                            return _y ;
                            break ;
                        case 2:
                            return _z ;
                            break ;
                    }
                return 0.0 ;
            }

            /** Accessing x,y,z with bracket operator for assignment */
            inline double& operator[](int i)
            {
                switch(i)
                    {
                        case 0:
                            return _x ;
                            break ;
                        case 1:
                            return _y ;
                            break ;
                        case 2:
                            return _z ;
                            break ;
                    }
                static double dummy(0.0) ;
                return dummy ;
            }

            /** Azimuthal angle - cylindrical and spherical */
            inline double phi() const
            {
                return _x == 0.0 && _y == 0.0 ? 0.0 : atan2(_y, _x);
            }

            /** Transversal component - cylindrical 'r' */
            inline double rho() const
            {
                return trans() ;
            }

            /** Transversal component */
            inline double trans() const
            {
                return sqrt(_x * _x + _y * _y) ;
            }

            /** Transversal component squared */
            inline double trans2() const
            {
                return  _x * _x + _y * _y  ;
            }

            /** Spherical r/magnitude */
            inline double r() const
            {
                return sqrt(_x * _x + _y * _y + _z * _z) ;
            }


            /** Spherical r/magnitude, squared */
            inline double r2() const
            {
                return  _x * _x + _y * _y + _z * _z  ;
            }

            /** Polar angle - spherical */
            inline double theta() const
            {
                return _x == 0.0 && _y == 0.0 && _z == 0.0 ? 0.0 : atan2(rho(), _z) ;
            }

            /** Scalar product */
            inline double dot(const Vector3D& v) const
            {
                return _x * v.x() + _y * v.y() + _z * v.z() ;
            }


            /** Vector product */
            inline Vector3D cross(const Vector3D& v) const
            {
                return Vector3D(_y * v.z() - _z * v.y() , _z * v.x() - _x * v.z() , _x * v.y() - _y * v.x())  ;
            }

            /** Parallel unit vector */
            inline Vector3D unit() const
            {
                double n = r() ;
                return Vector3D(_x / n , _y / n , _z / n) ;
            }


            /// direct access to data as const double*
            inline operator const double*() const
            {
                return &_x ;
            }

            /// direct access to data as double* - allows modification
            inline double* array()
            {
                return &_x ;
            }

            /** Explicit, unchecked conversion to anything that has a c'tor T(x,y,z).
             *  Example:
             *  CLHEP::Vector3D clhv = v.to< CLHEP::Vector3D>() ;
             *  @see operator T()
             */
            template <class T>
            inline T to() const
            {
                return T(_x, _y, _z) ;
            }

        protected:
            double _x, _y, _z ;

            // helper classes and function to allow
            // different c'tors selected at compile time
        public:

            struct Cartesian   { } ;
            struct Cylindrical { } ;
            struct Spherical   { } ;

            static Cartesian   cartesian()
            {
                return Cartesian() ;
            }
            static Cylindrical cylindrical()
            {
                return Cylindrical() ;
            }
            static Spherical   spherical()
            {
                return Spherical() ;
            }

    } ;

    /** Addition of two vectors */
    inline Vector3D operator+(const Vector3D& a, const Vector3D& b)
    {
        return Vector3D(a.x() + b.x()  , a.y() + b.y(), a.z() + b.z()) ;
    }

    /** Subtraction of two vectors */
    inline Vector3D operator-(const Vector3D& a, const Vector3D& b)
    {
        return Vector3D(a.x() - b.x()  , a.y() - b.y(), a.z() - b.z()) ;
    }

    /** Comparison of two vectors */
    inline bool operator==(const Vector3D& a, const Vector3D& b)
    {
        if(a.x() == b.x()  &&  a.y() == b.y() && a.z() == b.z())
            return true;
        else
            return false;
    }

    /** Multiplication with scalar */
    inline Vector3D operator*(double s , const Vector3D& v)
    {
        return Vector3D(s * v.x()  , s * v.y()  ,  s * v.z()) ;
    }

    /** Negative vector */
    inline Vector3D operator-(const Vector3D& v)
    {
        return Vector3D(-v.x(), - v.y(), - v.z()) ;
    }

    /// operator for scalar product
    inline double operator*(const Vector3D& v0, const Vector3D& v1)
    {
        return v0.dot(v1) ;
    }

    // template specializations for constructors of  different coordinate systems

    /** Cartesian c'tor  - example: <br>
     *  Vector3D  v( x, y, c , Vector3D::cartesian ) ;
     */
    template <>
    inline Vector3D::Vector3D(double x, double y, double z, Vector3D::Cartesian(&)()) :
        _x(x),  _y(y), _z(z)
    {}

    /** Cylindrical c'tor  - example: <br>
     *  Vector3D  v( rho, phi, z , Vector3D::cylindrical ) ;
     */
    template <>
    inline Vector3D::Vector3D(double rho, double phi, double z, Vector3D::Cylindrical(&)()) : _z(z)
    {
        _x = rho * cos(phi) ;
        _y = rho * sin(phi) ;
    }


    /** Spherical c'tor  - example: <br>
     *  Vector3D  v( r, phi, theta , Vector3D::spherical ) ;
     */
    template <>
    inline Vector3D::Vector3D(double r, double phi, double theta, Vector3D::Spherical(&)())
    {
        double rst =  r * sin(theta) ;
        _x = rst * cos(phi) ;
        _y = rst * sin(phi) ;
        _z = r * cos(theta) ;
    }

    /** Output operator */
    inline std::ostream & operator << (std::ostream & os, const Vector3D &v)
    {
        os << "( " << v[0] << ", " << v[1] << ", " << v[2] << " )" ;
        return os ;
    }
}
#endif // VECTOR3D_HH
#endif // AIDATT_USE_DD4HEP
