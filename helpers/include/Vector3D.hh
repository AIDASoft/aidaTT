#ifndef VECTOR3D_HH
#define VECTOR3D_HH

//#include <cmath>
#include <math.h>

namespace aidaTT
{
    /** Simple three dimensional vector providing the components for cartesian,
     *  cylindrical and spherical coordinate  systems - internal reperesentation is
     *  cartesian.
     *  @author F. Gaede, DESY
     */

    class Vector3D
    {

        public:

            /** Default c'tor - zero vector */
            Vector3D() : _x(0.0), _y(0.0), _z(0.0) {}


            /** Constructor for float array.*/
            Vector3D(const float* v) : _x(v[0]), _y(v[1]), _z(v[2]) {}

            /** Constructor for double array.*/
            Vector3D(const double* v) : _x(v[0]), _y(v[1]), _z(v[2]) {}

            /** Default corrdinate system for initialization is cartesian */
            Vector3D(double x, double y, double z) :
                _x(x),
                _y(y),
                _z(z)
            {
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
            inline double z()
            {
                return  _z ;
            }

            /** Assign to  cartesian x coordinate */
            inline double x()
            {
                return  _x ;
            }

            /**  Assign to  cartesian y coordinate */
            inline double y()
            {
                return  _y ;
            }

            /**  Assign to cartesian z coordinate */
            inline double z() const
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

                return Vector3D(_y * v.z() - _z * v.y() ,
                                _z * v.x() - _x * v.z() ,
                                _x * v.y() - _y * v.x())  ;
            }

            /** Parallel unit vector */
            inline Vector3D unit() const
            {
                double n = r() ;
                return Vector3D(_x / n , _y / n , _z / n) ;
            }


            /// /** Output operator */
            /// std::ostream & operator << (std::ostream & os, const Vector3D &v) ;

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


    };
}
#endif // VECTOR3D_HH
