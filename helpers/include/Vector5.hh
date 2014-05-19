#ifndef VECTOR5_HH
#define VECTOR5_HH

#include <ostream>
#include <gsl/gsl_vector.h>
#include <vector>

namespace aidaTT
{

    /** this is a simple abstraction of a five-element vector that will work with the matrix class
     *
     * it is a straightforward encapsulation of the gsl vector class with fixed size
     *
     */


    ///~ a five element vector that abstracts the gsl vector
    class Vector5
    {
            friend class fiveByFiveMatrix;
        public:
            /** the default construction, it initializes all entries to zero **/
            Vector5();

            /** copy constructor **/
            Vector5(const Vector5&);

            /** the construction with a vector
             **/
            Vector5(const std::vector<double>&);

            /** the construction from five doubles
             **/
            Vector5(double, double, double, double, double);

            /** the destructor **/
            ~Vector5();

            /** assignment operator **/
            Vector5& operator=(const Vector5&);

            /** direct read access to the individual elements by index**/
            double operator()(unsigned int index) const;

            /** direct write access to the individual elements by index**/
            double& operator()(unsigned int index);

        private:
            gsl_vector* _vector;
    };


    inline std::ostream & operator << (std::ostream & os, const Vector5& v)
    {
        os << " {" << v(0) << " , " << v(1) << " , " << v(2) << " , " << v(3) << " , " << v(4)  << "} "  ;
        return os ;
    }
}
#endif // VECTOR5_HH
