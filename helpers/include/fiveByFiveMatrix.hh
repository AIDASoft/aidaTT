#ifndef FIVEBYFIVEMATRIX_HH
#define FIVEBYFIVEMATRIX_HH

#include <gsl/gsl_matrix.h>
#include <vector>

#include <ostream>

#include "Vector5.hh"
#include "Vector3D.hh"

namespace aidaTT
{

    /** this is the main working horse class for a 5x5 matrix
     *
     * it is a straightforward encapsulation of the gsl matrix class with fixed size
     *
     * it controls resources, therefore assignment and copying needs to be done with care
     *
     */


    ///~ a five by five matrix that abstracts the gsl matrix (for now)
    class fiveByFiveMatrix
    {
            friend class IPropagation;

        public:
            /** the default construction, it initializes all entries to zero **/
            fiveByFiveMatrix();

            /** copy construction **/
            fiveByFiveMatrix(const fiveByFiveMatrix&);

            /** the construction with a vector, this is ROW wise;
             * the sixth element of the vector is the first in the second row!
             **/
            fiveByFiveMatrix(const std::vector<double>&);

            /** the destructor **/
            ~fiveByFiveMatrix();

            /** the assignment operator **/
            fiveByFiveMatrix& operator=(const fiveByFiveMatrix&);

            /** direct read access to the individual matrix elements by index**/
            double operator()(unsigned int row, unsigned int column) const;

            /** direct write access to the individual matrix elements by index**/
            double& operator()(unsigned int row, unsigned int column);

            /** define matrix-matrix multipliation **/
            virtual fiveByFiveMatrix operator*(const fiveByFiveMatrix&);

            /** define matrix-vector multiplication **/
            virtual Vector5 operator*(const Vector5&);

            /** make a unit matrix out of the given matrix **/
            inline void Unit()
            {
                gsl_matrix_set_identity(_matrix);
            };

            /** get array to construct other matrix representations */
            double* array() const;

        private:
            gsl_matrix* _matrix;
    };



    inline std::ostream & operator << (std::ostream & os, const fiveByFiveMatrix& m)
    {
        os << " [" << m(0, 0) << " , " << m(0, 1) << " , " << m(0, 2) << " , " << m(0, 3) << " , " << m(0, 4) << " / "
           << m(1, 0) << " , " << m(1, 1) << " , " << m(1, 2) << " , " << m(1, 3) << " , " << m(1, 4) << " / "
           << m(2, 0) << " , " << m(2, 1) << " , " << m(2, 2) << " , " << m(2, 3) << " , " << m(2, 4) << " / "
           << m(3, 0) << " , " << m(3, 1) << " , " << m(3, 2) << " , " << m(3, 3) << " , " << m(3, 4) << " / "
           << m(4, 0) << " , " << m(4, 1) << " , " << m(4, 2) << " , " << m(4, 3) << " , " << m(4, 4) << "]"  ;
        return os ;
    }



    typedef fiveByFiveMatrix fiveByFiveMatrixSym;
    typedef fiveByFiveMatrixSym fullCovariance;


    fiveByFiveMatrix curvilinearToPerigeeJacobian(const Vector5&, const Vector3D&);
    fiveByFiveMatrix perigeeToCurvilinearJacobian(const Vector5&, const Vector3D&);

    fiveByFiveMatrix perigeeToILDJacobian(const Vector5&);
    fiveByFiveMatrix ildToPerigeeJacobian(const Vector5&);

    fiveByFiveMatrix curvilinearToILDJacobian(const Vector5&, const Vector3D&);
    fiveByFiveMatrix ildToCurvilinearJacobian(const Vector5&, const Vector3D&);
}

#endif // FIVEBYFIVEMATRIX_HH
