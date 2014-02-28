#ifndef FIVEBYFIVEMATRIX_HH
#define FIVEBYFIVEMATRIX_HH

#include <gsl/gsl_matrix.h>
#include <vector>

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
        public:
            /** the default construction, it initializes all entries to zero **/
            fiveByFiveMatrix();
            
            /** the construction with a vector, this is ROW wise; 
             * the sixth element of the vector is the first in the second row!
             **/
            fiveByFiveMatrix(const std::vector<double>&);
            
            /** the destructor **/
            ~fiveByFiveMatrix();

            /** direct read access to the individual matrix elements by index**/
            double operator()(unsigned int row, unsigned int column) const;

            /** direct write access to the individual matrix elements by index**/
            double& operator()(unsigned int row, unsigned int column);

            /** make a unit matrix out of the given matrix **/
            inline void Unit() { gsl_matrix_set_identity(_matrix); };

        private:
            gsl_matrix* _matrix;
    };

    typedef fiveByFiveMatrix fiveByFiveMatrixSym;
    typedef fiveByFiveMatrixSym fullCovariance;
}

#endif // FIVEBYFIVEMATRIX_HH
