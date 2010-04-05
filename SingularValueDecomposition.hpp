   /** Singular Value Decomposition.
   <P>
   For an m-by-n matrix A with m >= n, the singular value decomposition is
   an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
   an n-by-n orthogonal matrix V so that A = U*S*V'.
   <P>
   The singular values, sigma[k] = S[k][k], are ordered so that
   sigma[0] >= sigma[1] >= ... >= sigma[n-1].
   <P>
   The singular value decompostion always exists, so the constructor will
   never fail.  The matrix condition number and the effective numerical
   rank can be computed from this decomposition.
   */

#ifndef _BOOST_UBLAS_SINGULARVALUEDECOMPOSITION_
#define _BOOST_UBLAS_SINGULARVALUEDECOMPOSITION_

#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/exception.hpp>

namespace boost { namespace numeric { namespace ublas {
            
class SingularValueDecomposition {

    typedef vector<double> Vector;
    typedef matrix<double> Matrix;

/* ------------------------
   Class variables
 * ------------------------ */

   /** Arrays for internal storage of U and V.
   @serial internal storage of U.
   @serial internal storage of V.
   */
    Matrix U, V;

   /** Array for internal storage of singular values.
   @serial internal storage of singular values.
   */
   Vector s;

   /** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
   */
   int m, n;

public:
/* ------------------------
   Constructor
 * ------------------------ */

   /** Construct the singular value decomposition
   @param A    Rectangular matrix
   @return     Structure to access U, S and V.
   */

    SingularValueDecomposition (const Matrix &Arg, bool wantu = true, bool wantv = true);
    
/* ------------------------
   Public Methods
 * ------------------------ */

   /** Return the left singular vectors
   @return     U
   */

   const Matrix& getU () const {
      return U;
   }

   /** Return the right singular vectors
   @return     V
   */

   const Matrix& getV () const {
      return V;
   }

   /** Return the one-dimensional array of singular values
   @return     diagonal of S.
   */

   const Vector& getSingularValues () const {
      return s;
   }

   /** Return the diagonal matrix of singular values
   @return     S
   */

   Matrix getS () const {
      Matrix S(n,n);
      S.clear();
      for (int i = 0; i < n; i++) {
         S(i,i) = s(i);
      }
      return S;
   }

   /** Two norm
   @return     max(S)
   */

   double norm2 () const {
      return s(0);
   }

   /** Two norm condition number
   @return     max(S)/min(S)
   */

   double cond () const {
      return s(0)/s[std::min(m,n)-1];
   }

   /** Effective numerical matrix rank
   @return     Number of nonnegligible singular values.
   */

    int rank () const;    
};

}}}

#endif
