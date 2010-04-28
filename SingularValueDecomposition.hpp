   /** Singular Value Decomposition.
   <P>
   For an m-by-n matrix A, the singular value decomposition is
   an m-by-(m or n) orthogonal matrix U, an (m or n)-by-n diagonal matrix S, and
   an n-by-n orthogonal matrix V so that A = U*S*V'.
   <P>
   The singular values, sigma[k] = S[k][k], are ordered so that
   sigma[0] >= sigma[1] >= ... >= sigma[n-1].
   <P>
   The singular value decompostion always exists, so the constructor will
   never fail.  The matrix condition number and the effective numerical
   rank can be computed from this decomposition.
   */

// This version includes modifications and fixes by Andreas Kyrmegalos
// explanation: http://cio.nist.gov/esd/emaildir/lists/jama/msg01430.html
// final version: http://cio.nist.gov/esd/emaildir/lists/jama/msg01431.html

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
   @serial U column dimension.
   */
    int m, n, ncu;

   /** Column specification of matrix U
   @serial U column dimension toggle
   */
   
   bool thin;
   
   /** Construct the singular value decomposition
   @param A    Rectangular matrix
   @param thin  If true U is economy sized
   @param wantu If true generate the U matrix
   @param wantv If true generate the V matrix
   @return     Structure to access U, S and V.
   */
   void init (const Matrix &Arg, bool thin, bool wantu, bool wantv);

public:
/* ------------------------
   Old Constructor
 * ------------------------ */
   /** Construct the singular value decomposition
   @param Arg  Rectangular matrix
   @return     Structure to access U, S and V.
   */
   
   SingularValueDecomposition (const Matrix &Arg) {
      init(Arg,true,true,true);
   }
   
/* ------------------------
   Constructor
 * ------------------------ */

   /** Construct the singular value decomposition
   @param A    Rectangular matrix
   @param thin  If true U is economy sized
   @param wantu If true generate the U matrix
   @param wantv If true generate the V matrix
   @return     Structure to access U, S and V.
   */

   SingularValueDecomposition (const Matrix &Arg, bool thin, bool wantu, bool wantv) {
      init(Arg,thin,wantu,wantv);
   }
    
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
      Matrix S(m>=n?(thin?n:ncu):ncu,n);
      S.clear();
      for (int i = std::min(m,n)-1; i >= 0; i--) {
         S(i,i) = s(i);
      }
      return S;
   }

   /** Return the diagonal matrix of the reciprocals of the singular values
   @return     S+
   */
   
   Matrix getreciprocalS () const {
      Matrix S(n,m>=n?(thin?n:ncu):ncu);
      S.clear();
      for (int i = std::min(m,n)-1; i>=0; i--)
         S(i,i) = s(i)==0.0?0.0:1.0/s(i);
      return S;
   }
   
   /** Return the Moore-Penrose (generalized) inverse
    *  Slightly modified version of Kim van der Linde's code
   @param omit if true tolerance based omitting of negligible singular values
   @return     A+
   */
   
   Matrix inverse(bool omit) const {
      Matrix inverse(n,m);
      if(rank()> 0) {
         Vector reciprocalS(s.size());
         if (omit) {
            double tol = std::max(m,n)*s(0)*std::pow(2.0,-52);
            for (int i = s.size()-1;i>=0;i--)
               reciprocalS(i) = std::abs(s(i))<tol?0.0:1.0/s(i);
         }
         else
            for (int i=s.size()-1;i>=0;i--)
               reciprocalS(i) = s(i)==0.0?0.0:1.0/s(i);
         int min = std::min(n, ncu);
         for (int i = n-1; i >= 0; i--)
            for (int j = m-1; j >= 0; j--)
               for (int k = min-1; k >= 0; k--)
                  inverse(i,j) += V(i,k) * reciprocalS(k) * U(j,k);
      } 
      return inverse;
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
      return s(0)/s(std::min(m,n)-1);
   }

   /** Effective numerical matrix rank
   @return     Number of nonnegligible singular values.
   */

    int rank () const;    
};

}}}

#endif
