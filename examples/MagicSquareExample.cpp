
/** Example of use of Matrix Class, featuring magic squares. **/

#include <iostream>
#include <string>
#include <iomanip>
#include <limits>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "QRDecomposition.hpp"
#include "LUDecomposition.hpp"
#include "SingularValueDecomposition.hpp"
#include "CholeskyDecomposition.hpp"
#include "EigenvalueDecomposition.hpp"

using namespace boost::numeric::ublas;
using std::cout;
using std::string;
using std::setw;
using std::setprecision;
using std::fixed;
using std::endl;

typedef matrix<double> Matrix;
typedef symmetric_matrix<double> SymmetricMatrix;
typedef identity_matrix<double> IdentityMatrix;
typedef scalar_matrix<double> ScalarMatrix;
typedef vector<double> Vector;
typedef vector<int> PivotVector;

   /** Generate magic square test matrix. **/

   static Matrix magic (int n) {

      Matrix M(n,n);

      // Odd order

      if ((n % 2) == 1) {
         int a = (n+1)/2;
         int b = (n+1);
         for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
               M(i,j) = n*((i+j+a) % n) + ((i+2*j+b) % n) + 1;
            }
         }

      // Doubly Even Order

      } else if ((n % 4) == 0) {
         for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
               if (((i+1)/2)%2 == ((j+1)/2)%2) {
                  M(i,j) = n*n-n*i-j;
               } else {
                  M(i,j) = n*i+j+1;
               }
            }
         }

      // Singly Even Order

      } else {
         int p = n/2;
         int k = (n-2)/4;
         Matrix A = magic(p);
         for (int j = 0; j < p; j++) {
            for (int i = 0; i < p; i++) {
               double aij = A(i,j);
               M(i,j) = aij;
               M(i,j+p) = aij + 2*p*p;
               M(i+p,j) = aij + 3*p*p;
               M(i+p,j+p) = aij + p*p;
            }
         }
         for (int i = 0; i < p; i++) {
            for (int j = 0; j < k; j++) {
               double t = M(i,j); M(i,j) = M(i+p,j); M(i+p,j) = t;
            }
            for (int j = n-k+1; j < n; j++) {
               double t = M(i,j); M(i,j) = M(i+p,j); M(i+p,j) = t;
            }
         }
         double t = M(k,0); M(k,0) = M(k+p,0); M(k+p,0) = t;
         t = M(k,k); M(k,k) = M(k+p,k); M(k+p,k) = t;
      }
      return M;
   }
   
   /** Shorten spelling of print. **/

   static void print (string s) {
       cout << s;
   }

   /** Format double with Fw.d. **/

   static string fixedWidthDoubletoString (double x, int w, int d) {
      std::ostringstream oss;
      oss << setw(w) << fixed << setprecision(d) << x;
      return oss.str();
   }

   /** Format integer with Iw. **/

   static string fixedWidthIntegertoString (int n, int w) {
      std::ostringstream oss;
      oss << setw(w) << n;
      return oss.str();
   }


   template<class matrix_type>
   static typename matrix_type::value_type trace(const matrix_type& M)
   {
       typename matrix_type::size_type n = std::min(M.size1(), M.size2());
       matrix_vector_range<matrix_type> diag(const_cast<matrix_type&>(M), range (0,n), range (0,n));
       return sum(diag);
   }

int main (int argc, char **argv) {

   /* 
    | Tests LU, QR, SVD and symmetric Eig decompositions.
    |
    |   n       = order of magic square.
    |   trace   = diagonal sum, should be the magic sum, (n^3 + n)/2.
    |   max_eig = maximum eigenvalue of (A + A')/2, should equal trace.
    |   rank    = linear algebraic rank,
    |             should equal n if n is odd, be less than n if n is even.
    |   cond    = L_2 condition number, ratio of singular values.
    |   lu_res  = test of LU factorization, norm1(L*U-A(p,:))/(n*eps).
    |   qr_res  = test of QR factorization, norm1(Q*R-A)/(n*eps).
    */

      print("\n    Test of Matrix Class, using magic squares.\n");
      print("    See MagicSquareExample.main() for an explanation.\n");
      print("\n      n     trace       max_eig   rank        cond      lu_res      qr_res\n\n");
 
      //Date start_time = new Date();
      double eps = std::numeric_limits<double>::epsilon();
      for (int n = 3; n <= 32; n++) {
         print(fixedWidthIntegertoString(n,7));

         Matrix M = magic(n);

         int t = (int) trace(M);
         print(fixedWidthIntegertoString(t,10));

         EigenvalueDecomposition<double> E(SymmetricMatrix(0.5*(M+trans(M))));
         Vector d = E.getRealEigenvalues();
         print(fixedWidthDoubletoString(d(n-1),14,3));

         SingularValueDecomposition<double> SVD(M);
         int r = SVD.rank();
         print(fixedWidthIntegertoString(r,7));

         double c = SVD.cond();
         print(c < 1/eps ? fixedWidthDoubletoString(c,12,3) :
            "         Inf");

         LUDecomposition LU(M);
         Matrix L = LU.getL();
         Matrix U = LU.getU();
         PivotVector p = LU.getPivot();
         // cout << "M=" << M << endl;
         // cout << "L=" << L << endl;
         // cout << "U=" << U << endl;
         // cout << "p=" << p << endl;
         Matrix R = prod(L,U);
         for(int i=0; i<(int)p.size(); i++)
             row(R,i) -= row(M,p(i));
         double res = norm_1(R)/(n*eps);
         print(fixedWidthDoubletoString(res,12,3));

         QRDecomposition QR(M);
         Matrix Q = QR.getQ();
         R = QR.getR();
         R = prod(Q,R) - M;
         res = norm_1(R)/(n*eps);
         print(fixedWidthDoubletoString(res,12,3));

         print("\n");
      }
      //Date stop_time = new Date();
      //double etime = (stop_time.getTime() - start_time.getTime())/1000.;
      //print("\nElapsed Time = " + 
      //   fixedWidthDoubletoString(etime,12,3) + " seconds\n");
      print("Adios\n");
   }
