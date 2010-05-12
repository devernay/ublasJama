// This version includes modifications and fixes by Andreas Kyrmegalos
// explanation: http://cio.nist.gov/esd/emaildir/lists/jama/msg01430.html
// final version: http://cio.nist.gov/esd/emaildir/lists/jama/msg01431.html

#include <cmath>
#include <limits>
#include <boost/math/special_functions/hypot.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include "SingularValueDecomposition.hpp"

//#define subcolumn(M,c,start,stop) subrange(column(M,c),start,stop)
//#define subrow(M,r,start,stop) subrange(row(M,r),start,stop)
#define subcolumn(M,c,start,stop) matrix_vector_slice<Matrix> ((M), slice((start),1,(stop)-(start)), slice((c),0,(stop)-(start)))
#define subrow(M,r,start,stop) matrix_vector_slice<Matrix> ((M), slice((r),0,(stop)-(start)), slice((start),1,(stop)-(start)))

namespace boost {
    namespace numeric {
        namespace ublas {

   /** sqrt(a^2 + b^2) without under/overflow. **/
   /*
   static inline double hypot(double a, double b) {
      double r;
      if (std::abs(a) > std::abs(b)) {
         r = b/a;
         r = std::abs(a)*std::sqrt(1+r*r);
      } else if (b != 0) {
         r = a/b;
         r = std::abs(b)*std::sqrt(1+r*r);
      } else {
         r = 0.0;
      }
      return r;
   }
   */

void SingularValueDecomposition::init (const Matrix &Arg, bool thin, bool wantu, bool wantv) {

      // Derived from LINPACK code.
      // Initialize.
      Matrix A = Arg;
      m = Arg.size1();
      n = Arg.size2();
      this->thin = thin;

      ncu = thin?std::min(m,n):m;
      s = Vector(std::min(m+1,n));
      if (wantu) {
          U = Matrix(m,ncu,0.0);
      }
      if (wantv) {
          V = Matrix(n,n,0.0);
      }
      Vector e(n);
      Vector work(m);

      // Reduce A to bidiagonal form, storing the diagonal elements
      // in s and the super-diagonal elements in e.

      int nct = std::min(m-1,n);
      int nrt = std::max(0,std::min(n-2,m));
      int lu = std::max(nct,nrt);
      for (int k = 0; k < lu; k++) {
         if (k < nct) {

            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            // s(k) = norm of elements k..m-1 of column k of A
            s(k) = norm_2(subcolumn(A,k,k,m));
            if (s(k) != 0.0) {
               if (A(k,k) < 0.0) {
                  s(k) = -s(k);
               }
               // divide elements k..m-1 of column k of A by s(k)
               subcolumn(A,k,k,m) /= s(k);
               A(k,k) += 1.0;
            }
            s(k) = -s(k);
         }
         for (int j = k+1; j < n; j++) {
            if ((k < nct) && (s(k) != 0.0))  {

            // Apply the transformation.

               // t = dot-product of elements k..m-1 of columns k and j of A
               double t = inner_prod(subcolumn(A,k,k,m),subcolumn(A,j,k,m));
               t = -t/A(k,k);
               // elements k..m-1 of column j of A +=  t*(elements k..m-1 of column k of A)
               subcolumn(A,j,k,m) += t*subcolumn(A,k,k,m);
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e(j) = A(k,j);
         }
         if (wantu && (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            // elements k..m-1 of column k of U = elements k..m-1 of column k of A
            subcolumn(U,k,k,m) = subcolumn(A,k,k,m);
         }
         if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            // e(k) = norm of elements k+1..n-1 of e
            e(k) = norm_2(subrange(e,k+1,n));
            if (e(k) != 0.0) {
               if (e(k+1) < 0.0) {
                  e(k) = -e(k);
               }
               // divide elements k+1..n-1 of e by e(k)
               subrange(e,k+1,n) /= e(k);
               e(k+1) += 1.0;
            }
            e(k) = -e(k);
            if ((k+1 < m) && (e(k) != 0.0)) {

            // Apply the transformation.

               for (int i = k+1; i < m; i++) {
                  work(i) = 0.0;
               }
               // elements k+1..m-1 of work = A.submatrix(k+1..m-1,k+1..n-1)*elements k+1..n-1 of e
               noalias(subrange(work, k+1, m)) = prod(subrange(A,k+1,m,k+1,n),subrange(e,k+1,n));
               for (int j = k+1; j < n; j++) {
                  double t = -e(j)/e(k+1);
                  // elements k+1..m-1 of column j of A += t*elements k+1..m-1 of work
                  subcolumn(A,j,k+1,m) += t*subrange(work,k+1,m);
               }
            }
            if (wantv) {

            // Place the transformation in V for subsequent
            // back multiplication.

               // elements k+1..n-1 of column k of V = elements k+1..n-1 of e
               subcolumn(V,k,k+1,n) = subrange(e,k+1,n);
            }
         }
      }

      // Set up the final bidiagonal matrix or order p.

      int p = std::min(n,m+1);
      if (nct < n) {
         s(nct) = A(nct,nct);
      }
      if (m < p) {
         s(p-1) = 0.0;
      }
      if (nrt+1 < p) {
         e(nrt) = A(nrt,p-1);
      }
      e(p-1) = 0.0;

      // If required, generate U.

      if (wantu) {
         for (int j = nct; j < ncu; j++) {
            // set column j of U to zero
            std::fill(column(U,j).begin(),column(U,j).end(),double(/*zero*/));
            U(j,j) = 1.0;
         }
         for (int k = nct-1; k >= 0; k--) {
            if (s(k) != 0.0) {
               for (int j = k+1; j < ncu; j++) {
                  // t = dot-product of elements k..m-1 of columns k and j of U
                  double t = inner_prod(subcolumn(U,k,k,m),subcolumn(U,j,k,m));
                  t /= -U(k,k);
                  // elements k..m-1 of column j of U +=  t*(elements k..m-1 of column k of U)
                  subcolumn(U,j,k,m) += t*subcolumn(U,k,k,m);
               }
               // elements k..m-1 of column k of U *= -1.
               subcolumn(U,k,k,m) *= -1.0;
               U(k,k) += 1.0;
               if(k-1 > 0) {
                  // set elements 0..k-2 of column k of U to zero.
                  for (int i = 0; i < k-1; i++) { 
                     U(i,k) = 0.0;
                  }
               }
            } else {
               // set column k of U to zero
               std::fill(column(U,k).begin(),column(U,k).end(),double(/*zero*/));
               U(k,k) = 1.0;
            }
         }
      }

      // If required, generate V.

      if (wantv) {
         for (int k = n-1; k >= 0; k--) {
            if ((k < nrt) && (e(k) != 0.0)) {
               for (int j = k+1; j < n; j++) {
                  // t = dot-product of elements k+1..n-1 of columns k and j of V
                  double t = inner_prod(subcolumn(V,k,k+1,n),subcolumn(V,j,k+1,n));
                  t /= -V(k+1,k);
                  // elements k+1..n-1 of column j of V +=  t*(elements k+1..n-1 of column k of V)
                  subcolumn(V,j,k+1,n) += t*subcolumn(V,k,k+1,n);
               }
            }
            // set column k of V to zero
            std::fill(column(V,k).begin(),column(V,k).end(),double(/*zero*/));
            V(k,k) = 1.0;
         }
      }

      // Main iteration loop for the singular values.

      int pp = p-1;
      int iter = 0;
      double eps = std::numeric_limits<double>::epsilon();
      double tiny = std::numeric_limits<double>::min();
      while (p > 0) {
         int k,kase;

         // Here is where a test for too many iterations would go.

         // This section of the program inspects for
         // negligible elements in the s and e arrays.  On
         // completion the variables kase and k are set as follows.

         // kase = 1     if s(p) and e[k-1] are negligible and k<p
         // kase = 2     if s(k) is negligible and k<p
         // kase = 3     if e[k-1] is negligible, k<p, and
         //              s(k), ..., s(p) are not negligible (qr step).
         // kase = 4     if e(p-1) is negligible (convergence).

         for (k = p-2; k >= -1; k--) {
            if (k == -1) {
               break;
            }
            if (std::abs(e(k)) <=
                  tiny + eps*(std::abs(s(k)) + std::abs(s(k+1)))) {
               e(k) = 0.0;
               break;
            }
         }
         if (k == p-2) {
            kase = 4;
         } else {
            int ks;
            for (ks = p-1; ks >= k; ks--) {
               if (ks == k) {
                  break;
               }
               double t = (ks != p ? std::abs(e(ks)) : 0.) + 
                          (ks != k+1 ? std::abs(e(ks-1)) : 0.);
               if (std::abs(s(ks)) <= tiny + eps*t)  {
                  s[ks] = 0.0;
                  break;
               }
            }
            if (ks == k) {
               kase = 3;
            } else if (ks == p-1) {
               kase = 1;
            } else {
               kase = 2;
               k = ks;
            }
         }
         k++;

         // Perform the task indicated by kase.

         switch (kase) {

            // Deflate negligible s(p).

            case 1: {
               double f = e(p-2);
               e[p-2] = 0.0;
               for (int j = p-2; j >= k; j--) {
                  double t = boost::math::hypot(s(j),f);
                  double cs = s(j)/t;
                  double sn = f/t;
                  s(j) = t;
                  if (j != k) {
                     f = -sn*e(j-1);
                     e(j-1) = cs*e(j-1);
                  }
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs*V(i,j) + sn*V(i,p-1);
                        V(i,p-1) = -sn*V(i,j) + cs*V(i,p-1);
                        V(i,j) = t;
                     }
                  }
               }
            }
            break;

            // Split at negligible s(k).

            case 2: {
               double f = e(k-1);
               e[k-1] = 0.0;
               for (int j = k; j < p; j++) {
                  double t = boost::math::hypot(s(j),f);
                  double cs = s(j)/t;
                  double sn = f/t;
                  s(j) = t;
                  f = -sn*e(j);
                  e(j) = cs*e(j);
                  if (wantu) {
                     for (int i = 0; i < m; i++) {
                        t = cs*U(i,j) + sn*U(i,k-1);
                        U(i,k-1) = -sn*U(i,j) + cs*U(i,k-1);
                        U(i,j) = t;
                     }
                  }
               }
            }
            break;

            // Perform one qr step.

            case 3: {

               // Calculate the shift.
   
               double scale = std::max(std::max(std::max(std::max(
                       std::abs(s(p-1)),std::abs(s(p-2))),std::abs(e(p-2))), 
                       std::abs(s(k))),std::abs(e(k)));
               double sp = s(p-1)/scale;
               double spm1 = s(p-2)/scale;
               double epm1 = e(p-2)/scale;
               double sk = s(k)/scale;
               double ek = e(k)/scale;
               double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
               double c = (sp*epm1)*(sp*epm1);
               double shift = 0.0;
               if ((b != 0.0) || (c != 0.0)) {
                  shift = std::sqrt(b*b + c);
                  if (b < 0.0) {
                     shift = -shift;
                  }
                  shift = c/(b + shift);
               }
               double f = (sk + sp)*(sk - sp) + shift;
               double g = sk*ek;
   
               // Chase zeros.
   
               for (int j = k; j < p-1; j++) {
                  double t = boost::math::hypot(f,g);
                  double cs = f/t;
                  double sn = g/t;
                  if (j != k) {
                     e(j-1) = t;
                  }
                  f = cs*s(j) + sn*e(j);
                  e(j) = cs*e(j) - sn*s(j);
                  g = sn*s(j+1);
                  s(j+1) = cs*s(j+1);
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs*V(i,j) + sn*V(i,j+1);
                        V(i,j+1) = -sn*V(i,j) + cs*V(i,j+1);
                        V(i,j) = t;
                     }
                  }
                  t = boost::math::hypot(f,g);
                  cs = f/t;
                  sn = g/t;
                  s(j) = t;
                  f = cs*e(j) + sn*s(j+1);
                  s(j+1) = -sn*e(j) + cs*s(j+1);
                  g = sn*e(j+1);
                  e(j+1) = cs*e(j+1);
                  if (wantu && (j < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = cs*U(i,j) + sn*U(i,j+1);
                        U(i,j+1) = -sn*U(i,j) + cs*U(i,j+1);
                        U(i,j) = t;
                     }
                  }
               }
               e(p-2) = f;
               iter++;
            }
            break;

            // Convergence.

            case 4: {
               // Make the singular values positive.
   
               if (s(k) <= 0.0) {
                  s(k) = (s(k) < 0.0 ? -s(k) : 0.0);
                  if (wantv) {
                     // multiply column k of V by -1
                     column(V,k) *= -1.0;
                  }
               }
   
               // Order the singular values.
   
               while (k < pp) {
                  if (s(k) >= s(k+1)) {
                     break;
                  }
                  double t = s(k);
                  s(k) = s(k+1);
                  s(k+1) = t;
                  if (wantv && (k < n-1)) {
                     // swap columns k and k+1 of V
                     column(V,k).swap(column(V,k+1));
                  }
                  if (wantu && (k < m-1)) {
                     // swap columns k and k+1 of U
                     column(U,k).swap(column(U,k+1));
                  }
                  k++;
               }
               iter = 0;
               p--;
            }
            break;
         }
      }
   }

   /** Effective numerical matrix rank
   @return     Number of nonnegligible singular values.
   */

unsigned int SingularValueDecomposition::rank () const {
      double tol = std::max(m,n)*s[0]*std::numeric_limits<double>::epsilon();
      int r = 0;
      for (unsigned i = 0; i < s.size(); i++) {
         if (s(i) > tol) {
            r++;
         }
      }
      return r;
   }

}}}
