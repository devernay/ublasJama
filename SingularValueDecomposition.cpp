// This version includes modifications and fixes by Andreas Kyrmegalos
// explanation: http://cio.nist.gov/esd/emaildir/lists/jama/msg01430.html
// final version: http://cio.nist.gov/esd/emaildir/lists/jama/msg01431.html

#include <cmath>
#include <boost/math/special_functions/hypot.hpp>
#include "SingularValueDecomposition.hpp"

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
          U = Matrix(m,ncu);
          U.clear();
      }
      if (wantv) {
          V = Matrix(n,n);
          V.clear();
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
            s(k) = 0;
            for (int i = k; i < m; i++) { // s(k) = norm of elements k..m-1 of column k of A
                s(k) = boost::math::hypot(s(k),A(i,k));
            }
            if (s(k) != 0.0) {
               if (A(k,k) < 0.0) {
                  s(k) = -s(k);
               }
               for (int i = k; i < m; i++) { // divide elements k..m-1 of column k of A by s(k)
                  A(i,k) /= s(k);
               }
               A(k,k) += 1.0;
            }
            s(k) = -s(k);
         }
         for (int j = k+1; j < n; j++) {
            if ((k < nct) && (s(k) != 0.0))  {

            // Apply the transformation.

               double t = 0;
               for (int i = k; i < m; i++) { // t = dot-product of elements k..m-1 of columns k and j of A
                  t += A(i,k)*A(i,j);
               }
               t = -t/A(k,k);
               for (int i = k; i < m; i++) { // elements k..m-1 of column j of A +=  t*(elements k..m-1 of column k of A)
                  A(i,j) += t*A(i,k);
               }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e(j) = A(k,j);
         }
         if (wantu && (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (int i = k; i < m; i++) { // elements k..m-1 of column k of U = elements k..m-1 of column k of A
               U(i,k) = A(i,k);
            }
         }
         if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e(k) = 0;
            for (int i = k+1; i < n; i++) { // e(k) = norm of elements k+1..n-1 of e
                e(k) = boost::math::hypot(e[k],e[i]);
            }
            if (e(k) != 0.0) {
               if (e(k+1) < 0.0) {
                  e(k) = -e(k);
               }
               for (int i = k+1; i < n; i++) { // divide elements k+1..n-1 of e by e(k)
                  e(i) /= e(k);
               }
               e(k+1) += 1.0;
            }
            e(k) = -e(k);
            if ((k+1 < m) && (e(k) != 0.0)) {

            // Apply the transformation.

               for (int i = k+1; i < m; i++) {
                  work(i) = 0.0;
               }
               for (int j = k+1; j < n; j++) { // elements k+1..n-1 of work = A.submatrix(k+1..n-1,k+1..m-1)*elements k+1..n-1 of e
                  for (int i = k+1; i < m; i++) {
                     work(i) += e(j)*A(i,j);
                  }
               }
               for (int j = k+1; j < n; j++) {
                  double t = -e(j)/e(k+1);
                  for (int i = k+1; i < m; i++) { // elements k+1..m-1 of column j of A += t*elements k+1..m-1 of work
                     A(i,j) += t*work(i);
                  }
               }
            }
            if (wantv) {

            // Place the transformation in V for subsequent
            // back multiplication.

               for (int i = k+1; i < n; i++) { // elements k+1..n-1 of column k of V = elements k+1..n-1 of e
                  V(i,k) = e(i);
               }
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
            for (int i = 0; i < m; i++) { // set column j of U to zero
               U(i,j) = 0.0;
            }
            U(j,j) = 1.0;
         }
         for (int k = nct-1; k >= 0; k--) {
            if (s(k) != 0.0) {
               for (int j = k+1; j < ncu; j++) {
                  double t = 0;
                  for (int i = k; i < m; i++) { // t = dot-product of elements k..m-1 of columns k and j of U
                     t += U(i,k)*U(i,j);
                  }
                  t = -t/U(k,k);
                  for (int i = k; i < m; i++) { // elements k..m-1 of column j of U +=  t*(elements k..m-1 of column k of U)
                     U(i,j) += t*U(i,k);
                  }
               }
               for (int i = k; i < m; i++ ) { // elements k..m-1 of column k of U *= -1.
                  U(i,k) = -U(i,k);
               }
               U(k,k) += 1.0;
               if(k-1 > 0) {
                  for (int i = 0; i < k-1; i++) { // set elements 0..k-2 of column k of U to zero.
                     U(i,k) = 0.0;
                  }
               }
            } else {
               for (int i = 0; i < m; i++) { // set column k of U to zero
                  U(i,k) = 0.0;
               }
               U(k,k) = 1.0;
            }
         }
      }

      // If required, generate V.

      if (wantv) {
         for (int k = n-1; k >= 0; k--) {
            if ((k < nrt) && (e(k) != 0.0)) {
               for (int j = k+1; j < n; j++) {
                  double t = 0;
                  for (int i = k+1; i < n; i++) { // t = dot-product of elements k+1..n-1 of columns k and j of V
                      t += V(i,k)*V(i,j);
                  }
                  t = -t/V(k+1,k);
                  for (int i = k+1; i < n; i++) { // elements k+1..n-1 of column j of V +=  t*(elements k+1..n-1 of column k of V)
                     V(i,j) += t*V(i,k);
                  }
               }
            }
            for (int i = 0; i < n; i++) { // set column k of V to zero
               V(i,k) = 0.0;
            }
            V(k,k) = 1.0;
         }
      }

      // Main iteration loop for the singular values.

      int pp = p-1;
      int iter = 0;
      double eps = std::pow(2.0,-52);
      double tiny = std::pow(2.0,-966);
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
                     for (int i = 0; i < n; i++) { // multiply column k of V by -1
                        V(i,k) = -V(i,k);
                     }
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
                     for (int i = 0; i < n; i++) { // swap columns k and k+1 of V
                        t = V(i,k+1); V(i,k+1) = V(i,k); V(i,k) = t;
                     }
                  }
                  if (wantu && (k < m-1)) {
                     for (int i = 0; i < m; i++) { // swap columns k and k+1 of U
                        t = U(i,k+1); U(i,k+1) = U(i,k); U(i,k) = t;
                     }
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

int SingularValueDecomposition::rank () const {
      double tol = std::max(m,n)*s[0]*std::pow(2.0,-52);
      int r = 0;
      for (unsigned i = 0; i < s.size(); i++) {
         if (s(i) > tol) {
            r++;
         }
      }
      return r;
   }

}}}
