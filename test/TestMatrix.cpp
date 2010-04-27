/** TestMatrix tests the functionality of the Jama Matrix class and associated decompositions.
<P>
Run the test from the command line using
<BLOCKQUOTE><PRE><CODE>
 java Jama.test.TestMatrix 
</CODE></PRE></BLOCKQUOTE>
Detailed output is provided indicating the functionality being tested
and whether the functionality is correctly implemented.   Exception handling
is also tested.  
<P>
The test is designed to run to completion and give a summary of any implementation errors
encountered. The final output should be:
<BLOCKQUOTE><PRE><CODE>
      TestMatrix completed.
      Total errors reported: n1
      Total warning reported: n2
</CODE></PRE></BLOCKQUOTE>
If the test does not run to completion, this indicates that there is a 
substantial problem within the implementation that was not anticipated in the test design.  
The stopping point should give an indication of where the problem exists.
**/
#include <iostream>
#include <iomanip>
#include <boost/numeric/ublas/io.hpp>
#include "QRDecomposition.hpp"
#include "LUDecomposition.hpp"
#include "SingularValueDecomposition.hpp"
#include "CholeskyDecomposition.hpp"
#include "EigenvalueDecomposition.hpp"

using namespace boost::numeric::ublas;
using std::cout;
using std::endl;
using std::string;
using std::setw;
using std::setprecision;
using std::fixed;

typedef matrix<double> Matrix;
typedef identity_matrix<double> IdentityMatrix;
typedef scalar_matrix<double> ScalarMatrix;
typedef vector<double> Vector;
typedef vector<int> PivotVector;

/** private utility routines **/

/** Check magnitude of difference of scalars. **/

static void check(double x, double y) {
    double eps = std::pow(2.0,-52.0);
    if (x == 0 && std::abs(y) < 10*eps) return;
    if (y == 0 && std::abs(x) < 10*eps) return;
    if (std::abs(x-y) > 10*eps*std::max(std::abs(x),std::abs(y))) {
        std::ostringstream oss;
        oss << "The difference x-y is too large: x = " << x << "  y = " << y;
        throw internal_logic(oss.str().c_str());
    }
}

/** Check norm of difference of "vectors". **/

#if 0
static void check(const Vector& x, const Vector& y) {
    BOOST_UBLAS_CHECK(x.size() == y.size(), bad_size("Attempt to compare vectors of different lengths"));
        for (unsigned i=0;i<x.size();i++) {
            check(x[i],y[i]);
        } 
}
#endif

/** Check norm of difference of Matrices. **/

static void check(const Matrix& X, const Matrix& Y) {
    double eps = std::pow(2.0,-52.0);
    if (norm_1(X) == 0. && norm_1(Y) < 10*eps) return;
    if (norm_1(Y) == 0. && norm_1(X) < 10*eps) return;
    if (norm_1(X-Y) > 1000*eps*std::max(norm_1(X),norm_1(Y))) {
        std::ostringstream oss;
        oss << "The norm of (X-Y) is too large: " << norm_1(X-Y);
        throw internal_logic(oss.str().c_str());
    }
}

/** Print appropriate messages for successful outcome try **/

static void try_success (string s,string e) {
    cout << ">    " << s << "success\n";
    if ( e != "" ) {
        cout << ">      Message: " << e << "\n";
    }
}
/** Print appropriate messages for unsuccessful outcome try **/

static int try_failure (int count, string s,string e) {
    cout << ">    " << s << "*** failure ***\n>      Message: " << e << "\n";
    return ++count;
}

/** Print appropriate messages for unsuccessful outcome try **/

#if 0
static int try_warning (int count, string s,string e) {
    cout << ">    " << s << "*** warning ***\n>      Message: " << e << "\n";
    return ++count;
}
#endif

/** Print a row vector. **/

#if 1
static void print(Vector x, int w, int d) {
    // Use format Fw.d for all elements.
    cout << "\n";
    for(unsigned i=0; i<x.size(); i++) {
        cout << setw(w) << fixed << setprecision(d) << x(i);
    }
    cout << "\n";
}

/** Print a matrix. **/

static void print(Matrix m, int w, int d) {
    // Use format Fw.d for all elements.
    cout << "[\n";
    for(unsigned i=0; i<m.size1(); i++) {
        print(matrix_row<Matrix>(m,i),w,d);
    }
    cout << "]\n";
}
#endif

int main (int argc, char **argv) {
      Matrix A,B,C,Z,O,I,R,S,X,SUB,M,T,SQ,DEF,SOL;
      int errorCount=0;
      int warningCount=0;
    
      double columnwise[12] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.};
      double rankdef[3][4] = {{1.,4.,7.,10.},{2.,5.,8.,11.},{3.,6.,9.,12.}};
      double subavals[2][3] = {{5.,8.,11.},{6.,9.,12.}};
      double pvals[3][3] = {{4.,1.,1.},{1.,2.,3.},{1.,3.,6.}};
      double evals[4][4] = 
         {{0.,1.,0.,0.},{1.,0.,2.e-7,0.},{0.,-2.e-7,0.,1.},{0.,0.,1.,0.}};
      double sqSolution[2][1] = {{13.},{15.}};
      double condmat[2][2] = {{1.,3.},{7.,9.}};

/**
      LA methods:
         transpose
         times
         cond
         rank
         det
         trace
         norm1
         norm2
         normF
         normInf
         solve
         solveTranspose
         inverse
         chol
         eig
         lu
         qr
         svd 
**/

      cout << "\nTesting linear algebra methods...\n";

      A = Matrix(4,3);
      for(unsigned i=0; i<A.size1(); i++) {
         for(unsigned j=0; j<A.size2(); j++) {
            A(i,j) = columnwise[i+j*4];
         }
      }
      QRDecomposition QR(A);
      R = QR.getR();
      try {
         check(A,prod(QR.getQ(),R));
         try_success("QRDecomposition...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"QRDecomposition...","incorrect QR decomposition calculation");
      }
      SingularValueDecomposition SVD(A);
      try {
         Matrix US = prod(SVD.getU(),SVD.getS());
         check(A,prod(US,trans(SVD.getV())));
         try_success("SingularValueDecomposition...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"SingularValueDecomposition...","incorrect singular value decomposition calculation");
      }
      SingularValueDecomposition SVDID(IdentityMatrix(3,3));
      try {
          cout << "U=";
          print(SVDID.getU(),6,2);
          cout << "S=";
          print(SVDID.getS(),6,2);
          cout << "V=";
          print(SVDID.getV(),6,2);
          Matrix US = prod(SVDID.getU(),SVDID.getS());
         check(IdentityMatrix(3,3),prod(US,trans(SVDID.getV())));
         try_success("SingularValueDecomposition(Identity33)...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"SingularValueDecomposition(Identity33)...","incorrect singular value decomposition calculation");
      }
      DEF = Matrix(3,4);
      for(unsigned i=0; i<DEF.size1(); i++) {
         for(unsigned j=0; j<DEF.size2(); j++) {
            DEF(i,j) = rankdef[i][j];
         }
      }
      SVD = SingularValueDecomposition(trans(DEF)); // SVD only works for m >= n 
      try {
         check(SVD.rank(),std::min(DEF.size1(),DEF.size2())-1);
         try_success("rank()...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"rank()...","incorrect rank calculation");
      }
      B = Matrix(2,2);
      for(unsigned i=0; i<B.size1(); i++) {
         for(unsigned j=0; j<B.size2(); j++) {
            B(i,j) = condmat[i][j];
         }
      }
      SVD = SingularValueDecomposition(B); 
      Vector singularvalues = SVD.getSingularValues();
      try {
         check(SVD.cond(),singularvalues(0)/singularvalues(std::min(B.size1(),B.size2())-1));
         try_success("cond()...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"cond()...","incorrect condition number calculation");
      }
      int n = A.size2();
      A = subrange(A,0,n,0,n);
      A(0,0) = 0.;
      LUDecomposition LU(A);
      try {
         // Compute the pivoted A
         B = Matrix(A.size1(),A.size2());
         PivotVector piv = LU.getPivot();
         for(unsigned i=0; i<piv.size(); i++) {
            row(B,i) = row(A,piv(i));
         }
         check(B,prod(LU.getL(),LU.getU()));
         try_success("LUDecomposition...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"LUDecomposition...","incorrect LU decomposition calculation");
      }
      
      QR = QRDecomposition(A);
      X = QR.inverse();
      try {
         check(prod(A,X),IdentityMatrix(3,3));
         try_success("inverse()...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"inverse()...","incorrect inverse calculation");
      }
      SUB = Matrix(2,3);
      for(unsigned i=0; i<SUB.size1(); i++) {
         for(unsigned j=0; j<SUB.size2(); j++) {
            SUB(i,j) = subavals[i][j];
         }
      }
      O = ScalarMatrix(SUB.size1(),1,1.0);
      SOL = Matrix(2,1);
      for(unsigned i=0; i<SOL.size1(); i++) {
        for(unsigned j=0; j<SOL.size2(); j++) {
           SOL(i,j) = sqSolution[i][j];
        }
      }
      SQ = subrange(SUB,0,SUB.size1(),0,SUB.size1());
      try {
         check(QRDecomposition(SQ).solve(SOL),O); 
         try_success("solve()...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"solve()...",e.what());
      }
      A = Matrix(3,3);
      for(unsigned i=0; i<A.size1(); i++) {
         for(unsigned j=0; j<A.size2(); j++) {
            A(i,j) = pvals[i][j];
         }
      }
      CholeskyDecomposition Chol(A); 
      Matrix L = Chol.getL();
      try {
         check(A,prod(L,trans(L)));
         try_success("CholeskyDecomposition...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"CholeskyDecomposition...","incorrect Cholesky decomposition calculation");
      }
      X = Chol.solve(IdentityMatrix(3,3));
      try {
         check(prod(A,X),IdentityMatrix(3,3));
         try_success("CholeskyDecomposition solve()...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"CholeskyDecomposition solve()...","incorrect Choleskydecomposition solve calculation");
      }
      EigenvalueDecomposition Eig(A);
      Matrix D = Eig.getD();
      Matrix V = Eig.getV();
      try {
         check(prod(A,V),prod(V,D));
         try_success("EigenvalueDecomposition (symmetric)...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"EigenvalueDecomposition (symmetric)...","incorrect symmetric Eigenvalue decomposition calculation");
      }
      A = Matrix(4,4);
      for(unsigned i=0; i<A.size1(); i++) {
         for(unsigned j=0; j<A.size2(); j++) {
            A(i,j) = evals[i][j];
         }
      }
      Eig = EigenvalueDecomposition(A);
      D = Eig.getD();
      V = Eig.getV();
      try {
         check(prod(A,V),prod(V,D));
         try_success("EigenvalueDecomposition (nonsymmetric)...","");
      } catch ( std::exception e ) {
         errorCount = try_failure(errorCount,"EigenvalueDecomposition (nonsymmetric)...","incorrect nonsymmetric Eigenvalue decomposition calculation");
      }

      cout << "\nTestMatrix completed.\n";
      cout << "Total errors reported: " << errorCount << "\n";
      cout << "Total warnings reported: " << warningCount << "\n";
   }

