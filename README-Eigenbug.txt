I found Creag's bug declaration at http://cio.nist.gov/esd/emaildir/lists/jama/msg01525.html , and indeed Jama gets stuck, as well as my C++ port of Jama (ublasJama).

Looking at where it got stuck, I found out that hqr2 was the problem, so I checked carefully the jama code against the EISPACK code http://www.netlib.org/eispack/hqr2.f

the part called "Double QR step involving rows l:n and columns m:n" in Jama (lines 160 to 260 in hqr2.f) seems to have a bug.
The lines
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0d0) go to 260
(260 is the end of the loop, so go to 260 really means "continue")

Are translated as:
               if (k != m) {
...
                  x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                  if (x != 0.0) {
                     p = p / x;
                     q = q / x;
                     r = r / x;
                  }
               } // end of if(k!=m)
               if (x == 0.0) {
                  break;
               }
which means that x is compared to zero even if it was not computed (k == m), and if it's zero, then we break (exit the loop) instead of continue (next iteration).

Here's what looks like a correct version:
            for (int k = m; k <= n-1; k++) {
               boolean notlast = (k != n-1);
               if (k != m) {
                  p = H[k][k-1];
                  q = H[k+1][k-1];
                  r = (notlast ? H[k+2][k-1] : 0.0);
                  x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                  if (x == 0.0) {
                      continue;
                  }
                  p = p / x;
                  q = q / x;
                  r = r / x;
               }
               s = Math.sqrt(p * p + q * q + r * r);
...


After this fix, Jama doesn't get stuck with the matrix mentioned by Creag, and the result is correct.

Frederic Devernay <frederic.devernay@m4x.org>
