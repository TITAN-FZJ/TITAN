module mod_dnsqe
  use mod_kind, only: sp,dp
  implicit none
  real(dp), private :: TOL1 = 1.e-15_dp
contains
      SUBROUTINE DNSQE (FCN, JAC, IOPT, N, X, FVEC, TOL, NPRINT, INFO, WA, LWA)
!***BEGIN PROLOGUE  DNSQE
!***PURPOSE  An easy-to-use code to find a ZERO of a system of N
!            nonlinear functions in N variables by a modification of
!            the Powell hybrid method.
!***LIBRARY   SLATEC
!***CATEGORY  F2A
!***TYPE      REAL(dp) :: (SNSQE-S, DNSQE-D)
!***KEYWORDS  EASY-TO-USE, NONLINEAR SQUARE SYSTEM,
!             POWELL HYBRID METHOD, ZEROS
!***AUTHOR  Hiebert, K. L. (SNLA)
!***DESCRIPTION
!
! 1. Purpose.
!
!       The purpose of DNSQE is to find a ZERO of a system of N
!       nonlinear functions in N variables by a modification of the
!       Powell hybrid method.  This is done by using the more general
!       nonlinear equation solver DNSQ.  The user must provide a
!       subroutine which calculates the functions.  The user has the
!       option of either to provide a subroutine which calculates the
!       Jacobian or to let the code calculate it by a forward-difference
!       approximation.  This code is the combination of the MINPACK
!       codes (Argonne) HYBRD1 and HYBRJ1.
!
! 2. Subroutine and Type Statements.
!
!       SUBROUTINE DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,
!      *                  WA,LWA)
!       INTEGER      :: IOPT,N,NPRINT,INFO,LWA
!       REAL(dp) :: TOL
!       REAL(dp) :: X(N),FVEC(N),WA(LWA)
!       EXTERNAL FCN,JAC
!
! 3. Parameters.
!
!       Parameters designated as input parameters must be specified on
!       entry to DNSQE and are not changed on exit, while parameters
!       designated as output parameters need not be specified on entry
!       and are set to appropriate values on exit from DNSQE.
!
!       FCN is the name of the user-supplied subroutine which calculates
!         the functions.  FCN must be declared in an external statement
!         in the user calling program, and should be written as follows.
!
!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER      :: N,IFLAG
!         REAL(dp) :: X(N),FVEC(N)
!         ----------
!         Calculate the functions at X and
!         return this vector in FVEC.
!         ----------
!         RETURN
!         END
!
!         The value of IFLAG should not be changed by FCN unless the
!         user wants to terminate execution of DNSQE.  In this case set
!         IFLAG to a negative integer.
!
!       JAC is the name of the user-supplied subroutine which calculates
!         the Jacobian.  If IOPT=1, then JAC must be declared in an
!         external statement in the user calling program, and should be
!         written as follows.
!
!         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
!         INTEGER      :: N,LDFJAC,IFLAG
!         REAL(dp) :: X(N),FVEC(N),FJAC(LDFJAC,N)
!         ----------
!         Calculate the Jacobian at X and return this
!         matrix in FJAC.  FVEC contains the function
!         values at X and should not be altered.
!         ----------
!         RETURN
!         END
!
!         The value of IFLAG should not be changed by JAC unless the
!         user wants to terminate execution of DNSQE. In this case set
!         IFLAG to a negative integer.
!
!         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
!
!       IOPT is an input variable which specifies how the Jacobian will
!         be calculated.  If IOPT=1, then the user must supply the
!         Jacobian through the subroutine JAC.  If IOPT=2, then the
!         code will approximate the Jacobian by forward-differencing.
!
!       N is a positive integer input variable set to the number of
!         functions and variables.
!
!       X is an array of length N.  On input X must contain an initial
!         estimate of the solution vector.  On output X contains the
!         final estimate of the solution vector.
!
!       FVEC is an output array of length N which contains the functions
!         evaluated at the output X.
!
!       TOL is a nonnegative input variable.  Termination occurs when
!         the algorithm estimates that the relative error between X and
!         the solution is at most TOL.  Section 4 contains more details
!         about TOL.
!
!       NPRINT is an integer input variable that enables controlled
!         printing of iterates if it is positive.  In this case, FCN is
!         called with IFLAG = 0 at the beginning of the first iteration
!         and every NPRINT iterations thereafter and immediately prior
!         to return, with X and FVEC available for printing. Appropriate
!         print statements must be added to FCN(see example).  If NPRINT
!         is not positive, no special calls of FCN with IFLAG = 0 are
!         made.
!
!       INFO is an integer output variable.  If the user has terminated
!         execution, INFO is set to the (negative) value of IFLAG.  See
!         description of FCN and JAC. Otherwise, INFO is set as follows.
!
!         INFO = 0  Improper input parameters.
!
!         INFO = 1  Algorithm estimates that the relative error between
!                   X and the solution is at most TOL.
!
!         INFO = 2  Number of calls to FCN has reached or exceeded
!                   100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2.
!
!         INFO = 3  TOL is too small.  No further improvement in the
!                   approximate solution X is possible.
!
!         INFO = 4  Iteration is not making good progress.
!
!         Sections 4 and 5 contain more details about INFO.
!
!       WA is a work array of length LWA.
!
!       LWA is a positive integer input variable not less than
!         (3*N**2+13*N))/2.
!
! 4. Successful Completion.
!
!       The accuracy of DNSQE is controlled by the convergence parameter
!       TOL.  This parameter is used in a test which makes a comparison
!       between the approximation X and a solution XSOL.  DNSQE
!       terminates when the test is satisfied.  If TOL is less than the
!       machine precision (as defined by the  function D1MACH(4)), then
!       DNSQE only attempts to satisfy the test defined by the machine
!       precision.  Further progress is not usually possible.  Unless
!       high precision solutions are required, the recommended value
!       for TOL is the square root of the machine precision.
!
!       The test assumes that the functions are reasonably well behaved,
!       and, if the Jacobian is supplied by the user, that the functions
!       and the Jacobian are coded consistently. If these conditions are
!       not satisfied, then DNSQE may incorrectly indicate convergence.
!       The coding of the Jacobian can be checked by the subroutine
!       DCKDER.  If the Jacobian is coded correctly or IOPT=2, then
!       the validity of the answer can be checked, for example, by
!       rerunning DNSQE with a tighter tolerance.
!
!       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
!         vector Z, then this test attempts to guarantee that
!
!               DENORM(X-XSOL) <= TOL*DENORM(XSOL).
!
!         If this condition is satisfied with TOL = 10**(-K), then the
!         larger components of X have K significant decimal digits and
!         INFO is set to 1.  There is a danger that the smaller
!         components of X may have large relative errors, but the fast
!         rate of convergence of DNSQE usually avoids this possibility.
!
! 5. Unsuccessful Completion.
!
!       Unsuccessful termination of DNSQE can be due to improper input
!       parameters, arithmetic interrupts, an excessive number of
!       function evaluations, errors in the functions, or lack of good
!       progress.
!
!       Improper Input Parameters.  INFO is set to 0 if IOPT < 1, or
!         IOPT > 2, or N <= 0, or TOL < 0.E0, or
!         LWA < (3*N**2+13*N)/2.
!
!       Arithmetic Interrupts.  If these interrupts occur in the FCN
!         subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of X by DNSQE.  In this
!         case, it may be possible to remedy the situation by not
!         evaluating the functions here, but instead setting the
!         components of FVEC to numbers that exceed those in the initial
!         FVEC.
!
!       Excessive Number of Function Evaluations.  If the number of
!         calls to FCN reaches 100*(N+1) for IOPT=1 or 200*(N+1) for
!         IOPT=2, then this indicates that the routine is converging
!         very slowly as measured by the progress of FVEC, and INFO is
!         set to 2.  This situation should be unusual because, as
!         indicated below, lack of good progress is usually diagnosed
!         earlier by DNSQE, causing termination with INFO = 4.
!
!       Errors In the Functions.  When IOPT=2, the choice of step length
!         in the forward-difference approximation to the Jacobian
!         assumes that the relative errors in the functions are of the
!         order of the machine precision.  If this is not the case,
!         DNSQE may fail (usually with INFO = 4).  The user should
!         then either use DNSQ and set the step length or use IOPT=1
!         and supply the Jacobian.
!
!       Lack of Good Progress.  DNSQE searches for a ZERO of the system
!         by minimizing the sum of the squares of the functions.  In so
!         doing, it can become trapped in a region where the minimum
!         does not correspond to a ZERO of the system and, in this
!         situation, the iteration eventually fails to make good
!         progress.  In particular, this will happen if the system does
!         not have a ZERO.  If the system has a ZERO, rerunning DNSQE
!         from a different starting point may be helpful.
!
! 6. Characteristics of The Algorithm.
!
!       DNSQE is a modification of the Powell Hybrid method.  Two of
!       its main characteristics involve the choice of the correction as
!       a convex combination of the Newton and scaled gradient
!       directions, and the updating of the Jacobian by the rank-1
!       method of Broyden.  The choice of the correction guarantees
!       (under reasonable conditions) global convergence for starting
!       points far from the solution and a fast rate of convergence.
!       The Jacobian is calculated at the starting point by either the
!       user-supplied subroutine or a forward-difference approximation,
!       but it is not recalculated until the rank-1 method fails to
!       produce satisfactory progress.
!
!       Timing.  The time required by DNSQE to solve a given problem
!         depends on N, the behavior of the functions, the accuracy
!         requested, and the starting point.  The number of arithmetic
!         operations needed by DNSQE is about 11.5*(N**2) to process
!         each evaluation of the functions (call to FCN) and 1.3*(N**3)
!         to process each evaluation of the Jacobian (call to JAC,
!         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
!         the timing of DNSQE will be strongly influenced by the time
!         spent in FCN and JAC.
!
!       Storage.  DNSQE requires (3*N**2 + 17*N)/2 single precision
!         storage locations, in addition to the storage required by the
!         program.  There are no internally declared storage arrays.
!
! *Long Description:
!
! 7. Example.
!
!       The problem is to determine the values of X(1), X(2), ..., X(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*X(1))*X(1)           -2*X(2)                   = -1
!               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
!                                   -X(8) + (3-2*X(9))*X(9) = -1
!
!       **********
!
!       PROGRAM TEST
! C
! C     DRIVER FOR DNSQE EXAMPLE.
! C
!       INTEGER      :: J,N,IOPT,NPRINT,INFO,LWA,NWRITE
!       REAL(dp) :: TOL,FNORM
!       REAL(dp) :: X(9),FVEC(9),WA(180)
!       REAL(dp) :: DENORM,D1MACH
!       EXTERNAL FCN
!       DATA NWRITE /6/
! C
!       IOPT = 2
!       N = 9
! C
! C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
! C
!       DO 10 J = 1, 9
!          X(J) = -1.E0
!    10    CONTINUE
!
!       LWA = 180
!       NPRINT = 0
! C
! C     SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
! C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
! C     THIS IS THE RECOMMENDED SETTING.
! C
!       TOL = SQRT(D1MACH(4))
! C
!       CALL DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
!       FNORM = DENORM(N,FVEC)
!       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
!       STOP
!  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
!      *        5X,' EXIT PARAMETER',16X,I10 //
!      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
!       END
!       SUBROUTINE FCN(N,X,FVEC,IFLAG)
!       INTEGER      :: N,IFLAG
!       REAL(dp) :: X(N),FVEC(N)
!       INTEGER      :: K
!       REAL(dp) :: ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
!       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
! C
!       DO 10 K = 1, N
!          TEMP = (THREE - TWO*X(K))*X(K)
!          TEMP1 = ZERO
!          IF (K /= 1) TEMP1 = X(K-1)
!          TEMP2 = ZERO
!          IF (K /= N) TEMP2 = X(K+1)
!          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
!    10    CONTINUE
!       RETURN
!       END
!
!       RESULTS OBTAINED WITH DIFFERENT COMPILERS OR MACHINES
!       MAY BE SLIGHTLY DIFFERENT.
!
!       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07
!
!       EXIT PARAMETER                         1
!
!       FINAL APPROXIMATE SOLUTION
!
!       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
!       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
!       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
!
!***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
!                 tions. In Numerical Methods for Nonlinear Algebraic
!                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
!                 1988.
!***ROUTINES CALLED  DNSQ, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNSQE
      INTEGER      :: INDEX, INFO, IOPT, J, LR, LWA, MAXFEV, ML, MODE, MU, N, &
           NFEV, NJEV, NPRINT
      REAL(dp) :: EPSFCN, FACTOR, FVEC(*), ONE, TOL, WA(*), &
           X(*), XTOL, ZERO
      EXTERNAL FCN, JAC
      SAVE FACTOR, ONE, ZERO
      DATA FACTOR,ONE,ZERO /1.0D2,1.0_dp,0.0_dp/
!     BEGIN BLOCK PERMITTING ...EXITS TO 20
!***FIRST EXECUTABLE STATEMENT  DNSQE
         INFO = 0
!
!        CHECK THE INPUT PARAMETERS FOR ERRORS.
!
!     ...EXIT
         IF (IOPT < 1 .OR. IOPT > 2 .OR. N <= 0 &
             .OR. TOL < ZERO .OR. LWA < (3*N**2 + 13*N)/2) &
            GO TO 20
!
!        CALL DNSQ.
!
         MAXFEV = 100*(N + 1)
         IF (IOPT == 2) MAXFEV = 2*MAXFEV
         XTOL = TOL
         ML = N - 1
         MU = N - 1
         EPSFCN = D1MACH(4)
!          EPSFCN = ZERO
         MODE = 2
         DO 10 J = 1, N
            WA(J) = ONE
   10    CONTINUE
         LR = (N*(N + 1))/2
         INDEX = 6*N + LR
         CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,WA(INDEX+1),N,XTOL,MAXFEV,ML, &
                   MU,EPSFCN,WA(1),MODE,FACTOR,NPRINT,INFO,NFEV,NJEV, &
                   WA(6*N+1),LR,WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1), &
                   WA(5*N+1))
         IF (INFO == 5) INFO = 4
   20 CONTINUE
      IF (INFO == 0) CALL XERMSG ('SLATEC', 'DNSQE', 'INVALID INPUT PARAMETER.', 2, 1)

!     LAST CARD OF SUBROUTINE DNSQE.
!
      END SUBROUTINE DNSQE


      SUBROUTINE DNSQ (FCN, JAC, IOPT, N, X, FVEC, FJAC, LDFJAC, XTOL, MAXFEV, ML, MU, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO, NFEV, NJEV, R, LR, QTF, WA1, WA2, WA3, WA4)
!***BEGIN PROLOGUE  DNSQ
!***PURPOSE  Find a ZERO of a system of a N nonlinear functions in N
!            variables by a modification of the Powell hybrid method.
!***LIBRARY   SLATEC
!***CATEGORY  F2A
!***TYPE      REAL(dp) :: (SNSQ-S, DNSQ-D)
!***KEYWORDS  NONLINEAR SQUARE SYSTEM, POWELL HYBRID METHOD, ZEROS
!***AUTHOR  Hiebert, K. L. (SNLA)
!***DESCRIPTION
!
! 1. Purpose.
!
!       The purpose of DNSQ is to find a ZERO of a system of N nonlinear
!       functions in N variables by a modification of the Powell
!       hybrid method.  The user must provide a subroutine which
!       calculates the functions.  The user has the option of either to
!       provide a subroutine which calculates the Jacobian or to let the
!       code calculate it by a forward-difference approximation.
!       This code is the combination of the MINPACK codes (Argonne)
!       HYBRD and HYBRDJ.
!
! 2. Subroutine and Type Statements.
!
!       SUBROUTINE DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,
!      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,
!      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4)
!       INTEGER      :: IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR
!       REAL(dp) :: XTOL,EPSFCN,FACTOR
!       REAL(dp)
!       X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),
!      *     WA1(N),WA2(N),WA3(N),WA4(N)
!       EXTERNAL FCN,JAC
!
! 3. Parameters.
!
!       Parameters designated as input parameters must be specified on
!       entry to DNSQ and are not changed on exit, while parameters
!       designated as output parameters need not be specified on entry
!       and are set to appropriate values on exit from DNSQ.
!
!       FCN is the name of the user-supplied subroutine which calculates
!         the functions.  FCN must be declared in an EXTERNAL statement
!         in the user calling program, and should be written as follows.
!
!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER      :: N,IFLAG
!         REAL(dp) :: X(N),FVEC(N)
!         ----------
!         CALCULATE THE FUNCTIONS AT X AND
!         RETURN THIS VECTOR IN FVEC.
!         ----------
!         RETURN
!         END
!
!         The value of IFLAG should not be changed by FCN unless the
!         user wants to terminate execution of DNSQ.  In this case set
!         IFLAG to a negative integer.
!
!       JAC is the name of the user-supplied subroutine which calculates
!         the Jacobian.  If IOPT=1, then JAC must be declared in an
!         EXTERNAL statement in the user calling program, and should be
!         written as follows.
!
!         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
!         INTEGER      :: N,LDFJAC,IFLAG
!         REAL(dp) :: X(N),FVEC(N),FJAC(LDFJAC,N)
!         ----------
!         Calculate the Jacobian at X and return this
!         matrix in FJAC.  FVEC contains the function
!         values at X and should not be altered.
!         ----------
!         RETURN
!         END
!
!         The value of IFLAG should not be changed by JAC unless the
!         user wants to terminate execution of DNSQ.  In this case set
!         IFLAG to a negative integer.
!
!         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
!
!       IOPT is an input variable which specifies how the Jacobian will
!         be calculated.  If IOPT=1, then the user must supply the
!         Jacobian through the subroutine JAC.  If IOPT=2, then the
!         code will approximate the Jacobian by forward-differencing.
!
!       N is a positive integer input variable set to the number of
!         functions and variables.
!
!       X is an array of length N.  On input X must contain an initial
!         estimate of the solution vector.  On output X contains the
!         final estimate of the solution vector.
!
!       FVEC is an output array of length N which contains the functions
!         evaluated at the output X.
!
!       FJAC is an output N by N array which contains the orthogonal
!         matrix Q produced by the QR factorization of the final
!         approximate Jacobian.
!
!       LDFJAC is a positive integer input variable not less than N
!         which specifies the leading dimension of the array FJAC.
!
!       XTOL is a nonnegative input variable.  Termination occurs when
!         the relative error between two consecutive iterates is at most
!         XTOL.  Therefore, XTOL measures the relative error desired in
!         the approximate solution.  Section 4 contains more details
!         about XTOL.
!
!       MAXFEV is a positive integer input variable.  Termination occurs
!         when the number of calls to FCN is at least MAXFEV by the end
!         of an iteration.
!
!       ML is a nonnegative integer input variable which specifies the
!         number of subdiagonals within the band of the Jacobian matrix.
!         If the Jacobian is not banded or IOPT=1, set ML to at
!         least N - 1.
!
!       MU is a nonnegative integer input variable which specifies the
!         number of superdiagonals within the band of the Jacobian
!         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at
!         least N - 1.
!
!       EPSFCN is an input variable used in determining a suitable step
!         for the forward-difference approximation.  This approximation
!         assumes that the relative errors in the functions are of the
!         order of EPSFCN.  If EPSFCN is less than the machine
!         precision, it is assumed that the relative errors in the
!         functions are of the order of the machine precision.  If
!         IOPT=1, then EPSFCN can be ignored (treat it as a dummy
!         argument).
!
!       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
!         internally set.  If MODE = 2, DIAG must contain positive
!         entries that serve as implicit (multiplicative) scale factors
!         for the variables.
!
!       MODE is an integer input variable.  If MODE = 1, the variables
!         will be scaled internally.  If MODE = 2, the scaling is
!         specified by the input DIAG.  Other values of MODE are
!         equivalent to MODE = 1.
!
!       FACTOR is a positive input variable used in determining the
!         initial step bound.  This bound is set to the product of
!         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else to
!         FACTOR itself.  In most cases FACTOR should lie in the
!         interval (.1,100.).  100. is a generally recommended value.
!
!       NPRINT is an integer input variable that enables controlled
!         printing of iterates if it is positive.  In this case, FCN is
!         called with IFLAG = 0 at the beginning of the first iteration
!         and every NPRINT iterations thereafter and immediately prior
!         to return, with X and FVEC available for printing. appropriate
!         print statements must be added to FCN(see example).  If NPRINT
!         is not positive, no special calls of FCN with IFLAG = 0 are
!         made.
!
!       INFO is an integer output variable.  If the user has terminated
!         execution, INFO is set to the (negative) value of IFLAG.  See
!         description of FCN and JAC. Otherwise, INFO is set as follows.
!
!         INFO = 0  Improper input parameters.
!
!         INFO = 1  Relative error between two consecutive iterates is
!                   at most XTOL.
!
!         INFO = 2  Number of calls to FCN has reached or exceeded
!                   MAXFEV.
!
!         INFO = 3  XTOL is too small.  No further improvement in the
!                   approximate solution X is possible.
!
!         INFO = 4  Iteration is not making good progress, as measured
!                   by the improvement from the last five Jacobian
!                   evaluations.
!
!         INFO = 5  Iteration is not making good progress, as measured
!                   by the improvement from the last ten iterations.
!
!         Sections 4 and 5 contain more details about INFO.
!
!       NFEV is an integer output variable set to the number of calls to
!         FCN.
!
!       NJEV is an integer output variable set to the number of calls to
!         JAC. (If IOPT=2, then NJEV is set to ZERO.)
!
!       R is an output array of length LR which contains the upper
!         triangular matrix produced by the QR factorization of the
!         final approximate Jacobian, stored rowwise.
!
!       LR is a positive integer input variable not less than
!         (N*(N+1))/2.
!
!       QTF is an output array of length N which contains the vector
!         (Q transpose)*FVEC.
!
!       WA1, WA2, WA3, and WA4 are work arrays of length N.
!
!
! 4. Successful completion.
!
!       The accuracy of DNSQ is controlled by the convergence parameter
!       XTOL.  This parameter is used in a test which makes a comparison
!       between the approximation X and a solution XSOL.  DNSQ
!       terminates when the test is satisfied.  If the convergence
!       parameter is less than the machine precision (as defined by the
!       function D1MACH(4)), then DNSQ only attempts to satisfy the test
!       defined by the machine precision.  Further progress is not
!       usually possible.
!
!       The test assumes that the functions are reasonably well behaved,
!       and, if the Jacobian is supplied by the user, that the functions
!       and the Jacobian are coded consistently.  If these conditions
!       are not satisfied, then DNSQ may incorrectly indicate
!       convergence.  The coding of the Jacobian can be checked by the
!       subroutine DCKDER. If the Jacobian is coded correctly or IOPT=2,
!       then the validity of the answer can be checked, for example, by
!       rerunning DNSQ with a tighter tolerance.
!
!       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
!         vector Z and D is the diagonal matrix whose entries are
!         defined by the array DIAG, then this test attempts to
!         guarantee that
!
!               DENORM(D*(X-XSOL)) <= XTOL*DENORM(D*XSOL).
!
!         If this condition is satisfied with XTOL = 10**(-K), then the
!         larger components of D*X have K significant decimal digits and
!         INFO is set to 1.  There is a danger that the smaller
!         components of D*X may have large relative errors, but the fast
!         rate of convergence of DNSQ usually avoids this possibility.
!         Unless high precision solutions are required, the recommended
!         value for XTOL is the square root of the machine precision.
!
!
! 5. Unsuccessful Completion.
!
!       Unsuccessful termination of DNSQ can be due to improper input
!       parameters, arithmetic interrupts, an excessive number of
!       function evaluations, or lack of good progress.
!
!       Improper Input Parameters.  INFO is set to 0 if IOPT .LT .1,
!         or IOPT > 2, or N <= 0, or LDFJAC < N, or
!         XTOL < 0.E0, or MAXFEV <= 0, or ML < 0, or MU < 0,
!         or FACTOR <= 0.E0, or LR < (N*(N+1))/2.
!
!       Arithmetic Interrupts.  If these interrupts occur in the FCN
!         subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of X by DNSQ.  In this
!         case, it may be possible to remedy the situation by rerunning
!         DNSQ with a smaller value of FACTOR.
!
!       Excessive Number of Function Evaluations.  A reasonable value
!         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2.
!         If the number of calls to FCN reaches MAXFEV, then this
!         indicates that the routine is converging very slowly as
!         measured by the progress of FVEC, and INFO is set to 2. This
!         situation should be unusual because, as indicated below, lack
!         of good progress is usually diagnosed earlier by DNSQ,
!         causing termination with info = 4 or INFO = 5.
!
!       Lack of Good Progress.  DNSQ searches for a ZERO of the system
!         by minimizing the sum of the squares of the functions.  In so
!         doing, it can become trapped in a region where the minimum
!         does not correspond to a ZERO of the system and, in this
!         situation, the iteration eventually fails to make good
!         progress.  In particular, this will happen if the system does
!         not have a ZERO.  If the system has a ZERO, rerunning DNSQ
!         from a different starting point may be helpful.
!
!
! 6. Characteristics of The Algorithm.
!
!       DNSQ is a modification of the Powell Hybrid method.  Two of its
!       main characteristics involve the choice of the correction as a
!       convex combination of the Newton and scaled gradient directions,
!       and the updating of the Jacobian by the rank-1 method of
!       Broyden.  The choice of the correction guarantees (under
!       reasonable conditions) global convergence for starting points
!       far from the solution and a fast rate of convergence.  The
!       Jacobian is calculated at the starting point by either the
!       user-supplied subroutine or a forward-difference approximation,
!       but it is not recalculated until the rank-1 method fails to
!       produce satisfactory progress.
!
!       Timing.  The time required by DNSQ to solve a given problem
!         depends on N, the behavior of the functions, the accuracy
!         requested, and the starting point.  The number of arithmetic
!         operations needed by DNSQ is about 11.5*(N**2) to process
!         each evaluation of the functions (call to FCN) and 1.3*(N**3)
!         to process each evaluation of the Jacobian (call to JAC,
!         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
!         the timing of DNSQ will be strongly influenced by the time
!         spent in FCN and JAC.
!
!       Storage.  DNSQ requires (3*N**2 + 17*N)/2 single precision
!         storage locations, in addition to the storage required by the
!         program.  There are no internally declared storage arrays.
!
! *Long Description:
!
! 7. Example.
!
!       The problem is to determine the values of X(1), X(2), ..., X(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*X(1))*X(1)           -2*X(2)                   = -1
!               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
!                                   -X(8) + (3-2*X(9))*X(9) = -1
! C     **********
!
!       PROGRAM TEST
! C
! C     Driver for DNSQ example.
! C
!       INTEGER      :: J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR,
!      *        NWRITE
!       REAL(dp) :: XTOL,EPSFCN,FACTOR,FNORM
!       REAL(dp) :: X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9),
!      *     WA1(9),WA2(9),WA3(9),WA4(9)
!       REAL(dp) :: DENORM,D1MACH
!       EXTERNAL FCN
!       DATA NWRITE /6/
! C
!       IOPT = 2
!       N = 9
! C
! C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
! C
!       DO 10 J = 1, 9
!          X(J) = -1.E0
!    10    CONTINUE
! C
!       LDFJAC = 9
!       LR = 45
! C
! C     SET XTOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
! C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
! C     THIS IS THE RECOMMENDED SETTING.
! C
!       XTOL = SQRT(D1MACH(4))
! C
!       MAXFEV = 2000
!       ML = 1
!       MU = 1
!       EPSFCN = 0.E0
!       MODE = 2
!       DO 20 J = 1, 9
!          DIAG(J) = 1.E0
!    20    CONTINUE
!       FACTOR = 1.E2
!       NPRINT = 0
! C
!       CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU,
!      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
!      *           R,LR,QTF,WA1,WA2,WA3,WA4)
!       FNORM = DENORM(N,FVEC)
!       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N)
!       STOP
!  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
!      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
!      *        5X,' EXIT PARAMETER',16X,I10 //
!      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
!       END
!       SUBROUTINE FCN(N,X,FVEC,IFLAG)
!       INTEGER      :: N,IFLAG
!       REAL(dp) :: X(N),FVEC(N)
!       INTEGER      :: K
!       REAL(dp) :: ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
!       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
! C
!       IF (IFLAG /= 0) GO TO 5
! C
! C     INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE.
! C
!       RETURN
!     5 CONTINUE
!       DO 10 K = 1, N
!          TEMP = (THREE - TWO*X(K))*X(K)
!          TEMP1 = ZERO
!          IF (K /= 1) TEMP1 = X(K-1)
!          TEMP2 = ZERO
!          IF (K /= N) TEMP2 = X(K+1)
!          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
!    10    CONTINUE
!       RETURN
!       END
!
!       Results obtained with different compilers or machines
!       may be slightly different.
!
!       Final L2 norm of the residuals  0.1192636E-07
!
!       Number of function evaluations        14
!
!       Exit parameter                         1
!
!       Final approximate solution
!
!       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
!       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
!       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
!
!***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
!                 tions. In Numerical Methods for Nonlinear Algebraic
!                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
!                 1988.
!***ROUTINES CALLED  D1MACH, D1MPYQ, D1UPDT, DDOGLG, DENORM, DFDJC1,
!                    DQFORM, DQRFAC, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNSQ
!       REAL(dp) :: D1MACH,DENORM
      INTEGER      :: I, IFLAG, INFO, IOPT, ITER, IWA(1), J, JM1, L, LDFJAC, &
           LR, MAXFEV, ML, MODE, MU, N, NCFAIL, NCSUC, NFEV, NJEV, &
           NPRINT, NSLOW1, NSLOW2
      REAL(dp) :: ACTRED, DELTA, DIAG(*), EPSFCN, EPSMCH, FACTOR, &
           FJAC(LDFJAC,*), FNORM, FNORM1, FVEC(*), ONE, P0001, P001, &
           P1, P5, PNORM, PRERED, QTF(*), R(*), RATIO, SUM, TEMP, &
           WA1(*), WA2(*), WA3(*), WA4(*), X(*), XNORM, XTOL, ZERO
      EXTERNAL FCN, JAC
      LOGICAL      :: JEVAL,SING
      SAVE ONE, P1, P5, P001, P0001, ZERO
      DATA ONE,P1,P5,P001,P0001,ZERO &
           /1.0_dp,1.0e-1_dp,5.0e-1_dp,1.0e-3_dp,1.0e-4_dp,0.0_dp/
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 320
!***FIRST EXECUTABLE STATEMENT  DNSQ
!          CALL D1MACHSUB(4,EPSMCH)
         EPSMCH = D1MACH(4)
!
         INFO = 0
         IFLAG = 0
         NFEV = 0
         NJEV = 0
!
!        CHECK THE INPUT PARAMETERS FOR ERRORS.
!
!     ...EXIT
         IF (IOPT < 1 .OR. IOPT > 2 .OR. N <= 0 &
             .OR. XTOL < ZERO .OR. MAXFEV <= 0 .OR. ML < 0 &
             .OR. MU < 0 .OR. FACTOR <= ZERO .OR. LDFJAC < N &
             .OR. LR < (N*(N + 1))/2) GO TO 320
         IF (MODE /= 2) GO TO 20
            DO 10 J = 1, N
!     .........EXIT
               IF (DIAG(J) <= ZERO) GO TO 320
   10       CONTINUE
   20    CONTINUE
!
!        EVALUATE THE FUNCTION AT THE STARTING POINT
!        AND CALCULATE ITS NORM.
!
         IFLAG = 1
         CALL FCN(N,X,FVEC,IFLAG)
         NFEV = 1
!     ...EXIT
         IF (IFLAG < 0) GO TO 320
!          CALL DENORMSUB(N,FVEC,FNORM)
         FNORM = DENORM(N,FVEC)
!
!        INITIALIZE ITERATION COUNTER AND MONITORS.
!
         ITER = 1
         NCSUC = 0
         NCFAIL = 0
         NSLOW1 = 0
         NSLOW2 = 0
!
!        BEGINNING OF THE OUTER LOOP.
!
   30    CONTINUE
!           BEGIN BLOCK PERMITTING ...EXITS TO 90
               JEVAL = .TRUE.
!
!              CALCULATE THE JACOBIAN MATRIX.
!
               IF (IOPT == 2) GO TO 40
!
!                 USER SUPPLIES JACOBIAN
!
                  CALL JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
                  NJEV = NJEV + 1
               GO TO 50
   40          CONTINUE
!
!                 CODE APPROXIMATES THE JACOBIAN
!
                  IFLAG = 2
                  CALL DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU, &
                              EPSFCN,WA1,WA2)
                  NFEV = NFEV + MIN(ML+MU+1,N)
   50          CONTINUE
!
!     .........EXIT
               IF (IFLAG < 0) GO TO 320
!
!              COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
!
               CALL DQRFAC(N,N,FJAC,LDFJAC,.FALSE.,IWA,WA1,WA2,WA3)
!
!              ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
!              TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
!
!           ...EXIT
               IF (ITER /= 1) GO TO 90
               IF (MODE == 2) GO TO 70
                  DO 60 J = 1, N
                     DIAG(J) = WA2(J)
                     IF (ABS(WA2(J)) < TOL1) DIAG(J) = ONE
   60             CONTINUE
   70          CONTINUE
!
!              ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED
!              X AND INITIALIZE THE STEP BOUND DELTA.
!
               DO 80 J = 1, N
                  WA3(J) = DIAG(J)*X(J)
   80          CONTINUE
!                CALL DENORMSUB(N,WA3,XNORM)
               XNORM = DENORM(N,WA3)
               DELTA = FACTOR*XNORM
               IF (ABS(DELTA) < TOL1) DELTA = FACTOR
   90       CONTINUE
!
!           FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
!
            DO 100 I = 1, N
               QTF(I) = FVEC(I)
  100       CONTINUE
            DO 140 J = 1, N
               IF (ABS(FJAC(J,J)) < TOL1) GO TO 130
                  SUM = ZERO
                  DO 110 I = J, N
                     SUM = SUM + FJAC(I,J)*QTF(I)
  110             CONTINUE
                  TEMP = -SUM/FJAC(J,J)
                  DO 120 I = J, N
                     QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  120             CONTINUE
  130          CONTINUE
  140       CONTINUE
!
!           COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
!
            SING = .FALSE.
            DO 170 J = 1, N
               L = J
               JM1 = J - 1
               IF (JM1 < 1) GO TO 160
               DO 150 I = 1, JM1
                  R(L) = FJAC(I,J)
                  L = L + N - I
  150          CONTINUE
  160          CONTINUE
               R(L) = WA1(J)
               IF (ABS(WA1(J)) < TOL1) SING = .TRUE.
  170       CONTINUE
!
!           ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
!
            CALL DQFORM(N,N,FJAC,LDFJAC,WA1)
!
!           RESCALE IF NECESSARY.
!
            IF (MODE == 2) GO TO 190
               DO 180 J = 1, N
                  DIAG(J) = MAX(DIAG(J),WA2(J))
  180          CONTINUE
  190       CONTINUE
!
!           BEGINNING OF THE INNER LOOP.
!
  200       CONTINUE
!
!              IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
!
               IF (NPRINT <= 0) GO TO 210
                  IFLAG = 0
                  IF (MOD(ITER-1,NPRINT) == 0) &
                     CALL FCN(N,X,FVEC,IFLAG)
!     ............EXIT
                  IF (IFLAG < 0) GO TO 320
  210          CONTINUE
!
!              DETERMINE THE DIRECTION P.
!
               CALL DDOGLG(N,R,DIAG,QTF,DELTA,WA1,WA2,WA3)
!
!              STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
!
               DO 220 J = 1, N
                  WA1(J) = -WA1(J)
                  WA2(J) = X(J) + WA1(J)
                  WA3(J) = DIAG(J)*WA1(J)
  220          CONTINUE
!                CALL DENORMSUB(N,WA3,PNORM)
               PNORM = DENORM(N,WA3)
!
!              ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
!
               IF (ITER == 1) DELTA = MIN(DELTA,PNORM)
!
!              EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
!
               IFLAG = 1
               CALL FCN(N,WA2,WA4,IFLAG)
               NFEV = NFEV + 1
!     .........EXIT
               IF (IFLAG < 0) GO TO 320
!                CALL DENORMSUB(N,WA4,FNORM1)
               FNORM1 = DENORM(N,WA4)
!
!              COMPUTE THE SCALED ACTUAL REDUCTION.
!
               ACTRED = -ONE
               IF (FNORM1 < FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
!
!              COMPUTE THE SCALED PREDICTED REDUCTION.
!
               L = 1
               DO 240 I = 1, N
                  SUM = ZERO
                  DO 230 J = I, N
                     SUM = SUM + R(L)*WA1(J)
                     L = L + 1
  230             CONTINUE
                  WA3(I) = QTF(I) + SUM
  240          CONTINUE
!                CALL DENORMSUB(N,WA3,TEMP)
               TEMP = DENORM(N,WA3)
               PRERED = ZERO
               IF (TEMP < FNORM) PRERED = ONE - (TEMP/FNORM)**2
!
!              COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
!              REDUCTION.
!
               RATIO = ZERO
               IF (PRERED > ZERO) RATIO = ACTRED/PRERED
!
!              UPDATE THE STEP BOUND.
!
               IF (RATIO >= P1) GO TO 250
                  NCSUC = 0
                  NCFAIL = NCFAIL + 1
                  DELTA = P5*DELTA
               GO TO 260
  250          CONTINUE
                  NCFAIL = 0
                  NCSUC = NCSUC + 1
                  IF (RATIO >= P5 .OR. NCSUC > 1) &
                     DELTA = MAX(DELTA,PNORM/P5)
                  IF (ABS(RATIO-ONE) <= P1) DELTA = PNORM/P5
  260          CONTINUE
!
!              TEST FOR SUCCESSFUL ITERATION.
!
               IF (RATIO < P0001) GO TO 280
!
!                 SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
!
                  DO 270 J = 1, N
                     X(J) = WA2(J)
                     WA2(J) = DIAG(J)*X(J)
                     FVEC(J) = WA4(J)
  270             CONTINUE
!                   CALL DENORMSUB(N,WA2,XNORM)
                  XNORM = DENORM(N,WA2)
                  FNORM = FNORM1
                  ITER = ITER + 1
  280          CONTINUE
!
!              DETERMINE THE PROGRESS OF THE ITERATION.
!
               NSLOW1 = NSLOW1 + 1
               IF (ACTRED >= P001) NSLOW1 = 0
               IF (JEVAL) NSLOW2 = NSLOW2 + 1
               IF (ACTRED >= P1) NSLOW2 = 0
!
!              TEST FOR CONVERGENCE.
!
               IF (DELTA <= XTOL*XNORM .OR. FNORM <= N*EPSMCH) INFO = 1
!     .........EXIT
               IF (INFO /= 0) GO TO 320
!
!              TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
!
               IF (NFEV >= MAXFEV) INFO = 2
               IF (P1*MAX(P1*DELTA,PNORM) <= EPSMCH*XNORM) INFO = 3
               IF (NSLOW2 == 5) INFO = 4
               IF (NSLOW1 == 10) INFO = 5
!     .........EXIT
               IF (INFO /= 0) GO TO 320
!
!              CRITERION FOR RECALCULATING JACOBIAN
!
!           ...EXIT
               IF (NCFAIL == 2) GO TO 310
!
!              CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
!              AND UPDATE QTF IF NECESSARY.
!
               DO 300 J = 1, N
                  SUM = ZERO
                  DO 290 I = 1, N
                     SUM = SUM + FJAC(I,J)*WA4(I)
  290             CONTINUE
                  WA2(J) = (SUM - WA3(J))/PNORM
                  WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
                  IF (RATIO >= P0001) QTF(J) = SUM
  300          CONTINUE
!
!              COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
!
               CALL D1UPDT(N,N,R,WA1,WA2,WA3,SING)
               CALL D1MPYQ(N,N,FJAC,LDFJAC,WA2,WA3)
               CALL D1MPYQ(1,N,QTF,1,WA2,WA3)
!
!              END OF THE INNER LOOP.
!
               JEVAL = .FALSE.
            GO TO 200
  310       CONTINUE
!
!           END OF THE OUTER LOOP.
!
         GO TO 30
  320 CONTINUE
!
!     TERMINATION, EITHER NORMAL OR USER IMPOSED.
!
      IF (IFLAG < 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT > 0) CALL FCN(N,X,FVEC,IFLAG)
      IF (INFO < 0) CALL XERMSG ('SLATEC', 'DNSQ', 'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
      IF (INFO == 0) CALL XERMSG ('SLATEC', 'DNSQ', 'INVALID INPUT PARAMETER.', 2, 1)
      IF (INFO == 2) CALL XERMSG ('SLATEC', 'DNSQ', 'TOO MANY FUNCTION EVALUATIONS.', 9, 1)
      IF (INFO == 3) CALL XERMSG ('SLATEC', 'DNSQ', 'XTOL TOO SMALL. NO FURTHER IMPROVEMENT POSSIBLE.', 3, 1)
      IF (INFO > 4) CALL XERMSG ('SLATEC', 'DNSQ', 'ITERATION NOT MAKING GOOD PROGRESS.', 1, 1)

!     LAST CARD OF SUBROUTINE DNSQ.
!
      END SUBROUTINE DNSQ


      FUNCTION D1MACH (I)
!
!     D1MACH can be used to obtain machine-dependent parameters for
!     double precision numbers.  The results for various values of I are:
!
!     D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!     D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!     D1MACH(3) = B**(-T), the smallest relative spacing.
!     D1MACH(4) = B**(1-T), the largest relative spacing.
!     D1MACH(5) = LOG10(B)
!
!     Assume double precision numbers are represented in the T-digit,
!     base-B form
!
!          sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!     where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
!
!     The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
!     I1MACH(10) = B, the base.
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!  ***REFERENCES
!
!    P. Fox, A. Hall and N. Schryer, Framework for a portable library,
!       ACM Transactions on Mathematical Software 4, 177-188 (1978).
!
!    David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray
!      Hill, NJ. submitted to ACM Transactions on Mathematical Software,
!      August 1996.
!
!  ***REVISION HISTORY  (YYMMDD)
!     790101  DATE WRITTEN
!     960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)
!   --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER      :: I
      REAL(dp) :: D1MACH, B, X = 1._dp

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1)
          D1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          D1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          STOP 'D1MACH -- input arg out of bounds'
      END SELECT

      END FUNCTION D1MACH

      SUBROUTINE D1MPYQ (M, N, A, LDA, V, W)
!***BEGIN PROLOGUE  D1MPYQ
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      REAL(dp) :: (R1MPYQ-S, D1MPYQ-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an M by N matrix A, this subroutine computes A*Q where
!     Q is the product of 2*(N - 1) transformations
!
!           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!     and GV(I), GW(I) are Givens rotations in the (I,N) plane which
!     eliminate elements in the I-th and N-th planes, respectively.
!     Q itself is not given, rather the information to recover the
!     GV, GW rotations is supplied.
!
!     The SUBROUTINE statement is
!
!       SUBROUTINE D1MPYQ(M,N,A,LDA,V,W)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A.
!
!       N IS a positive integer input variable set to the number
!         of columns of A.
!
!       A is an M by N array. On input A must contain the matrix
!         to be postmultiplied by the orthogonal matrix Q
!         described above. On output A*Q has replaced A.
!
!       LDA is a positive integer input variable not less than M
!         which specifies the leading dimension of the array A.
!
!       V is an input array of length N. V(I) must contain the
!         information necessary to recover the Givens rotation GV(I)
!         described above.
!
!       W is an input array of length N. W(I) must contain the
!         information necessary to recover the Givens rotation GW(I)
!         described above.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  D1MPYQ
      INTEGER      :: I, J, LDA, M, N, NM1, NMJ
      REAL(dp) :: A(LDA,*), COS, ONE, SIN, TEMP, V(*), W(*)
      SAVE ONE
      DATA ONE /1.0_dp/
!
!     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.
!
!***FIRST EXECUTABLE STATEMENT  D1MPYQ
      NM1 = N - 1
      IF (NM1 < 1) GO TO 50
      DO 20 NMJ = 1, NM1
         J = N - NMJ
         IF (ABS(V(J)) > ONE) COS = ONE/V(J)
         IF (ABS(V(J)) > ONE) SIN = SQRT(ONE-COS**2)
         IF (ABS(V(J)) <= ONE) SIN = V(J)
         IF (ABS(V(J)) <= ONE) COS = SQRT(ONE-SIN**2)
         DO 10 I = 1, M
            TEMP = COS*A(I,J) - SIN*A(I,N)
            A(I,N) = SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
!
!     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
!
      DO 40 J = 1, NM1
         IF (ABS(W(J)) > ONE) COS = ONE/W(J)
         IF (ABS(W(J)) > ONE) SIN = SQRT(ONE-COS**2)
         IF (ABS(W(J)) <= ONE) SIN = W(J)
         IF (ABS(W(J)) <= ONE) COS = SQRT(ONE-SIN**2)
         DO 30 I = 1, M
            TEMP = COS*A(I,J) + SIN*A(I,N)
            A(I,N) = -SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE

!     LAST CARD OF SUBROUTINE D1MPYQ.
!
      END SUBROUTINE D1MPYQ

      SUBROUTINE D1UPDT (M, N, S, U, V, W, SING)
!***BEGIN PROLOGUE  D1UPDT
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      REAL(dp) :: (R1UPDT-S, D1UPDT-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an M by N lower trapezoidal matrix S, an M-vector U,
!     and an N-vector V, the problem is to determine an
!     orthogonal matrix Q such that
!
!                   t
!           (S + U*V )*Q
!
!     is again lower trapezoidal.
!
!     This subroutine determines Q as the product of 2*(N - 1)
!     transformations
!
!           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!     where GV(I), GW(I) are Givens rotations in the (I,N) plane
!     which eliminate elements in the I-th and N-th planes,
!     respectively. Q itself is not accumulated, rather the
!     information to recover the GV, GW rotations is returned.
!
!     The SUBROUTINE statement is
!
!       SUBROUTINE D1UPDT(M,N,S,LS,U,V,W,SING)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of S.
!
!       N is a positive integer input variable set to the number
!         of columns of S. N must not exceed M.
!
!       S is an array of length LS. On input S must contain the lower
!         trapezoidal matrix S stored by columns. On output S contains
!         the lower trapezoidal matrix produced as described above.
!
!       LS is a positive integer input variable not less than
!         (N*(2*M-N+1))/2. (NOT USED)
!
!       U is an input array of length M which must contain the
!         vector U.
!
!       V is an array of length N. On input V must contain the vector
!         V. On output V(I) contains the information necessary to
!         recover the Givens rotation GV(I) described above.
!
!       W is an output array of length M. W(I) contains information
!         necessary to recover the Givens rotation GW(I) described
!         above.
!
!       SING is a LOGICAL      :: output variable. SING is set TRUE if any
!         of the diagonal elements of the output S are ZERO. Otherwise
!         SING is set FALSE.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  D1UPDT
!       REAL(dp) :: D1MACH
      INTEGER      :: I, J, JJ, L, M, N, NM1, NMJ
      REAL(dp) :: COS, COTAN, GIANT, ONE, P25, P5, S(*), &
           SIN, TAN, TAU, TEMP, U(*), V(*), W(*), ZERO
      LOGICAL      :: SING
      SAVE ONE, P5, P25, ZERO
      DATA ONE,P5,P25,ZERO /1.0_dp,5.0e-1_dp,2.5e-1_dp,0.0_dp/
!
!     GIANT IS THE LARGEST MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  D1UPDT
!       CALL D1MACHSUB(2,GIANT)
      GIANT = D1MACH(2)
!
!     INITIALIZE THE DIAGONAL ELEMENT POINTER.
!
      JJ = (N*(2*M - N + 1))/2 - (M - N)
!
!     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.
!
      L = JJ
      DO 10 I = N, M
         W(I) = S(L)
         L = L + 1
   10    CONTINUE
!
!     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
!     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.
!
      NM1 = N - 1
      IF (NM1 < 1) GO TO 70
      DO 60 NMJ = 1, NM1
         J = N - NMJ
         JJ = JJ - (M - J + 1)
         W(J) = ZERO
         IF (ABS(V(J)) < TOL1) GO TO 50
!
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!        J-TH ELEMENT OF V.
!
         IF (ABS(V(N)) >= ABS(V(J))) GO TO 20
            COTAN = V(N)/V(J)
            SIN = P5/SQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (ABS(COS)*GIANT > ONE) TAU = ONE/COS
            GO TO 30
   20    CONTINUE
            TAN = V(J)/V(N)
            COS = P5/SQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
   30    CONTINUE
!
!        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
!        NECESSARY TO RECOVER THE GIVENS ROTATION.
!
         V(N) = SIN*V(J) + COS*V(N)
         V(J) = TAU
!
!        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
!
         L = JJ
         DO 40 I = J, M
            TEMP = COS*S(L) - SIN*W(I)
            W(I) = SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.
!
      DO 80 I = 1, M
         W(I) = W(I) + V(N)*U(I)
   80    CONTINUE
!
!     ELIMINATE THE SPIKE.
!
      SING = .FALSE.
      IF (NM1 < 1) GO TO 140
      DO 130 J = 1, NM1
         IF (ABS(W(J)) < TOL1) GO TO 120
!
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!        J-TH ELEMENT OF THE SPIKE.
!
         IF (ABS(S(JJ)) >= ABS(W(J))) GO TO 90
            COTAN = S(JJ)/W(J)
            SIN = P5/SQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (ABS(COS)*GIANT > ONE) TAU = ONE/COS
            GO TO 100
   90    CONTINUE
            TAN = W(J)/S(JJ)
            COS = P5/SQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
  100    CONTINUE
!
!        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
!
         L = JJ
         DO 110 I = J, M
            TEMP = COS*S(L) + SIN*W(I)
            W(I) = -SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
  110       CONTINUE
!
!        STORE THE INFORMATION NECESSARY TO RECOVER THE
!        GIVENS ROTATION.
!
         W(J) = TAU
  120    CONTINUE
!
!        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
!
         IF (ABS(S(JJ)) < TOL1) SING = .TRUE.
         JJ = JJ + (M - J + 1)
  130    CONTINUE
  140 CONTINUE
!
!     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.
!
      L = JJ
      DO 150 I = N, M
         S(L) = W(I)
         L = L + 1
  150    CONTINUE
      IF (ABS(S(JJ)) < TOL1) SING = .TRUE.

!     LAST CARD OF SUBROUTINE D1UPDT.
!
      END SUBROUTINE D1UPDT

      SUBROUTINE DDOGLG (N, R, DIAG, QTB, DELTA, X, WA1, WA2)
!***BEGIN PROLOGUE  DDOGLG
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      REAL(dp) :: (DOGLEG-S, DDOGLG-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an M by N matrix A, an N by N nonsingular diagonal
!     matrix D, an M-vector B, and a positive number DELTA, the
!     problem is to determine the convex combination X of the
!     Gauss-Newton and scaled gradient directions that minimizes
!     (A*X - B) in the least squares sense, subject to the
!     restriction that the Euclidean norm of D*X be at most DELTA.
!
!     This subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     QR factorization of A. That is, if A = Q*R, where Q has
!     orthogonal columns and R is an upper triangular matrix,
!     then DDOGLG expects the full upper triangle of R and
!     the first N components of (Q transpose)*B.
!
!     The subroutine statement is
!
!       SUBROUTINE DDOGLG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
!
!     where
!
!       N is a positive integer input variable set to the order of R.
!
!       R is an input array of length LR which must contain the upper
!         triangular matrix R stored by rows.
!
!       LR is a positive integer input variable not less than
!         (N*(N+1))/2.
!
!       DIAG is an input array of length N which must contain the
!         diagonal elements of the matrix D.
!
!       QTB is an input array of length N which must contain the first
!         N elements of the vector (Q transpose)*B.
!
!       DELTA is a positive input variable which specifies an upper
!         bound on the Euclidean norm of D*X.
!
!       X is an output array of length N which contains the desired
!         convex combination of the Gauss-Newton direction and the
!         scaled gradient direction.
!
!       WA1 and WA2 are work arrays of length N.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH, DENORM
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DDOGLG
!       REAL(dp) :: D1MACH,DENORM
      INTEGER      :: I, J, JJ, JP1, K, L, N
      REAL(dp) :: ALPHA, BNORM, DELTA, DIAG(*), EPSMCH, GNORM, &
           ONE, QNORM, QTB(*), R(*), SGNORM, SUM, TEMP, WA1(*), &
           WA2(*), X(*), ZERO
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0_dp,0.0_dp/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
!***FIRST EXECUTABLE STATEMENT  DDOGLG
!       CALL D1MACHSUB(4,EPSMCH)
      EPSMCH = D1MACH(4)
!
!     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
!
      JJ = (N*(N + 1))/2 + 1
      DO 50 K = 1, N
         J = N - K + 1
         JP1 = J + 1
         JJ = JJ - K
         L = JJ + 1
         SUM = ZERO
         IF (N < JP1) GO TO 20
         DO 10 I = JP1, N
            SUM = SUM + R(L)*X(I)
            L = L + 1
   10       CONTINUE
   20    CONTINUE
         TEMP = R(JJ)
         IF (ABS(TEMP) > TOL1) GO TO 40
         L = J
         DO 30 I = 1, J
            TEMP = MAX(TEMP,ABS(R(L)))
            L = L + N - I
   30       CONTINUE
         TEMP = EPSMCH*TEMP
         IF (ABS(TEMP) < TOL1) TEMP = EPSMCH
   40    CONTINUE
         X(J) = (QTB(J) - SUM)/TEMP
   50    CONTINUE
!
!     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
!
      DO 60 J = 1, N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*X(J)
   60    CONTINUE
!       CALL DENORMSUB(N,WA2,QNORM)
      QNORM = DENORM(N,WA2)
      IF (QNORM <= DELTA) GO TO 140
!
!     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
!     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
!
      L = 1
      DO 80 J = 1, N
         TEMP = QTB(J)
         DO 70 I = J, N
            WA1(I) = WA1(I) + R(L)*TEMP
            L = L + 1
   70       CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   80    CONTINUE
!
!     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
!     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
!
!       CALL DENORMSUB(N,WA1,GNORM)
      GNORM = DENORM(N,WA1)
      SGNORM = ZERO
      ALPHA = DELTA/QNORM
      IF (ABS(GNORM) < TOL1) GO TO 120
!
!     CALCULATE THE POINT ALONG THE SCALED GRADIENT
!     AT WHICH THE QUADRATIC IS MINIMIZED.
!
      DO 90 J = 1, N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   90    CONTINUE
      L = 1
      DO 110 J = 1, N
         SUM = ZERO
         DO 100 I = J, N
            SUM = SUM + R(L)*WA1(I)
            L = L + 1
  100       CONTINUE
         WA2(J) = SUM
  110    CONTINUE
!       CALL DENORMSUB(N,WA2,TEMP)
      TEMP = DENORM(N,WA2)
      SGNORM = (GNORM/TEMP)/TEMP
!
!     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
!
      ALPHA = ZERO
      IF (SGNORM >= DELTA) GO TO 120
!
!     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
!     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
!     AT WHICH THE QUADRATIC IS MINIMIZED.
!
!       CALL DENORMSUB(N,QTB,BNORM)
      BNORM = DENORM(N,QTB)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2 &
             + SQRT((TEMP-(DELTA/QNORM))**2 &
                     +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP
  120 CONTINUE
!
!     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
!     DIRECTION AND THE SCALED GRADIENT DIRECTION.
!
      TEMP = (ONE - ALPHA)*MIN(SGNORM,DELTA)
      DO 130 J = 1, N
         X(J) = TEMP*WA1(J) + ALPHA*X(J)
  130    CONTINUE
  140 CONTINUE

!     LAST CARD OF SUBROUTINE DDOGLG.
!
      END SUBROUTINE DDOGLG

      FUNCTION DENORM (N, X)
!***BEGIN PROLOGUE  DENORM
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      REAL(dp) :: (ENORM-S, DENORM-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an N-vector X, this function calculates the
!     Euclidean norm of X.
!
!     The Euclidean norm is computed by accumulating the sum of
!     squares in three different sums. The sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. Non-destructive underflows are permitted. Underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     The definitions of small, intermediate and large components
!     depend on two constants, RDWARF and RGIANT. The main
!     restrictions on these constants are that RDWARF**2 not
!     underflow and RGIANT**2 not overflow. The constants
!     given here are suitable for every known computer.
!
!     The function statement is
!
!       REAL FUNCTION DENORM(N,X)
!
!     where
!
!       N is a positive integer input variable.
!
!       X is an input array of length N.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DENORM
      REAL(dp) :: DENORM
      INTEGER      :: I, N
      REAL(dp) :: AGIANT, FLOATN, ONE, RDWARF, RGIANT, S1, S2, S3, &
           X(*), X1MAX, X3MAX, XABS, ZERO
      SAVE ONE, ZERO, RDWARF, RGIANT
      DATA ONE,ZERO,RDWARF,RGIANT /1.0_dp,0.0_dp,3.834e-20_dp,1.304D19/
!***FIRST EXECUTABLE STATEMENT  DENORM
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = ABS(X(I))
         IF (XABS > RDWARF .AND. XABS < AGIANT) GO TO 70
            IF (XABS <= RDWARF) GO TO 30
!
!              SUM FOR LARGE COMPONENTS.
!
               IF (XABS <= X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
!
!              SUM FOR SMALL COMPONENTS.
!
               IF (XABS <= X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (ABS(XABS) > TOL1) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
!
!           SUM FOR INTERMEDIATE COMPONENTS.
!
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
!
!     CALCULATION OF NORM.
!
      IF (ABS(S1) < TOL1) GO TO 100
         DENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (ABS(S2) < TOL1) GO TO 110
            IF (S2 >= X3MAX) &
               DENORM = SQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 < X3MAX) &
               DENORM = SQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            DENORM = X3MAX*SQRT(S3)
  120    CONTINUE
  130 CONTINUE
!
!     LAST CARD OF FUNCTION DENORM.
!
      END FUNCTION DENORM

      SUBROUTINE DFDJC1 (FCN, N, X, FVEC, FJAC, LDFJAC, IFLAG, ML, MU, EPSFCN, WA1, WA2)
!***BEGIN PROLOGUE  DFDJC1
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      REAL(dp) :: (FDJAC1-S, DFDJC1-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine computes a forward-difference approximation
!     to the N by N Jacobian matrix associated with a specified
!     problem of N functions in N variables. If the Jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     The subroutine statement is
!
!       SUBROUTINE DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
!                         WA1,WA2)
!
!     where
!
!       FCN is the name of the user-supplied subroutine which
!         calculates the functions. FCN must be declared
!         in an EXTERNAL statement in the user calling
!         program, and should be written as follows.
!
!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER      :: N,IFLAG
!         REAL(dp) :: X(N),FVEC(N)
!         ----------
!         Calculate the functions at X and
!         return this vector in FVEC.
!         ----------
!         RETURN
!
!         The value of IFLAG should not be changed by FCN unless
!         the user wants to terminate execution of DFDJC1.
!         In this case set IFLAG to a negative integer.
!
!       N is a positive integer input variable set to the number
!         of functions and variables.
!
!       X is an input array of length N.
!
!       FVEC is an input array of length N which must contain the
!         functions evaluated at X.
!
!       FJAC is an output N by N array which contains the
!         approximation to the Jacobian matrix evaluated at X.
!
!       LDFJAC is a positive integer input variable not less than N
!         which specifies the leading dimension of the array FJAC.
!
!       IFLAG is an integer variable which can be used to terminate
!         the execution of DFDJC1. See description of FCN.
!
!       ML is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         Jacobian matrix. If the Jacobian is not banded, set
!         ML to at least N - 1.
!
!       EPSFCN is an input variable used in determining a suitable
!         step length for the forward-difference approximation. This
!         approximation assumes that the relative errors in the
!         functions are of the order of EPSFCN. If EPSFCN is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       MU is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         Jacobian matrix. If the Jacobian is not banded, set
!         MU to at least N - 1.
!
!       WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at
!         least N, then the Jacobian is considered dense, and WA2 is
!         not referenced.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DFDJC1
!       REAL(dp) :: D1MACH
      INTEGER      :: I, IFLAG, J, K, LDFJAC, ML, MSUM, MU, N
      REAL(dp) :: EPS, EPSFCN, EPSMCH, FJAC(LDFJAC,*), &
           FVEC(*), H, TEMP, WA1(*), WA2(*), X(*), ZERO
      EXTERNAL FCN
      SAVE ZERO
      DATA ZERO /0.0_dp/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
!***FIRST EXECUTABLE STATEMENT  DFDJC1
!       CALL D1MACHSUB(4,EPSMCH)
      EPSMCH = D1MACH(4)
!
      EPS = SQRT(MAX(EPSFCN,EPSMCH))
      MSUM = ML + MU + 1
      IF (MSUM < N) GO TO 40
!
!        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
!
         DO 20 J = 1, N
            TEMP = X(J)
            H = EPS*ABS(TEMP)
            IF (ABS(H) < TOL1) H = EPS
            X(J) = TEMP + H
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG < 0) GO TO 30
            X(J) = TEMP
            DO 10 I = 1, N
               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
         GO TO 110
   40 CONTINUE
!
!        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
!
         DO 90 K = 1, MSUM
            DO 60 J = K, N, MSUM
               WA2(J) = X(J)
               H = EPS*ABS(WA2(J))
               IF (ABS(H) < TOL1) H = EPS
               X(J) = WA2(J) + H
   60          CONTINUE
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG < 0) GO TO 100
            DO 80 J = K, N, MSUM
               X(J) = WA2(J)
               H = EPS*ABS(WA2(J))
               IF (ABS(H) < TOL1) H = EPS
               DO 70 I = 1, N
                  FJAC(I,J) = ZERO
                  IF (I >= J - MU .AND. I <= J + ML) &
                     FJAC(I,J) = (WA1(I) - FVEC(I))/H
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE

!     LAST CARD OF SUBROUTINE DFDJC1.
!
      END SUBROUTINE DFDJC1

      SUBROUTINE DQFORM (M, N, Q, LDQ, WA)
!***BEGIN PROLOGUE  DQFORM
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      REAL(dp) :: (QFORM-S, DQFORM-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine proceeds from the computed QR factorization of
!     an M by N matrix A to accumulate the M by M orthogonal matrix
!     Q from its factored form.
!
!     The subroutine statement is
!
!       SUBROUTINE DQFORM(M,N,Q,LDQ,WA)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A and the order of Q.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       Q is an M by M array. On input the full lower trapezoid in
!         the first MIN(M,N) columns of Q contains the factored form.
!         On output Q has been accumulated into a square matrix.
!
!       LDQ is a positive integer input variable not less than M
!         which specifies the leading dimension of the array Q.
!
!       WA is a work array of length M.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQFORM
      INTEGER      :: I, J, JM1, K, L, LDQ, M, MINMN, N, NP1
      REAL(dp) :: ONE, Q(LDQ,*), SUM, TEMP, WA(*), ZERO
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0_dp,0.0_dp/
!
!     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
!
!***FIRST EXECUTABLE STATEMENT  DQFORM
      MINMN = MIN(M,N)
      IF (MINMN < 2) GO TO 30
      DO 20 J = 2, MINMN
         JM1 = J - 1
         DO 10 I = 1, JM1
            Q(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
!
!     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
!
      NP1 = N + 1
      IF (M < NP1) GO TO 60
      DO 50 J = NP1, M
         DO 40 I = 1, M
            Q(I,J) = ZERO
   40       CONTINUE
         Q(J,J) = ONE
   50    CONTINUE
   60 CONTINUE
!
!     ACCUMULATE Q FROM ITS FACTORED FORM.
!
      DO 120 L = 1, MINMN
         K = MINMN - L + 1
         DO 70 I = K, M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   70       CONTINUE
         Q(K,K) = ONE
         IF (ABS(WA(K)) < TOL1) GO TO 110
         DO 100 J = K, M
            SUM = ZERO
            DO 80 I = K, M
               SUM = SUM + Q(I,J)*WA(I)
   80          CONTINUE
            TEMP = SUM/WA(K)
            DO 90 I = K, M
               Q(I,J) = Q(I,J) - TEMP*WA(I)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE

!     LAST CARD OF SUBROUTINE DQFORM.
!
      END SUBROUTINE DQFORM

      SUBROUTINE DQRFAC (M, N, A, LDA, PIVOT, IPVT, SIGMA, ACNORM, WA)
!***BEGIN PROLOGUE  DQRFAC
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNLS1, DNLS1E, DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      REAL(dp) :: (QRFAC-S, DQRFAC-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   **** Double Precision version of QRFAC ****
!
!     This subroutine uses Householder transformations with column
!     pivoting (optional) to compute a QR factorization of the
!     M by N matrix A. That is, DQRFAC determines an orthogonal
!     matrix Q, a permutation matrix P, and an upper trapezoidal
!     matrix R with diagonal elements of nonincreasing magnitude,
!     such that A*P = Q*R. The Householder transformation for
!     column K, K = 1,2,...,MIN(M,N), is of the form
!
!                           T
!           I - (1/U(K))*U*U
!
!     where U has zeros in the first K-1 positions. The form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding LINPACK subroutine.
!
!     The subroutine statement is
!
!       SUBROUTINE DQRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       A is an M by N array. On input A contains the matrix for
!         which the QR factorization is to be computed. On output
!         the strict upper trapezoidal part of A contains the strict
!         upper trapezoidal part of R, and the lower trapezoidal
!         part of A contains a factored form of Q (the non-trivial
!         elements of the U vectors described above).
!
!       LDA is a positive integer input variable not less than M
!         which specifies the leading dimension of the array A.
!
!       PIVOT is a LOGICAL      :: input variable. If pivot is set .TRUE.,
!         then column pivoting is enforced. If pivot is set .FALSE.,
!         then no column pivoting is done.
!
!       IPVT is an integer output array of length LIPVT. IPVT
!         defines the permutation matrix P such that A*P = Q*R.
!         Column J of P is column IPVT(J) of the identity matrix.
!         If pivot is .FALSE., IPVT is not referenced.
!
!       LIPVT is a positive integer input variable. If PIVOT is
!             .FALSE., then LIPVT may be as small as 1. If PIVOT is
!             .TRUE., then LIPVT must be at least N.
!
!       SIGMA is an output array of length N which contains the
!         diagonal elements of R.
!
!       ACNORM is an output array of length N which contains the
!         norms of the corresponding columns of the input matrix A.
!         If this information is not needed, then ACNORM can coincide
!         with SIGMA.
!
!       WA is a work array of length N. If pivot is .FALSE., then WA
!         can coincide with SIGMA.
!
!***SEE ALSO  DNLS1, DNLS1E, DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH, DENORM
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQRFAC
      INTEGER      :: M,N,LDA
      INTEGER      :: IPVT(*)
      LOGICAL      :: PIVOT
      SAVE ONE, P05, ZERO
      REAL(dp) :: A(LDA,*),SIGMA(*),ACNORM(*),WA(*)
      INTEGER      :: I,J,JP1,K,KMAX,MINMN
      REAL(dp) :: AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
!       REAL(dp) :: D1MACH,DENORM
      DATA ONE,P05,ZERO /1.0_dp,5.0e-2_dp,0.0_dp/
!***FIRST EXECUTABLE STATEMENT  DQRFAC
!       CALL D1MACHSUB(4,EPSMCH)
      EPSMCH = D1MACH(4)
!
!     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
!
      DO 10 J = 1, N
!          CALL DENORMSUB(M,A(1,J),ACNORM(J))
         ACNORM(J) = DENORM(M,A(1,J))
         SIGMA(J) = ACNORM(J)
         WA(J) = SIGMA(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
!
!     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
!
      MINMN = MIN(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
!
!        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
!
         KMAX = J
         DO 20 K = J, N
            IF (SIGMA(K) > SIGMA(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX == J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         SIGMA(KMAX) = SIGMA(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
!
!        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
!        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
!
!          CALL DENORMSUB(M-J+1,A(J,J),AJNORM)
         AJNORM = DENORM(M-J+1,A(J,J))
         IF (ABS(AJNORM) < TOL1) GO TO 100
         IF (A(J,J) < ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
!
!        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
!        AND UPDATE THE NORMS.
!
         JP1 = J + 1
         IF (N < JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. ABS(SIGMA(K)) < TOL1) GO TO 80
            TEMP = A(J,K)/SIGMA(K)
            SIGMA(K) = SIGMA(K)*SQRT(MAX(ZERO,ONE-TEMP**2))
            IF (P05*(SIGMA(K)/WA(K))**2 > EPSMCH) GO TO 80
!             CALL DENORMSUB(M-J,A(JP1,K),SIGMA(K))
            SIGMA(K) = DENORM(M-J,A(JP1,K))
            WA(K) = SIGMA(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         SIGMA(J) = -AJNORM
  110    CONTINUE

!     LAST CARD OF SUBROUTINE DQRFAC.
!
      END SUBROUTINE DQRFAC

      SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***PURPOSE  Symbolic dump (should be locally written).
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (FDUMP-A)
!***KEYWORDS  ERROR, XERMSG
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
      END SUBROUTINE FDUMP

      INTEGER FUNCTION I1MACH (I)
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  The results for various values of I are:
!
!   I/O unit numbers (**MAY REQUIRE USER CUSTOMIZATION**):
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit (obsolete, will cause error)
!     I1MACH( 4) = the standard error message unit.
!                  (the error message unit is usually 0 in UNIX systems)
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!                  (obsolete, will cause an error)
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!          sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!     where 0 <= X(I) < A for I=0,...,S-1.
!
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!     where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
!
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!***REFERENCES
!
!  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
!     ACM Transactions on Mathematical Software 4, 177-188 (1978).
!
!  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray
!    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
!    August 1996.
!
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   960411  Modified for Fortran 90 (BE after suggestions by Eric Grosse)
! --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER      :: I
      REAL(sp) :: X_single = 1.0_sp
      REAL(dp) :: X_double = 1._dp

      SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          STOP 'I1MACH: input arg = 3 is obsolete'
        CASE (4)
          I1MACH = 0 ! Error message unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          STOP 'I1MACH: input arg = 6 is obsolete'
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X_single)
        CASE (11)
          I1MACH = DIGITS(X_single)
        CASE (12)
          I1MACH = MINEXPONENT(X_single)
        CASE (13)
          I1MACH = MAXEXPONENT(X_single)
        CASE (14)
          I1MACH = DIGITS(X_double)
        CASE (15)
          I1MACH = MINEXPONENT(X_double)
        CASE (16)
          I1MACH = MAXEXPONENT(X_double)
        CASE DEFAULT
          STOP 'I1MACH: input argument out of bounds'
      END SELECT

      END FUNCTION I1MACH

      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
!***BEGIN PROLOGUE  J4SAVE
!***SUBSIDIARY
!***PURPOSE  Save or recall global variables needed by error
!            handling routines.
!***LIBRARY   SLATEC (XERROR)
!***TYPE      INTEGER      :: (J4SAVE-I)
!***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
!***SEE ALSO  XERMSG
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900205  Minor modifications to prologue.  (WRB)
!   900402  Added TYPE section.  (WRB)
!   910411  Added KEYWORDS section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  J4SAVE
      LOGICAL      :: ISET
      INTEGER      :: IWHICH, IVALUE, J4SAVE
      INTEGER      :: IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      END FUNCTION J4SAVE

      SUBROUTINE XERHLT (MESSG)
!***BEGIN PROLOGUE  XERHLT
!***SUBSIDIARY
!***PURPOSE  Abort program execution and print error message.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERHLT-A)
!***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        ***Note*** machine dependent routine
!        XERHLT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     Description of Parameters
!        MESSG is as in XERMSG.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to delete length of character
!           and changed routine name from XERABT to XERHLT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERHLT
      WRITE(*,"(A)") MESSG
      STOP
      END SUBROUTINE XERHLT

      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
!***BEGIN PROLOGUE  XERMSG
!***PURPOSE  Process error messages for SLATEC and other libraries.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERMSG-A)
!***KEYWORDS  ERROR MESSAGE, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
!   XERMSG processes a diagnostic message in a manner determined by the
!   value of LEVEL and the current value of the library error control
!   flag, KONTRL.  See subroutine XSETF for details.
!
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
!
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
!
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
!
!                   CALL XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
!
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
!
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
!
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to kount out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
!
!    NERR     An integer value that is chosen by the library routine's
!             author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer.
!
!    LEVEL    An integer value in the range 0 to 2 that indicates the
!             level (severity) of the error.  Their meanings are
!
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
!
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
!
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
!
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
!
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!***REVISION HISTORY  (YYMMDD)
!   880101  DATE WRITTEN
!   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!           THERE ARE TWO BASIC CHANGES.
!           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!               OF LOWER CASE.
!   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!           THE PRINCIPAL CHANGES ARE
!           1.  CLARIFY COMMENTS IN THE PROLOGUES
!           2.  RENAME XRPRNT TO XERPRN
!           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!               CHARACTER FOR NEW RECORDS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           CLEAN UP THE CODING.
!   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!           PREFIX.
!   891013  REVISED TO CORRECT COMMENTS.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!           NERR /= 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!           XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERMSG
      INTEGER NERR, LEVEL, LKNTRL, MKNTRL, MAXMES
      INTEGER KDUMMY, I, KOUNT, LTEMP
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*72  TEMP
!***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
!
!       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!          SHOULD BE PRINTED.
!
!       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!
      IF (NERR<-9999999 .OR. NERR>99999999 .OR. NERR==0 .OR. &
         LEVEL<-1 .OR. LEVEL>2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // &
            'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// &
            'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
!
!       RECORD THE MESSAGE.
!
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!       HANDLE PRINT-ONCE WARNING MESSAGES.
!
      IF (LEVEL==-1 .AND. KOUNT>1) RETURN

      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
!
!       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!       ZERO AND THE ERROR IS NOT FATAL.
!
      IF (LEVEL<2 .AND. LKNTRL==0) GO TO 30
      IF (LEVEL==0 .AND. KOUNT>MAXMES) GO TO 30
      IF (LEVEL==1 .AND. KOUNT>MAXMES .AND. MKNTRL==1) GO TO 30
      IF (LEVEL==2 .AND. KOUNT>MAX(1,MAXMES)) GO TO 30
!
!       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
!       IS NOT ZERO.
!
      IF (LKNTRL /= 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!       1.  LEVEL OF THE MESSAGE
!              'INFORMATIVE MESSAGE'
!              'POTENTIALLY RECOVERABLE ERROR'
!              'FATAL ERROR'
!       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!              'PROG CONTINUES'
!              'PROG ABORTED'
!       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!              'TRACEBACK REQUESTED'
!              'TRACEBACK NOT REQUESTED'
!       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!       EXCEED 74 CHARACTERS.
!       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
!
      IF (LKNTRL > 0) THEN
!
!       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
         IF (LEVEL <= 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL == 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
!
!       THEN WHETHER THE PROGRAM WILL CONTINUE.
!
         IF ((MKNTRL==2 .AND. LEVEL>=1) .OR. &
             (MKNTRL==1 .AND. LEVEL==2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
!
!       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
         IF (LKNTRL > 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!       NOW SEND OUT THE MESSAGE.
!
      CALL XERPRN (' *  ', -1, MESSG, 72)
!
!       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!          TRACEBACK.
!
      IF (LKNTRL > 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) /= ' ') GO TO 20
   10    CONTINUE
!
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
!
!       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
      IF (LKNTRL /= 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
!
!       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 IF (LEVEL<=0 .OR. (LEVEL==1 .AND. MKNTRL<=1)) RETURN
!
!       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
      IF (LKNTRL>0 .AND. KOUNT<MAX(1,MAXMES)) THEN
         IF (LEVEL == 1) THEN
            CALL XERPRN (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      END SUBROUTINE XERMSG

      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
!***BEGIN PROLOGUE  XERPRN
!***SUBSIDIARY
!***PURPOSE  Print error messages processed by XERMSG.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERPRN-A)
!***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
! This routine sends one or more lines to each of the (up to five)
! LOGICAL      :: units to which error messages are to be sent.  This routine
! is called several times by XERMSG, sometimes with a single line to
! print and sometimes with a (potentially very long) message that may
! wrap around into multiple lines.
!
! PREFIX  Input argument of type CHARACTER.  This argument contains
!         characters to be put at the beginning of each line before
!         the body of the message.  No more than 16 characters of
!         PREFIX will be used.
!
! NPREF   Input argument of type INTEGER.  This argument is the number
!         of characters to use from PREFIX.  If it is negative, the
!         intrinsic function LEN is used to determine its length.  If
!         it is ZERO, PREFIX is not used.  If it exceeds 16 or if
!         LEN(PREFIX) exceeds 16, only the first 16 characters will be
!         used.  If NPREF is positive and the length of PREFIX is less
!         than NPREF, a copy of PREFIX extended with blanks to length
!         NPREF will be used.
!
! MESSG   Input argument of type CHARACTER.  This is the text of a
!         message to be printed.  If it is a long message, it will be
!         broken into pieces for printing on multiple lines.  Each line
!         will start with the appropriate prefix and be followed by a
!         piece of the message.  NWRAP is the number of characters per
!         piece; that is, after each NWRAP characters, we break and
!         start a new line.  In addition the characters '$$' embedded
!         in MESSG are a sentinel for a new line.  The counting of
!         characters up to NWRAP starts over for each new line.  The
!         value of NWRAP typically used by XERMSG is 72 since many
!         older error messages in the SLATEC Library are laid out to
!         rely on wrap-around every 72 characters.
!
! NWRAP   Input argument of type INTEGER.  This gives the maximum size
!         piece into which to break MESSG for printing on multiple
!         lines.  An embedded '$$' ends a line, and the kount restarts
!         at the following character.  If a line break does not occur
!         on a blank (it would split a word) that word is moved to the
!         next line.  Values of NWRAP less than 16 will be treated as
!         16.  Values of NWRAP greater than 132 will be treated as 132.
!         The actual line length will be NPREF + NWRAP after NPREF has
!         been adjusted to fall between 0 and 16 and NWRAP has been
!         adjusted to fall between 16 and 132.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   880621  DATE WRITTEN
!   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
!           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
!           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
!           SLASH CHARACTER IN FORMAT STATEMENTS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
!           LINES TO BE PRINTED.
!   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
!           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
!   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Added code to break messages between words.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER      :: NPREF, NWRAP, N, I, LPREF, LWRAP, LENMSG, NEXTC, LPIECE
      INTEGER      :: IDELTA
      CHARACTER*148 CBUFF
      INTEGER      :: IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
!***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
!
!       A ZERO VALUE FOR A LOGICAL      :: UNIT NUMBER MEANS TO USE THE STANDARD
!       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
!       ERROR MESSAGE UNIT.
!
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) == 0) IU(I) = N
   10 CONTINUE
!
!       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
!       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
!       THE REST OF THIS ROUTINE.
!
      IF ( NPREF < 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF /= 0) CBUFF(1:LPREF) = PREFIX
!
!       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
!       TIME FROM MESSG TO PRINT ON ONE LINE.
!
      LWRAP = MAX(16, MIN(132, NWRAP))
!
!       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
!
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) /= ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
!
!       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
!
      IF (LENMSG == 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
!
!       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
!       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
!       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
!       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
!
!       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
!       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
!       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
!       OF THE SECOND ARGUMENT.
!
!       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
!       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
!       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
!       POSITION NEXTC.
!
!       LPIECE == 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
!                       REMAINDER OF THE CHARACTER STRING.  LPIECE
!                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
!                       WHICHEVER IS LESS.
!
!       LPIECE == 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
!                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
!                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
!                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
!                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
!                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
!                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
!                       SHOULD BE INCREMENTED BY 2.
!
!       LPIECE > LWRAP+1  REDUCE LPIECE TO LWRAP.
!
!       ELSE            THIS LAST CASE MEANS 2 <= LPIECE <= LWRAP+1
!                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
!                       PROPERLY HANDLES THE END CASE WHERE LPIECE ==
!                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
!                       AT THE END OF A LINE.
!
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE == 0) THEN
!
!       THERE WAS NO NEW LINE SENTINEL FOUND.
!
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE < LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) == ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE == 1) THEN
!
!       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
!       DON'T PRINT A BLANK LINE.
!
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE > LWRAP+1) THEN
!
!       LPIECE SHOULD BE SET DOWN TO LWRAP.
!
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) == ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
!
!       IF WE ARRIVE HERE, IT MEANS 2 <= LPIECE <= LWRAP+1.
!       WE SHOULD DECREMENT LPIECE BY ONE.
!
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
!
!       PRINT
!
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
!
      IF (NEXTC <= LENMSG) GO TO 50
      END SUBROUTINE XERPRN

      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!***BEGIN PROLOGUE  XERSVE
!***SUBSIDIARY
!***PURPOSE  Record that an error has occurred.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (XERSVE-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER      ::  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!
!        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!
! *Arguments:
!
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or ZERO if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
!
! *Description:
!
!   Record that this error occurred and possibly dump and clear the
!   tables.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   800319  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900413  Routine modified to remove reference to KFLAG.  (WRB)
!   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!           sequence, use IF-THEN-ELSE, make number of saved entries
!           easily changeable, changed routine name from XERSAV to
!           XERSVE.  (RWC)
!   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERSVE
      INTEGER LENTAB
      PARAMETER (LENTAB=10)
      INTEGER      :: LUN(5)
      INTEGER      :: KFLAG, NERR, LEVEL, ICOUNT, KOUNT, KOUNTX, NMSG, NUNIT
      INTEGER      :: KUNIT, IUNIT, I, NERTAB, LEVTAB
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
!***FIRST EXECUTABLE STATEMENT  XERSVE
!
      IF (KFLAG<=0) THEN
!
!        Dump the table.
!
         IF (NMSG==0) RETURN
!
!        Print to each unit.
!
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT==0) IUNIT = I1MACH(4)
!
!           Print the table header.
!
            WRITE (IUNIT,9000)
!
!           Print body of table.
!
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
!
!           Print number of other errors.
!
            IF (KOUNTX/=0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
!
!        Clear the error tables.
!
         IF (KFLAG==0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
!
!        PROCESS A MESSAGE...
!        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB==LIBTAB(I) .AND. SUB==SUBTAB(I) .AND. MES==MESTAB(I) .AND. NERR==NERTAB(I) .AND. LEVEL==LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
!
         IF (NMSG<LENTAB) THEN
!
!           Empty slot found for new message.
!
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
!
!           Table is full.
!
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
!
!     Formats.
!
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / ' LIBRARY    SUBROUTINE MESSAGE START             NERR', '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END SUBROUTINE XERSVE

      SUBROUTINE XGETUA (IUNITA, N)
!***BEGIN PROLOGUE  XGETUA
!***PURPOSE  Return unit number(s) to which error messages are being
!            sent.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XGETUA-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.
!
!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of ZERO refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N < 5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XGETUA
      INTEGER IUNITA, N, I, INDEX
      DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I==1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      END SUBROUTINE XGETUA
end module mod_dnsqe
