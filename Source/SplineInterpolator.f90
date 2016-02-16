!***************************************************************************************
!*                              MODULE SplineInterpolator
!***************************************************************************************
!*
!*        CONTAINS:           SplineType                              : Datatype for interpolation
!*
!*                            SetupSpline( SplineData, X, [Y], [Z], func) : Setup subroutine
!*                            GetSpline( SplineData, x, [y], [z])         : Compute spline 
!*                            DisposeSpline( SplineData )                 : Dispose module data
!*
!***************************************************************************************
!*
!*        AUTHOR(S):          1st: M.Bonfanti, January 2010.
!*
!***************************************************************************************
!*
!*        UPDATES(S):         20 Jan 2010 - 3D interpolation added
!*                            01 Feb 2010 - same datatype for any dimension (1,2,3) &
!*                                          1D interpolation interface
!*
!***************************************************************************************
!*
!*        NOTE(S):            The method implemented for the computation of 2D and 3D spline
!*                            are far from being efficient in case of the computation of the
!*                            interpolation function on a regular grid
!*                            On the other hand, this algorithm is general and can be used also
!*                            for computing the interpolating spline on scattered points
!*
!***************************************************************************************
!*
!*        REFS(S):            1)     "NUMERICAL RECIPES IN FORTRAN 90:
!*                                   The Art of PARALLEL Scientific Computing"
!*                                   Chapter B3. Interpolation and Extrapolation 
!*                                   ISBN 0-521-57439-0 
!*                                   Copyright (C) 1986-1996 by Cambridge University Press 
!*
!***************************************************************************************

MODULE SplineInterpolator

   USE NRUtility;               !  constants and utility modules
   IMPLICIT NONE

   TYPE SplineType
         INTEGER                          :: NDim   ! no of indipendent vars
         REAL, POINTER, DIMENSION(:)      :: x => NULL()
         REAL, POINTER, DIMENSION(:)      :: y => NULL()  ! Indipendent Var GRID
         REAL, POINTER, DIMENSION(:)      :: z => NULL()
         REAL, POINTER, DIMENSION(:,:,:)  :: fun => NULL() ! Values of the function in grid points
         REAL, POINTER, DIMENSION(:,:,:)  :: dd => NULL()  ! 2nd Deriv of the function in grid points
         LOGICAL               :: SplineSetUp = .FALSE.
   END TYPE SplineType

   INTERFACE SetupSpline
      MODULE PROCEDURE SetupSpline1D, SetupSpline2D, SetupSpline3D
   END INTERFACE

   INTERFACE GetSpline
      MODULE PROCEDURE GetSpline1D, GetSpline2D, GetSpline3D
   END INTERFACE

   INTERFACE DerivSpline
      MODULE PROCEDURE DerivSpline1D
   END INTERFACE

   CONTAINS

!----------------------------------------------------------------------------
!   FOLLOWING ARE SETUP SUBS
!----------------------------------------------------------------------------

!* * * * * * * * * * * * * * * * * SETUPSPLINE1D * * * * * * * * * * * * * * * *
!* Sets up the 1D Spline interpolation, setting the data in SplineData
!* The input arguments are: X (grid of the indipendent var) and
!* funct (the function on the grid points)
   SUBROUTINE SetupSpline1D( SplineData, X, fun )
   IMPLICIT NONE
   TYPE(SplineType), INTENT(INOUT) :: SplineData
   REAL, DIMENSION(:), INTENT(IN)  :: X
   REAL, DIMENSION(:), INTENT(IN)  :: fun
   INTEGER                         :: nx, status

      IF (SplineData%SplineSetUp) CALL nrerror("SplineData already setup in SetupSpline1D")
      nx = assert_eq( size(fun,1), size(X), "SetupSpline1D" )
      ALLOCATE(SplineData%x(nx), SplineData%fun(nx,1,1), &
                                 SplineData%dd(nx,1,1), STAT=status)
      SplineData%x          = X
      SplineData%fun(:,1,1) = fun

      CALL spline(SplineData%x,SplineData%fun(:,1,1),1.0e30,1.0e30,SplineData%dd(:,1,1)) 
      SplineData%SplineSetUp = .TRUE.
      SplineData%NDim = 1
   END SUBROUTINE SetupSpline1D

!* * * * * * * * * * * * * * * * * SETUPSPLINE2D * * * * * * * * * * * * * * * *
!* Sets up the 2D Spline interpolation, setting the data in SplineData
!* The input arguments are: X and Y the grid along the two dimensions and
!* funct the function on the grid points
   SUBROUTINE SetupSpline2D( SplineData, X, Y, fun )
   IMPLICIT NONE
   TYPE(SplineType), INTENT(INOUT)    :: SplineData
   REAL, DIMENSION(:), INTENT(IN)     :: X, Y
   REAL, DIMENSION(:,:), INTENT(IN)   :: fun
   INTEGER                            :: nx, ny, status

      IF (SplineData%SplineSetUp) CALL nrerror("SplineData already setup in SetupSpline2D")
      nx = assert_eq( size(fun,1), size(X), "SetupSpline2D" )
      ny = assert_eq( size(fun,2), size(Y), "SetupSpline2D" )

      ALLOCATE(SplineData%x(nx), SplineData%y(ny),  &
                             SplineData%fun(nx,ny,1), SplineData%dd(nx,ny,1), STAT=status)
      SplineData%x          = X
      SplineData%y          = Y
      SplineData%fun(:,:,1) = fun

      CALL splie2(SplineData%x,SplineData%y,SplineData%fun(:,:,1),SplineData%dd(:,:,1))
      SplineData%SplineSetUp = .TRUE.
      SplineData%NDim = 2

   END SUBROUTINE SetupSpline2D

!* * * * * * * * * * * * * * * * * SETUPSPLINE3D * * * * * * * * * * * * * * * *
!* Sets up the 3D Spline interpolation, setting the data in SplineData
!* The input arguments are: X,Y,Z the grid along the 3 dimensions and
!* funct the function on the grid points
   SUBROUTINE SetupSpline3D( SplineData, X, Y, Z, fun )
   IMPLICIT NONE
   TYPE(SplineType), INTENT(INOUT)      :: SplineData
   REAL, DIMENSION(:), INTENT(IN)       :: X, Y, Z
   REAL, DIMENSION(:,:,:), INTENT(IN)   :: fun
   INTEGER                              :: nx, ny, nz, status

      IF (SplineData%SplineSetUp) CALL nrerror("SplineData already setup in SetupSpline3D")
      nx = assert_eq( size(fun,1), size(X), "SetupSpline3D" )
      ny = assert_eq( size(fun,2), size(Y), "SetupSpline3D" )
      nz = assert_eq( size(fun,3), size(Z), "SetupSpline3D" )

      ALLOCATE(SplineData%x(nx), SplineData%y(ny), SplineData%z(nz), &
                             SplineData%fun(nx,ny,nz), SplineData%dd(nx,ny,nz), STAT=status)
      SplineData%x   = X
      SplineData%y   = Y
      SplineData%z   = Z
      SplineData%fun = fun

      CALL splie3(SplineData%x,SplineData%y,SplineData%z,SplineData%fun,SplineData%dd)
      SplineData%SplineSetUp = .TRUE.
      SplineData%NDim = 3

   END SUBROUTINE SetupSpline3D
   
!----------------------------------------------------------------------------
!   FOLLOWING IS THE DISPOSE SUB
!----------------------------------------------------------------------------

!* * * * * * * * * * * * * * * * * DISPOSESPLINE * * * * * * * * * * * * * * * *
!* Dispose the spline data type 
   SUBROUTINE DisposeSpline( SplineData )
   IMPLICIT NONE
   TYPE(SplineType), INTENT(INOUT) :: SplineData
   INTEGER  :: status

      IF (.NOT. SplineData%SplineSetUp) CALL nrerror("SplineData not setup in DisposeSpline")
      DEALLOCATE(SplineData%x, SplineData%fun, SplineData%dd, STAT=status)
      IF (ASSOCIATED(SplineData%y)) DEALLOCATE(SplineData%y)
      IF (ASSOCIATED(SplineData%z)) DEALLOCATE(SplineData%z)
      SplineData%SplineSetUp = .FALSE.
   END SUBROUTINE DisposeSpline

!----------------------------------------------------------------------------
!   AND FINALLY THE COMPUTE SUBS
!----------------------------------------------------------------------------

!* * * * * * * * * * * * * * * * * GETSPLINE1D * * * * * * * * * * * * * * * *
!* Get the value of the Spline interpolating function defined in SplineData
   FUNCTION GetSpline1D( SplineData, x )
   IMPLICIT NONE
   REAL                           :: GetSpline1D
   TYPE(SplineType), INTENT(IN)   :: SplineData
   REAL, INTENT(IN)               :: x

      IF (.NOT. SplineData%SplineSetUp) CALL nrerror("SplineData not SetUp in GetSpline1D")
      IF (SplineData%NDim /= 1)  CALL nrerror("Inconsistend dimensions in calling GetSpline1D")
      GetSpline1D = splint(SplineData%x,SplineData%fun(:,1,1),SplineData%dd(:,1,1), x )
   END FUNCTION GetSpline1D

!* * * * * * * * * * * * * * * * * DERIVSPLINE1D * * * * * * * * * * * * * * * *
!* Get the value of the derivative of the Spline interpolating function defined in SplineData
   FUNCTION DerivSpline1D( SplineData, x )
   IMPLICIT NONE
   REAL                           :: DerivSpline1D
   TYPE(SplineType), INTENT(IN)   :: SplineData
   REAL, INTENT(IN)               :: x

      IF (.NOT. SplineData%SplineSetUp) CALL nrerror("SplineData not SetUp in GetSpline1D")
      IF (SplineData%NDim /= 1)  CALL nrerror("Inconsistend dimensions in calling GetSpline1D")
      DerivSpline1D = splind(SplineData%x,SplineData%fun(:,1,1),SplineData%dd(:,1,1), x )
   END FUNCTION DerivSpline1D

!* * * * * * * * * * * * * * * * * GETSPLINE2D * * * * * * * * * * * * * * * *
!* Get the value of the Spline interpolating function defined in SplineData
   FUNCTION GetSpline2D( SplineData, x, y )
   IMPLICIT NONE
   REAL                           :: GetSpline2D
   TYPE(SplineType), INTENT(IN)   :: SplineData
   REAL, INTENT(IN)               :: x,y

      IF (.NOT. SplineData%SplineSetUp) CALL nrerror("SplineData not SetUp in GetSpline2D")
      IF (SplineData%NDim /= 2)  CALL nrerror("Inconsistend dimensions in calling GetSpline2D")
      GetSpline2D = splin2(SplineData%x,SplineData%y,SplineData%fun(:,:,1),SplineData%dd(:,:,1), x, y )
   END FUNCTION GetSpline2D

!* * * * * * * * * * * * * * * * * GETSPLINE3D * * * * * * * * * * * * * * * *
!* Get the value of the 3D Spline interpolating data defined in SplineData
   FUNCTION GetSpline3D( SplineData, x, y, z )
   IMPLICIT NONE
   REAL                           :: GetSpline3D
   TYPE(SplineType), INTENT(IN)   :: SplineData
   REAL, INTENT(IN)               :: x,y,z

      IF (.NOT. SplineData%SplineSetUp) CALL nrerror("SplineData not SetUp in GetSpline3D")
      IF (SplineData%NDim /= 3)  CALL nrerror("Inconsistend dimensions in calling GetSpline3D")
      GetSpline3D = splin3(SplineData%x,SplineData%y,SplineData%z,SplineData%fun,SplineData%dd, x, y, z)
   END FUNCTION GetSpline3D



!----------------------------------------------------------------------------
!   FOLLOWING SUB ARE ADAPTED FROM NR 
!----------------------------------------------------------------------------

!* * * * * * * * * * * * * * * * * SPLIND * * * * * * * * * * * * * * * * * * * * * 
!* Given the arrays xa and ya, which tabulate a function (with the xai  s in increasing 
!* or decreasing order), and given the array y2a, which is the output from spline above,
!* and given a value of x, this routine returns the derivative of a cubic-spline
!* interpolated value.
   FUNCTION splind(xa,ya,y2a,x) 
   IMPLICIT NONE 

   REAL, DIMENSION(:), INTENT(IN) :: xa,ya,y2a 
   REAL, INTENT(IN)               :: x 
   REAL                           :: splind
   INTEGER                        :: khi,klo,n 
   REAL                           :: a,b,h, aprim, bprim, cprim, dprim 

      n=assert_eq(size(xa),size(ya),size(y2a), "splint" ) 
      klo=max(min(locate(xa,x),n-1),1) 
      khi=klo+1 !klo and khi now bracket the input value of x
      h=xa(khi)-xa(klo)
      if (h == 0.0) call nrerror("bad xa input in splint") !The xa s must be distinct.
      a=(xa(khi)-x)/h !Cubic spline polynomial is now evaluated.
      b=(x-xa(klo))/h 
      aprim = -1.0/h
      bprim = +1.0/h
      cprim = ( aprim*( 3.0*a**2-1.0 )*h**2 ) / 6.0
      dprim = ( bprim*( 3.0*b**2-1.0 )*h**2 ) / 6.0
      splind= aprim*ya(klo) + bprim*ya(khi) + cprim*y2a(klo) + dprim*y2a(khi)

   END FUNCTION splind

!* * * * * * * * * * * * * * * * * SPLIE3 * * * * * * * * * * * * * * * * * * * * * 
!* Given an M × N x O tabulated function ya, and tabulated independent variables (x1a,x2a)
!* this routine constructs one-dimensional natural cubic splines for 1d cut of ya and
!* returns the second derivatives in the M × N x O array y2a. (The arrays x1a and x2a
!* are included in the argument list merely for consistency with routine splin3.)
   SUBROUTINE splie3(x1a,x2a,x3a,ya,y2a) 
   IMPLICIT NONE

   REAL, DIMENSION(:),   INTENT(IN)    :: x1a,x2a,x3a
   REAL, DIMENSION(:,:,:), INTENT(IN)  :: ya
   REAL, DIMENSION(:,:,:), INTENT(OUT) :: y2a
   INTEGER                             :: j,m,n,odum

      m    = assert_eq(size(x1a),size(ya,1),size(y2a,1), "splie2: m" )
      n    = assert_eq(size(x2a),size(ya,2),size(y2a,2), "splie2: n" )
      odum = assert_eq(size(x3a),size(ya,3),size(y2a,3), "splie2: odum" )

      DO j=1,odum
            CALL splie2(x1a,x2a,ya(:,:,j),y2a(:,:,j))
      END DO

   END SUBROUTINE splie3

!* * * * * * * * * * * * * * * * * SPLIN3 * * * * * * * * * * * * * * * * * * * * * 
!* Given x1a, x2a, ya as described in splie2 and y2a as produced by that routine; 
!* and given a desired interpolating point x1,x2; this routine returns an interpolated
!*  function value by bicubic spline interpolation. 
   FUNCTION splin3(x1a,x2a,x3a,ya,y2a,x1,x2,x3)
   IMPLICIT NONE

   REAL, DIMENSION(:), INTENT(IN)      :: x1a,x2a,x3a
   REAL, DIMENSION(:,:,:), INTENT(IN)  :: ya,y2a
   REAL, INTENT(IN)                    :: x1,x2,x3
   REAL                                :: splin3
   INTEGER                             :: j,m,n,odum 
   REAL, DIMENSION(size(x3a))          :: yytmp,y2tmp2 

      m    = assert_eq(size(x1a),size(ya,1),size(y2a,1), "splin3: m" )
      n    = assert_eq(size(x2a),size(ya,2),size(y2a,2), "splin3: n" )
      odum = assert_eq(size(x3a),size(ya,3),size(y2a,3), "splin3: odum" )

      do j=1,odum
         yytmp(j)=splin2(x1a,x2a,ya(:,:,j),y2a(:,:,j),x1,x2)
         ! Perform evaluations of the 2d splines constructed by splie3,
         ! using the 2-dimensional spline evaluator splin2
      end do 

      call spline(x3a,yytmp,1.0e30,1.0e30,y2tmp2) 
      ! Construct the one-dimensional column spline and evaluate it. 
      splin3=splint(x3a,yytmp,y2tmp2,x3)

   END FUNCTION splin3



!----------------------------------------------------------------------------
!   HERE ARE NR SUBROUTINE (TAKEN FROM NR for FORTRAN90)
!----------------------------------------------------------------------------

!* * * * * * * * * * * * * * * * * SPLIE2 * * * * * * * * * * * * * * * * * * * * * 
!* Given an M × N tabulated function ya, and N tabulated independent variables x2a, 
!* this routine constructs one-dimensional natural cubic splines of the rows of ya and
!* returns the second derivatives in the M × N array y2a. (The array x1a is included
!*  in the argument list merely for consistency with routine splin2.)
   SUBROUTINE splie2(x1a,x2a,ya,y2a) 
   IMPLICIT NONE

   REAL, DIMENSION(:),   INTENT(IN)  :: x1a,x2a
   REAL, DIMENSION(:,:), INTENT(IN)  :: ya
   REAL, DIMENSION(:,:), INTENT(OUT) :: y2a
   INTEGER                           :: j,m,ndum

      m=assert_eq(size(x1a),size(ya,1),size(y2a,1), "splie2: m" )
      ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2), "splie2: ndum" )
      do j=1,m
         call spline(x2a,ya(j,:),1.0e30,1.0e30,y2a(j,:))
        !  Values 1 × 1030 signal a natural spline.
      end do

   END SUBROUTINE splie2

!* * * * * * * * * * * * * * * * * SPLIN2 * * * * * * * * * * * * * * * * * * * * * 
!* Given x1a, x2a, ya as described in splie2 and y2a as produced by that routine; 
!* and given a desired interpolating point x1,x2; this routine returns an interpolated
!*  function value by bicubic spline interpolation. 
   FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
   IMPLICIT NONE

   REAL, DIMENSION(:), INTENT(IN)   :: x1a,x2a
   REAL, DIMENSION(:,:), INTENT(IN) :: ya,y2a
   REAL, INTENT(IN)                 :: x1,x2 
   REAL                             :: splin2
   INTEGER                          :: j,m,ndum 
   REAL, DIMENSION(size(x1a))       :: yytmp,y2tmp2 

      m=assert_eq(size(x1a),size(ya,1),size(y2a,1), "splin2: m" ) 
      ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2), "splin2: ndum" ) 

      do j=1,m 
         yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2) 
         ! Performm evaluations of the row splines constructed by splie2,
         ! using the one-dimensional spline evaluator splint
      end do 

      call spline(x1a,yytmp,1.0e30,1.0e30,y2tmp2) 
      ! Construct the one-dimensional column spline and evaluate it. 
      splin2=splint(x1a,yytmp,y2tmp2,x1) 

   END FUNCTION splin2

!* * * * * * * * * * * * * * * * * SPLINE * * * * * * * * * * * * * * * * * * * * * 
!* Given arrays x and y of length N containing a tabulated function, i.e., yi = f(xi),
!* with x1 < x2 < ... < xN, and given values yp1 and ypn for the  rst derivative of the 
!* interpolating function at points 1 and N, respectively, this routine returns an array
!* y2 of length N that contains the second derivatives of the interpolating function at 
!* the tabulated points xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine
!* is signaled to set the corresponding boundary condition for a natural spline, with zero 
!* second derivative on that boundary.
   SUBROUTINE spline(x,y,yp1,ypn,y2)
   IMPLICIT NONE

   REAL, DIMENSION(:), INTENT(IN) :: x,y
   REAL, INTENT(IN) :: yp1,ypn
   REAL, DIMENSION(:), INTENT(OUT) :: y2
   INTEGER :: n
   REAL, DIMENSION(size(x)) :: a,b,c,r

      n=assert_eq(size(x),size(y),size(y2), "spline" )

      c(1:n-1)=x(2:n)-x(1:n-1) !Set up the tridiagonal equations
      r(1:n-1)=6.0*((y(2:n)-y(1:n-1))/c(1:n-1))
      r(2:n-1)=r(2:n-1)-r(1:n-2)
      a(2:n-1)=c(1:n-2)
      b(2:n-1)=2.0*(c(2:n-1)+a(2:n-1))
      b(1)=1.0
      b(n)=1.0

      if (yp1 > 0.99e30) then !The lower boundary condition is set either to be natural
         r(1)=0.0
         c(1)=0.0
      else   !or else to have a specified 1st derivative
         r(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
         c(1)=0.5
      end if

      if (ypn > 0.99e30) then !The upper boundary condition is set either to be natural
         r(n)=0.0
         a(n)=0.0
      else !or else to have a specified 1st derivative. 
         r(n)=(-3.0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn) 
         a(n)=0.5 
      end if 

      call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n)) 

   END SUBROUTINE spline

!* * * * * * * * * * * * * * * * * SPLINT * * * * * * * * * * * * * * * * * * * * * 
!* Given the arrays xa and ya, which tabulate a function (with the xai  s in increasing 
!* or decreasing order), and given the array y2a, which is the output from spline above,
!* and given a value of x, this routine returns a cubic-spline interpolated value. 
!* The arrays xa, ya and y2a are all of the same size.
!* We will find the right place in the table by means of locate s bisection algorithm. 
!* This is optimal if sequential calls to this routine are at random values of x. 
!* If sequential calls are in order, and closely spaced, one would do better to store 
!* previous values of klo and khi and test if they remain appropriate on the next call. 
   FUNCTION splint(xa,ya,y2a,x) 
   IMPLICIT NONE 

   REAL, DIMENSION(:), INTENT(IN) :: xa,ya,y2a 
   REAL, INTENT(IN)               :: x 
   REAL                           :: splint
   INTEGER                        :: khi,klo,n 
   REAL                           :: a,b,h 

      n=assert_eq(size(xa),size(ya),size(y2a), "splint" ) 
      klo=max(min(locate(xa,x),n-1),1) 
      khi=klo+1 !klo and khi now bracket the input value of x
      h=xa(khi)-xa(klo)
      if (h == 0.0) call nrerror("bad xa input in splint") !The xa s must be distinct. 
      a=(xa(khi)-x)/h !Cubic spline polynomial is now evaluated. 
      b=(x-xa(klo))/h 
      splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0 

   END FUNCTION splint

!* * * * * * * * * * * * * * * * * LOCATE * * * * * * * * * * * * * * * * * * * * * 
!* Given an array xx(1:N), and given a value x, returns a value j such that x is 
!* between xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing.
!*  j = 0 or j = N is returned to indicate that x is out of range. 
   FUNCTION locate(xx,x) 
   IMPLICIT NONE 

   REAL, DIMENSION(:), INTENT(IN) :: xx 
   REAL, INTENT(IN)               :: x 
   INTEGER                        :: locate 
   INTEGER                        :: n,jl,jm,ju 
   LOGICAL                        :: ascnd 

      n=size(xx) 
      ascnd = (xx(n) >= xx(1)) ! True if ascending order of table, false otherwise. 
      jl=0                     ! Initialize lower 
      ju=n+1                   ! and upper limits. 

      do 
         if (ju-jl <= 1) exit  ! Repeat until this condition is satis ed. 
         jm=(ju+jl)/2          ! Compute a midpoint, 
         if (ascnd .eqv. (x >= xx(jm))) then 
            jl=jm              ! and replace either the lower limit 
         else 
            ju=jm              ! or the upper limit, as appropriate. 
         end if 
      end do 

      if (x == xx(1)) then !Then set the output, being careful with the endpoints. 
         locate=1 
      else if (x == xx(n)) then 
         locate=n-1 
      else 
         locate=jl 
      end if 

   END FUNCTION locate

!* * * * * * * * * * * * * * * * * TRIDAG * * * * * * * * * * * * * * * * * * * * * 
!* Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1)
!* using a serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides)
!* have size N, while a and c (off-diagonal elements) are size N-1
   SUBROUTINE tridag(a,b,c,r,u) 
   IMPLICIT NONE

   REAL, DIMENSION(:), INTENT(IN)  :: a,b,c,r
   REAL, DIMENSION(:), INTENT(OUT) :: u
   REAL, DIMENSION(size(b))        :: gam   !One vector of workspace, gam is needed.
   INTEGER                         :: n,j
   REAL                            :: bet

      n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/), "tridag_ser" )
      bet=b(1)
      if (bet == 0.0) call nrerror( "tridag_ser: Error at code stage 1" )
      !If this happens then you should rewrite your equations as a set of order N-1, with u2 trivially eliminated.
      u(1)=r(1)/bet

      do j=2,n !Decomposition and forward substitution.
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j-1)*gam(j)
         if (bet == 0.0) &  !Algorithm fails; see below routine in Vol. 1.
                  call nrerror( "tridag_ser: Error at code stage 2" )
         u(j)=(r(j)-a(j-1)*u(j-1))/bet
      end do

      do j=n-1,1,-1  !Backsubstitution.
         u(j)=u(j)-gam(j+1)*u(j+1)
      end do

   END SUBROUTINE tridag

END MODULE SplineInterpolator
