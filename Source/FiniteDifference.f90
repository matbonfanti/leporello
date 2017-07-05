!******************************************************************************
!*                      MODULE FiniteDifference
!******************************************************************************
!
!>  \brief     Compute derivatives with finite difference
!>  \details   Subroutines which compute derivatives by means of \n
!>             finite difference formulas: \n
!>             1) GetGradient compute gradient vector from potential
!>             2) GetHessian  compute hessian matrix from potential
!>             3) GetHessianFromForces computes Hessian matrix as
!>                  first derivatives of the analytic forces
!>             4) TestForces  test implemented analytic forces by
!>                  comparing them with finite difference ones
!
!******************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             20 February 2017
!>
!******************************************************************************
!
!>  \par Updates
!>  \arg DATE : brief description of the change
!>
!>  \todo   introduce a initialization mechanism to define the type of
!>               finite difference formula to use
!>  \todo   include some automatic testing procedures which might be
!>               used to identify a proper interval for the displacement
!>               e.g. compute difference (choosing some appropriate matrix norm)
!>               of the hessian with respect to a reference hessian for
!>               difference values of smalldelta
!>
!******************************************************************************

MODULE FiniteDifference
#include "preprocessoptions.cpp"
   USE RandomNumberGenerator
   IMPLICIT NONE

      PRIVATE
      PUBLIC :: GetGradient, GetHessian, GetHessianFromForces
      PUBLIC :: TestAnalyticForces

      ! Finite difference - 4 points formula for first derivative
      !> Displacements in units of delta for 4pts first derivative formula
!       REAL, DIMENSION(4), PARAMETER :: DeltasI = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
      !> Coefficients for 4pts first derivative formula
!       REAL, DIMENSION(4), PARAMETER :: CoeffsI = (/ +1./12., -8./12., +8./12., -1./12. /)

      ! Finite difference - 6 points formula for first derivative
      !> Displacements in units of delta for 6pts first derivative formula
      REAL, DIMENSION(6), PARAMETER :: DeltasI = (/ -3.0,    -2.0,    -1.0,    +1.0,     +2.0,    +3.0  /)
      !> Coefficients for 6pts first derivative formula
      REAL, DIMENSION(6), PARAMETER :: CoeffsI = (/ -1./60., 3./20.,  -3./4.,  3./4., -3./20.,   1./60. /)

      ! Finite difference - 9 points formula for second derivative
      !> Displacements in units of delta for 9pts second derivative formula
      REAL, DIMENSION(9), PARAMETER :: DeltasII = &
        (/ -4.0    , -3.0   , -2.0  , -1.0 ,  0.0     , +1.0 , +2.0  , +3.0   , +4.0     /)
      !> Coefficients for 9pts second derivative formula
      REAL, DIMENSION(9), PARAMETER :: CoeffsII = &
        (/ -1./560., 8./315., -1./5., 8./5., -205./72., 8./5., -1./5., 8./315., -1./560. /)

      ! Default value for delta for numerical first derivatives
      REAL, PARAMETER :: SmallDeltaI = 0.01
      ! Default value for delta for numerical second derivatives
      REAL, PARAMETER :: SmallDeltaII = 0.005

      !> Number of random generated points for the derivative testing
      INTEGER, PARAMETER :: NPointsTest = 10.**3

#if defined(LOG_FILE)
    CHARACTER(17), SAVE :: LogStr = " FiniteDifference |"
#endif

!       ! use preprocession option to write debug information to log file...
! #if defined(LOG_FILE)
!       __OPEN_LOG_FILE
!       WRITE(__LOG_UNIT,*)  LogStr," Pippo Pluto Paperino"
!       __CLOSE_LOG_FILE
! #endif
!
!       ! ...or to print them as output
! #if defined(VERBOSE_OUTPUT)
!       WRITE(*,*) "Pippo pluto paperino"
! #endif

   CONTAINS

!==============================================================================
!                                 SUBROUTINES
!==============================================================================


!******************************************************************************
!> Compute the gradient of a potential given as external subroutine
!> When DeltaInp is not given, default value of the module is taken.
!>
!> @param AtPoint       Input vector with the coords where to compute gradient
!> @param GetPotential  Function to evaluate the potential
!> @param DeltaInp      Optional, magnitude of the finite coords displacements
!> @returns             Array with the first derivatives of the pot in AtPoint
!******************************************************************************
FUNCTION GetGradient( AtPoint, GetPotential, DeltaInp ) RESULT( Grad )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)   :: AtPoint
      REAL, DIMENSION(size(AtPoint))   :: Grad
      REAL, OPTIONAL                   :: DeltaInp

      INTERFACE
         REAL FUNCTION GetPotential( X )
            REAL, DIMENSION(:), INTENT(IN)  :: X
         END FUNCTION GetPotential
      END INTERFACE

      REAL, DIMENSION(size(AtPoint)) :: Coordinates
      REAL :: Potential, SmallDelta
      INTEGER :: i,k

      ! Define the displacement length
      IF (PRESENT(DeltaInp)) THEN
         SmallDelta = DeltaInp
      ELSE
         SmallDelta = SmallDeltaI
      ENDIF

      ! Initialize gradient to 0
      Grad(:) = 0.0

      DO i = 1, size(AtPoint)       ! Cycle over the number of coordinates
         DO k = 1, size(DeltasI)         !  Cycle over the finite displacements

               ! Define small displacement from the point where compute the derivative
               Coordinates(:) = AtPoint(:)
               Coordinates(i) = Coordinates(i) + DeltasI(k)*SmallDelta

               ! Compute potential and forces in the displaced coordinate
               Potential = GetPotential( Coordinates )
               ! Increment numerical derivative of the analytical derivative
               Grad(i) = Grad(i) + CoeffsI(k)*Potential

         END DO
      END DO

      ! Divide by the finite displament length
      Grad(:) = Grad(:)/SmallDelta

   END FUNCTION GetGradient

!==============================================================================


!******************************************************************************
!> Compute the Hessian of a potential given as external subroutine
!> When DeltaInp is not given, default value of the module is taken.
!>
!> @param AtPoint       Input vector with the coords where to compute Hessian
!> @param GetPotential  Function to evaluate the potential
!> @param DeltaInp      Optional, magnitude of the finite coords displacements
!> @returns             Matrix of 2nd derivatives of the potential in AtPoint
!******************************************************************************
   FUNCTION GetHessian( AtPoint, GetPotential, DeltaInp ) RESULT( Hessian )
      REAL, DIMENSION(:), INTENT(IN)                 :: AtPoint
      REAL, DIMENSION(size(AtPoint),size(AtPoint))   :: Hessian
      REAL, OPTIONAL                                 :: DeltaInp

      INTERFACE
         REAL FUNCTION GetPotential( X )
            REAL, DIMENSION(:), INTENT(IN)  :: X
         END FUNCTION GetPotential
      END INTERFACE

      REAL, DIMENSION(size(AtPoint)) :: Coordinates
      REAL :: Potential, SmallDelta
      INTEGER :: i,j, m,n

      ! Define the displacement length
      IF (PRESENT(DeltaInp)) THEN
         SmallDelta = DeltaInp
      ELSE
         SmallDelta = SmallDeltaII
      ENDIF

      ! Initialize hessian to 0
      Hessian(:,:) = 0.0

      ! Diagonal elements
      DO i = 1, size(AtPoint)
         DO n = 1, size(DeltasII)
            ! Define small displacement from the point where compute the derivative
            Coordinates(:) = AtPoint(:)
            Coordinates(i) = Coordinates(i) + DeltasII(n)*SmallDelta
            ! Compute potential and forces in the displaced coordinate
            Potential = GetPotential( Coordinates )
            ! Increment numerical derivative
            Hessian(i,i) = Hessian(i,i) + CoeffsII(n)*Potential
         END DO
      END DO

      DO i = 2, size(AtPoint)
         DO j = 1, i-1
            ! Off-diagonal elements, lower triangle
            DO m = 1, size(DeltasI)
               DO n = 1, size(DeltasI)
                  ! Define small displacement from the point where compute the derivative
                  Coordinates(:) = AtPoint(:)
                  Coordinates(i) = Coordinates(i) + DeltasI(m)*SmallDelta
                  Coordinates(j) = Coordinates(j) + DeltasI(n)*SmallDelta
                  ! Compute potential and forces in the displaced coordinate
                  Potential = GetPotential( Coordinates )
                  ! Increment numerical derivative
                  Hessian(i,j) = Hessian(i,j) + CoeffsI(m)*CoeffsI(n)*Potential
               END DO
            END DO
            ! Off-diagonal elements, upper triangle
            Hessian(j,i) = Hessian(i,j)
         END DO
      END DO

      ! Divide by the squared finite displament length
      Hessian(:,:) = Hessian(:,:) / SmallDelta**2

   END FUNCTION GetHessian

!==============================================================================


!******************************************************************************
!> Compute the Hessian of the potential from its Forces (Forces = -Gradient)
!>
!> @param AtPoint       Input vector with the coords where to compute Hessian
!> @param GetPotential  Function to evaluate the potential and forces
!> @param DeltaInp      Optional, magnitude of the finite coords displacements
!> @returns             Matrix of 2nd derivatives of the potential in AtPoint
!**************************************************************************************
   FUNCTION GetHessianFromForces( AtPoint, GetPotAndForces, DeltaInp ) RESULT( Hessian )
      REAL, DIMENSION(:), INTENT(IN)                 :: AtPoint
      REAL, DIMENSION(size(AtPoint),size(AtPoint))   :: Hessian
      REAL, OPTIONAL                                 :: DeltaInp

      INTERFACE
         REAL FUNCTION GetPotAndForces( X, Force )
            REAL, DIMENSION(:), INTENT(IN)  :: X
            REAL, DIMENSION(:), INTENT(OUT) :: Force
         END FUNCTION GetPotAndForces
      END INTERFACE

      REAL, DIMENSION(size(AtPoint)) :: Coordinates, FirstDerivative
      REAL :: Potential, SmallDelta
      INTEGER :: i,k

      ! Define the displacement length
      IF (PRESENT(DeltaInp)) THEN
         SmallDelta = DeltaInp
      ELSE
         SmallDelta = SmallDeltaI
      ENDIF

      ! Initialize hessian to 0
      Hessian(:,:) = 0.0

      DO i = 1, size(AtPoint)       ! Cycle over the number of coordinates
         DO k = 1, size(DeltasI)         !  Cycle over the finite displacements

            ! Define small displacement from the point where compute the derivative
            Coordinates(:) = AtPoint(:)
            Coordinates(i) = Coordinates(i) + DeltasI(k)*SmallDelta

            ! Compute potential and forces in the displaced coordinate
            Potential = GetPotAndForces( Coordinates, FirstDerivative )
            FirstDerivative = - FirstDerivative

            ! Increment numerical derivative of the analytical derivative
            Hessian(i,:) = Hessian(i,:) + CoeffsI(k)*FirstDerivative(:)

         END DO
      END DO

      ! Divide by the finite displament length
      Hessian(:,:) = Hessian(:,:)/SmallDelta
!       CALL TheOneWithMatrixPrintedLineAfterLine( SecDeriv )

   END FUNCTION GetHessianFromForces

!==============================================================================


!******************************************************************************
!> Subroutine to test the computation of the forces by comparing
!> analytically and numerically (finite difference) derivatives
!> Points are sampled randomly (uniform distribution) in the interval
!> between the values defined by CoordMin and CoordMax
!>
!> @param GetPotential  Function to evaluate the potential and forces
!> @param CoordMin      Array with the NDim min values of the coordinates
!> @param CoordMax      Array with the NDim max values of the coordinates
!*******************************************************************************
   SUBROUTINE TestAnalyticForces( GetPotAndForces, CoordMin, CoordMax, DeltaInp )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: CoordMin, CoordMax
      REAL, OPTIONAL, INTENT(IN)     :: DeltaInp

      INTERFACE
         REAL FUNCTION GetPotAndForces( X, Force )
            REAL, DIMENSION(:), INTENT(IN)  :: X
            REAL, DIMENSION(:), INTENT(OUT) :: Force
         END FUNCTION GetPotAndForces
      END INTERFACE

      TYPE(RNGInternalState) :: Random

      REAL, DIMENSION(size(CoordMin))  :: AtPoint, Coordinates, AnalyticalDerivs, NumericalDerivs
      REAL, DIMENSION(size(CoordMin))  :: Average, Deviation, Dummy
      REAL    :: V, SmallDelta
      INTEGER :: iPnt, iCoord, iDispl, NSelected, NDim

      ! Define the displacement length
      IF (PRESENT(DeltaInp)) THEN
         SmallDelta = DeltaInp
      ELSE
         SmallDelta = SmallDeltaI
      ENDIF

      ! Check the number of degree of freedom
      NDim = size(CoordMin)
      CALL ERROR( size(CoordMax) /= NDim, "PotentialModule.TestForces: array dimension mismatch" )

      ! Initialize random number generator
      CALL SetSeed( Random, -512 )

      Average = 0.0; Deviation = 0.0
      NSelected = 0
      DO iPnt = 1, NPointsTest

         ! generate random numbers for the coordinates
         DO iCoord = 1, NDim
            AtPoint(iCoord) = CoordMin(iCoord) + (CoordMax(iCoord) - CoordMin(iCoord)) * UniformRandomNr(Random)
         END DO
         ! Compute analytical derivatives
         V = GetPotAndForces( AtPoint, AnalyticalDerivs )
         AnalyticalDerivs = - AnalyticalDerivs

         IF ( V > -0.18d0 .AND. V < 0.05d0 ) THEN
            NSelected = NSelected + 1
            ! Compute numerical derivatives
            NumericalDerivs(:) = 0.0
            DO iCoord = 1, NDim
               DO iDispl = 1, size(DeltasI)
                  ! Define small displacement from the point where compute the derivative
                  Coordinates(:) = AtPoint(:)
                  Coordinates(iCoord) = Coordinates(iCoord) + DeltasI(iDispl)*SmallDelta
                  ! Compute potential in the displaced coordinate
                  V = GetPotAndForces( Coordinates, Dummy )
                  ! Increment numerical derivative
                  NumericalDerivs(iCoord) = NumericalDerivs(iCoord) + CoeffsI(iDispl)*V
               END DO
            END DO
            NumericalDerivs(:) = NumericalDerivs(:) / SmallDelta

            WRITE(888,*) " "
            WRITE(888,*) V
            DO iCoord = 1, NDim
               WRITE(888,*) NumericalDerivs(iCoord), AnalyticalDerivs(iCoord)
            END DO

            ! Accumulate to compute average and root mean squared deviations
            Average = Average + ( NumericalDerivs - AnalyticalDerivs )
            Deviation = Deviation + ( NumericalDerivs - AnalyticalDerivs )**2
         END IF
      END DO

      ! Normalize averages
      Average = Average / NSelected
      Deviation = SQRT(Deviation / NSelected)

      ! Print results to screen
      PRINT "(2/,A)",    " ***************************************************"
      PRINT "(A,F10.5)", "           TESTING POTENTIAL DERIVATIVES"
      PRINT "(A,/)" ,    " ***************************************************"

      WRITE(*,300) NPointsTest
      WRITE(*,302) NSelected

      DO iCoord = 1, NDim
         WRITE(*,301) iCoord, CoordMin(iCoord), "au", CoordMax(iCoord), "au",&
                              Average(iCoord), "au", Deviation(iCoord), "au"
      END DO

      300 FORMAT(" * Number of points where dV's are evaluated:   ",I10 )
      302 FORMAT(" * Number of points kept for the averages:      ",I10,2/)
      301 FORMAT(" * Coordinate ",I3,/,&
                  "      min for the sampling:                     "F10.4,1X,A,/,&
                  "      max for the sampling:                     "F10.4,1X,A,/,&
                  "      average deviation:                        "F10.4,1X,A,/,&
                  "      RMS deviation:                            "F10.4,1X,A,2/ )

   END SUBROUTINE TestAnalyticForces


!==============================================================================
!                              END OF MODULE
!==============================================================================
END MODULE FiniteDifference

!    REAL FUNCTION DirectionalSecondDerivative( AtPoint, Direction, NormalMode )
!       IMPLICIT NONE
!       REAL, DIMENSION(NDim), INTENT(IN) :: AtPoint, Direction, NormalMode
!
!       REAL, DIMENSION(4), PARAMETER :: Deltas = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
!       REAL, DIMENSION(4), PARAMETER :: Coeffs = (/ +1./12., -8./12., +8./12., -1./12. /)
!
!       REAL, DIMENSION(NDim) :: Coordinates, FirstDerivative
!       REAL :: Potential, Norm
!       INTEGER :: k
!
!       Norm = 0.0
!       DO k = 1, size(Direction)
!          Norm = Norm + Direction(k)**2
!       END DO
!       CALL ERROR( ABS(SQRT(Norm)-1.0) > 1E-12 , " DirectionalSecondDerivative: Direction is not normalized" )
!       Norm = 0.0
!       DO k = 1, size(NormalMode)
!          Norm = Norm + NormalMode(k)**2
!       END DO
!       CALL ERROR( ABS(SQRT(Norm)-1.0) > 1E-12 , " DirectionalSecondDerivative: NormalMode is not normalized" )
!
!       DirectionalSecondDerivative = 0.0
!
!       DO k = 1, size(Deltas)
!          Coordinates(:) = AtPoint(:) + SmallDelta*Deltas(k)*NormalMode(:)
!          Potential = VHSticking( Coordinates, FirstDerivative )
!          DirectionalSecondDerivative = DirectionalSecondDerivative - &
!                           Coeffs(k)* TheOneWithVectorDotVector(Direction,FirstDerivative) /SmallDelta
!       END DO
!
!    END FUNCTION DirectionalSecondDerivative
!
! !*************************************************************************************************
