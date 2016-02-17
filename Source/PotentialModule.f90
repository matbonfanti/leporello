!***************************************************************************************
!*                      MODULE PotentialModule
!***************************************************************************************
!
!>  \brief      System potential subroutines
!>  \details    This class define the necessary subroutines to define and compute
!>              the system potential plus some other potential related subroutines
!
!***************************************************************************************
!
!>  \author     Matteo Bonfanti
!>  \version    1.0
!>  \date       8 February 2016
!>
!***************************************************************************************
!
!>  \pre        To use the class the potential needs to be setup with the 
!>              subroutine SetupPotential 
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!
!***************************************************************************************

MODULE PotentialModule
#include "preprocessoptions.cpp"
   USE RandomNumberGenerator
   USE UnitConversion

   PRIVATE
   PUBLIC :: SetupPotential                                        !< setup subroutine
   PUBLIC :: GetXLabel                                             !< info subroutines
   PUBLIC :: GetPotential, GetPotAndForces, GetSecondDerivatives   !< get potential and forces, and second derivatives
   PUBLIC :: SteepLocator, NewtonLocator          !< optimization and stationary points

   !> Setup variable for the potential
   LOGICAL, SAVE :: PotentialModuleIsSetup = .FALSE.

   !> Number of dimensions of the potential
   INTEGER, SAVE :: NDim

   !> Labels of the potential coordinates
   CHARACTER(5), DIMENSION(:), ALLOCATABLE, SAVE :: CoordLabels

   !> Number of random generated points for the derivative testing
   INTEGER, PARAMETER :: NPointsTest = 10.**3
   !> Step for computing the numerical derivatives with finite differences
   REAL, PARAMETER :: SmallDelta = 0.00001

   CONTAINS

!===============================================================================================================================

!**************************************************************************************
!> Setup subroutine for the system potential. 
!> Currently nothing needs to be set, but the subroutine is included for 
!> future developments. The input optional variable DerivTest triggers the 
!> testing procedure for the derivatives of the potential. Note that the subroutine 
!> terminates the program after the test!
!>
!> @param DerivTest     Logical optional variable to set the derivative test
!**************************************************************************************
      SUBROUTINE SetupPotential( DerivTest )
         IMPLICIT NONE
         LOGICAL, INTENT(IN), OPTIONAL :: DerivTest
         REAL, DIMENSION(:), ALLOCATABLE :: CoordMin, CoordMax         !< Coordinate intervals where to check the derivatives

         ! exit if module is setup
         IF ( PotentialModuleIsSetup ) RETURN
 
         ! Potential module is set up
         PotentialModuleIsSetup = .TRUE.
 
         ! Set the number of dimensions
         NDim = 3

         ! Set the labels of the coordinates
         ALLOCATE( CoordLabels(NDim) )
         CoordLabels(1) = "Hi_z"
         CoordLabels(2) = "Ht_z"
         CoordLabels(3) = "C_z"

#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " System potential has been setup"
      WRITE(__LOG_UNIT,*) " ... (details) ... "
      __CLOSE_LOG_FILE
#endif

         ! Tests on the computation of the forces
         IF ( PRESENT(DerivTest) ) THEN
            IF ( DerivTest ) THEN

               ! Set the intervals where to evaluate the derivatives
               ALLOCATE( CoordMin(NDim), CoordMax(NDim) )
               CoordMin(1) = 1.0/MyConsts_Bohr2Ang ; CoordMax(1) = 6.0/MyConsts_Bohr2Ang
               CoordMin(2) = 0.5/MyConsts_Bohr2Ang ; CoordMax(2) = 3.0/MyConsts_Bohr2Ang
               CoordMin(3) = -1.0/MyConsts_Bohr2Ang ; CoordMax(3) = 1.0/MyConsts_Bohr2Ang

               ! Call the test subroutine
               CALL TestForces( CoordMin, CoordMax )

               ! Deallocate memory
               DEALLOCATE( CoordMin, CoordMax )
               ! Terminate the program
               CALL AbortWithError( " Terminating the program after derivatives tests " )

            END IF
         END IF

      END SUBROUTINE SetupPotential

!===============================================================================================================================

!**************************************************************************************
!> Function that returns the label of the coordinates of the potential.
!>
!> @param iCoord     Integer ordinal number of the coordinate
!> @returns Label    Label of the iCoord-th coordinate
!**************************************************************************************
      CHARACTER(5) FUNCTION GetXLabel( iCoord ) RESULT( Label )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: iCoord

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetXLabel : Module not Setup" )
         ! Potential module is set up
         PotentialModuleIsSetup = .TRUE.

         CALL ERROR( iCoord > NDim .OR. iCoord < 1, "PotentialModule.GetXLabel : wrong coordinate number" )
         Label = CoordLabels(iCoord)
 
      END FUNCTION GetXLabel

!===============================================================================================================================

!**************************************************************************************
!> Subrotine to compute potential for the 3D system model, i.e.
!> collinear H + H + C, input and output are in internal coordinates
!>
!> @param Positions    Array with NDim coordinates, 1) z of the incident
!>                     H atom, 2) z of the target H atom, 3) z of the carbon 
!> @returns vv         Output potential
!**************************************************************************************
      REAL FUNCTION GetPotential( Positions ) RESULT(V) 
         IMPLICIT NONE
         REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
         REAL, DIMENSION(SIZE(Positions))        :: Dummy

         INTERFACE
            SUBROUTINE ER_3D (zi, zt, zc, vv, dv_zi, dv_zt, dv_zc)
               REAL*8, INTENT (IN) :: zi, zt, zc
               REAL*8, INTENT(OUT) :: vv, dv_zi, dv_zt, dv_zc
            END SUBROUTINE
          END INTERFACE

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetPotential : Module not Setup" )
         ! Check the number of degree of freedom
         CALL ERROR( size(Positions) /= NDim, "PotentialModule.GetPotential: Positions array dimension mismatch" )

         ! Compute energy
         CALL ER_3D( Positions(1), Positions(2), Positions(3), V, Dummy(1), Dummy(2), Dummy(3) )

      END FUNCTION GetPotential

!===============================================================================================================================

!**************************************************************************************
!>  Subrotine to compute potential and forces, for the 3D system model, i.e.
!>  collinear H + H + C, input and output are in internal coordinates
!>
!>  @param Positions    Array with 3 cartesian Z coordinates, 1) z of the incident
!>                     H atom, 2) z of the target H atom, 3) z of the carbon 
!>  @param Forces       Output 3D array with the derivatives of the potential
!>  @returns vv         Output potential
!**************************************************************************************
      REAL FUNCTION GetPotAndForces( Positions, Forces ) RESULT(V) 
         IMPLICIT NONE
         REAL, DIMENSION(:), TARGET, INTENT(IN)                :: Positions     
         REAL, DIMENSION(SIZE(Positions)), TARGET, INTENT(OUT) :: Forces 

         INTERFACE
            SUBROUTINE ER_3D (zi, zt, zc, vv, dv_zi, dv_zt, dv_zc)
               REAL*8, INTENT (IN) :: zi, zt, zc
               REAL*8, INTENT(OUT) :: vv, dv_zi, dv_zt, dv_zc
            END SUBROUTINE
          END INTERFACE

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetPotAndForces : module not set" )
         ! Check the number of degree of freedom
         CALL ERROR( size(Positions) /= NDim, "PotentialModule.GetPotAndForces: input array dimension mismatch" )

         ! Compute energy
         CALL ER_3D( Positions(1), Positions(2), Positions(3), V, Forces(1), Forces(2), Forces(3) )
         Forces(:) = - Forces(:)

      END FUNCTION GetPotAndForces

!===============================================================================================================================

!**************************************************************************************
!> Compute the second derivatives of the potential
!>
!> @param AtPoint   Input vector with the coordinates where to compute H
!> @returns         Matrix of 2nd derivatives of the potential in AtPoint
!**************************************************************************************
      FUNCTION GetSecondDerivatives( AtPoint ) RESULT( SecDeriv )
         REAL, DIMENSION(:), INTENT(IN)                 :: AtPoint
         REAL, DIMENSION(size(AtPoint),size(AtPoint))   :: SecDeriv

         REAL, DIMENSION(4), PARAMETER :: Deltas = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
         REAL, DIMENSION(4), PARAMETER :: Coeffs = (/ +1./12., -8./12., +8./12., -1./12. /) 

         REAL, DIMENSION(size(AtPoint)) :: Coordinates, FirstDerivative
         REAL :: Potential
         INTEGER :: i,k

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetSecondDerivatives : module not set" )
         ! Check the number of degree of freedom
         CALL ERROR( size(AtPoint) /= NDim, "PotentialModule.GetSecondDerivatives: input array dimension mismatch" )

         DO i = 1, NDim
            DO k = 1, size(Deltas)

               ! Define small displacement from the point where compute the derivative
               Coordinates(:) = AtPoint(:)
               Coordinates(i) = Coordinates(i) + Deltas(k)*SmallDelta

               ! Compute potential and forces in the displaced coordinate
               Potential = GetPotAndForces( Coordinates, FirstDerivative )
               FirstDerivative = - FirstDerivative

               ! Increment numerical derivative of the analytical derivative
               SecDeriv(i,:) = SecDeriv(i,:) + Coeffs(k)*FirstDerivative(:)

            END DO
         END DO

         SecDeriv(:,:) = SecDeriv(:,:)/SmallDelta
!        CALL TheOneWithMatrixPrintedLineAfterLine( SecDeriv )

      END FUNCTION GetSecondDerivatives
 
!===============================================================================================================================

!**************************************************************************************
!> Project the Hessian to impose a set of constraints defined by a logical mask.
!>
!> @param Hessian   Input matrix with the Hessian
!> @param Mask      Input vector, with the logical mask of the constraints
!> @returns         Constrained Hessian matrix of dimension COUNT(Maks)
!**************************************************************************************
      FUNCTION ConstrainedHessian( Hessian, Mask )
         IMPLICIT NONE
         REAL, DIMENSION(:,:), INTENT(IN) :: Hessian
         LOGICAL, DIMENSION(SIZE(Hessian,1)), INTENT(IN) :: Mask
         REAL, DIMENSION(COUNT(Mask),COUNT(Mask)) :: ConstrainedHessian
         INTEGER :: i,j, n1,n2
         
         n1 = 0
         DO i = 1, SIZE(Hessian,2)
            IF ( .NOT. Mask(i) ) CYCLE
            n1 = n1+1
            n2 = 0
            DO j = 1, SIZE(Hessian,1)
               IF ( .NOT. Mask(j) ) CYCLE
               n2 = n2+1
               ConstrainedHessian(n2,n1) = Hessian(j,i)
            END DO
         END DO

      END FUNCTION ConstrainedHessian

!===============================================================================================================================

!**************************************************************************************
!> Project a vector to impose a set of constraints defined by a logical mask.
!>
!> @param Vector    Input vector to project
!> @param Mask      Input vector, with the logical mask of the constraints
!> @returns         Constrained vector of dimension COUNT(Maks)
!**************************************************************************************
      FUNCTION ConstrainedVector( Vector, Mask )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(IN) :: Vector
         LOGICAL, DIMENSION(SIZE(Vector)), INTENT(IN) :: Mask
         REAL, DIMENSION(COUNT(Mask)) :: ConstrainedVector
         INTEGER :: i, n

         n = 0
         DO i = 1, size(Vector)
            IF ( Mask(i) ) THEN
               n = n + 1
               ConstrainedVector(n) = Vector(i)
            END IF
         END DO
         
      END FUNCTION ConstrainedVector

!===============================================================================================================================

!**************************************************************************************
!>  Minimize potential with a steepest descent algorithm.
!>
!> @param StartX        Input vector with the starting point of the minimization
!> @param NMaxIter      Max number of iterations 
!> @param GradThresh    Threshold on the gradient 
!> @returns             Vector with the result of the steepest descent search
!**************************************************************************************
      FUNCTION SteepLocator( StartX, NMaxIter, GradThresh ) RESULT( StationaryPoint )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(IN)   :: StartX
         INTEGER, INTENT(IN)              :: NMaxIter
         REAL, INTENT(IN)                 :: GradThresh
         REAL, DIMENSION(size(StartX))    :: StationaryPoint

         REAL, DIMENSION(size(StartX))              :: CurrentX
         REAL, DIMENSION(size(StartX),size(StartX)) :: EigenVectors, Hessian
         REAL, DIMENSION(size(StartX))              :: Forces, EigenValues
         REAL :: V, GradNorm
         INTEGER :: NIter

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.SteepLocator : module not set" )

         ! Check the number of degree of freedom
         CALL ERROR( size(StartX) /= NDim, "PotentialModule.SteepLocator: input array dimension mismatch" )

#if defined(LOG_FILE)
         __OPEN_LOG_FILE
         ! Print info to screen
         WRITE(__LOG_UNIT,"(/,A,/)") " SteepLocator | Stationary point locator with steepest descent method "
         WRITE(__LOG_UNIT,"(A,/)")   " SteepLocator | Looking for a minimum of the potential"
#endif

         ! Start at initial position
         CurrentX = StartX

#if defined(LOG_FILE)
         WRITE(__LOG_UNIT,"(A12,A20)") "SteepLocator | N Iteration",  "Gradient Norm"
         WRITE(__LOG_UNIT,*)           "SteepLocator | ---------------------------------------------"
#endif

         Iterations: DO NIter = 1, NMaxIter

            ! Compute Hessian and forces at current position
            V = GetPotAndForces( CurrentX, Forces )
            ! Compute norm of the gradient
            GradNorm  = SQRT(TheOneWithVectorDotVector(Forces, Forces))
            ! Move to new coordinates
            CurrentX(:) = CurrentX(:) + Forces(:)

#if defined(LOG_FILE)
            ! Print info to screen
            IF ( MOD(NIter-1,NMaxIter/20) == 0 ) WRITE(__LOG_UNIT,"(A,I12,E20.6,E20.6)") "SteepLocator | ", NIter, GradNorm
#endif

            ! Check convergence criteria
            IF ( GradNorm < GradThresh ) THEN
               EXIT Iterations
            END IF

         END DO Iterations

#if defined(LOG_FILE)
         WRITE(__LOG_UNIT,"(A,I12,E20.6,E20.6)") "SteepLocator | ",NIter, GradNorm
#endif

         ! Check max number of iterations
         CALL WARN( NIter == NMaxIter+1, " PotentialModule.SteepLocator: max number of iterations reached " )
         ! Store final point
         StationaryPoint = CurrentX

#if defined(LOG_FILE)
         ! Check the number of imaginary frequencies
         Hessian = GetSecondDerivatives( StationaryPoint )
         CALL TheOneWithDiagonalization( Hessian, EigenVectors, EigenValues )
         WRITE(__LOG_UNIT,"(/,A,I3,A,/)") " SteepLocator | Final stationary point has ", COUNT( EigenValues < 0.0 ),   &
                                 " imaginary frequency/ies "
         __CLOSE_LOG_FILE
#endif

      END FUNCTION SteepLocator

!===============================================================================================================================

!**************************************************************************************
!>  Find stationary point with Newton's algorithm.
!>  see J. Chem. Phys., Vol. 121, No. 20, 22 November 2004
!>
!> @param StartX            Input vector with the starting point of the minimization
!> @param NMaxIter          Max number of iterations 
!> @param GradThresh        Threshold on the gradient 
!> @param DisplThresh       Threshold on the displacement
!> @param Mask              Mask to set constraints on the search
!> @param TransitionState   Logical variable to set TS search instead of minimum
!> @returns                 Vector with the result of the Newton search
!**************************************************************************************
      FUNCTION NewtonLocator( StartX, NMaxIter, GradThresh, DisplThresh, Mask, TransitionState ) RESULT( StationaryPoint )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(IN)    :: StartX
         INTEGER, INTENT(IN)               :: NMaxIter
         REAL, INTENT(IN)                  :: GradThresh, DisplThresh
         LOGICAL, INTENT(IN), OPTIONAL     :: TransitionState
         LOGICAL, DIMENSION(size(StartX)), INTENT(IN), OPTIONAL :: Mask
         REAL, DIMENSION(size(StartX))     :: StationaryPoint

         REAL, DIMENSION(size(StartX)) :: CurrentX, Forces
         REAL, DIMENSION(size(StartX),size(StartX)) :: Hessian

         REAL, DIMENSION(:), ALLOCATABLE :: EigenValues, Factors, WrkForces
         REAL, DIMENSION(:,:), ALLOCATABLE :: EigenVectors, WrkHessian

         REAL :: V, DisplNorm, GradNorm, Factor
         INTEGER :: NIter, i, n, NOpt
         LOGICAL :: TSCheck, SteepestDescent

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.NewtonLocator : module not set" )

         ! Check the number of degree of freedom
         CALL ERROR( size(StartX) /= NDim, "PotentialModule.NewtonLocator: input array dimension mismatch" )

         ! Decide whether to look for a TS
         IF (PRESENT( TransitionState )) THEN
            TSCheck = TransitionState
         ELSE
            TSCheck = .FALSE.
         ENDIF

#if defined(LOG_FILE)
         __OPEN_LOG_FILE
         ! Print info to screen
         WRITE(__LOG_UNIT,"(/,A,/)") " NewtonLocator | Stationary point locator with Newton's method "
         IF (.NOT. TSCheck) THEN
            WRITE(__LOG_UNIT,"(A)") " NewtonLocator | Looking for a minimum of the potential"
         ELSE
            WRITE(__LOG_UNIT,"(A)") " NewtonLocator | Looking for a first order saddle point"
         END IF
#endif

         ! Number of non constrained variables
         IF ( PRESENT(Mask) ) THEN
            NOpt = COUNT(Mask)
#if defined(LOG_FILE)
            WRITE(__LOG_UNIT,"(A,I3,A,/)") " NewtonLocator | Optimization of ",NOpt," variables "
#endif
         ELSE
            NOpt = NDim
#if defined(LOG_FILE)
            WRITE(__LOG_UNIT,"(A,/)") " NewtonLocator | All the variables will be optimized "
#endif
         END IF

         ! Allocate memory
         ALLOCATE( EigenValues(NOpt), Factors(NOpt), WrkForces(NOpt) )
         ALLOCATE( EigenVectors(NOpt, NOpt), WrkHessian(NOpt, NOpt) )

         ! Start at initial position
         CurrentX = StartX
    
#if defined(LOG_FILE)
         WRITE(__LOG_UNIT,"(A12,A20,A20)") " NewtonLocator | N Iteration", "Displacement Norm", "Gradient Norm"
         WRITE(__LOG_UNIT,*)               " NewtonLocator |-------------------------------------------------------------"
#endif

         Iterations: DO NIter = 1, NMaxIter

            ! Get potential, 1st and 2nd derivatives
            V = GetPotAndForces( CurrentX, Forces )
            Hessian = GetSecondDerivatives( CurrentX )

            ! Apply constraints
            IF ( PRESENT(Mask) ) THEN
               WrkForces = ConstrainedVector( Forces, Mask )
               WrkHessian = ConstrainedHessian( Hessian, Mask )
            ELSE 
               WrkForces = Forces
               WrkHessian = Hessian
            ENDIF

            ! Compute norm of the gradient
            GradNorm  = SQRT(TheOneWithVectorDotVector(WrkForces, WrkForces) / NOpt )

            ! When the gradient is large, switch off newton and use gradient only
            IF ( GradNorm > 1.E-1 ) THEN
               SteepestDescent = .TRUE.
            ELSE
               SteepestDescent = .FALSE.
            END IF

            ! Diagonalize Hessian to transform coords to normal modes
            CALL TheOneWithDiagonalization( WrkHessian, EigenVectors, EigenValues )
            ! Transform forces to normal modes
            WrkForces = TheOneWithMatrixVectorProduct( TheOneWithTransposeMatrix(EigenVectors), WrkForces )

            ! Weigths forces with eigenvalues (hence obtain displacements)
            IF (.NOT.SteepestDescent ) THEN
               DO i = 1, NOpt
                  IF ( TSCheck ) THEN 
                     Factors(i) = 0.5 * ( ABS(EigenValues(i)) + SQRT( EigenValues(i)**2 + 4.0 * WrkForces(i)**2  ) )
                  ELSE
                     Factors(i) = EigenValues(i)
                  END IF
               END DO
               WrkForces(:) =  WrkForces(:) / Factors(:)
            END IF

            ! In case of TS search, change the sign of the step in the direction of the eigenvector with lowest eigenvalue 
            IF ( TSCheck ) THEN
               i = MINLOC( EigenValues, 1 )
               WrkForces(i) = - WrkForces(i) 
            END IF

            ! Tranform displacements back to original coordinates
            WrkForces = TheOneWithMatrixVectorProduct( EigenVectors, WrkForces )

            ! Compute norm of the displacement
            DisplNorm = SQRT(TheOneWithVectorDotVector(WrkForces, WrkForces) / NOpt )

            ! Move to new coordinates
            IF ( PRESENT(Mask) ) THEN
               n = 0
               DO i = 1, size(CurrentX)
                  IF ( Mask(i) ) THEN
                     n = n + 1 
                     CurrentX(i) = CurrentX(i) + WrkForces(n)
                  END IF
               END DO
            ELSE
               CurrentX(:) = CurrentX(:) + WrkForces(:)
            END IF

#if defined(LOG_FILE)
            ! Print info to screen
            IF ( MOD(NIter-1,NMaxIter/20) == 0 ) WRITE(__LOG_UNIT,"(A,I12,E20.6,E20.6)") &
                               " NewtonLocator |", NIter, DisplNorm, GradNorm
#endif

            ! Check convergence criteria
            IF ( GradNorm < GradThresh .AND. DisplNorm < DisplThresh ) THEN
               EXIT Iterations
            END IF

         END DO Iterations

#if defined(LOG_FILE)
         WRITE(__LOG_UNIT,"(A,I12,E20.6,E20.6)") " NewtonLocator |",NIter, DisplNorm, GradNorm
#endif
         ! Check max number of iterations
         CALL WARN( NIter == NMaxIter+1, " NewtonLocator: max number of iterations reached " )

         ! Store final point
         StationaryPoint = CurrentX


#if defined(LOG_FILE)
         ! Check the number of imaginary frequencies
         IF ( PRESENT(Mask) ) THEN
            ! Compute constrained Hessian at current position
            Hessian = GetSecondDerivatives( CurrentX )
            WrkHessian = ConstrainedHessian( Hessian, Mask )
         ELSE
            ! compute Hessian at current position
            WrkHessian = GetSecondDerivatives( CurrentX )
         END IF
         ! Diagonalize Hessian to transform coords to normal modes
         CALL TheOneWithDiagonalization( WrkHessian, EigenVectors, EigenValues )
         WRITE(__LOG_UNIT,"(/,A,I3,A,/)") " NewtonLocator | Final stationary point has ", COUNT( EigenValues < 0.0 ), &
                                   " imaginary frequency/ies "
         __CLOSE_LOG_FILE
#endif

         DEALLOCATE( EigenValues, Factors, EigenVectors, WrkHessian, WrkForces )

      END FUNCTION NewtonLocator

!===============================================================================================================================

!**************************************************************************************
!> Subrotine to test the computation of the forces by comparing
!> analytically and numerically (finite difference) derivatives
!> Points are sampled randomly (uniform distribution) in the interval 
!> between the values defined by CoordMin and CoordMax
!>
!> @param CoordMin    Array with the NDim min values of the coordinates
!> @param CoordMax    Array with the NDim max values of the coordinates
!**************************************************************************************
      SUBROUTINE TestForces( CoordMin, CoordMax ) 
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(IN) :: CoordMin, CoordMax
! 
!          REAL, DIMENSION(4), PARAMETER :: Deltas = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
!          REAL, DIMENSION(4), PARAMETER :: Coeffs = (/ +1./12., -8./12., +8./12., -1./12. /) 

         REAL, DIMENSION(6), PARAMETER :: Deltas = (/ -3.0,    -2.0,    -1.0,    +1.0,     +2.0,    +3.0  /)
         REAL, DIMENSION(6), PARAMETER :: Coeffs = (/ -1./60., 3./20.,  -3./4.,  3./4., -3./20.,   1./60. /) 


         TYPE(RNGInternalState) :: Random

         REAL, DIMENSION(:), ALLOCATABLE  :: AtPoint, Coordinates, AnalyticalDerivs, NumericalDerivs
         REAL, DIMENSION(:), ALLOCATABLE  :: Average, Deviation

         REAL    :: V
         INTEGER :: iPnt, iCoord, iDispl, NSelected
         CHARACTER(16) :: ForceUnit

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.TestForces : Module not Setup" )
         ! Check the number of degree of freedom
         CALL ERROR( size(CoordMin) /= NDim, "PotentialModule.TestForces: array dimension mismatch" )
         CALL ERROR( size(CoordMax) /= NDim, "PotentialModule.TestForces: array dimension mismatch" )

         ! Allocate memory
         ALLOCATE( Coordinates(NDim), AtPoint(NDim), AnalyticalDerivs(NDim), NumericalDerivs(NDim) )
         ALLOCATE( Average(NDim), Deviation(NDim) )

         ! Initialize random number generator
         CALL SetSeed( Random, -512 )

         Average = 0.0
         Deviation = 0.0
         NSelected = 0
         DO iPnt = 1, NPointsTest
            
            ! generate random numbers for the coordinates
            DO iCoord = 1, NDim 
               AtPoint(iCoord) = CoordMin(iCoord) + (CoordMax(iCoord) - CoordMin(iCoord)) * UniformRandomNr(Random)
            END DO

            ! Compute analytical derivatives
            V = GetPotAndForces( AtPoint, AnalyticalDerivs )

            IF ( V < 0.5 ) THEN
               NSelected = NSelected + 1

               ! Compute numerical derivatives
               NumericalDerivs(:) = 0.0
               DO iCoord = 1, NDim 
                  DO iDispl = 1, size(Deltas)
                     ! Define small displacement from the point where compute the derivative
                     Coordinates(:) = AtPoint(:)
                     Coordinates(iCoord) = Coordinates(iCoord) + Deltas(iDispl)*SmallDelta
                     ! Compute potential in the displaced coordinate
                     V = GetPotential( Coordinates )
                     ! Increment numerical derivative
                     NumericalDerivs(iCoord) = NumericalDerivs(iCoord) + Coeffs(iDispl)*V
                  END DO
               END DO
               NumericalDerivs(:) = NumericalDerivs(:) / SmallDelta

               ! Accumulate to compute average and root mean squared deviations
               Average = Average + ( NumericalDerivs - AnalyticalDerivs )
               Deviation = Deviation + ( NumericalDerivs - AnalyticalDerivs )**2

!                IF ( ANY( ABS(AnalyticalDerivs(:) - NumericalDerivs(:)) > 1.E-3 ) ) THEN 
!                PRINT*, " "
!                PRINT*, V*EnergyConversion(InternalUnits,InputUnits)
!                PRINT*, AtPoint(:)*LengthConversion(InternalUnits,InputUnits)
!                PRINT*, (AnalyticalDerivs(:) - NumericalDerivs(:)) *ForceConversion(InternalUnits,InputUnits)
!                PRINT*, SUM(NumericalDerivs(:)) *ForceConversion(InternalUnits,InputUnits)
!                PRINT*, SUM(AnalyticalDerivs(:)) *ForceConversion(InternalUnits,InputUnits)
!                END IF

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
            
         ForceUnit = TRIM(EnergyUnit(InputUnits))//"/"//TRIM(LengthUnit(InputUnits))
         DO iCoord = 1, NDim 
            WRITE(*,301) GetXLabel(iCoord), &
                         CoordMin(iCoord)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                         CoordMax(iCoord)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                         Average(iCoord)*ForceConversion(InternalUnits,InputUnits), ForceUnit,                &
                         Deviation(iCoord)*ForceConversion(InternalUnits,InputUnits), ForceUnit

         END DO

         300 FORMAT(" * Number of points where dV's are evaluated:   ",I10 )
         302 FORMAT(" * Number of points kept for the averages:      ",I10,2/)
         301 FORMAT(" * Coordinate ",A,/,&
                    "      min for the sampling:                     "F10.4,1X,A,/,& 
                    "      max for the sampling:                     "F10.4,1X,A,/,& 
                    "      average deviation:                        "F10.4,1X,A,/,& 
                    "      RMS deviation:                            "F10.4,1X,A,2/ )


         ! Deallocate memory
         DEALLOCATE( Coordinates, AtPoint, AnalyticalDerivs, NumericalDerivs )
         DEALLOCATE( Average, Deviation )

      END SUBROUTINE TestForces
      
!===============================================================================================================================

END MODULE PotentialModule

! 
! ! ************************************************************************************
! 
!       ! Setup initial conditions for the H atom + C slab for 
!       ! a simulation of vibrational relaxation
!       ! The slab is fixed in the equilibrium position with no momentum ( classical 0K )
!       ! The initial position and momenta of C and H are randomly chosen among a set of 
!       ! conditions which are given as input
!       ! data are initialized in ATOMIC UNITS
!       SUBROUTINE ZeroKelvinSlabConditions( Positions, Velocities, CHInitialConditions, RandomNr )
!          IMPLICIT NONE
!          REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
!          REAL, DIMENSION(:,:), INTENT(IN) :: CHInitialConditions
!          TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
!          INTEGER :: NDoF, iBath, NRandom, NInit
!          REAL :: Value
! 
!          ! Check the number of non frozen degree of freedom
!          NDoF = size( Positions )
!          CALL ERROR( size(Velocities) /= NDoF, "PotentialModule.ZeroKelvinSlabConditions: array dimension mismatch" )
! 
!          ! Check if the nr of dimension is compatible with the slab maximum size
!          CALL ERROR( (NDoF > 124) .OR. (NDoF < 4), "PotentialModule.ZeroKelvinSlabConditions: wrong number of DoFs" )
! 
!          ! Check the nr of starting conditions given ( there should be 8 coordinates: 4 positions and 4 momenta )
!          NRandom = size( CHInitialConditions, 1 )
!          CALL ERROR( size( CHInitialConditions, 2 ) /= 8, "PotentialModule.ZeroKelvinSlabConditions: wrong number of coords " )
!       
!          ! Set the velocities to zero
!          Velocities(:) = 0.0
! 
!          ! Set the slab in the equilibrium geometry
!          Positions(5:NDoF) = MinSlab(1:NDoF-4)
! 
!          ! Choose a random initial set of coordinates
!          NInit = CEILING( UniformRandomNr(RandomNr)*real(NRandom)  )
! 
!          ! Accordingly set position and velocity
!          Positions(1:4) = CHInitialConditions( NInit, 1:4 )
!          Velocities(1:4) = CHInitialConditions( NInit, 5:8 )
! 
!       END SUBROUTINE ZeroKelvinSlabConditions
! 
! 
!       ! Setup initial conditions for the H atom + C slab
!       ! data are initialized in ATOMIC UNITS
!       SUBROUTINE ThermalEquilibriumConditions( Positions, Velocities, Temperature, MassHydro, MassCarb, RandomNr )
!          IMPLICIT NONE
! 
!          REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
!          REAL, INTENT(IN)  :: Temperature, MassCarb, MassHydro
!          TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
!          INTEGER           :: nCarbon, NDoF
!          REAL              :: SigmaCarbonVelocity, SigmaHydroVelocity
! 
!          ! All the atoms are initially at the equilibrium position for stable chemisorption 
!          ! Value for the puckering are taken from J. Phys. Chem. B, 2006, 110, 18811-18817
!          ! Equilibrium position of zH obtained instead from plot of the PES
!          ! Velocities are sampled according to a Maxwell-Boltzmann distribution at temperature T
! 
!          ! Check the number of non frozen degree of freedom
!          NDoF = size( Positions )
!          CALL ERROR( size(Velocities) /= NDoF, "PotentialModule.ThermalEquilibriumConditions: array dimension mismatch" )
! 
!          ! Check if the nr of dimension is compatible with the slab maximum size
!          CALL ERROR( (NDoF > 124) .OR. (NDoF < 4), "PotentialModule.ThermalEquilibriumConditions: wrong number of DoFs" )
!             
!          ! Equilibrium position of H atom
!          Positions(1) = 0.0000
!          Positions(2) = 0.0000
!          Positions(3) = 1.483 / MyConsts_Bohr2Ang
! 
!          ! Equilibrium position of C1 atom
!          Positions(4) = C1Puckering
! 
!          ! Equilibrium position of the other carbon atoms 
!          DO nCarbon = 5,NDoF
!             Positions(nCarbon)   = 0.0
!          END DO
! 
!          ! Compute st deviation of Maxwell-Boltzmann distribution ( for the VELOCITY, not momenta!)
!          SigmaCarbonVelocity = sqrt( Temperature / MassCarb )
!          SigmaHydroVelocity  = sqrt( Temperature / MassHydro )
! 
!          ! Random velocities according to Maxwell-Boltzmann
!          IF ( CollinearPES ) THEN
!                Velocities(1) = 0.0
!                Velocities(2) = 0.0
!          ELSE 
!                Velocities(1) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
!                Velocities(2) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
!          END IF
!          Velocities(3) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
!          DO nCarbon = 4,NDoF
!             Velocities(nCarbon) = GaussianRandomNr( RandomNr ) * SigmaCarbonVelocity 
!          END DO
! !         Velocities(4) = 0.0      ! TO FIX EVEN C1 ATOM
! 
! !          DO nCarbon = 1, Size( EdgeCarbons ) 
! !             Velocities(EdgeCarbons(nCarbon)+3) = 0.0
! !          END DO
! 
!       END SUBROUTINE ThermalEquilibriumConditions
! 
! ! ************************************************************************************
! 
 
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
