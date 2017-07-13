!******************************************************************************
!*                      MODULE Optimize
!******************************************************************************
!
!>  \brief     Optimization of the potential
!>  \details   A given potential (defined with a function that computes \n
!>             potential and forces) is optimized with two algorithms:  \n
!>             a simple steepest descent and a Newton algorithm using   \n
!>             numerical Hessian. This Newton algorithm can also be     \n
!>             used to optimize transition states (first order saddle   \n
!>             points of the potential energy surface).
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
!******************************************************************************

MODULE Optimize
#include "preprocessoptions.cpp"
   USE FiniteDifference
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SteepLocator, NewtonLocator        !< optimization and stationary points

#if defined(LOG_FILE)
    CHARACTER(17), SAVE :: LogStr = " Optimize |"
#endif

   !> Wrapper for the optimization with/without hessian function from input
   INTERFACE NewtonLocator
      MODULE PROCEDURE NewtonLocator_WithHessian, NewtonLocator_WithoutHessian
   END INTERFACE

   CONTAINS

!==============================================================================
!                                 SUBROUTINES
!==============================================================================


!*******************************************************************************
!> Project the Hessian to impose a set of constraints defined by a logical mask.
!>
!> @param Hessian   Input matrix with the Hessian
!> @param Mask      Input vector, with the logical mask of the constraints
!> @returns         Constrained Hessian matrix of dimension COUNT(Maks)
!*******************************************************************************
   FUNCTION ConstrainedHessian( Hessian, Mask, NMask )
      IMPLICIT NONE
      REAL, DIMENSION(:,:), INTENT(IN) :: Hessian
      LOGICAL, DIMENSION(SIZE(Hessian,1)), INTENT(IN) :: Mask
      INTEGER, INTENT(IN) :: NMask
      REAL, DIMENSION(NMask,NMask) :: ConstrainedHessian
      INTEGER :: i,j, n1,n2

      CALL ERROR( NMask /= COUNT(Mask), &
         " PotentialModule.ConstrainedHessian: wrong input dimension" )

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

!==============================================================================


!******************************************************************************
!> Project a vector to impose a set of constraints defined by a logical mask.
!>
!> @param Vector    Input vector to project
!> @param Mask      Input vector, with the logical mask of the constraints
!> @param N         Input integer with the number of true values of Mask
!> @returns         Constrained vector of dimension COUNT(Maks)
!******************************************************************************
      FUNCTION ConstrainedVector( Vector, Mask, NMask )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(IN) :: Vector
         LOGICAL, DIMENSION(SIZE(Vector)), INTENT(IN) :: Mask
         INTEGER, INTENT(IN) :: NMask
         REAL, DIMENSION(NMask) :: ConstrainedVector
         INTEGER :: i, n

         CALL ERROR( NMask /= COUNT(Mask), &
            " PotentialModule.ConstrainedVector: wrong input dimension" )

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
!> @param SmallDelta    Displacement to compute numerical Hessian (freq at the last step)
!> @returns             Vector with the result of the steepest descent search
!**************************************************************************************
   FUNCTION SteepLocator( GetPotAndForces, StartX, NMaxIter, GradThresh, SmallDelta ) RESULT( StationaryPoint )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)   :: StartX
      INTEGER, INTENT(IN)              :: NMaxIter
      REAL, INTENT(IN)                 :: GradThresh
      REAL, INTENT(IN)                 :: SmallDelta
      REAL, DIMENSION(size(StartX))    :: StationaryPoint

      INTERFACE
         REAL FUNCTION GetPotAndForces( X, Force )
            REAL, DIMENSION(:), INTENT(IN)  :: X
            REAL, DIMENSION(:), INTENT(OUT) :: Force
         END FUNCTION GetPotAndForces
      END INTERFACE

      REAL, DIMENSION(size(StartX))              :: CurrentX
      REAL, DIMENSION(size(StartX),size(StartX)) :: EigenVectors, Hessian
      REAL, DIMENSION(size(StartX))              :: Forces, EigenValues
      REAL :: V, GradNorm
      INTEGER :: NIter

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
      Hessian = GetHessianFromForces( StationaryPoint, GetPotAndForces, DeltaInp=SmallDelta )
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
!> @param SmallDelta    Displacement to compute numerical Hessian (freq at the last step)
!> @param Mask              Mask to set constraints on the search
!> @param TransitionState   Logical variable to set TS search instead of minimum
!> @returns                 Vector with the result of the Newton search
!**************************************************************************************
   FUNCTION NewtonLocator_WithoutHessian( GetPotAndForces, StartX, NMaxIter, GradThresh, DisplThresh, &
                              SmallDelta, Mask, TransitionState ) RESULT( StationaryPoint )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)    :: StartX
      INTEGER, INTENT(IN)               :: NMaxIter
      REAL, INTENT(IN)                  :: GradThresh, DisplThresh
      REAL, INTENT(IN)                  :: SmallDelta
      LOGICAL, INTENT(IN), OPTIONAL     :: TransitionState
      LOGICAL, DIMENSION(size(StartX)), INTENT(IN), OPTIONAL :: Mask
      REAL, DIMENSION(size(StartX))     :: StationaryPoint

      INTERFACE
         REAL FUNCTION GetPotAndForces( X, Force )
            REAL, DIMENSION(:), INTENT(IN)  :: X
            REAL, DIMENSION(:), INTENT(OUT) :: Force
         END FUNCTION GetPotAndForces
      END INTERFACE

      REAL, DIMENSION(size(StartX)) :: CurrentX, Forces
      REAL, DIMENSION(size(StartX),size(StartX)) :: Hessian

      REAL, DIMENSION(:), ALLOCATABLE :: EigenValues, Factors, WrkForces
      REAL, DIMENSION(:,:), ALLOCATABLE :: EigenVectors, WrkHessian

      REAL :: V, DisplNorm, GradNorm, Factor
      INTEGER :: NIter, i, n, NOpt
      LOGICAL :: TSCheck, SteepestDescent

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
         NOpt = size(StartX)
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
      WRITE(__LOG_UNIT,"(A28,A20,A20)") " NewtonLocator | N Iteration", "Displacement Norm", "Gradient Norm"
      WRITE(__LOG_UNIT,"(A)")           " NewtonLocator |-------------------------------------------------------------"
#endif

      Iterations: DO NIter = 1, NMaxIter

         ! Get potential, 1st and 2nd derivatives
         V = GetPotAndForces( CurrentX, Forces )
         Hessian = GetHessianFromForces( CurrentX, GetPotAndForces, DeltaInp=SmallDelta )

         ! Apply constraints
         IF ( PRESENT(Mask) ) THEN
            WrkForces = ConstrainedVector( Forces, Mask, COUNT(Mask) )
            WrkHessian = ConstrainedHessian( Hessian, Mask, COUNT(Mask) )
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
               Factors(i) = 0.5 * ( ABS(EigenValues(i)) + SQRT( EigenValues(i)**2 + 4.0 * WrkForces(i)**2  ) )
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
         Hessian = GetHessianFromForces( CurrentX, GetPotAndForces, DeltaInp=SmallDelta )
         WrkHessian = ConstrainedHessian( Hessian, Mask, COUNT(Mask) )
      ELSE
         ! compute Hessian at current position
         WrkHessian = GetHessianFromForces( CurrentX, GetPotAndForces, DeltaInp=SmallDelta )
      END IF
      ! Diagonalize Hessian to transform coords to normal modes
      CALL TheOneWithDiagonalization( WrkHessian, EigenVectors, EigenValues )
      WRITE(__LOG_UNIT,"(/,A,I3,A,/)") " NewtonLocator | Final stationary point has ", COUNT( EigenValues < 0.0 ), &
                                 " imaginary frequency/ies "
      __CLOSE_LOG_FILE
#endif

      DEALLOCATE( EigenValues, Factors, EigenVectors, WrkHessian, WrkForces )

   END FUNCTION NewtonLocator_WithoutHessian

!===============================================================================================================================

   FUNCTION NewtonLocator_WithHessian( GetPotAndForces, GetHessian, StartX, NMaxIter, GradThresh, DisplThresh, &
                              SmallDelta, Mask, TransitionState ) RESULT( StationaryPoint )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)    :: StartX
      INTEGER, INTENT(IN)               :: NMaxIter
      REAL, INTENT(IN)                  :: GradThresh, DisplThresh
      REAL, INTENT(IN)                  :: SmallDelta
      LOGICAL, INTENT(IN), OPTIONAL     :: TransitionState
      LOGICAL, DIMENSION(size(StartX)), INTENT(IN), OPTIONAL :: Mask
      REAL, DIMENSION(size(StartX))     :: StationaryPoint

      INTERFACE
         REAL FUNCTION GetPotAndForces( X, Force )
            REAL, DIMENSION(:), INTENT(IN)  :: X
            REAL, DIMENSION(:), INTENT(OUT) :: Force
         END FUNCTION GetPotAndForces
      END INTERFACE

      INTERFACE
         FUNCTION GetHessian( X )
            REAL, DIMENSION(:), INTENT(IN)    :: X
            REAL, DIMENSION(SIZE(X),SIZE(X))  :: GetHessian
         END FUNCTION GetHessian
      END INTERFACE

      REAL, DIMENSION(size(StartX)) :: CurrentX, Forces
      REAL, DIMENSION(size(StartX),size(StartX)) :: Hessian

      REAL, DIMENSION(:), ALLOCATABLE :: EigenValues, Factors, WrkForces
      REAL, DIMENSION(:,:), ALLOCATABLE :: EigenVectors, WrkHessian

      REAL :: V, DisplNorm, GradNorm, Factor
      INTEGER :: NIter, i, n, NOpt
      LOGICAL :: TSCheck, SteepestDescent

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
         NOpt = size(StartX)
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
      WRITE(__LOG_UNIT,"(A28,A20,A20)") " NewtonLocator | N Iteration", "Displacement Norm", "Gradient Norm"
      WRITE(__LOG_UNIT,"(A)")           " NewtonLocator |-------------------------------------------------------------"
#endif

      Iterations: DO NIter = 1, NMaxIter

         ! Get potential, 1st and 2nd derivatives
         V = GetPotAndForces( CurrentX, Forces )
         Hessian = GetHessian( CurrentX )

         ! Apply constraints
         IF ( PRESENT(Mask) ) THEN
            WrkForces = ConstrainedVector( Forces, Mask, COUNT(Mask) )
            WrkHessian = ConstrainedHessian( Hessian, Mask, COUNT(Mask) )
         ELSE
            WrkForces = Forces
            WrkHessian = Hessian
         ENDIF

         ! Compute norm of the gradient
         GradNorm  = SQRT(TheOneWithVectorDotVector(WrkForces, WrkForces) / NOpt )

         ! When the gradient is large, switch off newton and use gradient only
         IF ( GradNorm > 1.E-6 ) THEN
            SteepestDescent = .TRUE.
         ELSE
            SteepestDescent = .FALSE.
         END IF

         IF ( .NOT. SteepestDescent .OR. TSCheck ) THEN

            ! Diagonalize Hessian to transform coords to normal modes
            CALL TheOneWithDiagonalization( WrkHessian, EigenVectors, EigenValues )
            ! Transform forces to normal modes
            WrkForces = TheOneWithMatrixVectorProduct( TheOneWithTransposeMatrix(EigenVectors), WrkForces )

            ! Weigths forces with eigenvalues (hence obtain displacements)
            IF (.NOT.SteepestDescent ) THEN
               DO i = 1, NOpt
                  Factors(i) = 0.5 * ( ABS(EigenValues(i)) + SQRT( EigenValues(i)**2 + 4.0 * WrkForces(i)**2  ) )
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

         ENDIF

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
         Hessian = GetHessian( CurrentX )
         WrkHessian = ConstrainedHessian( Hessian, Mask, COUNT(Mask) )
      ELSE
         ! compute Hessian at current position
         WrkHessian = GetHessian( CurrentX )
      END IF
      ! Diagonalize Hessian to transform coords to normal modes
      CALL TheOneWithDiagonalization( WrkHessian, EigenVectors, EigenValues )
      WRITE(__LOG_UNIT,"(/,A,I3,A,/)") " NewtonLocator | Final stationary point has ", COUNT( EigenValues < 0.0 ), &
                                 " imaginary frequency/ies "
      __CLOSE_LOG_FILE
#endif

      DEALLOCATE( EigenValues, Factors, EigenVectors, WrkHessian, WrkForces )

   END FUNCTION NewtonLocator_WithHessian

!==============================================================================
!                              END OF MODULE
!==============================================================================
END MODULE Optimize