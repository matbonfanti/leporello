!***************************************************************************************
!*                              MODULE PotentialAnalysis
!***************************************************************************************
!
!>  \brief     Subroutines for an Potential Analysis of the 3D system PES
!>  \details   This module contains subroutines to analyze the potential  \n
!>             for the system: energies of the relevant asymptotic geometries, \n
!>             plots of the potential, computing the MEP \n
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             15 January 2015
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!>                 
!***************************************************************************************
MODULE PotentialAnalysis
#include "preprocessoptions.cpp"
   USE SharedData
   USE InputField
   USE UnitConversion
   USE PotentialModule
   USE PrintTools

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: PotentialAnalysis_ReadInput, PotentialAnalysis_Initialize, PotentialAnalysis_Run, PotentialAnalysis_Dispose

   !> Dimensions of the potential
   INTEGER :: NDim

   !> Input variables for the 3D potential plot 
   REAL :: ZHIncMin, ZHIncMax, ZHTarMin, ZHTarMax, ZCMin, ZCMax           !< boundaries of the grid
   INTEGER :: NpointZHInc, NpointZHTar, NpointZC                          !< nr of points of the grid
   REAL :: GridSpacing                                                    !< grid spacing

   ! Parameters for optimizations
   INTEGER :: MaxOptSteps  !< Maximum number of steps for optimizations
   REAL :: OptThreshold    !< Convergence criterium for optimizations

   ! Parameters for Minimum Energy Path computation
   INTEGER :: MaxMEPNrSteps  !< Max number of MEP steps
   REAL    :: MEPStep        !< Size of each MEP step


      CONTAINS

!===============================================================================================================================

!**************************************************************************************
!> Read from input unit the variable which are specific for the
!> analysis of the potential
!>
!> @param InputData     Datatype with an already setup input unit
!**************************************************************************************
   SUBROUTINE PotentialAnalysis_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

      ! parameters for the plot 
      CALL SetFieldFromInput( InputData, "Plot_GridSpacing", GridSpacing )
      GridSpacing = GridSpacing * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_ZHIncMin", ZHIncMin )
      ZHIncMin = ZHIncMin * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_ZHIncMax", ZHIncMax )
      ZHIncMax = ZHIncMax * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_ZHTarMin", ZHTarMin )
      ZHTarMin = ZHTarMin * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_ZHTarMax", ZHTarMax )
      ZHTarMax = ZHTarMax * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_ZCMin", ZCMin )
      ZCMin = ZCMin * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_ZCMax", ZCMax )
      ZCMax = ZCMax * LengthConversion(InputUnits, InternalUnits)

      ! Set variables for minima and TS optimization
      CALL SetFieldFromInput( InputData, "MaxOptSteps", MaxOptSteps, 10**6 )
      CALL SetFieldFromInput( InputData, "OptThreshold", OptThreshold, 1.E-6 )

      ! Set variables for Minimum Energy Path computation
      CALL SetFieldFromInput( InputData, "MaxMEPNrSteps", MaxMEPNrSteps, 1000 )
      CALL SetFieldFromInput( InputData, "MEPStep", MEPStep, 0.001 )
      MEPStep = MEPStep * LengthConversion(InputUnits, InternalUnits)

   END SUBROUTINE PotentialAnalysis_ReadInput

!===============================================================================================================================

!**************************************************************************************
!> Initialization of the data for the analysis of the Potential:
!> memory allocation, variable initialization and data type setup.
!>
!**************************************************************************************
   SUBROUTINE PotentialAnalysis_Initialize()
      IMPLICIT NONE

      NDim = 3

      ! Allocate memory and initialize vectors for positions, forces and masses
      ALLOCATE( X(NDim), A(NDim), MassVector(NDim) )
      MassVector = (/ MassHInc, MassHTar, MassC /)

   END SUBROUTINE PotentialAnalysis_Initialize

!===============================================================================================================================

!**************************************************************************************
!> Run the Potential Analysis.
!>
!**************************************************************************************
   SUBROUTINE PotentialAnalysis_Run()
      IMPLICIT NONE

      !> Memory to store data for the 3D potential plot 
      TYPE(VTKInfo), SAVE :: PotentialPlot                                   !< VTK data type
      REAL, DIMENSION(:), ALLOCATABLE :: PotentialArray                      !< array to store 2D potential
      REAL, DIMENSION(:), ALLOCATABLE :: ZCArray, ZHIncArray, ZHTarArray     !< array with coordinates grid

      !> Variables for the MEP analysis
      REAL, DIMENSION(:), ALLOCATABLE   :: XStart
      REAL, DIMENSION(:), ALLOCATABLE   :: EigenFreq
      REAL, DIMENSION(:,:), ALLOCATABLE :: EigenModes, Hessian
      REAL    :: EStart, E, GradNorm

      INTEGER :: i, j, k, nPoint, MaxStep
      LOGICAL :: Check

      PRINT "(2/,A)",    " ***************************************************"
      PRINT "(A,F10.5)", "               PES ANALYSIS "
      PRINT "(A,/)" ,    " ***************************************************"


      ! =========================================================
      !                (1) 3D pes in VTK format
      ! =========================================================

      PRINT "(/,A)",    " **** Write 3D cut of the PES to output VTR ****"

      ! Set grid dimensions
      NpointZHInc = INT((ZHIncMax-ZHIncMin)/GridSpacing) + 1
      NpointZHTar = INT((ZHTarMax-ZHTarMin)/GridSpacing) + 1
      NpointZC = INT((ZCMax-ZCMin)/GridSpacing) + 1

      ! Allocate temporary array to store potential data and coord grids
      ALLOCATE( PotentialArray( NpointZHInc * NpointZHTar * NpointZC ),   &
                ZCArray( NpointZC ), ZHIncArray( NpointZHInc ), ZHTarArray( NpointZHTar ) )
    
      ! Define coordinate grids
      ZCArray    = (/ ( ZCMin    + GridSpacing*(i-1), i=1,NpointZC)    /)
      ZHIncArray = (/ ( ZHIncMin + GridSpacing*(i-1), i=1,NpointZHInc) /)
      ZHTarArray = (/ ( ZHTarMin + GridSpacing*(i-1), i=1,NpointZHTar) /)

      ! Open VTK file
      CALL VTK_NewRectilinearSnapshot(PotentialPlot, FileName="V3D_Plot", X=ZHIncArray*LengthConversion(InternalUnits,InputUnits),& 
                   Y=ZHTarArray*LengthConversion(InternalUnits,InputUnits), Z=ZCArray*LengthConversion(InternalUnits,InputUnits) )

      nPoint = 0
      ! Cycle over the ZC coordinate values
      DO k = 1, NpointZC
         ! Cycle over the ZHTar coordinate values
         DO j = 1, NpointZHTar
            ! Cycle over the ZHTar coordinate values
            DO i = 1, NpointZHInc
               nPoint = nPoint + 1
               ! Compute potential 
               PotentialArray(nPoint) = GetPotential( (/ ZHIncArray(i), ZHTarArray(j), ZCArray(k) /) )
            END DO
         END DO
      END DO
      ! Print the potential to vtk file
      CALL VTK_AddScalarField( PotentialPlot, Name="Potential", Field=PotentialArray*EnergyConversion(InternalUnits,InputUnits), &
                                                                                                            LetFileOpen=.FALSE. )

      WRITE(*,"(/,A)") " * 3D PES as a func of zH_inc, zH_tar and zC written as VTR to file V3D_Plot.vtr"

      ! Deallocate memory
      DEALLOCATE( PotentialArray, ZCArray, ZHIncArray, ZHTarArray )

      ! =========================================================
      !                 (2) minimum energy path
      ! =========================================================

      PRINT "(/,A)",    " **** Starting geometry of the MEP ****"

      ! Allocate arrays for this section
      ALLOCATE( XStart(NDim), Hessian(NDim, NDim), EigenFreq(NDim), EigenModes(NDim,NDim) )

      ! guess reasonable coordinates of the minimum of the PES
      X(1) = 10.0/MyConsts_Bohr2Ang
      X(2) = (0.35+1.1)/MyConsts_Bohr2Ang
      X(3) = 0.35/MyConsts_Bohr2Ang

      ! Find minimum by Newton's optimization
      XStart = NewtonLocator( X, MaxOptSteps, OptThreshold, OptThreshold, (/ .FALSE., .TRUE., .TRUE. /) )
      ! Computing the energy at this geometry
      EStart = GetPotAndForces( XStart, A )

      ! Compute normal modes 
      ! Numerical hessian of the system potential
      Hessian = GetHessian( XStart )
      ! Diagonalize the hessian
      CALL TheOneWithDiagonalization( Hessian, EigenModes, EigenFreq )

      WRITE(*,"(/,A)") " Starting geometry and energy of the MEP " 
      DO i = 1, NDim 
         WRITE(*,501) GetXLabel(i), X(i)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits)
      END DO
      WRITE(*,502) EStart*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

      WRITE(*,"(/,A)") " Potential normal modes at the starting geometry" 
      DO i = 1, NDim 
         WRITE(*,503) i, SQRT(EigenFreq(i))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                      TRIM(LengthUnit(InternalUnits)), EigenModes(:,i)
      END DO

      PRINT "(/,A)",    " **** Minimum energy path ****"

      ! Write header lines screen 
      WRITE(*,700) TRIM(LengthUnit(InputUnits)), TRIM(LengthUnit(InputUnits)), TRIM(EnergyUnit(InputUnits))

      ! Start from XStart
      X(:) = XStart; E = GetPotential( X )
      ! Write to output file and to screen the starting point
      WRITE(*,601) -1, X(1)*LengthConversion(InternalUnits, InputUnits), &
         X(2)*LengthConversion(InternalUnits, InputUnits), E*EnergyConversion(InternalUnits, InputUnits)

      ! First follow the incident direction until significant gradient is found
      DO i = 1, 1000
         ! Move along the entrance channel
         X(1) = X(1) - 0.01
         ! Compute norm of the gradient
         E = GetPotAndForces( X, A ); GradNorm = SQRT(TheOneWithVectorDotVector(A, A))
         ! Check when the gradient is large
         IF ( GradNorm > 0.001 ) EXIT
      END DO

      ! Write to output file and to screen the starting point
      WRITE(*,601) 0, X(1)*LengthConversion(InternalUnits, InputUnits), &
         X(2)*LengthConversion(InternalUnits, InputUnits), E*EnergyConversion(InternalUnits, InputUnits)

      DO i = 2, MaxMEPNrSteps
         ! Following steps
         Check = FollowGradient( X, MEPStep )
         IF ( .NOT. Check )  EXIT
      END DO
      MaxStep = i - 1

      ! Write to screen final step
      E = GetPotential( X )
      WRITE(*,601) MaxStep, X(1)*LengthConversion(InternalUnits, InputUnits), &
         X(2)*LengthConversion(InternalUnits, InputUnits), E*EnergyConversion(InternalUnits, InputUnits)

      WRITE(*,"(/,A)") " * Last step reached... MEP written to file ________"

      501 FORMAT( " * ",A5,23X,1F15.6,1X,A,/ )
      502 FORMAT( " * Energy",22X,1F15.6,1X,A,/ )
      503 FORMAT( " Normal Mode ",I5," - frequency: ",1F15.2,1X,A, /, &
                  "    mass-scaled coords of the normal mode / ",A," : ",4F12.6, / )

      600 FORMAT ( A12,A20,A20,A20 )
      601 FORMAT ( I12,F20.6,F20.6,F20.6 )
      602 FORMAT ( F20.6,F20.6,F20.6,F20.6 )

      700 FORMAT ( "#  N Step   ", "   zH Coord / ", A6, "   zC Coord / ", A6, "     Energy / ", A6, /,  &
                   "#---------------------------------------------------------------------------------" )
      701 FORMAT ( "#---------------------------------------------------------------------------------" )

      ! Deallocate memory
      DEALLOCATE( XStart, Hessian, EigenFreq, EigenModes )

   END SUBROUTINE PotentialAnalysis_Run

!===============================================================================================================================

!**************************************************************************************
!> Free the memory which has been used for the simulation.
!>
!**************************************************************************************
   SUBROUTINE PotentialAnalysis_Dispose()
      IMPLICIT NONE

      DEALLOCATE( X, A, MassVector )

   END SUBROUTINE PotentialAnalysis_Dispose

!===============================================================================================================================

!**************************************************************************************
!> Compute the Hessian of the potential in multiplied by the square root 
!> of the masses, as it is needed for the normal modes analysis  
!>
!> @param AtPoint   Input vector with the coordinates where to compute H
!> @returns         Hessian matrix of the potential in AtPoint
!**************************************************************************************
      FUNCTION GetHessian( AtPoint ) RESULT( Hessian )
         REAL, DIMENSION(:), INTENT(IN)                 :: AtPoint
         REAL, DIMENSION(size(AtPoint),size(AtPoint))   :: Hessian
         INTEGER :: i,j

         ! Check the number of degree of freedom
         CALL ERROR( size(AtPoint) /= NDim, "PotentialModule.GetHessian: input array dimension mismatch" )

         ! Use subroutine GetSecondDerivatives to compute the matrix of the 2nd derivatives
         Hessian(:,:) = GetSecondDerivatives( AtPoint ) 
         ! Numerical hessian of the potential in mass weighted coordinates
         DO j = 1, NDim
            DO i = 1, NDim
               Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
            END DO
         END DO

      END FUNCTION GetHessian

!===============================================================================================================================

!**************************************************************************************
!> Follow gradient integrating the MEP differential equation, 
!> solved with a runge-kutta 4th order step.
!>
!> @param X            Current position
!> @param StepLength   Length of the runge-kutta step
!> @returns            Logical mask to set constraints on the position
!**************************************************************************************
   LOGICAL FUNCTION FollowGradient( X, StepLength, Mask )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(INOUT)  :: X
      REAL, INTENT(IN)                   :: StepLength
      LOGICAL, DIMENSION(SIZE(X)), INTENT(IN), OPTIONAL :: Mask

      REAL, DIMENSION(SIZE(X)) :: Forces
      REAL, DIMENSION(SIZE(X)) :: Kappa1, Kappa2, Kappa3, Kappa4, Step

      REAL :: V, FMax
      INTEGER :: i

      ! compute forces at current position and the maximum component
      V = GetPotAndForces( X, Forces )
      IF ( PRESENT(Mask) ) THEN
         FMax = MAXVAL( ABS(Forces), Mask )
      ELSE 
         FMax = MAXVAL( ABS(Forces) )
      END IF

      ! Check if forces would result in a negligible step
      IF ( FMax < 1.E-4 )  THEN
         ! Do not take any step
         FollowGradient = .FALSE.
      ELSE
         ! Make a runge kutta 4th order step
         Kappa1 = NormalizedForces( Forces, Mask )
         V = GetPotAndForces( X + 0.5*ABS(StepLength)*Kappa1, Forces )
         Kappa2 = NormalizedForces( Forces, Mask )
         V = GetPotAndForces( X + 0.5*ABS(StepLength)*Kappa2, Forces )
         Kappa3 = NormalizedForces( Forces, Mask )
         V = GetPotAndForces( X + ABS(StepLength)*Kappa3, Forces )
         Kappa4 = NormalizedForces( Forces, Mask )
         Step = ( Kappa1/6.0 + Kappa2/3.0 + Kappa3/3.0 + Kappa4/6.0 )

         ! move uncostrained coordinates in the direction of the forces
         IF ( PRESENT(Mask) ) THEN
            DO i = 1, size(X)
               IF ( Mask(i) ) X(i) = X(i) + ABS(StepLength) * Step(i)
            END DO
         ELSE
            X(:) = X(:) + ABS(StepLength) * Step(:)
         END IF
         ! Return TRUE value
         FollowGradient = .TRUE.
      END IF

   END FUNCTION FollowGradient

   FUNCTION NormalizedForces( Forces, Mask )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: Forces
      LOGICAL, DIMENSION(size(X)), INTENT(IN) :: Mask
      REAL, DIMENSION(size(X)) :: NormalizedForces
      REAL :: VectorNorm

      VectorNorm = SQRT(TheOneWithVectorDotVector( ConstrainedVector(Forces,Mask), &
                     ConstrainedVector(Forces,Mask) ))
      NormalizedForces(:) = Forces(:) / VectorNorm
      NormalizedForces(:) = NormalizedForces(:) / SQRT(MassVector(:))
   END FUNCTION NormalizedForces

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

END MODULE PotentialAnalysis
