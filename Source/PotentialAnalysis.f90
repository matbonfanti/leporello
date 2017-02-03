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
!>  \todo N.A.
!>                 
!***************************************************************************************
MODULE PotentialAnalysis
#include "preprocessoptions.cpp"
   USE SharedData
   USE InputField
   USE UnitConversion
   USE PotentialModule
   USE PrintTools
   USE FiniteDifference

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: PotentialAnalysis_ReadInput, PotentialAnalysis_Initialize, PotentialAnalysis_Run, PotentialAnalysis_Dispose

   !> Dimensions of the potential
   INTEGER, SAVE :: NDim

   !> Input variables for the 3D potential plot 
   REAL :: ZHIncMin, ZHIncMax, ZHTarMin, ZHTarMax, ZCMin, ZCMax           !< boundaries of the grid
   INTEGER :: NpointZHInc, NpointZHTar, NpointZC                          !< nr of points of the grid
   REAL :: GridSpacing                                                    !< grid spacing

   ! Parameters for optimizations
   INTEGER :: MaxOptSteps  !< Maximum number of steps for optimizations
   REAL :: OptThreshold    !< Convergence criterium for optimizations

   ! Parameters for Minimum Energy Path computation
   INTEGER :: MaxMEPNrSteps      !< Max number of MEP steps
   INTEGER :: PrintMEPNrSteps    !< Max number of MEP steps printed to file output
   REAL    :: MEPStep            !< Size of each MEP step
   REAL    :: StartMEPThreshold  !< gradient threshold to decide when to start to follow the gradient

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
      REAL :: LengthConv

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
      IF ( GetSystemDimension( ) == 3 ) THEN
         CALL SetFieldFromInput( InputData, "Plot_ZCMin", ZCMin )
         ZCMin = ZCMin * LengthConversion(InputUnits, InternalUnits)
         CALL SetFieldFromInput( InputData, "Plot_ZCMax", ZCMax )
         ZCMax = ZCMax * LengthConversion(InputUnits, InternalUnits)
      ENDIF

      ! Set variables for minima and TS optimization
      CALL SetFieldFromInput( InputData, "MaxOptSteps", MaxOptSteps, 10**6 )
      CALL SetFieldFromInput( InputData, "OptThreshold", OptThreshold, 1.E-6 )
      OptThreshold = OptThreshold * ForceConversion(InputUnits, InternalUnits)
      
      ! Set variables for Minimum Energy Path computation
      CALL SetFieldFromInput( InputData, "MaxMEPNrSteps", MaxMEPNrSteps, 1000 )
      CALL SetFieldFromInput( InputData, "MEPStep", MEPStep, 0.001 )
      MEPStep = MEPStep * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "PrintMEPNrSteps", PrintMEPNrSteps, MaxMEPNrSteps )
      CALL SetFieldFromInput( InputData, "StartMEPThreshold", StartMEPThreshold, 0.05 )
      StartMEPThreshold = StartMEPThreshold * ForceConversion(InputUnits, InternalUnits)

      ! SCREEN LOG OF THE INPUT VARIABLES
      
      LengthConv = LengthConversion(InternalUnits,InputUnits)

      WRITE(*, 900) GridSpacing*LengthConv, LengthUnit(InputUnits)
      WRITE(*, 903) GetXLabel(1), ZHIncMin*LengthConv, ZHIncMax*LengthConv, LengthUnit(InputUnits)
      WRITE(*, 903) GetXLabel(2), ZHTarMin*LengthConv, ZHTarMax*LengthConv, LengthUnit(InputUnits)
      IF ( GetSystemDimension( ) == 3 ) WRITE(*, 903) GetXLabel(3), ZCMin*LengthConv, ZCMax*LengthConv, LengthUnit(InputUnits)

      WRITE(*, 901) MaxOptSteps, OptThreshold*ForceConversion(InternalUnits,InputUnits), &
                    TRIM(EnergyUnit(InputUnits))//"/"//TRIM(LengthUnit(InputUnits))
                    
      WRITE(*, 902) StartMEPThreshold*ForceConversion(InternalUnits,InputUnits), &
                    TRIM(EnergyUnit(InputUnits))//"/"//TRIM(LengthUnit(InputUnits)), &
                    MEPStep*LengthConv, LengthUnit(InputUnits), MaxMEPNrSteps, PrintMEPNrSteps

      900 FORMAT(" * Plot of a 3D cut of the potential in VTK format ",         /,&
                 " * Grid spacing :                                ",F10.4,1X,A   )
      903 FORMAT(" * Interval along ",A5, ":                       ",F6.2," / ",F6.2,1X,A)    

      901 FORMAT(" * Stationary states optimization with Newton's method ",     /,&
                 " * Max nr of steps of the optimization:          ",I10,       /,&
                 " * Threshold on the gradient:                    ",E10.2,1X,A,/)

      902 FORMAT(" * Minimum energy path computation ",                         /,&
                 " * Gradient thresh., start MEP in in-channel:    ",F10.4,1X,A,/,&
                 " * Length of the 4th order RK integration step:  ",F10.4,1X,A,/,&
                 " * Max nr of integration steps:                  ",I10,       /,&
                 " * Max nr of steps printed to output:            ",I10,       /)
      
   END SUBROUTINE PotentialAnalysis_ReadInput

!===============================================================================================================================

!**************************************************************************************
!> Initialization of the data for the analysis of the Potential:
!> memory allocation, variable initialization and data type setup.
!>
!**************************************************************************************
   SUBROUTINE PotentialAnalysis_Initialize()
      IMPLICIT NONE

      NDim = GetSystemDimension()

      ! Allocate memory and initialize vectors for positions, forces and masses
      ALLOCATE( X(NDim), A(NDim), MassVector(NDim) )
      IF ( GetSystemDimension( ) == 3 ) THEN
         MassVector = (/ MassHInc, MassHTar, MassC /)
      ELSE
         MassVector = (/ MassHInc, MassHTar /)
      ENDIF

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
      REAL, DIMENSION(:), ALLOCATABLE    :: XStart
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LogMask
      REAL, DIMENSION(:), ALLOCATABLE    :: EigenFreq
      REAL, DIMENSION(:,:), ALLOCATABLE  :: EigenModes, Hessian
      REAL, DIMENSION(:,:), ALLOCATABLE  :: StoreMEP, StoreApproach
      REAL    :: EStart, E, GradNorm
      INTEGER :: MEPFullUnit

      !> Data for printing MEP in VTK format
      TYPE(VTKInfo), SAVE :: MEPTraj

      
      INTEGER :: i, j, k, nPoint, MaxStep, NPrint
      LOGICAL :: Check

      PRINT "(2/,A)",    " ***************************************************"
      PRINT "(A,F10.5)", "               PES ANALYSIS "
      PRINT "(A,/)" ,    " ***************************************************"


      ! =========================================================
      !                (1) 3D pes in VTK format
      ! =========================================================

      PRINT "(/,A)",    " **** Write cut of the PES to output VTR ****"

      ! Set grid dimensions
      NpointZHInc = INT((ZHIncMax-ZHIncMin)/GridSpacing) + 1
      NpointZHTar = INT((ZHTarMax-ZHTarMin)/GridSpacing) + 1
      IF ( GetSystemDimension( ) == 3 ) THEN
         NpointZC = INT((ZCMax-ZCMin)/GridSpacing) + 1
      ELSE 
         NpointZC = 1
      END IF

      ! Allocate temporary array to store potential data and coord grids
      ALLOCATE( PotentialArray( NpointZHInc * NpointZHTar * NpointZC ),   &
                ZCArray( NpointZC ), ZHIncArray( NpointZHInc ), ZHTarArray( NpointZHTar ) )
    
      ! Define coordinate grids
      ZHIncArray = (/ ( ZHIncMin + GridSpacing*(i-1), i=1,NpointZHInc) /)
      ZHTarArray = (/ ( ZHTarMin + GridSpacing*(i-1), i=1,NpointZHTar) /)
      IF ( GetSystemDimension( ) == 3 ) THEN
         ZCArray    = (/ ( ZCMin    + GridSpacing*(i-1), i=1,NpointZC)    /)
      ELSE 
         ZCArray    = (/ 0.0 /)
      END IF
      
      ! Open VTK file
      CALL VTK_NewRectilinearSnapshot(PotentialPlot,FileName="V3D_Plot",X=ZHIncArray*LengthConversion(InternalUnits,InputUnits),& 
                   Y=ZHTarArray*LengthConversion(InternalUnits,InputUnits), Z=ZCArray*LengthConversion(InternalUnits,InputUnits))

      nPoint = 0
      ! Cycle over the ZC coordinate values
      DO k = 1, NpointZC
         ! Cycle over the ZHTar coordinate values
         DO j = 1, NpointZHTar
            ! Cycle over the ZHTar coordinate values
            DO i = 1, NpointZHInc
               nPoint = nPoint + 1
               ! Compute potential 
               IF ( GetSystemDimension( ) == 3 ) THEN
                  PotentialArray(nPoint) = GetPotential( (/ ZHIncArray(i), ZHTarArray(j), ZCArray(k) /) )
               ELSE 
                  PotentialArray(nPoint) = GetPotential( (/ ZHIncArray(i), ZHTarArray(j) /) )
               END IF
            END DO
         END DO
      END DO
      ! Print the potential to vtk file
      CALL VTK_AddScalarField( PotentialPlot,Name="Potential",Field=PotentialArray*EnergyConversion(InternalUnits,InputUnits), &
                                                                                                            LetFileOpen=.FALSE.)
      IF ( GetSystemDimension( ) == 3 ) THEN
         WRITE(*,"(/,A)") " * 3D PES as a func of zH_inc, zH_tar and zC written as VTR to file V3D_Plot.vtr"
      ELSE 
         WRITE(*,"(/,A)") " * 2D PES as a func of zH_inc, zH_tar written as VTR to file V3D_Plot.vtr"
      END IF

      ! Deallocate memory
      DEALLOCATE( PotentialArray, ZCArray, ZHIncArray, ZHTarArray )

      ! =========================================================
      !                 (2) asymptotic initial geometry
      ! =========================================================

      PRINT "(/,A)",    " **** Starting geometry of the MEP ****"

      ! Allocate arrays for this section
      ALLOCATE( XStart(NDim), Hessian(NDim, NDim), EigenFreq(NDim), EigenModes(NDim,NDim), LogMask(NDim) )
      ALLOCATE( StoreApproach(NDim+2,0:1000), StoreMEP(NDim+2,0:MaxMEPNrSteps) )
      
      ! guess reasonable coordinates of the minimum of the PES
      X(1) = 10.0/MyConsts_Bohr2Ang
      X(2) = (0.35+1.1)/MyConsts_Bohr2Ang
      IF ( GetSystemDimension( ) == 3 ) X(3) = 0.35/MyConsts_Bohr2Ang

      ! Find minimum by Newton's optimization
      LogMask(:) = .TRUE.
      LogMask(1) = .FALSE.
      
      XStart = NewtonLocator( X, MaxOptSteps, OptThreshold, 1.0, LogMask )
      ! Computing the energy at this geometry
      EStart = GetPotAndForces( XStart, A )

      ! Compute normal modes 
      ! Numerical hessian of the system potential
      Hessian = GetMassScaledHessian( XStart )
      ! Diagonalize the hessian
      CALL TheOneWithDiagonalization( Hessian, EigenModes, EigenFreq )

      WRITE(*,"(/,A)") " Starting geometry and energy of the MEP " 
      DO i = 1, NDim 
         WRITE(*,501) GetXLabel(i), XStart(i)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits)
      END DO
      WRITE(*,502) EStart*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

      WRITE(*,"(/,A)") " Potential normal modes at the starting geometry" 
      DO i = 1, NDim 
         IF ( EigenFreq(i) > 0.0 ) THEN
            WRITE(*,503) i, SQRT(EigenFreq(i))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                      TRIM(LengthUnit(InternalUnits)), EigenModes(:,i)
         ELSE
            WRITE(*,504) i, SQRT(-EigenFreq(i))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                      TRIM(LengthUnit(InternalUnits)), EigenModes(:,i)
         END IF
      END DO

      ! =========================================================
      !                 (3) minimum energy path
      ! =========================================================
      
      PRINT "(/,A,/)",    " **** Minimum energy path ****"
      WRITE(*,701)

      ! Write header lines screen 
      WRITE(*,700) TRIM(LengthUnit(InputUnits)),TRIM(LengthUnit(InputUnits)),&
                   TRIM(LengthUnit(InputUnits)),TRIM(EnergyUnit(InputUnits))

      ! Start from XStart
      X(:) = XStart; E = GetPotential( X )
      ! Write to output file and to screen the starting point
      WRITE(*,601) -1., E*EnergyConversion(InternalUnits, InputUnits), X(:)*LengthConversion(InternalUnits, InputUnits)
      ! Store step in the data array
      StoreApproach(:,0) = (/ 0., E, X(1:NDim) /)

      ! First follow the incident direction until significant gradient is found
      DO i = 1, 1000
         ! Move along the entrance channel
         X(1) = X(1) - 0.05
         ! Compute norm of the gradient
         E = GetPotAndForces( X, A ); GradNorm = SQRT(TheOneWithVectorDotVector(A, A))
         ! Check when the gradient is large
         IF ( GradNorm > StartMEPThreshold ) EXIT
         StoreApproach(:,i) = (/ i*0.05, E, X(1:NDim) /)
      END DO
      NPrint = i - 1

      ! Open output file and write path to file
      MEPFullUnit = LookForFreeUnit()
      OPEN( FILE="MEP_Coord_Energy.dat", UNIT=MEPFullUnit )
      ! Header of the output file
      WRITE(MEPFullUnit,700) TRIM(EnergyUnit(InputUnits)), TRIM(LengthUnit(InputUnits)), &
                             TRIM(LengthUnit(InputUnits)), TRIM(LengthUnit(InputUnits))
      ! Write approach along the entrance channel
      DO i = 1, NPrint
         WRITE(MEPFullUnit,601) (StoreApproach(1,i)-StoreApproach(1,NPrint))*LengthConversion(InternalUnits, InputUnits), &         
                                StoreApproach(2,i)*EnergyConversion(InternalUnits, InputUnits), &
                                StoreApproach(3:2+NDim,i)*LengthConversion(InternalUnits, InputUnits)
      END DO

      ! Write to output file and to screen the starting point
      WRITE(*,601) 0., E*EnergyConversion(InternalUnits, InputUnits), X(:)*LengthConversion(InternalUnits, InputUnits)
      ! Store step in the data array
      StoreMEP(:,0) = (/ 0., E, X(1:NDim) /)

      NPrint = 0
      DO i = 1, MaxMEPNrSteps-1
         ! Following steps
         Check = FollowGradient( X, MEPStep )
         E = GetPotential( X )
         IF ( .NOT. Check )  EXIT
         IF ( MOD(i,MaxMEPNrSteps/20) == 0 ) THEN
            WRITE(*,601) i*MEPStep*LengthConversion(InternalUnits, InputUnits), E*EnergyConversion(InternalUnits, InputUnits), &
                         X(:)*LengthConversion(InternalUnits, InputUnits)
         END IF
         ! Store step in the data array (only for a number of steps equal to PrintMEPNrSteps )
         IF ( MOD(i, MaxMEPNrSteps/PrintMEPNrSteps) == 0 ) THEN
            NPrint = NPrint + 1
            StoreMEP(:,NPrint) = (/ i*MEPStep, E, X(1:NDim) /)
         END IF
      END DO

      MaxStep = i - 1
      ! Write to screen final step
      WRITE(*,701)
      E = GetPotential( X )
      WRITE(*,601) MaxStep*MEPStep*LengthConversion(InternalUnits, InputUnits), E*EnergyConversion(InternalUnits, InputUnits), &
                         X(:)*LengthConversion(InternalUnits, InputUnits)
      WRITE(*,701)
      ! Store final step in the data array 
      NPrint = NPrint + 1
      StoreMEP(:,NPrint) = (/ MaxStep*MEPStep, E, X(1:NDim) /)

      ! Write the MEP
      DO i = 1, NPrint
         WRITE(MEPFullUnit,601) StoreMEP(1,i)*LengthConversion(InternalUnits, InputUnits), &         
                                StoreMEP(2,i)*EnergyConversion(InternalUnits, InputUnits), &
                                StoreMEP(3:2+NDim,i)*LengthConversion(InternalUnits, InputUnits)
      END DO
      ! Close file
      CLOSE(MEPFullUnit)
      
      ! Write the MEP to VTV files
      CALL VTK_WriteTrajectory( MEPTraj,StoreMEP(3:2+NDim,1:NPrint)*LengthConversion(InternalUnits, InputUnits),"MEP_Trajectory")

      
      WRITE(*,"(/,A)") " * Last step reached... MEP written to file ________"

      ! =========================================================
      !                 (4) asymptotic final geometry
      ! =========================================================

      PRINT "(/,A)",    " **** Final asymptotic geometry ****"

      ! guess reasonable coordinates of the minimum of the PES in the asymptotic out channel
      X(1) = 20.0/MyConsts_Bohr2Ang
      X(2) = X(1)-0.35/MyConsts_Bohr2Ang
      IF ( GetSystemDimension( ) == 3 ) X(3) = 0.0/MyConsts_Bohr2Ang

      ! Find minimum by Newton's optimization
      LogMask(:) = .TRUE.
      X = NewtonLocator( X, MaxOptSteps, OptThreshold, 1.0, LogMask )
      ! Computing the energy at this geometry
      E = GetPotAndForces( X, A )

      ! Compute normal modes 
      ! Numerical hessian of the system potential
      Hessian = GetMassScaledHessian( X )
      ! Diagonalize the hessian
      CALL TheOneWithDiagonalization( Hessian, EigenModes, EigenFreq )

      WRITE(*,"(/,A)") " Final asymptotic geometry and energy of the MEP " 
      DO i = 1, NDim 
         WRITE(*,501) GetXLabel(i), X(i)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits)
      END DO
      WRITE(*,502) E*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

      WRITE(*,"(/,A)") " Potential normal modes at the final asymptotic geometry" 
      DO i = 1, NDim 
         IF ( EigenFreq(i) > 0.0 ) THEN
            WRITE(*,503) i, SQRT(EigenFreq(i))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                      TRIM(LengthUnit(InternalUnits)), EigenModes(:,i)
         ELSE
            WRITE(*,504) i, SQRT(-EigenFreq(i))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                      TRIM(LengthUnit(InternalUnits)), EigenModes(:,i)
         END IF
      END DO      
      
      501 FORMAT( " * ",A5,23X,1F15.6,1X,A )
      502 FORMAT( " * Energy",22X,1F15.6,1X,A )
      503 FORMAT( " Normal Mode ",I5," - real frequency: ",1F15.2,1X,A, /, &
                  "    mass-scaled coords of the normal mode / ",A," : ",4F12.6, / )
      504 FORMAT( " Normal Mode ",I5," - imag frequency: ",1F15.2,1X,A, /, &
                  "    mass-scaled coords of the normal mode / ",A," : ",4F12.6, / )

      601 FORMAT ( F12.4,F20.6,F20.6,F20.6,F20.6 )
      700 FORMAT ( "#       Step", 10X, "    Energy", 10X, "    zH inc", 10X, "    zH tar", 10X, "        zC", /,  &
                   "#           ", A20              , A20              , A20              , A20              , /,  &
                   "#---------------------------------------------------------------------------------------------" )
      701 FORMAT ( "#---------------------------------------------------------------------------------------------" )

      ! Deallocate memory
      DEALLOCATE( XStart, Hessian, EigenFreq, EigenModes, LogMask )
      DEALLOCATE( StoreMEP )

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
!> Compute the Hessian of the potential in mass square root scaled coordinates
!> as it is needed for the normal modes analysis  
!>
!> @param AtPoint   Input vector with the coordinates where to compute H
!> @returns         Hessian matrix of the potential in AtPoint
!**************************************************************************************
      FUNCTION GetMassScaledHessian( AtPoint ) RESULT( Hessian )
         REAL, DIMENSION(:), INTENT(IN)                 :: AtPoint
         REAL, DIMENSION(size(AtPoint),size(AtPoint))   :: Hessian
         INTEGER :: i,j

         ! Check the number of degree of freedom
         CALL ERROR( size(AtPoint) /= NDim, "PotentialModule.GetMassScaledHessian: input array dimension mismatch" )

         ! Use subroutine GetHessianFromForces to compute the matrix of the 2nd derivatives
         Hessian(:,:) = GetHessianFromForces( AtPoint, GetPotAndForces, 0.00001 ) 
         ! Numerical hessian of the potential in mass weighted coordinates
         DO j = 1, NDim
            DO i = 1, NDim
               Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
            END DO
         END DO

      END FUNCTION GetMassScaledHessian

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
      LOGICAL, DIMENSION(SIZE(X)) :: LogMask
      REAL, DIMENSION(SIZE(X)) :: Kappa1, Kappa2, Kappa3, Kappa4, Step

      REAL :: V, FMax
      INTEGER :: i

      IF ( PRESENT(Mask) ) THEN
         LogMask = Mask
      ELSE 
         LogMask = .TRUE.
      END IF

      ! compute forces at current position and the maximum component
      V = GetPotAndForces( X, Forces )
      FMax = MAXVAL( ABS(Forces), LogMask )

      ! Check if forces would result in a negligible step
      IF ( FMax < 1.E-4 )  THEN
         ! Do not take any step
         FollowGradient = .FALSE.
      ELSE
         ! Make a runge kutta 4th order step
         Kappa1 = NormalizedForces( Forces, LogMask )
         V = GetPotAndForces( X + 0.5*ABS(StepLength)*Kappa1, Forces )
         Kappa2 = NormalizedForces( Forces, LogMask )
         V = GetPotAndForces( X + 0.5*ABS(StepLength)*Kappa2, Forces )
         Kappa3 = NormalizedForces( Forces, LogMask )
         V = GetPotAndForces( X + ABS(StepLength)*Kappa3, Forces )
         Kappa4 = NormalizedForces( Forces, LogMask )
         Step = ( Kappa1/6.0 + Kappa2/3.0 + Kappa3/3.0 + Kappa4/6.0 )

         ! move uncostrained coordinates in the direction of the forces
         DO i = 1, size(X)
            IF ( LogMask(i) ) X(i) = X(i) + ABS(StepLength) * Step(i)
         END DO
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

      VectorNorm = SQRT(TheOneWithVectorDotVector( ConstrainedVector(Forces,Mask,COUNT(Mask)), &
                     ConstrainedVector(Forces,Mask,COUNT(Mask)) ))
      NormalizedForces(:) = Forces(:) / VectorNorm
      NormalizedForces(:) = NormalizedForces(:) / SQRT(MassVector(:))
   END FUNCTION NormalizedForces

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

END MODULE PotentialAnalysis
