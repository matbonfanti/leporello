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
   USE Optimize

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: PotentialAnalysis_ReadInput, PotentialAnalysis_Initialize, PotentialAnalysis_Run, PotentialAnalysis_Dispose

   !> Dimensions of the potential
   INTEGER, SAVE :: NDim

   !> Input variables for the 3D potential plot
   INTEGER :: IndexQ1, IndexQ2, IndexQ3                      !< indices of the three coords as x,y,z in the 3D plot
   REAL    :: Q1Min, Q1Max, Q2Min, Q2Max, Q3Min, Q3Max       !< boundaries of the grid
   INTEGER :: NpointQ1, NpointQ2, NpointQ3                   !< nr of points of the grid
   REAL    :: GridSpacing                                    !< grid spacing
   CHARACTER(100) :: QReferenceFile                          !< file with the reference coords for the dof's which are not scanned

   ! Parameters for optimizations
   INTEGER :: MaxOptSteps  !< Maximum number of steps for optimizations
   REAL :: OptThreshold    !< Convergence criterium for optimizations
   REAL :: DeltaFiniteDiff !< Magnitude of the displacement used to compute hessian

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
      CALL SetFieldFromInput( InputData, "Plot_IndexQ1",     IndexQ1 )
      CALL SetFieldFromInput( InputData, "Plot_IndexQ2",     IndexQ2 )
      CALL SetFieldFromInput( InputData, "Plot_IndexQ3",     IndexQ3,  0 ) ! NOTE: for IndexQ3 == 0, a 2D plot is required!
      CALL SetFieldFromInput( InputData, "Plot_Q1Min", Q1Min )
      Q1Min = Q1Min * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_Q1Max", Q1Max )
      Q1Max = Q1Max * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_Q2Min", Q2Min )
      Q2Min = Q2Min * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "Plot_Q2Max", Q2Max )
      Q2Max = Q2Max * LengthConversion(InputUnits, InternalUnits)
      IF ( IndexQ3 /= 0 ) THEN
         CALL SetFieldFromInput( InputData, "Plot_Q3Min", Q3Min )
         Q3Min = Q3Min * LengthConversion(InputUnits, InternalUnits)
         CALL SetFieldFromInput( InputData, "Plot_Q3Max", Q3Max )
         Q3Max = Q3Max * LengthConversion(InputUnits, InternalUnits)
      ENDIF
      CALL SetFieldFromInput( InputData, "Plot_QReferenceFile", QReferenceFile )

      ! Set variables for minima and TS optimization
      CALL SetFieldFromInput( InputData, "MaxOptSteps", MaxOptSteps, 10**6 )
      CALL SetFieldFromInput( InputData, "OptThreshold", OptThreshold, 1.E-6 )
      OptThreshold = OptThreshold * ForceConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "DeltaFiniteDiff", DeltaFiniteDiff, 1.E-5 )
      DeltaFiniteDiff = DeltaFiniteDiff * LengthConversion(InputUnits, InternalUnits)

      ! Set variables for Minimum Energy Path computation
      CALL SetFieldFromInput( InputData, "MaxMEPNrSteps", MaxMEPNrSteps, 1000 )
      CALL SetFieldFromInput( InputData, "MEPStep", MEPStep, 0.001 )
      MEPStep = MEPStep * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "PrintMEPNrSteps", PrintMEPNrSteps, MaxMEPNrSteps )
      CALL SetFieldFromInput( InputData, "StartMEPThreshold", StartMEPThreshold, 0.05 )
      StartMEPThreshold = StartMEPThreshold * ForceConversion(InputUnits, InternalUnits)

      ! SCREEN LOG OF THE INPUT VARIABLES

      LengthConv = LengthConversion(InternalUnits,InputUnits)

      IF ( IndexQ3 /= 0 )  WRITE(*, 900) GridSpacing*LengthConv, LengthUnit(InputUnits)
      IF ( IndexQ3 == 0 )  WRITE(*, 904) GridSpacing*LengthConv, LengthUnit(InputUnits)
      WRITE(*, 903) GetXLabel(IndexQ1), Q1Min*LengthConv, Q1Max*LengthConv, LengthUnit(InputUnits)
      WRITE(*, 903) GetXLabel(IndexQ2), Q2Min*LengthConv, Q2Max*LengthConv, LengthUnit(InputUnits)
      IF ( IndexQ3 /= 0 ) WRITE(*, 903) GetXLabel(IndexQ3), Q3Min*LengthConv, Q3Max*LengthConv, LengthUnit(InputUnits)

      WRITE(*, 901) MaxOptSteps, OptThreshold*ForceConversion(InternalUnits,InputUnits), &
                    TRIM(EnergyUnit(InputUnits))//"/"//TRIM(LengthUnit(InputUnits))

      WRITE(*, 902) StartMEPThreshold*ForceConversion(InternalUnits,InputUnits), &
                    TRIM(EnergyUnit(InputUnits))//"/"//TRIM(LengthUnit(InputUnits)), &
                    MEPStep*LengthConv, LengthUnit(InputUnits), MaxMEPNrSteps, PrintMEPNrSteps

      900 FORMAT(" * Plot of a 3D cut of the potential in VTK format ",         /,&
                 " * Grid spacing :                                ",F10.4,1X,A   )
      904 FORMAT(" * Plot of a 2D cut of the potential in VTK format ",         /,&
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

      ! Get the number of dimension of the potential
      NDim = GetSystemDimension()

      ! Allocate memory and initialize vectors for positions, forces and masses
      ALLOCATE( X(NDim), A(NDim), MassVector(NDim) )

      ! Define vector of the masses
      SELECT CASE( GetPotentialID() )
         CASE( ELEYRIDEAL_3D )
            IF ( NDim == 3 ) THEN
               MassVector(1:3) = (/ MassHInc, MassHTar, MassC /)
            ELSE IF ( NDim == 2 ) THEN
               MassVector(1:2) = (/ MassHInc, MassHTar /)
            END IF
         CASE( ELEYRIDEAL_7D )
            MassVector(1:3) = MassHInc; MassVector(4:6) = MassHTar; MassVector(7) = MassC
      END SELECT

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
      REAL, DIMENSION(:), ALLOCATABLE :: Q3Array, Q1Array, Q2Array           !< array with coordinates grid

      !> Variables for the MEP analysis
      REAL, DIMENSION(:), ALLOCATABLE    :: XStart, Dummy, XRef
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LogMask
      REAL, DIMENSION(:), ALLOCATABLE    :: EigenFreq
      REAL, DIMENSION(:,:), ALLOCATABLE  :: EigenModes, Hessian, BlockHessian
      REAL, DIMENSION(:,:), ALLOCATABLE  :: StoreMEP, StoreApproach
      INTEGER, DIMENSION(:), ALLOCATABLE :: Indices
      REAL    :: EStart, E, GradNorm
      INTEGER :: MEPFullUnit, InputGeomUnit

      !> Data for printing MEP in VTK format
      TYPE(VTKInfo), SAVE :: MEPTraj


      INTEGER :: i, j, k, nPoint, MaxStep, NPrint, NScatter
      LOGICAL :: Check

      PRINT "(2/,A)",    " ***************************************************"
      PRINT "(A,F10.5)", "               PES ANALYSIS "
      PRINT "(A,/)" ,    " ***************************************************"


      ! =========================================================
      !                (1) 3D pes in VTK format
      ! =========================================================

      PRINT "(/,A)",    " **** Write cut of the PES to output VTR ****"

      ! Read from file the value to use for the coordinates which are not scanned
      InputGeomUnit = LookForFreeUnit()
      OPEN( FILE= TRIM(ADJUSTL(QReferenceFile)), UNIT=InputGeomUnit )
      ! Read coordinates
      ALLOCATE( XRef(NDim) )
      DO i = 1, NDim
          READ(InputGeomUnit,*) XRef(i)
      END DO
      XRef = XRef*LengthConversion(InputUnits, InternalUnits)
      ! Close input file
      CLOSE( InputGeomUnit )

      ! Set grid dimensions
      NpointQ1 = INT((Q1Max-Q1Min)/GridSpacing) + 1
      NpointQ2 = INT((Q2Max-Q2Min)/GridSpacing) + 1
      NpointQ3 = 1; IF ( IndexQ3 /= 0 )  NpointQ3 = INT((Q3Max-Q3Min)/GridSpacing) + 1

      ! Allocate temporary array to store potential data and coord grids
      ALLOCATE( PotentialArray( NpointQ1 * NpointQ2 * NpointQ3 ),   &
                Q3Array( NpointQ3 ), Q1Array( NpointQ1 ), Q2Array( NpointQ2 ) )

      ! Define coordinate grids
      Q1Array = (/ ( Q1Min + GridSpacing*(i-1), i=1,NpointQ1) /)
      Q2Array = (/ ( Q2Min + GridSpacing*(i-1), i=1,NpointQ2) /)
      IF ( IndexQ3 /= 0 )  Q3Array    = (/ ( Q3Min    + GridSpacing*(i-1), i=1,NpointQ3)    /)

      ! Open VTK file
      CALL VTK_NewRectilinearSnapshot(PotentialPlot,FileName="V_Plot",X=Q1Array*LengthConversion(InternalUnits,InputUnits),&
                   Y=Q2Array*LengthConversion(InternalUnits,InputUnits), Z=Q3Array*LengthConversion(InternalUnits,InputUnits))

      nPoint = 0
      ! Cycle over the Q3 coordinate values
      DO k = 1, NpointQ3
         ! Cycle over the Q2 coordinate values
         DO j = 1, NpointQ2
            ! Cycle over the Q1 coordinate values
            DO i = 1, NpointQ1
               nPoint = nPoint + 1
               ! Define coordinates of the point
               XRef(IndexQ1) = Q1Array(i); XRef(IndexQ2) = Q2Array(j)
               IF ( IndexQ3 /= 0 ) XRef(IndexQ3) = Q3Array(k)
!                PRINT*, XRef
               ! Compute potential
               PotentialArray(nPoint) = GetPotential( XRef )
            END DO
         END DO
      END DO
      ! Print the potential to vtk file
      CALL VTK_AddScalarField( PotentialPlot,Name="Potential",Field=PotentialArray*EnergyConversion(InternalUnits,InputUnits), &
                                                                                                            LetFileOpen=.FALSE.)
      IF ( IndexQ3 /= 0  ) THEN
         WRITE(*,"(/,' * 3D PES as a func of ',A,', ',A,' and ',A,' written as VTR to file V_Plot.vtr')") &
                  TRIM(ADJUSTL(GetXLabel(IndexQ1))), TRIM(ADJUSTL(GetXLabel(IndexQ2))), TRIM(ADJUSTL(GetXLabel(IndexQ3)))
      ELSE
         WRITE(*,"(/,' * 2D PES as a func of ',A,' and ',A,' written as VTR to file V_Plot.vtr')") &
                  TRIM(ADJUSTL(GetXLabel(IndexQ1))), TRIM(ADJUSTL(GetXLabel(IndexQ2)))
      END IF

      ! Deallocate memory
      DEALLOCATE( PotentialArray, Q3Array, Q1Array, Q2Array, XRef )

      ! =========================================================
      !                 (2) asymptotic initial geometry
      ! =========================================================

      PRINT "(/,A)",    " **** Starting geometry of the MEP ****"

      ! Allocate arrays for geometry optimization
      ALLOCATE( XStart(NDim), LogMask(NDim), Dummy(NDim) )

      ! guess reasonable coordinates of the minimum of the PES
      CALL StartSystemForScattering( X, Dummy, MassVector, 100.0/MyConsts_Bohr2Ang, 0.0, 0.0, Task=2 )
      LogMask(:) =  GetInitialAsymptoteMask( )
      XStart = NewtonLocator( GetPotAndForces, X, MaxOptSteps, OptThreshold, 1.0, DeltaFiniteDiff, LogMask )
      ! Computing the energy at this geometry
      EStart = GetPotAndForces( XStart, A )

      ! Arrays for hessian diagonalization
      NScatter = GetScatterDimension( )
      ALLOCATE( Hessian(NDim, NDim), Indices(NDim-NScatter), BlockHessian(NDim-NScatter,NDim-NScatter), &
                EigenFreq(NDim-NScatter), EigenModes(NDim-NScatter,NDim-NScatter))

      ! compute the hessian in mass scaled coordinates
      Hessian(:,:) = GetMassScaledHessian( XStart )
      ! remove unbound coordinates to avoid numerical problems with singular matrix diagonalization
      Indices = GetInitialBoundIndices()
      BlockHessian = Hessian(Indices,Indices)
      ! and diagonalize it
      CALL TheOneWithDiagonalization(BlockHessian, EigenModes, EigenFreq )

      WRITE(*,"(/,A)") " Starting geometry and energy of the MEP "
      DO i = 1, NDim
         WRITE(*,501) GetXLabel(i), XStart(i)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits)
      END DO
      WRITE(*,502) EStart*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

      WRITE(*,"(/,A)") " Potential normal modes at the starting geometry"
      DO i = 1, NDim-NScatter
         IF ( EigenFreq(i) > 0.0 ) THEN
            WRITE(*,503) i, SQRT(EigenFreq(i))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                      TRIM(LengthUnit(InternalUnits)), EigenModes(:,i)
         ELSE
            WRITE(*,504) i, SQRT(-EigenFreq(i))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                      TRIM(LengthUnit(InternalUnits)), EigenModes(:,i)
         END IF
      END DO

      DEALLOCATE( Hessian, Indices, BlockHessian, EigenFreq, EigenModes )

      ! =========================================================
      !                 (3) minimum energy path
      ! =========================================================

      ! Allocate memory for MEP storage
      ALLOCATE( StoreApproach(NDim+2,0:1000), StoreMEP(NDim+2,0:MaxMEPNrSteps) )

      PRINT "(/,A,/)",    " **** Minimum energy path ****"
      WRITE(*,701)

      ! Write header lines screen
      WRITE(*,698, advance='no')
      DO i = 1, NDim
         WRITE(*,"(A10,10X)", advance='no') GetXLabel( i )
      END DO
      WRITE(*,699, advance='no') TRIM(EnergyUnit(InputUnits))
      DO i = 1, NDim
         WRITE(*,"(A20)", advance='no') TRIM(LengthUnit(InputUnits))
      END DO
      WRITE(*,700)

      ! Start from XStart
      X(:) = XStart; E = GetPotential( X )
      ! Write to output file and to screen the starting point
      WRITE(*,601) -1., E*EnergyConversion(InternalUnits, InputUnits), X(:)*LengthConversion(InternalUnits, InputUnits)
      ! Store step in the data array
      StoreApproach(:,0) = (/ 0., E, X(1:NDim) /)

      ! First follow the incident direction until significant gradient is found
      DO i = 1, 1000
         ! Move along the entrance channel
         IF ( GetSystemDimension( ) == 2 .OR. GetSystemDimension( ) == 3 ) THEN
            X(1) = X(1) - 0.05
         ELSE IF ( GetSystemDimension( ) == 7 ) THEN
            X(3) = X(3) - 0.05
         END IF
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
      WRITE(MEPFullUnit,698, advance='no')
      DO i = 1, NDim
         WRITE(MEPFullUnit,"(A10,10X)", advance='no') GetXLabel( i )
      END DO
      WRITE(MEPFullUnit,699, advance='no') TRIM(EnergyUnit(InputUnits))
      DO i = 1, NDim
         WRITE(MEPFullUnit,"(A20)", advance='no') TRIM(LengthUnit(InputUnits))
      END DO
      WRITE(MEPFullUnit,700)
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
      IF ( GetSystemDimension( ) == 2  ) THEN
         CALL VTK_WriteTrajectory(MEPTraj,StoreMEP((/3,4/),1:NPrint)*LengthConversion(InternalUnits,InputUnits),"MEP_Trajectory")
      ELSE IF ( GetSystemDimension( ) == 3 ) THEN
         CALL VTK_WriteTrajectory(MEPTraj,StoreMEP((/3,4,5/),1:NPrint)*LengthConversion(InternalUnits,InputUnits),"MEP_Trajectory")
      ELSE IF  ( GetSystemDimension( ) == 7 ) THEN
         CALL VTK_WriteTrajectory(MEPTraj,StoreMEP((/5,8,9/),1:NPrint)*LengthConversion(InternalUnits, InputUnits),"MEP_Trajectory")
      END IF

      WRITE(*,"(/,A)") " * Last step reached... MEP written to file ________"

      ! =========================================================
      !                 (4) asymptotic final geometry
      ! =========================================================

      PRINT "(/,A)",    " **** Final asymptotic geometry ****"

      ! guess reasonable coordinates of the minimum of the PES in the asymptotic out channel
      CALL StartSystemForScattering( X, Dummy, MassVector, 20.0/MyConsts_Bohr2Ang, 0.0, 0.0, Task=3 )
      LogMask(:) = .TRUE.
      ! Find minimum by Newton's optimization
      X = NewtonLocator( GetPotAndForces, X, MaxOptSteps, OptThreshold, 1.0, DeltaFiniteDiff, LogMask )
      ! Computing the energy at this geometry
      E = GetPotAndForces( X, A )

      ALLOCATE( Hessian(NDim, NDim), EigenFreq(NDim), EigenModes(NDim,NDim))

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
                  "    mass-scaled coords of the normal mode / ",A," : ",1000F12.6, / )
      504 FORMAT( " Normal Mode ",I5," - imag frequency: ",1F15.2,1X,A, /, &
                  "    mass-scaled coords of the normal mode / ",A," : ",1000F12.6, / )

      601 FORMAT (F12.4,1001F20.6 )
      698 FORMAT (   "#       Step", 10X, "    Energy", 11X )
      699 FORMAT (/, "#           ", A20 )
      700 FORMAT (/,"#------------------------------------------------------------------------------------------------------------")
      701 FORMAT (  "#------------------------------------------------------------------------------------------------------------")

      ! Deallocate memory
      DEALLOCATE( XStart, LogMask, Dummy )
      DEALLOCATE( StoreMEP, StoreApproach )
      DEALLOCATE( Hessian, EigenFreq, EigenModes )

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
         Hessian(:,:) = GetHessianFromForces( AtPoint, GetPotAndForces, DeltaFiniteDiff )
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
