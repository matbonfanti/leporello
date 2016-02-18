!***************************************************************************************
!*                              MODULE Scattering
!***************************************************************************************
!
!>  \brief     Subroutines for a scattering simulation of a atom surface collision
!>  \details   This module contains subroutines to compute averages  \n
!>             in an scattering simulation of a atom-surface system  \n
!>             at a given temperature of the surface and a given incident  \n
!>             energy of the projectile atom.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             18 February 2016
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!                
!>  \todo         N.A.
!
!***************************************************************************************
MODULE ScatteringSimulation
#include "preprocessoptions.cpp"
   USE SharedData
   USE InputField
   USE UnitConversion
   USE ClassicalEqMotion
   USE PotentialModule
   USE IndependentOscillatorsModel
   USE RandomNumberGenerator
!    USE PrintTools

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Scattering_ReadInput, Scattering_Initialize, Scattering_Run, Scattering_Dispose

   ! Initial conditions of the scatterer
   REAL    :: InitEKin                  !< Initial kinetic energy of the scattering H atom
   REAL    :: InitDistance              !< Initial Z position of the scattering H atom
!    INTEGER :: NRhoMax                   !< Nr of rho values to sample
!    REAL    :: DeltaRho                  !< Grid spacing in rho
  
   ! Initial equilibration of the slab
   REAL    :: Temperature               !< Temperature of the simulation
   REAL    :: EquilTStep                !< Time step for integrating equilibration dynamics
   REAL    :: EquilGamma                !< Friction parameter of the Langevin equation
   INTEGER :: NrEquilibSteps            !< Nr of time step of the equilibration
   INTEGER :: EquilibrationStepInterval !< Nr of tsteps between each print of equil info (=NrEquilibSteps/NrOfPrintSteps )

   ! Variables of the propagation
   INTEGER :: NrTrajs                   !< Nr of trajectories per impact parameter value
   INTEGER :: NrOfSteps                 !< Total nr of time step per trajectory
   REAL    :: TimeStep                  !< Time step for the integration of the classical EOM
   INTEGER :: NrOfPrintSteps            !< Nr of analysis steps
   INTEGER :: PrintStepInterval         !< Nr of tsteps between each print of propag info (=NrOfSteps/NrOfPrintSteps)

   ! Dimensions of the problem
   INTEGER :: NDim                      !< integer nr of dofs of the system+bath hamiltonian
   INTEGER :: NSys                      !< integer nr of dofs of the system only
   
   ! Time evolution dataset
   TYPE(Evolution),SAVE :: MolecularDynamics     !< Propagate in micro/macrocanonical ensamble to extract results
   TYPE(Evolution),SAVE :: Equilibration         !< Propagate in macrocanonical ensamble at given T to generate init conditions

!    ! Averages computed during propagation
!    REAL, DIMENSION(:), ALLOCATABLE      :: AverageZHydro      !< Average position of zH vs time
!    REAL, DIMENSION(:), ALLOCATABLE      :: AverageVHydro      !< Average velocity of zH vs time
!    REAL, DIMENSION(:), ALLOCATABLE      :: AverageZCarbon     !< Average position of zC vs time
!    REAL, DIMENSION(:), ALLOCATABLE      :: AverageVCarbon     !< Average velocity of zC vs time
!    REAL, DIMENSION(:,:), ALLOCATABLE    :: TrappingProb       !< trapping prob vs time and impact parameter

   ! Random number generator internal status
   TYPE(RNGInternalState), SAVE :: RandomNr    !< spectial type dat to store internal status of random nr generator

!    ! VTK file format
!    TYPE(VTKInfo), SAVE :: TrappingFile         !< derived datatype to print trapping probability

   CONTAINS

!===============================================================================================================================

!**************************************************************************************
!> Read from input unit the variable which are specific for a scattering
!> simulation at given T and Einc.
!>
!> @param InputData     Datatype with an already setup input unit
!**************************************************************************************
   SUBROUTINE Scattering_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

       ! PARAMETERS FOR THE INITIAL CONDITIONS OF H 
 
      ! Initial kinetic energy of the scattering H atom
      CALL SetFieldFromInput( InputData, "InitEKin", InitEKin )
      InitEKin = InitEKin * EnergyConversion(InputUnits, InternalUnits)
      ! Initial distance of the scatterer
      CALL SetFieldFromInput( InputData, "InitDistance", InitDistance )
      InitDistance = InitDistance * LengthConversion(InputUnits, InternalUnits)

!       ! Nr of rho values to sample
!       CALL SetFieldFromInput( InputData, "NRhoMax", NRhoMax )
!       ! Grid spacing in rho
!       IF ( NRhoMax > 0 ) THEN
!          CALL SetFieldFromInput( InputData, "DeltaRho", DeltaRho )
!          DeltaRho = DeltaRho * LengthConversion(InputUnits, InternalUnits)
!       ELSE  ! for collinear calculation, no deltarho is needed
!          DeltaRho = 0.0
!       ENDIF

      ! READ THE VARIABLES FOR THE TIME PROPAGATION 

      ! Nr of the trajectories of the simulation
      CALL SetFieldFromInput( InputData, "NrTrajs", NrTrajs )
      ! Timestep of the propagation
      CALL SetFieldFromInput( InputData, "TimeStep", TimeStep )
      TimeStep = TimeStep * TimeConversion(InputUnits, InternalUnits) 
      ! Nr of steps of the propagation
      CALL SetFieldFromInput( InputData, "NrOfSteps", NrOfSteps )
      ! Nr of steps in which the output is written 
      CALL SetFieldFromInput( InputData, "NrOfPrintSteps", NrOfPrintSteps )
      ! Accordingly set the interval between each printing step
      PrintStepInterval = NrOfSteps / NrOfPrintSteps

      ! READ THE VARIABLES FOR EQUILIBRATION OF THE BATH

      ! Temperature of the equilibrium simulation
      CALL SetFieldFromInput( InputData, "Temperature", Temperature, 0.0 )
      Temperature = Temperature * TemperatureConversion(InputUnits, InternalUnits)
!       ! Set temperature to zero in case of a quasiclassical simulation
!       IF ( (BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH) .AND. ZPECorrection ) Temperature = 0.0
      ! Set gamma of the equilibration Langevin dynamics
      CALL SetFieldFromInput( InputData, "EquilRelaxTime", EquilGamma)
      EquilGamma = 1. / ( EquilGamma * TimeConversion(InputUnits, InternalUnits) )
      ! Set the time step of the equilibration
      CALL SetFieldFromInput( InputData, "EquilTStep",  EquilTStep, TimeStep/MyConsts_fs2AU )
      EquilTStep = EquilTStep * TimeConversion(InputUnits, InternalUnits)
      ! Set nr of steps of the equilibration
      CALL SetFieldFromInput( InputData, "NrEquilibrSteps", NrEquilibSteps, int(10.0*(1.0/EquilGamma)/EquilTStep) )
      EquilibrationStepInterval = CEILING( real(NrEquilibSteps) / real(NrOfPrintSteps) )

      ! PRINT TO SCREEN INFORMATION ON THE INPUT VARIABLES

      IF ( .NOT. ZPECorrection ) THEN 
         WRITE(*, 900) Temperature*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits), NrTrajs
      ELSE IF ( ZPECorrection ) THEN 
         WRITE(*, 901) NrTrajs
      END IF
      
      WRITE(*, 902) InitDistance*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                    InitEKin*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

      WRITE(*, 903) TimeStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    TimeStep*NrOfSteps*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    NrOfSteps, NrOfPrintSteps

      WRITE(*, 904) 1./EquilGamma*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                 EquilTStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                 EquilTStep*NrEquilibSteps*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                 NrEquilibSteps

   900 FORMAT(" * Molecular dynamics with substrate at T = ",F5.1,1X,A,/, &
              " * Nr of trajectories:                          ",I10,/)

   901 FORMAT(" * Quasiclassical molecular dynamics with substrate at T = 0 K", &
              " * Nr of trajectories:                          ",I10,/) 

   902 FORMAT(" * Initial conditions of the atom-surface system ", /,&
              " * Initial distance of the scatterer:           ",F10.4,1X,A,/,&  
              " * Initial energy of the scatterer:             ",F10.4,1X,A,/ ) 
!               " * Nr of impact parameter values:               ",I10,       /,& 
!               " * Max value of impact parameter:               ",F10.4,1X,A,/ )

   903 FORMAT(" * Dynamical simulation variables               ",           /,&
              " * Propagation time step:                       ",F10.4,1X,A,/,& 
              " * Propagation total time:                      ",F10.4,1X,A,/,&
              " * Nr of time steps of each trajectory:         ",I10,       /,& 
              " * Nr of print steps of each trajectory:        ",I10,       / )

   904 FORMAT(" * Bath equilibration variables                 ",           /,&
              " * Relaxation time of the Langevin dynamics:    ",F10.4,1X,A,/,& 
              " * Equilibration time step:                     ",F10.4,1X,A,/,& 
              " * Equilibration total time:                    ",F10.4,1X,A,/,&
              " * Nr of equilibration steps:                   ",I10,       / )

   END SUBROUTINE Scattering_ReadInput

!===============================================================================================================================

!**************************************************************************************
!> Initialization of the data for the scattering simulation:
!> memory allocation, variable initialization and data type setup.
!>
!**************************************************************************************
   SUBROUTINE Scattering_Initialize()
      IMPLICIT NONE
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LangevinSwitchOn
      INTEGER :: iCoord

      ! Define the relevant problem dimensions
      NSys = GetSystemDimension()
      IF ( BathType == LANGEVIN_DYN ) THEN
         NDim = NSys
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         NDim = NSys + NBath
      END IF

      ! Allocate memory and initialize vectors for trajectory, acceleration and masses
      ALLOCATE( X(NDim), V(NDim), A(NDim), MassVector(NDim), LangevinSwitchOn(NDim) )

      ! Define vector of the masses
      IF ( BathType == LANGEVIN_DYN ) THEN
         MassVector = (/ MassHInc, MassHTar, MassC /)
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         MassVector = (/ MassHInc, MassHTar, MassC, (MassBath, iCoord=1,NBath) /)
      END IF

      ! Set variables for EOM integration in the microcanonical ensamble (system + bath)
      CALL EvolutionSetup( MolecularDynamics, NDim, MassVector, TimeStep )

      ! Set canonical dynamics at the end of the oscillator chain
      IF ( BathType == CHAIN_BATH .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn( NDim ) = .TRUE.
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! In case of Langevin relaxation, switch on gamma at the carbon atom
      IF ( BathType == LANGEVIN_DYN .AND. DynamicsGamma /= 0.0 ) THEN
         LangevinSwitchOn = (/ .FALSE. , .FALSE. , .TRUE. /)
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! Set variables for EOM integration with Langevin thermostat, during initial equilibration
      CALL EvolutionSetup( Equilibration, NDim, MassVector, EquilTStep )
      LangevinSwitchOn = (/ .FALSE., (.TRUE., iCoord=2,NDim ) /)
      CALL SetupThermostat( Equilibration, EquilGamma, Temperature, LangevinSwitchOn )

      DEALLOCATE( LangevinSwitchOn )

      ! Initialize random number seed
      CALL SetSeed( RandomNr, -1 )

       ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

!       ! Allocate and initialize the variables for the trajectory averages
!       ALLOCATE( AverageVHydro(0:NrOfPrintSteps), AverageZHydro(0:NrOfPrintSteps) )
!       ALLOCATE( AverageVCarbon(0:NrOfPrintSteps), AverageZCarbon(0:NrOfPrintSteps) )
! 
!       ! Allocate and initialize trapping probability
!       ALLOCATE( TrappingProb(0:NRhoMax,NrOfPrintSteps) )
!       TrappingProb = 0.0
! 
!       ! Check that a non collinear V is considered when running a non collinear simulation
!       IF (NRhoMax > 0) THEN
!          CALL ERROR( Collinear, " Scattering_Initialize: a non collinear potential is needed!" )
!       END IF 
! 
!       ! if XYZ files of the trajectories are required, allocate memory to store the traj
!       IF ( PrintType >= FULL  .AND. ( RunType == SCATTERING .OR. RunType == EQUILIBRIUM ) ) THEN
!            ALLOCATE( Trajectory( 16, ntime ) )
!       END IF

   END SUBROUTINE Scattering_Initialize

!===============================================================================================================================

!**************************************************************************************
!> Run the scattering simulation and compute trapping probability over time.
!>
!**************************************************************************************
   SUBROUTINE Scattering_Run()
      IMPLICIT NONE
!       INTEGER  ::  AvHydroOutputUnit, AvCarbonOutputUnit       ! UNITs FOR OUTPUT AND DEBUG
!       INTEGER  ::  CollinearTrapUnit, CrossSectionUnit, OpacityUnit
!       INTEGER  ::  jRho, iTraj, iStep, kStep, iCoord
!       INTEGER  ::  DebugUnitEn, DebugUnitCoord, DebugUnitVel
!       REAL     ::  ImpactPar, Time, CrossSection
!       REAL     ::  RndValue, CarbonFreq
!       REAL, DIMENSION(1,8)  :: AsymptoticCH
!       REAL     ::  TotEnergy, PotEnergy, KinEnergy, KinHydro, KinCarbon
!       REAL     ::  TempAverage, TempVariance, IstTemperature
!       REAL, DIMENSION(NRhoMax+1) :: ImpactParameterGrid
!       REAL, DIMENSION(NrOfPrintSteps) :: TimeGrid
!       CHARACTER(100) :: OutFileName
! 
!       IF ( PrintType >= FULL ) THEN
!          AvHydroOutputUnit = LookForFreeUnit()
!          OPEN( FILE="AverageHydro.dat", UNIT=AvHydroOutputUnit )
!          WRITE(AvHydroOutputUnit, "(A,I6,A,/)") "# average zH coord and vel vs time (fs | ang) - ",NrTrajs, " trajs"
! 
!          AvCarbonOutputUnit = LookForFreeUnit()
!          OPEN( FILE="AverageCarbon.dat", UNIT=AvCarbonOutputUnit )
!          WRITE(AvCarbonOutputUnit, "(A,I6,A,/)") "# average zC coord and vel vs time (fs | ang) - ",NrTrajs," trajs"
!       ENDIF
! 
!       PRINT "(2/,A)",    "***************************************************"
!       PRINT "(A,F10.5)", "         H-GRAPHITE SCATTERING SIMULATION"
!       PRINT "(A,/)" ,    "***************************************************"
! 
!       PRINT "(A,I5,A,I5,A)"," Running ", NrTrajs, " trajectories per ",NRhoMax+1," impact parameters ... "
! 
!       ! Initialize variables to print average results at a given impact parameter value
!       AverageVCarbon(:)         = 0.0
!       AverageZCarbon(:)         = 0.0
!       AverageVHydro(:)          = 0.0
!       AverageZHydro(:)          = 0.0
! 
!       ! scan over the impact parameter... 
!       ImpactParameter: DO jRho = 0,NRhoMax
! 
!          ! Set impact parameter and print message
!          ImpactPar = float(jRho)*DeltaRho
! 
!          IF ( NRhoMax == 0 ) THEN 
!             PRINT "(/,A,I5,A)"," Collinear scattering simulation ... "
!          ELSE
!             PRINT "(2/,A)",    "***************************************************"
!             PRINT "(A,F10.5)", "       IMPACT PARAMETER = ", ImpactPar*LengthConversion(InternalUnits,InputUnits)
!             PRINT "(A,/)" ,    "***************************************************"
!          END IF
! 
!          PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "
! 
!          ! run NTrajs number of trajectories at the current impact parameter
!          Trajectory: DO iTraj = 1, NrTrajs
! 
!             PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 
! 
!             !*************************************************************
!             ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
!             !*************************************************************
! 
!             ! Put zero point energy in each oscillator of the bath and in the carbon
!             IF ( (BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH) .AND. ZPECorrection ) THEN
! 
!                ! Set initial coordinates of H atom
!                X(1) = ImpactPar
!                X(2) = 0.0
!                X(3) = ZHInit
! 
!                ! Set initial velocities of H atom
!                V(1) = 0.0
!                V(2) = 0.0
!                V(3) = -sqrt( 2.0 * EKinZH / MassH )
!  
!                ! Set initial coord and vel of C atom
!                RndValue = UniformRandomNr( RandomNr )*2.0*MyConsts_PI
!                CarbonFreq = SQRT( CarbonForceConstant() / MassC )
!                X(4) = COS(RndValue)/SQRT(MassC*CarbonFreq)
!                V(4) = SIN(RndValue)*SQRT(CarbonFreq/MassC)
! 
!                ! Set initial conditions of the bath
!                CALL ZeroKelvinBathConditions( Bath, X(5:), V(5:), ZPECorrection, RandomNr )
! 
!                ! Compute starting potential and forces
!                A(:) = 0.0
!                PotEnergy = ScatteringPotential( X, A )
!                A(:) = A(:) / MassVector(:)
! 
!                ! compute kinetic energy and total energy
!                KinEnergy = EOM_KineticEnergy(Equilibration, V )
!                TotEnergy = PotEnergy + KinEnergy
!                IstTemperature = 2.0*KinEnergy/size(X)
! 
!             ELSE  ! Thermalize bath + carbon
! 
!                ! Set initial coordinates of C and H atoms
!                AsymptoticCH(1,1) = ImpactPar
!                AsymptoticCH(1,2) = 0.0
!                AsymptoticCH(1,3) = ZHInit
!                AsymptoticCH(1,4) = 0.0
! 
!                ! Set initial velocities of C and H atoms
!                AsymptoticCH(1,5) = 0.0
!                AsymptoticCH(1,6) = 0.0
!                AsymptoticCH(1,7) = -sqrt( 2.0 * EKinZH / MassH )
!                AsymptoticCH(1,8) = 0.0
! 
!                ! Set initial conditions
!                IF ( BathType == SLAB_POTENTIAL ) THEN 
!                   CALL ZeroKelvinSlabConditions( X, V, AsymptoticCH, RandomNr ) 
!                ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
!                   CALL ThermalEquilibriumBathConditions( Bath, X(5:), V(5:), Temperature, RandomNr )
!                ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
!                   CALL ThermalEquilibriumBathConditions( DblBath(1), X(5:NBath+4), V(5:NBath+4), Temperature, RandomNr )
!                   CALL ThermalEquilibriumBathConditions( DblBath(2), X(NBath+5:2*NBath+4), V(NBath+5:2*NBath+4), &
!                                                                               Temperature, RandomNr )
!                ELSE IF ( BathType == LANGEVIN_DYN ) THEN
!                   ! nothing to do
!                END IF
! 
!                ! During equilibration, fix H in asymptotic position with null velocity
!                ! and initialize carbon atom in the equilibrium position with null velocity
!                X(1:4) = AsymptoticCH( 1, 1:4 )
!                V(1:4) = 0.0
! 
!                IF ( BathType /= LANGEVIN_DYN ) THEN 
! 
!                   PRINT "(/,A,F6.1,1X,A)",   " Equilibrating the initial conditions at T = ", &
!                       Temperature*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)
! 
!                   ! Initialize temperature average and variance
!                   TempAverage = 0.0
!                   TempVariance = 0.0
! 
!                   ! Compute starting potential and forces
!                   A(:) = 0.0
!                   PotEnergy = ScatteringPotential( X, A )
!                   A(:) = A(:) / MassVector(:)
! 
!                   ! compute kinetic energy and total energy
!                   KinEnergy = EOM_KineticEnergy(Equilibration, V )
!                   TotEnergy = PotEnergy + KinEnergy
!                   IstTemperature = 2.0*KinEnergy/size(X)
! 
!                   ! Do an equilibration run
!                   EquilibrationCycle: DO iStep = 1, NrEquilibSteps
! 
!                      ! PROPAGATION for ONE TIME STEP
!                      CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, ScatteringPotential, PotEnergy, RandomNr )
! 
!                      ! compute kinetic energy and total energy
!                      KinEnergy = EOM_KineticEnergy(Equilibration, V )
!                      TotEnergy = PotEnergy + KinEnergy
!                      IstTemperature = 2.0*KinEnergy/size(X)
! 
!                      ! store temperature average and variance
!                      TempAverage = TempAverage + IstTemperature
!                      TempVariance = TempVariance + IstTemperature**2
! 
!                   END DO EquilibrationCycle
! 
!                   PRINT "(A)", " Equilibration completed! "
! 
!                   ! Compute average and standard deviation
!                   TempAverage = TempAverage / NrEquilibSteps 
!                   TempVariance = (TempVariance/NrEquilibSteps) - TempAverage**2
!                   ! output message with average values
!                   WRITE(*,500)  TempAverage*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits), &
!                                 sqrt(TempVariance)*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)
! 
!                END IF 
! 
!                ! for the slab potential translate C2-C3-C4 plane at z=0
!                IF ( BathType == SLAB_POTENTIAL )   X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
!                ! give right initial conditions to the H atom
!                X(1:3) = AsymptoticCH( 1, 1:3 )
!                V(1:3) = AsymptoticCH( 1, 5:7 )
! 
!             END IF
! 
! 
! !           !*************************************************************
! !           ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
! !           !*************************************************************
! 
!             ! Compute kinetic energy and total energy
!             KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
!             KinHydro  = EOM_KineticEnergy( MolecularDynamics, V, 3 )
!             KinCarbon = KinEnergy - KinHydro
!             TotEnergy = PotEnergy + KinEnergy
! 
!             ! PRINT INITIAL CONDITIONS of THE TRAJECTORY
!             WRITE(*,600)  PotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
!                           KinHydro*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),  &
!                           KinCarbon*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
!                           KinEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
!                           TotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)
! 
!             ! Increment averages at starting conditions
!             AverageVCarbon(0) = AverageVCarbon(0) + V(4)
!             AverageZCarbon(0) = AverageZCarbon(0) + X(4)
!             AverageVHydro(0)  = AverageVHydro(0)  + V(3)
!             AverageZHydro(0)  = AverageZHydro(0)  + X(3)
! 
!             ! Open unit for massive output, with detailed info on trajectories
!             IF ( PrintType == DEBUG ) THEN
!                WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Energy.dat"
!                DebugUnitEn = LookForFreeUnit()
!                OPEN( Unit=DebugUnitEn, File=OutFileName )
! 
!                WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Coord.dat"
!                DebugUnitCoord = LookForFreeUnit()
!                OPEN( Unit=DebugUnitCoord, File=OutFileName )
! 
!                WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Vel.dat"
!                DebugUnitVel = LookForFreeUnit()
!                OPEN( Unit=DebugUnitVel, File=OutFileName )
! 
!                ! Write initial values
!                WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | V, KH, KC, Kin, E / Eh "
!                WRITE(DebugUnitEn,800) 0.0,  PotEnergy, KinHydro, KinCarbon, KinEnergy, TotEnergy
! 
!                WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
!                WRITE(DebugUnitCoord,800) 0.0, X(1), X(2), X(3), X(4), (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, NDim-4) /)
! 
!                WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
!                WRITE(DebugUnitVel,800) 0.0, V(:)
!              ENDIF
! 
!             !*************************************************************
!             !         TIME EVOLUTION OF THE TRAJECTORY
!             !*************************************************************
!     
!             ! initialize counter for printing steps
!             kStep = 0
! 
!             PRINT "(/,A)", " Propagating the H-Graphene system in time... "
! 
!             ! cycle over nstep velocity verlet iterations
!             Propagation: DO iStep = 1, NrOfSteps
! 
!                ! Propagate for one timestep
!                CALL EOM_LangevinSecondOrder( MolecularDynamics, X, V, A, ScatteringPotential, PotEnergy, RandomNr )
! 
!                ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
!                IF ( BathType == SLAB_POTENTIAL ) THEN 
!                   X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
!                ENDIF
! 
!                ! output to write every nprint steps 
!                IF ( mod(iStep,PrintStepInterval) == 0 ) THEN
! 
!                   ! increment counter for printing steps
!                   kStep = kStep+1
!                   IF ( kStep > NrOfPrintSteps ) CYCLE 
! 
!                   ! increment ch coords averages
!                   AverageVCarbon(kStep) = AverageVCarbon(kStep) + V(4)
!                   AverageZCarbon(kStep) = AverageZCarbon(kStep) + X(4)
!                   AverageVHydro(kStep)  = AverageVHydro(kStep)  + V(3)
!                   AverageZHydro(kStep)  = AverageZHydro(kStep)  + X(3)
!                   
!                   ! Compute kinetic energy and total energy
!                   KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
!                   KinHydro  = EOM_KineticEnergy( MolecularDynamics, V, 3 )
!                   KinCarbon = KinEnergy - KinHydro
!                   TotEnergy = PotEnergy + KinEnergy
! 
!                   ! check if H is in the trapping region
!                   IF ( X(3) <= AdsorpLimit ) TrappingProb(jRho,kStep) = TrappingProb(jRho,kStep)+1.0
! 
!    !                ! Store the trajectory for XYZ printing
!    !                IF ( PrintType >= FULL  .AND. RunType == EQUILIBRIUM ) THEN
!    !                      Trajectory( :, kstep ) = 0.0
!    !                      Trajectory( 1:min(16,3+nevo) , kstep ) = X( 1:min(16,3+nevo) ) 
!    !                      NrOfTrajSteps = kstep
!    !                END IF
! 
!                   ! If massive level of output, print traj information to std out
!                   IF ( PrintType == DEBUG ) THEN
!                      WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, &
!                                            PotEnergy, KinHydro, KinCarbon, KinEnergy, TotEnergy
!                      WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(1:4), &
!                                              (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, NDim-4) /)
!                      WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
!                   END IF
! 
!                END IF 
!     
!             END DO Propagation
! 
! !             ! At the end of the propagation, write the xyz file of the trajectory, if requested
! !             IF ( PrintType >= FULL ) THEN
! ! 
! !                ! print the full trajectory as output file
! !                WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Traj_",jRho,"_",iTraj,".xyz"
! !                CALL WriteTrajectoryXYZ( Trajectory(:,1:NrOfTrajSteps)*MyConsts_Bohr2Ang, OutFileName, &
! !                                              GraphiteLatticeConstant()*MyConsts_Bohr2Ang )
! ! 
! 
!             PRINT "(A)", " Time propagation completed! "
! 
!             ! print the final values of this trajectory 
!             WRITE(*,700) X(3)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
!                          X(4)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
!                          TotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)
! 
!             IF ( PrintType == DEBUG ) THEN
!                   CLOSE( Unit=DebugUnitEn )
!                   CLOSE( Unit=DebugUnitCoord )
!                   CLOSE( Unit=DebugUnitVel )
!             ENDIF
! 
!          END DO Trajectory
! 
!          PRINT "(A)"," Done! "                                    
! 
!          ! Normalize coordinates averages 
!          AverageVCarbon(:) = AverageVCarbon(:)/NrTrajs
!          AverageZCarbon(:) = AverageZCarbon(:)/NrTrajs 
!          AverageVHydro(:)  = AverageVHydro(:)/NrTrajs
!          AverageZHydro(:)  = AverageZHydro(:)/NrTrajs
! 
!          IF ( PrintType >= FULL ) THEN
!             ! write impact parameter header in output files
!             WRITE(AvHydroOutputUnit,"(/,A,F10.5)")  &
!                        "# Impact parameter : ",ImpactPar*LengthConversion(InternalUnits,InputUnits)
!             WRITE(AvCarbonOutputUnit,"(/,A,F10.5)")  &
!                        "# Impact parameter : ",ImpactPar*LengthConversion(InternalUnits,InputUnits)
! 
!             ! Print time resolved data to output files
!             DO iStep = 1, NrOfPrintSteps
!                Time = FLOAT(PrintStepInterval*iStep)*TimeStep*TimeConversion(InternalUnits,InputUnits)
!                WRITE(AvHydroOutputUnit,800) Time, AverageZHydro(iStep)*LengthConversion(InternalUnits,InputUnits),  &
!                                                 AverageVHydro(iStep)*VelocityConversion(InternalUnits,InputUnits)
!                WRITE(AvCarbonOutputUnit,800) Time, AverageZCarbon(iStep)*LengthConversion(InternalUnits,InputUnits), &     
!                                                 AverageVCarbon(iStep)*VelocityConversion(InternalUnits,InputUnits)
!             END DO
!          END IF
! 
!       END DO ImpactParameter
! 
!       ! Normalize trapping probability 
!       TrappingProb(:,:) = TrappingProb(:,:) / NrTrajs
! 
!       ! Prepare time grid and impact parameter grid
!       TimeGrid = (/( float(PrintStepInterval*iStep)*TimeStep*TimeConversion(InternalUnits,InputUnits), &
!                                       iStep= 1,NrOfPrintSteps  )/)
!       ImpactParameterGrid = (/( float(jRho)*DeltaRho*LengthConversion(InternalUnits,InputUnits), jRho = 0,NRhoMax  )/)
! 
!       IF ( NRhoMax > 0 ) THEN
!          ! Print trapping to vtk file
!          CALL VTK_NewRectilinearSnapshot( TrappingFile, X=ImpactParameterGrid, Y=TimeGrid, FileName="Trapping" )
!          CALL VTK_AddScalarField( TrappingFile, "trapping_p", RESHAPE( TrappingProb(:,:), (/NrOfPrintSteps*(NRhoMax+1) /))  )
!       END IF
! 
!       ! Print to file collinear trapping
!       CollinearTrapUnit = LookForFreeUnit()
!       OPEN( FILE="CollinearTrapping.dat", UNIT=CollinearTrapUnit )
!       WRITE(CollinearTrapUnit, "(A,I6,A,/)") "# collinear trapping probability (" // TRIM(TimeUnit(InputUnits)) // &
!           " | adim) -",NrTrajs, " trajs"
!       DO iStep= 1,NrOfPrintSteps
!          WRITE(CollinearTrapUnit,800) TimeGrid(iStep), TrappingProb(0,iStep)
!       END DO
! 
!       ! Close open units
!       CLOSE( CollinearTrapUnit )
!       IF ( PrintType >= FULL ) CLOSE( AvCarbonOutputUnit )
!       IF ( PrintType >= FULL ) CLOSE( AvHydroOutputUnit )
! 
!       ! write cross section data 
!       IF ( NRhoMax > 0 ) THEN
! 
!          ! open file and write header
!          CrossSectionUnit = LookForFreeUnit()
!          OPEN( FILE="CrossSection.dat", UNIT=CrossSectionUnit )
!          WRITE(CrossSectionUnit, "(A,I6,A,/)") "# cross section vs time (" // TRIM(TimeUnit(InputUnits)) // &
!              " | " // TRIM(LengthUnit(InputUnits)) //"^2) -",NrTrajs, " trajs"
! 
!          ! loop over time steps, compute and write trapping cross section
!          DO iStep= 1,NrOfPrintSteps
! 
!             ! integrate over impact parameter
!             CrossSection=0.0
!             DO jRho = 0, NRhoMax
!                CrossSection = CrossSection + TrappingProb(jRho,iStep) * float(jRho)*DeltaRho
!             END DO
!             CrossSection = 2.0*MyConsts_PI * DeltaRho * CrossSection * LengthConversion(InternalUnits,InputUnits)**2
! 
!             ! print trapping cross section to file
!             WRITE(CrossSectionUnit,800) TimeGrid(iStep), CrossSection
! 
!          END DO
! 
!          ! Print to file final opacity function
!          OpacityUnit = LookForFreeUnit()
!          OPEN( FILE="OpacityFunction.dat", UNIT=OpacityUnit )
!          WRITE(OpacityUnit, "(A,F8.2,A,/)") "# opacity function @ time ", &
!             float(PrintStepInterval*iStep)*TimeStep*TimeConversion(InternalUnits,InputUnits), " "//TRIM(TimeUnit(InputUnits))
!          DO jRho = 0,NRhoMax
!            WRITE(OpacityUnit,800) ImpactParameterGrid(jRho+1), TrappingProb(jRho,NrOfPrintSteps) 
!          END DO
!          CLOSE(OpacityUnit)
! 
!       END IF
! 
!    500 FORMAT (/, " Equilibration averages                 ",/     &
!                   " * Average temperature          ",1F10.4,1X,A,/ &
!                   "   Standard deviation           ",1F10.4,1X,A,/ ) 
! 
!    600 FORMAT (/, " Initial condition of the MD trajectory ",/     &
!                   " * Potential Energy             ",1F10.4,1X,A,/ &
!                   " * Kinetic Energy of H atom     ",1F10.4,1X,A,/ &
!                   " * Kinetic Energy of C atoms    ",1F10.4,1X,A,/ &
!                   " * Total Kinetic Energy         ",1F10.4,1X,A,/ &
!                   " * Total Energy                 ",1F10.4,1X,A,/ ) 
! 
!    700 FORMAT (/, " * Final zH coord               ",1F10.4,1X,A,/ &
!                   " * Final zC coord               ",1F10.4,1X,A,/ &
!                   " * Final energy                 ",1F10.4,1X,A,/ ) 
! 
!    800 FORMAT(F12.5,1000F15.8)
! 
   END SUBROUTINE Scattering_Run

!===============================================================================================================================

!**************************************************************************************
!> Free the memory which has been used for the simulation.
!>
!**************************************************************************************
   SUBROUTINE Scattering_Dispose()
      IMPLICIT NONE

      ! Deallocate memory
      DEALLOCATE( X, V, A, MassVector )
!       DEALLOCATE( AverageVHydro, AverageZHydro, AverageVCarbon, AverageZCarbon, TrappingProb )

      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( Equilibration )

   END SUBROUTINE Scattering_Dispose

!===============================================================================================================================

!    REAL FUNCTION ScatteringPotential( Positions, Forces )
!       REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
!       REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
!       INTEGER :: NrDOF, i
!       REAL :: CouplingFs, DTimesX
! 
!       ! Check the number degrees of freedom
!       NrDOF = size( Positions )
!       CALL ERROR( size(Forces) /= NrDOF, "ScatteringSimulation.ScatteringPotential: array dimension mismatch" )
!       CALL ERROR( NrDOF /= NDim, "ScatteringSimulation.ScatteringPotential: wrong number of DoFs" )
! 
!       ! Initialize forces and potential
!       ScatteringPotential = 0.0
!       Forces(:)          = 0.0
! 
!       IF ( BathType == SLAB_POTENTIAL ) THEN 
!          ! Compute potential using the potential subroutine
!          ScatteringPotential = VHSticking( Positions, Forces )
! 
!       ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
!          ! Compute potential and forces of the system
!          ScatteringPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
!          ! Coupling function of the system
! #if defined(__ZH_DEPENDENT_COUPLING)
!          CouplingFs = Positions(4)-GetSpline1D(ZCofZHSpline, Positions(3))
! #else
!          CouplingFs = Positions(4)
! #endif
!          ! Add potential and forces of the bath and the coupling
!          CALL BathPotentialAndForces( Bath, CouplingFs, Positions(5:), &
!                                                 ScatteringPotential, Forces(4), Forces(5:), DTimesEffMode=DTimesX ) 
! #if defined(__ZH_DEPENDENT_COUPLING)
!          ! Add forces on ZH coming from ZH-dependent distortion correction
!          Forces(3) = Forces(3) + GetDistorsionForce(Bath)*CouplingFs*DerivSpline(ZCofZHSpline, Positions(3))
!          Forces(3) = Forces(3) - DTimesX * DerivSpline(ZCofZHSpline, Positions(3))
! #endif
! 
!       ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
!          ! Compute potential and forces of the system
!          ScatteringPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
!          ! Add potential and forces of the bath and the coupling
!          CALL BathPotentialAndForces( DblBath(1), Positions(4)-C1Puckering, Positions(5:NBath+4), ScatteringPotential, &
!                                                                                Forces(4), Forces(5:NBath+4) ) 
!          CALL BathPotentialAndForces( DblBath(2), Positions(4)-C1Puckering, Positions(NBath+5:2*NBath+4), ScatteringPotential, &
!                                                                                Forces(4), Forces(NBath+5:2*NBath+4) ) 
! 
!       ELSE IF ( BathType == LANGEVIN_DYN ) THEN
!          ! Compute potential and forces of the system
!          ScatteringPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
! 
!       END IF
! 
!    END FUNCTION ScatteringPotential
! 
!===============================================================================================================================

END MODULE ScatteringSimulation
