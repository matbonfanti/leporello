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
!>  \todo         implement general output for a non-collinear calculation
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
   USE PrintTools
   USE FiniteDifference
   USE Optimize

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Scattering_ReadInput, Scattering_Initialize, Scattering_Run, Scattering_Dispose

   ! Initial conditions of the scatterer
   REAL    :: InitEKin                  !< Initial kinetic energy of the scattering H atom
   REAL    :: InitDistance              !< Initial Z position of the scattering H atom
   INTEGER :: NRhoMax                   !< Nr of rho values to sample
   REAL    :: DeltaRho                  !< Grid spacing in rho
  
   ! Initial equilibration / initial conditions of the slab
   REAL    :: Temperature               !< Temperature of the simulation
   LOGICAL :: ZPECorrection             !< ZeroPointEnergy correction in the initial conditions of the bath (at 0 K)
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
   INTEGER :: NCoupled                  !< integer index of the dof which is coupled to the bath

   ! Time evolution dataset
   TYPE(Evolution),SAVE :: MolecularDynamics     !< Propagate in micro/macrocanonical ensamble to extract results
   TYPE(Evolution),SAVE :: Equilibration         !< Propagate in macrocanonical ensamble at given T to generate init conditions

   ! Averages computed during propagation
   INTEGER  :: NrEnAverages                                     !< nr of average values which are computed and stored
   REAL, DIMENSION(:,:,:), ALLOCATABLE    :: TrajOutcome        !< prob of the different channels vs time and impact parameter
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: EnergyAver         !< energy average values

   ! Random number generator internal status
   TYPE(RNGInternalState), SAVE :: RandomNr    !< spectial type dat to store internal status of random nr generator

   ! Arrays used for defining initial conditions of the trajectories  
   REAL, ALLOCATABLE, DIMENSION(:)     ::  XEquil             !< Optimized geometry at the asymptotic init channel
   REAL, ALLOCATABLE, DIMENSION(:,:)   ::  Hessian            !< Hessian matrix at XEquil
   REAL, ALLOCATABLE, DIMENSION(:,:)   ::  NormalModesVec     !< Normal modes vectors at XEquil
   REAL, ALLOCATABLE, DIMENSION(:)     ::  NormalModesVal     !< Normal modes squared frequencies at XEquil

  
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

      ! Nr of rho values to sample
      CALL SetFieldFromInput( InputData, "NRhoMax", NRhoMax, 0 )
      ! Grid spacing in rho
      IF ( NRhoMax > 0 ) THEN
         CALL SetFieldFromInput( InputData, "DeltaRho", DeltaRho )
         DeltaRho = DeltaRho * LengthConversion(InputUnits, InternalUnits)
      ELSE  ! for collinear calculation, no deltarho is needed
         DeltaRho = 0.0
      ENDIF

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

      ! READ THE VARIABLES FOR INITIAL CONDITIONS/EQUILIBRATION OF THE SUBSTRATE

      ! Temperature of the equilibrium simulation
      CALL SetFieldFromInput( InputData, "Temperature", Temperature, 0.0 )
      Temperature = Temperature * TemperatureConversion(InputUnits, InternalUnits)
      ! Quasi-classical correction of the initial conditions of the bath (ZPE), relevant only for 0K
      CALL SetFieldFromInput( InputData, "ZPECorrection", ZPECorrection, .FALSE. )
      ! Set temperature to zero in case of a quasiclassical simulation
      IF ( (BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH) .AND. ZPECorrection ) Temperature = 0.0

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
      NCoupled = NSys

      ! Allocate memory and initialize vectors for trajectory, acceleration and masses
      ALLOCATE( X(NDim), V(NDim), A(NDim), MassVector(NDim), LangevinSwitchOn(NDim) )

      ! Define vector of the masses
      IF ( BathType == LANGEVIN_DYN ) THEN
         IF ( NSys == 3 ) THEN 
            MassVector = (/ MassHInc, MassHTar, MassC /)
         ELSE IF ( NSys == 2 ) THEN 
            MassVector = (/ MassHInc, MassHTar /)
         END IF
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
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn(NCoupled) = .TRUE.
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

      ! Define the number of average values to be stored
      NrEnAverages = NSys+1+3+1+3  ! NSys+1 kin en, 3 pot en and 3 pot partitions
      
      ! Allocate and initialize the variables for the trajectory averages
      IF ( PrintType >= FULL ) THEN 
         ALLOCATE( EnergyAver(NrEnAverages,0:GetNrChannels()-1,0:NRhoMax,NrOfPrintSteps) )
         EnergyAver = 0.0
      END IF
      
      ! Allocate and initialize channel resolved probability array
      ALLOCATE( TrajOutcome(0:GetNrChannels()-1,0:NRhoMax,NrOfPrintSteps) )
      TrajOutcome = 0.0

      ! Check that a non collinear V is considered when running a non collinear simulation
      IF (NRhoMax > 0) THEN
         CALL ERROR( PESIsCollinear(), " Scattering_Initialize: a non collinear potential is needed!" )
      END IF 

      ! Allocate memory for geometry optimization and normal mode analysis
      ALLOCATE( Hessian(NDim,NDim), NormalModesVec(NDim,NDim), NormalModesVal(NDim), XEquil(NDim) )
      
   END SUBROUTINE Scattering_Initialize

!===============================================================================================================================

!**************************************************************************************
!> Run the scattering simulation and compute trapping probability over time.
!>
!**************************************************************************************
   SUBROUTINE Scattering_Run()
      IMPLICIT NONE
      ! Output units
      INTEGER :: TEquilUnit, DebugUnitEn, DebugUnitCoord, DebugUnitVel       ! DEBUG
      INTEGER  ::  CollinearProbUnit, EnergyAverUnit, NrmlMdsUnit!, CrossSectionUnit, OpacityUnit
      CHARACTER(100) :: OutFileName                                          ! string to define output file name
      ! Integer indices
      INTEGER  ::  jRho, iTraj, iStep, kStep, iCoord, iChan
      ! Real variables 
      REAL     ::  ImpactPar, Time, CrossSection
      REAL     ::  TotEnergy, PotEnergy, KinEnergy, KinScatter, KinSubstrate  ! Energy values
      REAL     ::  TempAverage, TempVariance, IstTemperature                  ! Temperature values
      REAL, DIMENSION(NrEnAverages) :: EnergyExpect                           ! Array to store energy values
      ! arrays to store grids useful for results output
      REAL, DIMENSION(NRhoMax+1) :: ImpactParameterGrid
      REAL, DIMENSION(NrOfPrintSteps) :: TimeGrid
      ! VTK file format
      TYPE(VTKInfo), SAVE :: PrintChannelP                                    !< derived datatype to print resolved probability
      ! temporary subarrays for diagonalization
      INTEGER, DIMENSION(NDim-1) :: Indices
      REAL, DIMENSION(NDim-1,NDim-1) :: BlockHessian, BlockEigenVec
      REAL, DIMENSION(NDim-1) :: BlockEigenVal
      REAL, DIMENSION(NDim) :: FinalEigenVal
      
      
      IF ( PrintType == DEBUG .AND. .NOT. ZPECorrection ) THEN
         ! Open output file to print the brownian realizations of ZH vs time
         TEquilUnit = LookForFreeUnit()
         OPEN( FILE="EquilibrationTemp.dat", UNIT=TEquilUnit )
         WRITE(TEquilUnit, "(A,I6,A,/)") "# ", NrTrajs, " temperature average and variance for each equilibration (fs | K)"
      END IF

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "            SCATTERING SIMULATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A,I5,A)"," Running ", NrTrajs, " trajectories per ",NRhoMax+1," impact parameters ... "

      ! Compute the normal modes at the minimum of the asymptotic out configuration
      
      ! Put system in the final scattering equilibrium geometry...
      CALL StartSystemForScattering( X(1:NSys), V(1:NSys), MassVector(1:NSys), 2000., 0.0, 0.0, Task=3 )
      ! ... and the rest of the degrees of freedom in a good guess for initial minimum
      IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         X(NSys+1:) = 0.0
      END IF
      
      ! Optimize the potential and find the position of the local minimum 
      XEquil = NewtonLocator( ScatteringPotential, X, 10**3, 1.E-6, 1.E-6, 1.E-3 )

      ! compute the hessian in mass scaled coordinates
      Hessian(:,:) = GetMassScaledHessian( XEquil ) 
      ! and diagonalize it
      Indices = (/ GetInitialBoundIndices(), ( NSys+iCoord, iCoord=1,NBath ) /)  
      BlockHessian = Hessian(Indices,Indices)
      CALL TheOneWithDiagonalization(Hessian, NormalModesVec, FinalEigenVal )

      ! Compute the normal modes at the minimum of the asymptotic incoming configuration
      
      ! Put initial scatterer asymptotically...
      CALL StartSystemForScattering( X(1:NSys), V(1:NSys), MassVector(1:NSys), 2000., 0.0, 0.0, Task=2 )
      ! ... and the rest of the degrees of freedom in a good guess for initial minimum
      IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         X(NSys+1:) = 0.0
      END IF

      ! Optimize the potential and find the position of the local minimum 
      XEquil = NewtonLocator( ScatteringPotential, X, 10**3, 1.E-6, 1.E-6, 1.E-3, &
             (/ GetInitialAsymptoteMask(), (.TRUE., iCoord=1,NBath ) /) )
      
      ! compute the hessian in mass scaled coordinates
      Hessian(:,:) = GetMassScaledHessian( XEquil ) 
      ! and diagonalize it
      Indices = (/ GetInitialBoundIndices(), ( NSys+iCoord, iCoord=1,NBath ) /)  
      BlockHessian = Hessian(Indices,Indices)
      CALL TheOneWithDiagonalization(BlockHessian, BlockEigenVec, BlockEigenVal )
      NormalModesVec = 0.0; NormalModesVec(1,1) = 1.0; NormalModesVal(1) = 0.0
      NormalModesVec(Indices,Indices) = BlockEigenVec
      NormalModesVal(Indices) = BlockEigenVal
      
      IF ( .NOT. ZPECorrection ) THEN 
         PRINT "(/,A)" , " Initial conditions for the substrate will be sampled in the normal modes frame of reference"
         PRINT "(A)" ,   " with a classical Boltzmann distribution of harmonic oscillators. "
         PRINT "(A,/)" , " Initial sampling will be then improved with a Langevin equilibration run. " 
      ELSE IF ( ZPECorrection ) THEN
         PRINT "(/,A)" , " Initial conditions for the substrate will be sampled in the normal modes frame of reference"
         PRINT "(A)" ,   " with microcanonical distribution of coordinates and momenta for a harmonic oscillator"
         PRINT "(A,/)" , " in the classical orbit of its zero point energy"
      END IF
        
      PotEnergy = ScatteringPotential( XEquil, A )
      PRINT "(A,F10.4,A)", " Asymptotic minimum of the potential has V equal to ",& 
                  PotEnergy*EnergyConversion(InternalUnits, InputUnits)," "//TRIM(EnergyUnit(InputUnits))
 
      ! Open output file and print the normal modes of the asymptotic initial geometry
      IF ( PrintType >= FULL  ) THEN
         PRINT "(A)" , " Writing normal modes frequencies to file NormalModes.dat"

         NrmlMdsUnit = LookForFreeUnit()
         OPEN( FILE="NormalModes.dat", UNIT=NrmlMdsUnit )
         
         WRITE(NrmlMdsUnit, "(A,/)")  &
               "# Normal modes frequencies for the initial asymptotic geometry ("//TRIM(FreqUnit(InputUnits))//")"
         DO iCoord = 1, NDim
            IF ( NormalModesVal(iCoord) >= 0. ) THEN
               WRITE(NrmlMdsUnit,*) iCoord, SQRT(NormalModesVal(iCoord))*FreqConversion(InternalUnits, InputUnits)
            ELSE
               WRITE(NrmlMdsUnit,*) iCoord, SQRT(-NormalModesVal(iCoord))*FreqConversion(InternalUnits, InputUnits), " *i"
            ENDIF
         END DO
         WRITE(NrmlMdsUnit, "(/,A,/)")  &
               "# Normal modes frequencies for the final asymptotic geometry ("//TRIM(FreqUnit(InputUnits))//")"
         DO iCoord = 1, NDim
            IF ( NormalModesVal(iCoord) >= 0. ) THEN
               WRITE(NrmlMdsUnit,*) iCoord, SQRT(FinalEigenVal(iCoord))*FreqConversion(InternalUnits, InputUnits)
            ELSE
               WRITE(NrmlMdsUnit,*) iCoord, SQRT(-FinalEigenVal(iCoord))*FreqConversion(InternalUnits, InputUnits), " *i"
            ENDIF
         END DO
         
         CLOSE( NrmlMdsUnit )
      END IF   
      
      
      ! scan over the impact parameter... 
      ImpactParameter: DO jRho = 0,NRhoMax
! 
         ! Set impact parameter and print message
         ImpactPar = float(jRho)*DeltaRho

         IF ( NRhoMax == 0 ) THEN 
            PRINT "(/,A,I5,A)"," Collinear scattering simulation ... "
         ELSE
            PRINT "(2/,A)",    "***************************************************"
            PRINT "(A,F10.5)", "       IMPACT PARAMETER = ", ImpactPar*LengthConversion(InternalUnits,InputUnits)
            PRINT "(A,/)" ,    "***************************************************"
         END IF

         PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

         ! run NTrajs number of trajectories at the current impact parameter
         Trajectory: DO iTraj = 1, NrTrajs

            IF ( PrintType == DEBUG .AND. .NOT. ZPECorrection ) THEN
               ! Temperature profile during equilibration
               WRITE(TEquilUnit,"(/,/,A,I5,/)" ) "# Trajectory nr. ", iTraj
            END IF

            PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

            !*************************************************************
            ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
            !*************************************************************

            ! Set initial conditions for the system
            ! in normal modes: classical sampling of maxwell-boltzmann or 
            ! microcanonical ensable at zero point energy
            CALL SubstrateInitialConditions( X(1:NDim), V(1:NDim) )

            ! Compute starting potential and forces
            A(:) = 0.0
            PotEnergy = ScatteringPotential( X, A )
            A(:) = A(:) / MassVector(:)

            ! compute kinetic energy and total energy
            KinEnergy = EOM_KineticEnergy(Equilibration, V )
            TotEnergy = PotEnergy + KinEnergy
            IstTemperature = 2.0*KinEnergy/(NDim-1)

            ! When the initial conditions are classical, do a Langevin equilibration
            IF ( .NOT. ZPECorrection ) THEN

               PRINT "(/,A,F6.1,1X,A)",   " Equilibrating the initial conditions at T = ", &
                   Temperature*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)

               ! Initialize temperature average and variance
               TempAverage = 0.0; TempVariance = 0.0

               ! Do an equilibration run
               EquilibrationCycle: DO iStep = 1, NrEquilibSteps

                  A(1) = 0.0; V(1) = 0.0

                  ! PROPAGATION for ONE TIME STEP
                  CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, ScatteringPotential, PotEnergy, RandomNr )
                  
                  ! compute kinetic energy and total energy
                  KinEnergy = EOM_KineticEnergy(Equilibration, V )
                  TotEnergy = PotEnergy + KinEnergy
                  IstTemperature = 2.0*KinEnergy/(NDim-1)

                  ! store temperature average and variance
                  TempAverage = TempAverage + IstTemperature
                  TempVariance = TempVariance + IstTemperature**2

                  ! every PrintStepInterval steps, write debug output
                  IF ( PrintType == DEBUG .AND. mod(iStep,EquilibrationStepInterval) == 0 ) THEN
                     ! Temperature profile during equilibration
                     WRITE(TEquilUnit,851)  real(iStep)*EquilTStep/MyConsts_fs2AU, &
                            TempAverage/iStep*TemperatureConversion(InternalUnits,InputUnits), &
                         sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)*TemperatureConversion(InternalUnits,InputUnits)
                  END IF
                  851 FORMAT( F12.5, 2F13.6 )

               END DO EquilibrationCycle

               PRINT "(A)", " Equilibration completed! "

               ! Compute average and standard deviation
               TempAverage = TempAverage / NrEquilibSteps 
               TempVariance = (TempVariance/NrEquilibSteps) - TempAverage**2
               ! output message with average values
               WRITE(*,500)  TempAverage*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits), &
                             sqrt(TempVariance)*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)

            END IF

            ! Set initial coordinates for the scattering, task = 2 fix the scatterer at initial conditions
            CALL StartSystemForScattering( X(1:NSys), V(1:NSys), MassVector(1:NSys), InitDistance, InitEKin, ImpactPar, Task=1 )
                        
            !*************************************************************
            ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
            !*************************************************************

            ! Compute energy expectation values 
            EnergyExpect = ExpectationValues( X, V )
            KinScatter   = EnergyExpect(1)
            KinSubstrate = SUM(EnergyExpect(2:NSys+1))
            KinEnergy    = KinScatter + KinSubstrate
            PotEnergy    = SUM(EnergyExpect(NSys+2:NSys+4))
            TotEnergy    = PotEnergy + KinEnergy

            ! PRINT INITIAL CONDITIONS of THE TRAJECTORY 
            WRITE(*,600)  PotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),    &
                          KinScatter*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),   &
                          KinSubstrate*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                          KinEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),    &
                          TotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

            ! Open unit for massive output, with detailed info on trajectories
            IF ( PrintType == DEBUG ) THEN
               WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Energy.dat"
               DebugUnitEn = LookForFreeUnit()
               OPEN( Unit=DebugUnitEn, File=OutFileName )

               WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Coord.dat"
               DebugUnitCoord = LookForFreeUnit()
               OPEN( Unit=DebugUnitCoord, File=OutFileName )

               WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Vel.dat"
               DebugUnitVel = LookForFreeUnit()
               OPEN( Unit=DebugUnitVel, File=OutFileName )

               ! Write initial values
               WRITE( DebugUnitEn,"(/,A)") "# TRAJECTORY ENERGY: time/fs | K_sub(NSys),K_bath,V_sub,V_bath,V_coup,V_p(...)/Eh"
               WRITE(DebugUnitEn,800) 0.0,  EnergyExpect(:)

               WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
               WRITE(DebugUnitCoord,800) 0.0, (/(X(iCoord), iCoord=1,NSys)/), (/(X(iCoord)+0.05*iCoord, iCoord=NSys+1,NDim)/)

               WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
               WRITE(DebugUnitVel,800) 0.0, V(:)
            ENDIF

            !*************************************************************
            !         TIME EVOLUTION OF THE TRAJECTORY
            !*************************************************************
    
            ! initialize counter for printing steps
            kStep = 0

            PRINT "(/,A)", " Propagating the system+bath in time... "

            ! cycle over nstep velocity verlet iterations
            Propagation: DO iStep = 1, NrOfSteps

               ! Propagate for one timestep
               CALL EOM_LangevinSecondOrder( MolecularDynamics, X, V, A, ScatteringPotential, PotEnergy, RandomNr )

               ! output to write every nprint steps 
               IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

                  ! increment counter for printing steps
                  kStep = kStep+1
                  IF ( kStep > NrOfPrintSteps ) CYCLE 

                  ! Compute energy expectation values 
                  EnergyExpect = ExpectationValues( X, V )
                  KinScatter   = EnergyExpect(1)
                  KinSubstrate = SUM(EnergyExpect(2:NSys+1))
                  KinEnergy    = KinScatter + KinSubstrate
                  PotEnergy    = SUM(EnergyExpect(NSys+2:NSys+4))
                  TotEnergy    = PotEnergy + KinEnergy

                  ! check the channel which corresponds to the current coordinates of the trajectory
                  iChan = GetCurrentChannel(X(1:NSys))
                  ! increment the array to compute channel probability
                  TrajOutcome(iChan,jRho,kStep) = TrajOutcome(iChan,jRho,kStep) + 1.0
                  ! increment the array to compute energy averages
                  IF ( PrintType >= FULL ) &
                     EnergyAver(:,iChan,jRho,kStep) = EnergyAver(:,iChan,jRho,kStep) + EnergyExpect
                  
                  ! If massive level of output, print traj information to std out
                  IF ( PrintType == DEBUG ) THEN
                     WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, EnergyExpect(:)
                     WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(1:4), &
                                             (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, NDim-4) /)
                     WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
                  END IF

               END IF 
    
            END DO Propagation

            PRINT "(A)", " Time propagation completed! "

            ! print the final values of this trajectory 
            WRITE(*,700) TRIM(GetChannelLabel(GetCurrentChannel(X(1:NSys))))
            DO iCoord = 1,NSys
               WRITE(*,702) iCoord, X(iCoord)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits)
            END DO
            WRITE(*,703) TotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

            IF ( PrintType == DEBUG ) THEN
                  CLOSE( Unit=DebugUnitEn )
                  CLOSE( Unit=DebugUnitCoord )
                  CLOSE( Unit=DebugUnitVel )
            ENDIF

         END DO Trajectory

         PRINT "(A)"," Done! "                                    

      END DO ImpactParameter

      ! Normalize energy average
      IF ( PrintType >= FULL ) THEN 
         DO iStep= 1, NrOfPrintSteps
            DO jRho = 0, NRhoMax
               DO iChan = 0, GetNrChannels( )-1
                  ! When the channel at given time and rho is populated, average is well defined
                  IF ( TrajOutcome(iChan,jRho,iStep) /= 0.0 ) &
                     EnergyAver(:,iChan,jRho,iStep) = &
                           EnergyAver(:,iChan,jRho,iStep) / TrajOutcome(iChan,jRho,iStep)
                  ! When the channel is not populated, average is left equal to zero
               END DO
            END DO
         END DO
      END IF
      
      ! Normalize probability 
      TrajOutcome(:,:,:) = TrajOutcome(:,:,:) / NrTrajs

      ! Prepare time grid and impact parameter grid
      TimeGrid = (/( float(PrintStepInterval*iStep)*TimeStep*TimeConversion(InternalUnits,InputUnits), &
                                      iStep= 1,NrOfPrintSteps  )/)
      ImpactParameterGrid = (/( float(jRho)*DeltaRho*LengthConversion(InternalUnits,InputUnits), jRho = 0,NRhoMax  )/)

      IF ( NRhoMax > 0 ) THEN
         ! Print trapping to vtk file
         CALL VTK_NewRectilinearSnapshot( PrintChannelP, X=ImpactParameterGrid, Y=TimeGrid, FileName="ChannelProb" )
         DO iChan = 0, GetNrChannels( )-2
            CALL VTK_AddScalarField( PrintChannelP, TRIM( GetChannelLabel(iChan) ), &
                     RESHAPE(TrajOutcome(iChan,:,:), (/ NrOfPrintSteps*(NRhoMax+1) /)), LetFileOpen=.TRUE. ) 
         END DO
         CALL VTK_AddScalarField( PrintChannelP, TRIM( GetChannelLabel(GetNrChannels( )-1) ), &
                  RESHAPE(TrajOutcome(GetNrChannels( )-1,:,:), (/ NrOfPrintSteps*(NRhoMax+1) /)) ) 
      END IF

      ! Print to file collinear trapping
      CollinearProbUnit = LookForFreeUnit()
      OPEN( FILE="CollinearProbabilities.dat", UNIT=CollinearProbUnit )
      WRITE(CollinearProbUnit, "(100A)") "# probability vs t (" // TRIM(TimeUnit(InputUnits)) // &
          ") -",(/ (GetChannelLabel( iChan ), iChan = 0, GetNrChannels( )-1) /)
      DO iStep= 1,NrOfPrintSteps
         WRITE(CollinearProbUnit,800) TimeGrid(iStep), TrajOutcome(:,0,iStep)
      END DO
      CLOSE(CollinearProbUnit)

      ! Print to file the average energies
      IF ( PrintType >= FULL ) THEN 
         EnergyAverUnit = LookForFreeUnit()
         ! Loop over the channels (one file per channel)
         DO iChan = 0, GetNrChannels( )-1
            ! Set file name and open unit
            WRITE(OutFileName,"(A,I3.3,A)") "AvEnergy_Channel_",iChan,".dat"
            OPEN( FILE=TRIM(OutFileName), UNIT=EnergyAverUnit )
            WRITE(EnergyAverUnit,"(A,/)")  "# Energy averages for channel: "//GetChannelLabel(iChan)
            WRITE(EnergyAverUnit,"(A,2/)") "# time / "// TRIM(TimeUnit(InputUnits)) // &
               " | K_sub(NSys), K_bath, V_sub, V_bath, V_coup, V_part(...) / " // TRIM(EnergyUnit(InputUnits))
            ! Loop over the impact parameter (section per each impact values are appended one after the other)
            DO jRho = 0, NRhoMax
               WRITE(EnergyAverUnit,"(A,F15.6,A/)")  "# Impact parameters ",ImpactParameterGrid(jRho+1),LengthUnit(InputUnits)
               DO iStep= 1, NrOfPrintSteps
                  WRITE(EnergyAverUnit,800) TimeGrid(iStep), &
                     EnergyAver(:,iChan,jRho,iStep)*EnergyConversion(InternalUnits,InputUnits)
               END DO
            END DO
            ! Close unit
            CLOSE(EnergyAverUnit)
         END DO
      END IF
      
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



   500 FORMAT (/, " Equilibration averages                 ",/     &
                  " * Average temperature              ",1F10.4,1X,A,/ &
                  "   Standard deviation               ",1F10.4,1X,A,/ ) 

   600 FORMAT (/, " Initial condition of the MD trajectory ",/     &
                  " * Potential Energy                 ",1F10.4,1X,A,/ &
                  " * Kinetic Energy of the projectile ",1F10.4,1X,A,/ &
                  " * Kinetic Energy of the substrate  ",1F10.4,1X,A,/ &
                  " * Total Kinetic Energy             ",1F10.4,1X,A,/ &
                  " * Total Energy                     ",1F10.4,1X,A,/ ) 

   700 FORMAT (/, " Trajectory in the ",A," channel " )
   702 FORMAT (   " * Final system coordinates ",I4,"  ",1F10.4,1X,A )
   703 FORMAT (   " * Final energy                     ",1F10.4,1X,A,/ ) 


   800 FORMAT(F12.5,1000F15.8)

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

      ! Deallocate memory for geometry optimization and normal mode analysis
      DEALLOCATE( Hessian, NormalModesVec, NormalModesVal, XEquil )

      ! Deallocate arrays for averages and probabilities
      IF (ALLOCATED(TrajOutcome)) DEALLOCATE(TrajOutcome)
      IF (ALLOCATED(EnergyAver)) DEALLOCATE(EnergyAver)
      
      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( Equilibration )

   END SUBROUTINE Scattering_Dispose

!===============================================================================================================================

   REAL FUNCTION ScatteringPotential( Positions, Forces )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)  :: Positions
      REAL, DIMENSION(:), INTENT(OUT) :: Forces
      REAL :: CouplingFs, DTimesX

      ! Check the number degrees of freedom
      CALL ERROR( size(Positions) /= NDim, "ScatteringSimulation.ScatteringPotential: array dimension mismatch (1)" )
      CALL ERROR( size(Forces) /= NDim, "ScatteringSimulation.ScatteringPotential: array dimension mismatch (2)" )

      ! Initialize forces and potential
      ScatteringPotential = 0.0
      Forces(:)           = 0.0

      IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! Compute potential and forces of the system
         ScatteringPotential = GetPotAndForces( Positions(1:NSys), Forces(1:NSys) )
         ! Coupling function of the system
         CouplingFs = Positions(NCoupled)
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( Bath, CouplingFs, Positions(NSys+1:), &
                           ScatteringPotential, Forces(NCoupled), Forces(NSys+1:), DTimesEffMode=DTimesX ) 
 
      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! Compute potential and forces of the system
         ScatteringPotential = GetPotAndForces( Positions(1:NSys), Forces(1:NSys) )
       END IF

   END FUNCTION ScatteringPotential
   
!===============================================================================================================================

!**************************************************************************************
!> Compute the Hessian of the system-bath potential in mass square root scaled 
!> coordinates as it is needed for the normal modes analysis  
!>
!> @param AtPoint   Input vector with the coordinates where to compute Hessian
!> @returns         Hessian matrix of the potential in AtPoint
!**************************************************************************************
   FUNCTION GetMassScaledHessian( AtPoint ) RESULT( Hessian )
      REAL, DIMENSION(NDim), INTENT(IN)  :: AtPoint
      REAL, DIMENSION(NDim,NDim)         :: Hessian
      
      REAL, DIMENSION(NBath)   :: CouplingsCoeffs
      REAL                     :: DistortionCoeff
      INTEGER                  :: i,j

      ! Initialize Hessian equal to 0 (necessary for 0 elements in system+bath matrix)
      Hessian(:,:) = 0.0

      ! Use subroutine GetHessianFromForces to compute the numerical hessian for the subsystem
      Hessian(1:NSys,1:NSys) = GetHessianFromForces( AtPoint(1:NSys), GetPotAndForces, 1.E-3 ) 
      ! and transform it to mass weighted coordinates
      DO j = 1, NSys
         DO i = 1, NSys
            Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
         END DO
      END DO
      
      IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! compute bath hessian, which is computed by the same subroutine, regardless bath type
         CALL  HessianOfTheBath( Bath, Hessian(NSys+1:NDim, NSys+1:NDim) ) 
         CALL  CouplingAndDistortionHessian( Bath, CouplingsCoeffs, DistortionCoeff )
         
         ! add the contribution from DISTORTION CORRECTION, equal in normal bath and chain bath 
         Hessian(NCoupled,NCoupled) = Hessian(NCoupled,NCoupled) + DistortionCoeff / MassVector(NCoupled)
      END IF
      
      ! Add contribution coming from COUPLING, which is different in Normal and Chain bath
      IF ( BathType == NORMAL_BATH ) THEN
         DO i = NSys+1, NDim
            Hessian(NCoupled,i) = Hessian(NCoupled,i) + CouplingsCoeffs(i-NSys)/SQRT(MassVector(NCoupled)*MassVector(i) ) 
            Hessian(i,NCoupled) = Hessian(i,NCoupled) + CouplingsCoeffs(i-NSys)/SQRT(MassVector(NCoupled)*MassVector(i) ) 
         END DO
      ELSE IF ( BathType == CHAIN_BATH ) THEN
         Hessian(NCoupled,NSys+1) = Hessian(NCoupled,NSys+1) + CouplingsCoeffs(1)/SQRT(MassVector(NCoupled)*MassVector(NSys+1)) 
         Hessian(NSys+1,NCoupled) = Hessian(NSys+1,NCoupled) + CouplingsCoeffs(1)/SQRT(MassVector(NCoupled)*MassVector(NSys+1)) 
      END IF
         
   END FUNCTION GetMassScaledHessian
   

! ===============================================================================================================================

!**************************************************************************************
!> Return current values of the energy expectation values at given time 
!> to be averaged over the trajectories ensamble.
!>
!> @param Positions        Array with current positions
!> @param Velocities       Array with current velocities
!> @returns Expectations   Array with the energy expectation values
!**************************************************************************************
   FUNCTION ExpectationValues( X, V ) RESULT( Expectations )
      IMPLICIT NONE
      REAL, DIMENSION(NrEnAverages) :: Expectations   
      REAL, DIMENSION(:), INTENT(IN)  :: X
      REAL, DIMENSION(:), INTENT(IN)  :: V

      REAL :: VCoupling, VBath, ReducedMass
      INTEGER :: i
      
      ! 1 -- NSys) kinetic energies of the subsystem coordinates
      DO i = 1, NSys 
         Expectations(i) = EOM_KineticEnergy(MolecularDynamics, V, i, i)
      END DO
      
      ! NSys+1 ) bath kinetic energy
      IF ( BathType == LANGEVIN_DYN ) THEN
         Expectations(NSys+1) = 0.0
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         Expectations(NSys+1) = EOM_KineticEnergy(MolecularDynamics, V, NDim, NSys+1)
      END IF
      
      ! NSys+2 ) potential energy of the subsystem
      Expectations(NSys+2) = GetPotential( X(1:NSys) )
      ! NSys+3 ) potential energy of the bath
      ! NSys+4 ) system-bath coupling energy 
      IF ( BathType == LANGEVIN_DYN ) THEN
         Expectations(NSys+3:NSys+4) = 0.0
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         CALL EnergyOfTheBath( Bath, X(NCoupled), X(NSys+1:), VCoupling, VBath )
         Expectations(NSys+3) = VBath
         Expectations(NSys+4) = VCoupling
      END IF
      
      ! Compute internal energy corresponding to Hinc-Htar fragment
      ReducedMass = 1.0/(1./MassVector(1)+1./MassVector(2))
      Expectations(NSys+5) = 0.5 * ReducedMass * (V(1)-V(2))**2
      
      ! Compute potential energy partitioning
      ! NSys+4) potential energy of the carbon for H1 and H2 far from surface at eq position
      ! NSys+5) potential energy of H-H, for non interacting planar carbon
      ! NSys+6) potential energy of the C-H, with other H non interacting
      Expectations(NSys+6:NSys+8) = GetVPartitions( X(1:NSys) )
      
   END FUNCTION ExpectationValues

! ===============================================================================================================================

!**************************************************************************************
!> Setup initial condition for the bath and those degrees of freedom which are 
!> not scattering coordinates. The scatterer is fixed far from the interaction
!> region and can be subsequently fixed in conditions which are suitable for 
!> relaxation or scattering. 
!> The initial state can be chosen to be either a finite temperature classical 
!> state or a quasi-classical T=0 K state.
!>
!> @param X              Output vector with the initial positions of the system
!> @param V              Output vector with the initial Velocities of the system
!**************************************************************************************
   SUBROUTINE SubstrateInitialConditions( X, V )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(OUT)  :: X
      REAL, DIMENSION(:), INTENT(OUT)  :: V
      REAL                             ::  Value, CarbonFreq, SigmaV
      INTEGER                          ::  i
      
      ! Check the number degrees of freedom
      CALL ERROR( size(X) /= NDim, "ScatteringSimulation.SubstrateInitialConditions: array dimension mismatch (1)" )
      CALL ERROR( size(V) /= NDim, "ScatteringSimulation.SubstrateInitialConditions: array dimension mismatch (2)" )

      DO i = 1, NDim            ! cycle over the coordinates
         ! Set sigma of the maxwell boltzmann distribution
         SigmaV = sqrt( Temperature )

         ! the coordinate is bound
         IF ( NormalModesVal(i) > 0. ) THEN
            ! set frequency of the normal mode
            CarbonFreq = SQRT( NormalModesVal(i) )

            ! Quasiclassical distribution
            IF ( ZPECorrection ) THEN
               Value = UniformRandomNr( RandomNr )*2.0*MyConsts_PI
               X(i) = COS(Value) / SQRT(CarbonFreq)
               V(i) = SIN(Value) * SQRT(CarbonFreq)
               
            ! Classical distribution
            ELSE IF ( .NOT. ZPECorrection ) THEN
               X(i) = GaussianRandomNr(RandomNr) * SigmaV / CarbonFreq
               V(i) = GaussianRandomNr(RandomNr) * SigmaV
            END IF
            
         ELSE
            ! unbound coordinate (set velocity according to maxwell-boltzman) 
            X(i) = 0.0
            V(i) = 0.0 !GaussianRandomNr(RandomNr) * SigmaV

         ENDIF
         
      END DO

      ! TRASFORM BACK TO ORIGINAL FRAME AND TO NOT-MASS-WEIGHTED COORDINATES
      X = TheOneWithMatrixVectorProduct( NormalModesVec, X )
      V = TheOneWithMatrixVectorProduct( NormalModesVec, V )
      DO i = 1, NDim
         X(i) = X(i) / SQRT( MassVector(i) )
         V(i) = V(i) / SQRT( MassVector(i) )
      END DO
      X = X + XEquil
      
   END SUBROUTINE SubstrateInitialConditions
   
      
END MODULE ScatteringSimulation
