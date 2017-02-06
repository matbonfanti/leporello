!***************************************************************************************
!*                              PROGRAM leporello
!***************************************************************************************
!>  \mainpage      Program leporello 
!>
!>  Classical simulations of H + H + C(graphene) + dissipative bath       \n
!>  * Model: 1D for Hinc + 1D for Htar + 1D Z coordinate for carbon atom  \n
!>  * Bath:  normal bath, chain bath                                      \n
!>  * Propagation: Velocity-Verlet in the microcanonical ensamble         \n
!>                 symplectic in the canonical ensamble                   \n
!>
!>  \author        Matteo Bonfanti 
!>  \version       1.0
!>  \date          9 February 2013
!>
!>
!***************************************************************************************
PROGRAM leporello
#include "preprocessoptions.cpp"
   USE SharedData
   USE InputField
   USE UnitConversion
   USE PotentialAnalysis
!    USE MinimumEnergyPath
!    USE VibrationalRelax
!    USE ThermalEquilibrium
   USE ScatteringSimulation
   USE IndependentOscillatorsModel
   USE PotentialModule

   IMPLICIT NONE

   ! Variable to handle the command line
   INTEGER :: NArgs
   LOGICAL :: Help = .FALSE.
   
   ! Input file name, set from command line arguments
   CHARACTER(120) :: InputFileName

   ! Derived type to handle input data
   TYPE(InputFile) :: InputData
   ! Units of input data, defined from the input file
   INTEGER     :: InputLength, InputEnergy, InputMass, InputTime, InputTemp, InputFreq

   ! keep track of initial time and final time
   INTEGER, DIMENSION(8)    :: Time1, Time2

#if defined(__PRINT_SPECTRAL_DENSITY)
   ! Unit to write bath parameters to output
   INTEGER :: SpectralDensityUnit
#endif

   __TIME_STRING__
   
   PRINT "(/,     '                    ==============================')"
   PRINT "(       '                               leporello          ')"
   PRINT "(       '                    ==============================',/)"
   PRINT "(       '                       Author: Matteo Bonfanti'      )"
   PRINT "(       '                       Release: ',A)", VERSIONTAG
   PRINT "(       '                       Compilation: ',A,1X,A,/)", __DATE__, __TIME__ 

   PRINT "(       '                       << Notte e giorno faticar      ')"
   PRINT "(       '                         per chi nulla sa gradir      ')"
   PRINT "(       '                         piova e vento sopportar,     ')"
   PRINT "(       '                        mangiar male e mal dormir...  ')"
   PRINT "(       '                         Voglio far il gentiluomo,    ')"
   PRINT "(       '                        e non voglio più servir... >> '/)"
   PRINT "(       '                      [Don Giovanni act I scene 1]  ',2/)"

#if defined(LOG_FILE) 
   __INIT_LOG_FILE
#endif
   
   CALL date_and_time (values=Time1)

   !*************************************************************
   !         COMMAND LINE ARGUMENT
   !*************************************************************

   ! Check and read from command line the input file name
   NArgs = COMMAND_ARGUMENT_COUNT()
   IF (NArgs<1) THEN
      Help = .TRUE.
   ELSE
      CALL GET_COMMAND_ARGUMENT( 1, InputFileName )
      IF ( trim(InputFileName) == "help" ) Help = .TRUE.
   ENDIF
   IF (Help) THEN ! Call help
      PRINT*, ' Launch this program as:'
      PRINT*, ' % leporello "InputFileName" '
      STOP
   ENDIF

   !*************************************************************
   !         INPUT SECTION 
   !*************************************************************

   ! Open and read from input file the input parameters of the calculation
   CALL OpenFile( InputData, InputFileName )

   ! read input units ( or set them to default value )

   !      ***** DEFAULT VALUES ******
   !      distance    - Angstrom
   !      energy      - electronVolt
   !      mass        - AMU
   !      time        - femtosecond
   !      temperature - Kelvin
   !      frequency   - cm-1

   CALL SetFieldFromInput( InputData, "InputLength", InputLength,  1 )
   CALL SetFieldFromInput( InputData, "InputEnergy", InputEnergy,  3 )
   CALL SetFieldFromInput( InputData, "InputMass",   InputMass,    8 )
   CALL SetFieldFromInput( InputData, "InputTime",   InputTime,   13 )
   CALL SetFieldFromInput( InputData, "InputTemp",   InputTemp,   16 )
   CALL SetFieldFromInput( InputData, "InputFreq",   InputFreq,   18 )
   CALL Initialize_UnitConversion( InputUnits, InputLength, InputEnergy, InputMass, 11, &
                                   InputTime, InputTemp, InputFreq )

   ! Define the kind of simulation to do
   CALL SetFieldFromInput( InputData, "RunType", RunType, 1 )
   CALL CheckRunType( RunType )
   ! Set the representation of the bath
   CALL SetFieldFromInput( InputData, "BathType", BathType )
   CALL CheckBathType( BathType )
   ! Set the level of output
   CALL SetFieldFromInput( InputData, "PrintType", PrintType )
   CALL CheckPrintType( PrintType )

   ! Hydrogen and carbon masses
   CALL SetFieldFromInput( InputData, "MassHTar", MassHTar )
   MassHTar = MassHTar * MassConversion(InputUnits, InternalUnits)
   CALL SetFieldFromInput( InputData, "MassHInc", MassHInc )
   MassHInc = MassHInc * MassConversion(InputUnits, InternalUnits)
   CALL SetFieldFromInput( InputData, "MassC", MassC )
   MassC = MassC * MassConversion(InputUnits, InternalUnits)

   ! SET THE BATH-REPRESENTATION DEPENDENT VARIABLES 

   IF ( BathType == NORMAL_BATH ) THEN

      ! No langevin oscillators in the normal bath
      DynamicsGamma = 0.0
      ! Mass of the bath oscillator
      CALL SetFieldFromInput( InputData, "MassBath", MassBath )
      MassBath = MassBath * MassConversion(InputUnits, InternalUnits)
      ! Nr of bath degrees of freedom
      CALL SetFieldFromInput( InputData, "NBath",  NBath ) 
      ! Read ohmic spectral density, when zero read spectral density from file
      CALL SetFieldFromInput( InputData, "OhmicGammaTimesMass", OhmicGammaTimesMass, 0.0 )
      IF ( OhmicGammaTimesMass /= 0 ) THEN
         OhmicGammaTimesMass = OhmicGammaTimesMass
         ! For an ohmic SD, cutoff freq is compulsory
         CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq )
      ELSE IF ( OhmicGammaTimesMass == 0.0 ) THEN
         ! Read file with spectral density
         CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )
         ! Read cutoff frequency of the bath, if BathCutOffFreq is not present, it is set to zero
         CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq, 0.0 )
      END IF
      BathCutOffFreq = BathCutOffFreq * FreqConversion(InputUnits, InternalUnits) 
      ! Set lower boundary for oscillator bath frequency
      CALL SetFieldFromInput( InputData, "BathLowerCutOffFreq", BathLowerCutOffFreq, 0.0 )
      BathLowerCutOffFreq = BathLowerCutOffFreq * FreqConversion(InputUnits, InternalUnits) 
      ! Set non linear system-bath coupling
      CALL SetFieldFromInput( InputData, "NonLinearCoupling", NonLinearCoupling, .FALSE. )
      IF ( NonLinearCoupling ) THEN
         CALL SetFieldFromInput( InputData, "AlphaCoupling", AlphaCoupling )
         AlphaCoupling = AlphaCoupling / LengthConversion(InputUnits, InternalUnits)
      END IF

   ELSE IF ( BathType == CHAIN_BATH ) THEN

      ! Langevin relaxation at the end of the chain
      CALL SetFieldFromInput( InputData, "RelaxAtChainEnd",  DynamicsGamma, 0.0 ) 
      IF ( DynamicsGamma /= 0 ) &
         DynamicsGamma = 1. / ( DynamicsGamma * TimeConversion(InputUnits, InternalUnits) )
      ! Mass of the bath oscillator
      CALL SetFieldFromInput( InputData, "MassBath", MassBath )
      MassBath = MassBath * MassConversion(InputUnits, InternalUnits)
      ! Nr of bath degrees of freedom
      CALL SetFieldFromInput( InputData, "NBath",  NBath )
      ! Read cutoff frequency of the bath, if BathCutOffFreq is not present, it is set to zero
      CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq, 0.0 )
      BathCutOffFreq = BathCutOffFreq * FreqConversion(InputUnits, InternalUnits) 
      ! Read file with normal modes freq and couplings
      CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )

   ELSE IF ( BathType == LANGEVIN_DYN ) THEN

      ! Langevin relaxation of the system (at the carbon atom)
      CALL SetFieldFromInput( InputData, "RelaxAtCarbon",  DynamicsGamma, 0.0 ) 
      IF ( DynamicsGamma /= 0 ) &
         DynamicsGamma = 1. / ( DynamicsGamma * TimeConversion(InputUnits, InternalUnits) )

   END IF

   ! Variable to set the potential derivatives testing
   CALL SetFieldFromInput( InputData, "DerivTesting", DerivTesting, .FALSE. )
   ! Variables to set sudden and adiabatic approximations
   CALL SetFieldFromInput( InputData, "AdiabaticV", AdiabaticV, .FALSE. )
   CALL SetFieldFromInput( InputData, "SuddenV", SuddenV, .FALSE. )
   CALL ERROR( SuddenV .AND. AdiabaticV, " impossible to choose adiabatic and sudden pot at the same time ")
   CALL ERROR( (SuddenV .OR. AdiabaticV) .AND. (BathType /= LANGEVIN_DYN .OR. DynamicsGamma /= 0), &
             " adiabatic and sudden pot available only for system-only microcanonical dynamics" )

   !*************************************************************
   !       PRINT OF THE INPUT DATA TO STD OUT
   !*************************************************************

   ! Write info about the kind of calculation
   SELECT CASE( RunType )
!       CASE( EQUILIBRIUM )
!          WRITE(*,"(/,A)") " * Atom-surface equilibrium simulation "
!       CASE( RELAXATION )
!          WRITE(*,"(/,A)") " * Atom-surface vibrational relaxation simulation "
      CASE( SCATTERING )
         WRITE(*,"(/,A)") " * Atom-surface sticking simulation "
      CASE( POTENTIAL )
         WRITE(*,"(/,A)") " * Analysis of the potential energy surfaces "
   END SELECT

   WRITE(*,898) MassHInc*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits), &
                MassHTar*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits), &
                MassC*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits)

   ! Write info about the bath representation
   SELECT CASE( BathType )

      CASE( NORMAL_BATH ) 
         IF ( OhmicGammaTimesMass == 0.0 ) THEN
            WRITE(*,900) NBath, MassBath*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits), &
                         BathCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),  &
                         BathLowerCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),  &
                         trim(adjustl(SpectralDensityFile))
         ELSE
            WRITE(*,910) NBath, MassBath*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits), &
                         BathCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),  &
                         OhmicGammaTimesMass
         END IF

      CASE( CHAIN_BATH )
         IF (DynamicsGamma /= 0. ) THEN
            WRITE(*,901) NBath, MassBath*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits),   &
                        BathCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),    &
                        1.0/DynamicsGamma*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits), &
                        trim(adjustl(SpectralDensityFile))
         ELSE
            WRITE(*,801) NBath, MassBath*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits),   &
                        BathCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),    &
                        trim(adjustl(SpectralDensityFile))
         END IF

         
      CASE( LANGEVIN_DYN )
         IF (DynamicsGamma /= 0. ) THEN
            WRITE(*,902) 1.0/DynamicsGamma*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits)
         ELSE
            WRITE(*,802) 
         END IF

   END SELECT

   ! Write info about the kind of output
   SELECT CASE( PrintType )
      CASE( MINIMAL )
         WRITE(*,"(A,/)") " * Minimal output will be written "
      CASE( FULL )
         WRITE(*,"(A,/)") " * All the averages will be written to output files "
      CASE( DEBUG )
         WRITE(*,"(A,/)") " * Detailed information on each trajectory will be printed "
   END SELECT

   898 FORMAT(" * Mass of the incident H atom:                 ",F10.4,1X,A,/,&
              " * Mass of the target H atom:                   ",F10.4,1X,A,/,&
              " * Mass of the C atom:                          ",F10.4,1X,A,/ )

   900 FORMAT(" * Bath is a set of independent HO coupled to the system ",/,&
              " * Nr of bath oscillators:                      ",I10,  /,& 
              " * Mass of the bath oscillator:                 ",F10.4,1X,A,/,& 
              " * Cutoff frequency of the bath:                ",F10.1,1X,A,/,& 
              " * Lower cutoff frequency:                      ",F10.1,1X,A,/,& 
              " * File with the spectral density:  "            ,A22,/ )

   910 FORMAT(" * Bath is a set of independent HO coupled to the system, ohmic SD ",/,&
              " * Nr of bath oscillators:                      ",I10,  /,& 
              " * Mass of the bath oscillator:                 ",F10.4,1X,A,/,& 
              " * Cutoff frequency of the bath:                ",F10.1,1X,A,/,& 
              " * Relaxation time of the ohmic SD              ",F10.4,1X,"au",/ )

   901 FORMAT(" * Bath is is a linear chain of harmonic oscillators ", /,&
              " * Nr of bath oscillators:                      ",I10,  /,& 
              " * Mass of the bath oscillator:                 ",F10.4,1X,A,/,& 
              " * Cutoff frequency of the bath:                ",F10.1,1X,A,/,& 
              " * Langevin relax time at the end of the chain: ",F10.4,1X,A,/,&
              " * File with the spectral density:  "            ,A22,  / )
   801 FORMAT(" * Bath is is a linear chain of harmonic oscillators ", /,&
              " * Nr of bath oscillators:                      ",I10,  /,& 
              " * Mass of the bath oscillator:                 ",F10.4,1X,A,/,& 
              " * Cutoff frequency of the bath:                ",F10.1,1X,A,/,& 
              " * Infinite relax time at the end of the chain  ",           /,&
              " * File with the spectral density:  "            ,A22,  / )

   902 FORMAT(" * Bath is effectively represented by Langevin dynamics ", /,&
              " * Relaxation time of Langevin dynamics:        ",F10.4,1X,A,/ )
   802 FORMAT(" * Bath is effectively represented by Langevin dynamics ", /,&
              " * Infinite relaxation time                             ", / )

   !*************************************************************
   !       POTENTIAL SETUP 
   !*************************************************************

   ! Setup potential energy surface
   IF ( AdiabaticV ) THEN
      CALL SetupPotential( ADIABATIC, DerivTesting  )
   ELSE IF ( SuddenV ) THEN 
      CALL SetupPotential( SUDDEN,    DerivTesting  )
   ELSE
      CALL SetupPotential( FULLPOT,   DerivTesting  )
   END IF

   !*************************************************************
   !       IO BATH SETUP 
   !*************************************************************
   
   ! If needed setup bath frequencies and coupling for oscillator bath models
   IF (  BathType == NORMAL_BATH ) THEN
         IF ( OhmicGammaTimesMass == 0.0 ) THEN
            CALL SetupIndepOscillatorsModel( Bath, NBath, 0, SpectralDensityFile, MassBath, BathCutOffFreq, BathLowerCutOffFreq )
         ELSE
            CALL SetupOhmicIndepOscillatorsModel( Bath, NBath, 0, OhmicGammaTimesMass, MassBath, BathCutOffFreq )
         END IF
         IF ( NonLinearCoupling ) CALL SetNonLinearCoupling( Bath, AlphaCoupling )
   ELSE IF (  BathType == CHAIN_BATH ) THEN
         CALL SetupIndepOscillatorsModel( Bath, NBath, 1, SpectralDensityFile, MassBath, BathCutOffFreq )
   END IF
   
   IF  ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) &
      PRINT "(/,A,F10.6,A,/)"," * Bath distorsion force constant:              ", GetDistorsionForce( Bath ), " au"

#if defined(__PRINT_SPECTRAL_DENSITY)
      SpectralDensityUnit = LookForFreeUnit()
      OPEN( FILE="ReadSpectralDensity.dat", UNIT=SpectralDensityUnit )

      IF ( OhmicGammaTimesMass == 0.0 ) THEN
         WRITE(SpectralDensityUnit, "(A,A)") "# Spectral Density read from file: ", TRIM(ADJUSTL(FileName))
      ELSE
         WRITE(SpectralDensityUnit, "(A,1F15.6)") "# Ohmic Spectral Density with mass*gamma = ", SysMassTimeGamma
      ENDIF

      IF ( BathType == CHAIN_BATH ) THEN
         WRITE(SpectralDensityUnit, "(A,A,/)") "# Bath in linear chain form "
      ELSE IF ( BathType == STANDARD_BATH ) THEN
         WRITE(SpectralDensityUnit, "(A,A,/)") "# Bath in normal form "
      END IF

      DO iBath = 1, Bath%BathSize
         WRITE(SpectralDensityUnit,"(2F20.12)") Bath%Frequencies(iBath)*FreqConversion(InternalUnits,InputUnits),&
                        Bath%Couplings(iBath) 
      END DO
      WRITE(SpectralDensityUnit,"(/)") 
#endif

   !*************************************************************
   !       SPECIFIC INPUT SECTION 
   !*************************************************************

   SELECT CASE( RunType )
!       CASE( EQUILIBRIUM )
!          CALL ThermalEquilibrium_ReadInput( InputData )
!       CASE( RELAXATION )
!          CALL VibrationalRelax_ReadInput( InputData )
      CASE( SCATTERING )
         CALL Scattering_ReadInput( InputData )
      CASE( POTENTIAL )
         CALL PotentialAnalysis_ReadInput( InputData )
   END SELECT

   CALL CloseFile( InputData )

   !*************************************************************
   !       INITIALIZE AND RUN CALCULATION 
   !*************************************************************

   SELECT CASE( RunType )
!       CASE( EQUILIBRIUM )
!          CALL ThermalEquilibrium_Initialize( )
!          CALL ThermalEquilibrium_Run()
!       CASE( RELAXATION )
!          CALL VibrationalRelax_Initialize( )
!          CALL VibrationalRelax_Run()
      CASE( SCATTERING )
         CALL Scattering_Initialize()
         CALL Scattering_Run()
      CASE( POTENTIAL )
         CALL PotentialAnalysis_Initialize()
         CALL PotentialAnalysis_Run()
   END SELECT

   !*************************************************************
   !       DISPOSE MEMORY AND TERMINATE EXECUTION
   !*************************************************************

   SELECT CASE( RunType )
!       CASE( EQUILIBRIUM )
!          CALL ThermalEquilibrium_Dispose()
!       CASE( RELAXATION )
!          CALL VibrationalRelax_Dispose()
      CASE( SCATTERING )
         CALL Scattering_Dispose()
      CASE( POTENTIAL )
         CALL PotentialAnalysis_Dispose()
   END SELECT

   IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) &
      CALL DisposeIndepOscillatorsModel( Bath )

   CALL date_and_time (values=Time2)
   
   WRITE(*,*)
   WRITE(*,"(A,F10.1,A)") " Execution Time : ",TimeDifference( Time2, Time1 )/1000.0, " / s "
   
#if defined(LOG_FILE) 
   __END_LOG_FILE
#endif
   
      CONTAINS
   
   FUNCTION TimeDifference( Time1, Time2 )
      INTEGER, DIMENSION(8), INTENT(IN)    :: Time1, Time2         
      REAL :: TimeDifference
   
      TimeDifference =  (Time1(3)-Time2(3))*24.0*60.0*60.0*1000. + &
                        (Time1(5)-Time2(5))*60.0*60.0*1000. + &
                        (Time1(6)-Time2(6))*60.0*1000. + &
                        (Time1(7)-Time2(7))*1000. + &
                        (Time1(8)-Time2(8))
         
   END FUNCTION TimeDifference
   
END PROGRAM leporello


