!***************************************************************************************
!*                      MODULE IndependentOscillatorsModel
!***************************************************************************************
!
!>  \brief      Independent oscillator model
!>  \details    This class define the necessary subroutines to define and compute
!>              the potential and the derivatives of a discrete independent oscillator
!>              model, both in the standard and in the chain representation
!
!***************************************************************************************
!
!>  \author     Matteo Bonfanti
!>  \version    2.0
!>  \date       8 March 2013
!>
!***************************************************************************************
!
!>  \pre        To use the class the potential needs to be setup with the 
!>              couplings and the frequencies of the oscillators
!>              This data is read from input file by the setup subroutine
!>              and then stored in the BathData object
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg 22 October 2013: to implement the possibility of having multiple baths, 
!>             a bath datatype has been implemented and all the subroutine have been adapted accordingly
!>  \arg 8 Novembre 2013: implemented debug output of frequencies and couplings
!>  \arg 28 January 2015: non linear coupling implemented for normal bath
!>  \arg 3 February 2015: sub BathPotentialAndForces has a new optional arg for effective mode
!>  \arg 8 February 2016: section of the __PRINT_SPECTRAL_DENSITY directive has been moved at an higher
!>                        level: no formatted output is written here except for the LOG file
!>  \arg 6 February 2017: in previous revisions the Hessian for the bath was erroneously computed!!!!
!>                        this wrong 0.5 factor has been removed in this version
!
!>  \todo    Initial thermal conditions of a chain bath can be refined by taking into account the coupling
!>  \todo    Quasi-classical 0K conditions of a chain bath need to be implemented
!
!***************************************************************************************

MODULE IndependentOscillatorsModel
#include "preprocessoptions.cpp"
   USE SplineInterpolator
   USE RandomNumberGenerator
   USE FFTWrapper

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetupIndepOscillatorsModel, SetupOhmicIndepOscillatorsModel, SetNonLinearCoupling
   PUBLIC :: BathPotentialAndForces, DisposeIndepOscillatorsModel
   PUBLIC :: ThermalEquilibriumBathConditions, ZeroKelvinBathConditions, BathOfRingsThermalConditions
   PUBLIC :: EnergyOfTheBath, GetDistorsionForce, CouplingAndDistortionHessian, HessianOfTheBath
   PUBLIC :: GetFirstEffectiveMode

   PUBLIC :: BathData

   INTEGER, PARAMETER :: STANDARD_BATH = 0        ! < normal bath, in which all the oscillators are coupled to the system
   INTEGER, PARAMETER :: CHAIN_BATH    = 1        ! < chain bath, in which all the oscillators are coupled in chain
    
   !> Derived data type to store all the relevant information of a single bath
   TYPE BathData
!       PRIVATE
      INTEGER :: BathType                             !< Integer parameter to store the bath type (chain, standard ... )
      INTEGER :: BathSize                             !< Number of oscillators of the bath
      REAL, DIMENSION(:), POINTER :: Frequencies      !< Harmonic frequencies of the bath (stored in AU)
      REAL, DIMENSION(:), POINTER :: Couplings        !< Coupling of the bath oscillators (stored in AU)
      REAL :: CutOff                                  !< Cutoff frequency of the normal bath (stored in AU)
      REAL :: LowerCutOff                             !< Lower cutoff frequency of the normal bath (stored in AU)
      REAL :: DeltaOmega                              !< Frequency spacing of the normal bath (stored in AU)
      REAL :: OscillatorsMass                         !< Mass of the oscillators
      REAL :: DistorsionForce                         !< Force constant of the distorsion correction
      LOGICAL :: NonLinearCoup = .FALSE.              !< Logical variable to set exponential system-bath coupling
      REAL :: NonLinCoup_Alpha                        !< Alpha coefficient of non linear coupling
      REAL :: NonLinCoup_D0                           !< D coefficient of non linear coupling
      REAL :: NonLinCoup_OmegaBar                     !< Omega Bar coefficient of non linear coupling
      REAL :: NonLinCoup_Beta                         !< Beta coefficient of non linear coupling
      LOGICAL :: BathIsSetup = .FALSE.                !< Logical variable to set status of the class
   END TYPE BathData

   !> Max nr of iterations for potential optimization
   INTEGER, PARAMETER :: MaxIter = 10000
   !> Threshold for conjugate gradient convergence
   REAL, PARAMETER :: GradEps = 1.0E-4

CONTAINS   

!===============================================================================================================


!*******************************************************************************
!                          SetupIndepOscillatorsModel
!*******************************************************************************
!>  This subroutine can be used to setup the parameters of the oscillator bath
!>  reading the couplings and the frequencies from input file
!> 
!> @param    FileName       Name of the file with the input parameters.
!*******************************************************************************
   SUBROUTINE SetupIndepOscillatorsModel( Bath, N, SetBathType, FileName, Mass, CutOffFreq, LowerCutOffFreq )
      IMPLICIT NONE
      TYPE(BathData)             :: Bath
      INTEGER, INTENT(IN)        :: N
      INTEGER, INTENT(IN)        :: SetBathType
      CHARACTER(*), INTENT(IN)   :: FileName
      REAL, INTENT(IN)           :: Mass
      REAL, INTENT(IN)           :: CutOffFreq
      REAL, OPTIONAL, INTENT(IN) :: LowerCutOffFreq

      INTEGER :: InpUnit, RdStatus
      INTEGER :: iBath, NData
      LOGICAL :: FileIsPresent
      REAL, DIMENSION(:), ALLOCATABLE :: RdFreq, RdSpectralDens
      TYPE(SplineType) :: SpectralDensitySpline
      REAL :: D0

      ! If data already setup give a warning and deallocate memory
      CALL WARN( Bath%BathIsSetup, "IndependentOscillatorsModel.SetupIndepOscillatorsModel: overwriting bath data" )
      IF (Bath%BathIsSetup) CALL DisposeIndepOscillatorsModel( Bath )

      ! Store number of bath degrees of freedom
      Bath%BathSize = N
      ! Set the type of bath
      Bath%BathType = SetBathType
      ! Set the mass of the oscillators
      Bath%OscillatorsMass = Mass

      ! Allocate memory
      ALLOCATE( Bath%Frequencies(Bath%BathSize), Bath%Couplings(Bath%BathSize) )

      ! Check if spectral density file exists
      INQUIRE( File = TRIM(ADJUSTL(FileName)), EXIST=FileIsPresent ) 
      CALL ERROR( .NOT. FileIsPresent,            &
                "IndependentOscillatorsModel.SetupIndepOscillatorsModel: spectral density file does not exists" )

      IF ( PRESENT( LowerCutOffFreq ) ) THEN
         Bath%LowerCutOff = LowerCutOffFreq
      ELSE
         Bath%LowerCutOff = 0.0
      END IF

      IF ( Bath%BathType == CHAIN_BATH ) THEN

            ! Open input file
            InpUnit = LookForFreeUnit()
            OPEN( File = TRIM(ADJUSTL(FileName)), Unit = InpUnit )

            ! Read force constant of the distorsion
            READ(InpUnit,*) Bath%DistorsionForce

            ! Read frequencies and couplings
            DO iBath = 1, Bath%BathSize
               READ(InpUnit,*, IOSTAT=RdStatus) Bath%Frequencies(iBath), Bath%Couplings(iBath)
               ! if data is missing, a rubin continuation of the chain is assumed 
               IF ( RdStatus /= 0 ) THEN
                  IF ( CutOffFreq > 0.0 ) THEN
                     Bath%Frequencies(iBath) = CutOffFreq / SQRT(2.)
                     Bath%Couplings(iBath) = 0.4999 * Bath%Frequencies(iBath)**2
                  ELSE 
                     Bath%Frequencies(iBath) = Bath%Frequencies(iBath-1)
                     Bath%Couplings(iBath) = Bath%Couplings(iBath-1)
                  ENDIF
               END IF
            END DO

            ! Scale with the masses
            Bath%Couplings(1) = Bath%Couplings(1) * sqrt( Bath%OscillatorsMass )
            DO iBath = 2, Bath%BathSize
               Bath%Couplings(iBath) = Bath%Couplings(iBath) * Bath%OscillatorsMass
            END DO
            D0 = Bath%Couplings(1)

            ! Close input file
            CLOSE( InpUnit )

      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN

            ! Count the number of data in the input file
            NData = CountLinesInFile( FileName )

            ! Allocate temporary arrays to store the input spectral density
            ALLOCATE( RdFreq(NData), RdSpectralDens(NData) )

            ! Open input file
            InpUnit = LookForFreeUnit()
            OPEN( File = TRIM(ADJUSTL(FileName)), Unit = InpUnit )

            ! Read frequencies and couplings
            DO iBath = 1, NData
               READ(InpUnit,*) RdFreq(iBath), RdSpectralDens(iBath)
               ! Transform frequencies from cm-1 to au
               RdFreq(iBath) = RdFreq(iBath) * MyConsts_cmmin1toAU
            END DO

            ! Close input file
            CLOSE( InpUnit )

            ! Set cutoff frequency
            IF ( CutOffFreq > 0.0 ) THEN
               Bath%CutOff = CutOffFreq               ! if is given from input
            ELSE 
               Bath%CutOff = RdSpectralDens( NData )  ! otherwise is the last freq in the input file
            ENDIF
            Bath%DeltaOmega = (Bath%CutOff - Bath%LowerCutOff) / real( Bath%BathSize)

            ! Spline interpolation of the input spectral density
            CALL SetupSpline( SpectralDensitySpline, RdFreq, RdSpectralDens)

            ! initialize distortion constant
            Bath%DistorsionForce = 0.0

            ! Compute frequencies and couplings 
            D0 = 0.0
            DO iBath = 1,  Bath%BathSize
               Bath%Frequencies(iBath) = Bath%LowerCutOff + iBath * Bath%DeltaOmega
               Bath%Couplings(iBath) = SQRT( 2.0 * Bath%OscillatorsMass * Bath%Frequencies(iBath) * Bath%DeltaOmega *      & 
                    ABS(GetSpline( SpectralDensitySpline, Bath%Frequencies(iBath) )) / MyConsts_PI )
               ! Compute force constant of the distorsion
               Bath%DistorsionForce = Bath%DistorsionForce + Bath%Couplings(iBath)**2 / &
                                                                          ( Bath%OscillatorsMass * Bath%Frequencies(iBath)**2 ) 
               D0 = D0 + Bath%Couplings(iBath)**2
            ENDDO
            D0 = sqrt(D0)

            ! Deallocate memory for spectral density interpolation
            CALL DisposeSpline ( SpectralDensitySpline )
            DEALLOCATE( RdFreq, RdSpectralDens )

      ENDIF

      ! Default coupling is linear
      Bath%NonLinearCoup = .FALSE.

      ! Module is setup
      Bath%BathIsSetup = .TRUE.

#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " Independent oscillator model potential has been setup"
      WRITE(__LOG_UNIT,*) " ... (details) ... "
      WRITE(__LOG_UNIT,*) " Distorsion frequency coefficient (atomic units): ", Bath%DistorsionForce
      WRITE(__LOG_UNIT,*) " D0 (atomic units): ", D0
      WRITE(__LOG_UNIT,*) " Oscillators mass (atomic units): ", Bath%OscillatorsMass
      __CLOSE_LOG_FILE
#endif

   END SUBROUTINE SetupIndepOscillatorsModel


!===============================================================================================================

!*******************************************************************************
!                       SetupOhmicIndepOscillatorsModel
!*******************************************************************************
!>  This subroutine can be used to setup the parameters of the oscillator bath
!>  set the frequency and couplings for an ohmic spectral density with given gamma
!> 
!> @param    
!*******************************************************************************
   SUBROUTINE SetupOhmicIndepOscillatorsModel( Bath, N, SetBathType, SysMassTimeGamma, BathMass, CutOffFreq )
      IMPLICIT NONE
      TYPE(BathData)             :: Bath
      INTEGER, INTENT(IN)        :: N
      INTEGER, INTENT(IN)        :: SetBathType
      REAL, INTENT(IN)           :: SysMassTimeGamma
      REAL, INTENT(IN)           :: BathMass
      REAL, INTENT(IN)           :: CutOffFreq

      REAL :: SpectralDens, D0
      INTEGER :: iBath

      ! If data already setup give a warning and deallocate memory
      CALL WARN( Bath%BathIsSetup, "IndependentOscillatorsModel.SetupOhmicIndepOscillatorsModel: overwriting bath data" )
      IF (Bath%BathIsSetup) CALL DisposeIndepOscillatorsModel( Bath )

      ! Store number of bath degrees of freedom
      Bath%BathSize = N
      ! Set the type of bath
      Bath%BathType = SetBathType
      ! Set the mass of the oscillators
      Bath%OscillatorsMass = BathMass

      ! Allocate memory
      ALLOCATE( Bath%Frequencies(Bath%BathSize), Bath%Couplings(Bath%BathSize) )

      IF ( Bath%BathType == CHAIN_BATH ) THEN

            CALL AbortWithError( "SetupOhmicIndepOscillatorsModel: chain bath not implemented" )

      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN

            ! Set cutoff frequency and frequency spacing
            Bath%CutOff = CutOffFreq
            Bath%DeltaOmega = Bath%CutOff / real( Bath%BathSize)

            ! initialize distortion constant
            Bath%DistorsionForce = 0.0

            ! Compute frequencies and couplings 
            D0 = 0.0
            DO iBath = 1,  Bath%BathSize
               Bath%Frequencies(iBath) = iBath * Bath%DeltaOmega
               SpectralDens = SysMassTimeGamma * Bath%Frequencies(iBath)
               Bath%Couplings(iBath) = SQRT( 2.0 * Bath%OscillatorsMass * Bath%Frequencies(iBath) * Bath%DeltaOmega *      & 
                    SpectralDens / MyConsts_PI )
               ! Compute force constant of the distorsion
               Bath%DistorsionForce = Bath%DistorsionForce + Bath%Couplings(iBath)**2 / &
                                                                          ( Bath%OscillatorsMass * Bath%Frequencies(iBath)**2 ) 
               D0 = D0 + Bath%Couplings(iBath)**2
            ENDDO
            D0 = sqrt(D0)

      ENDIF

      ! Module is setup
      Bath%BathIsSetup = .TRUE.



#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " Independent oscillator model potential has been setup"
      WRITE(__LOG_UNIT,*) " ... (details) ... "
      WRITE(__LOG_UNIT,*) " Distorsion frequency coefficient (atomic units): ", Bath%DistorsionForce
      WRITE(__LOG_UNIT,*) " D0 (atomic units): ", D0
      WRITE(__LOG_UNIT,*) " Oscillators mass (atomic units): ", Bath%OscillatorsMass
      __CLOSE_LOG_FILE
#endif

   END SUBROUTINE SetupOhmicIndepOscillatorsModel

!===============================================================================================================

!*******************************************************************************
!                       SetNonLinearCoupling
!*******************************************************************************
!>  Setup non linear coupling between the system and the bath.
!>  Hcoup = D * (exp(-alpha X)-1) / alpha * s
!>  where D = sqrt(sum c_k^2)/sqrt(mass) and X = sqrt(mass)*sum c_k q_k
!> 
!> @param    
!*******************************************************************************
   SUBROUTINE SetNonLinearCoupling( Bath, Alpha )
      IMPLICIT NONE
      TYPE(BathData)   :: Bath
      REAL, INTENT(IN) :: Alpha
      INTEGER :: iBath

      CALL ERROR( Bath%BathType == CHAIN_BATH, &
            "IndependentOscillatorsModel.SetNonLinearCoupling: exp coupling only with normal bath" )

      ! Compute D0 and OmegaBar
      Bath%NonLinCoup_D0 = 0.0
      Bath%NonLinCoup_OmegaBar = 0.0
      DO iBath = 1, Bath%BathSize
         Bath%NonLinCoup_D0 = Bath%NonLinCoup_D0 + Bath%Couplings(iBath)**2
         Bath%NonLinCoup_OmegaBar = Bath%NonLinCoup_OmegaBar + Bath%Couplings(iBath)**2 / Bath%Frequencies(iBath)**2
      END DO
      Bath%NonLinCoup_D0 = SQRT( Bath%NonLinCoup_D0 )
      Bath%NonLinCoup_OmegaBar = Bath%NonLinCoup_D0 / SQRT(Bath%NonLinCoup_OmegaBar)

      ! Store alpha
      Bath%NonLinCoup_Alpha = Alpha

      ! Compute beta coefficient
      Bath%NonLinCoup_Beta = Bath%NonLinCoup_D0 / Bath%OscillatorsMass / Bath%NonLinCoup_OmegaBar**2

      ! Set non linear coupling
      Bath%NonLinearCoup = .TRUE.

#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " Non-linear coupling in the bath coords has been setup"
      WRITE(__LOG_UNIT,*) " Alpha (atomic units):     ", Bath%NonLinCoup_Alpha
      WRITE(__LOG_UNIT,*) " D0 (atomic units):        ", Bath%NonLinCoup_D0
      WRITE(__LOG_UNIT,*) " Omega bar (atomic units): ", Bath%NonLinCoup_OmegaBar
      WRITE(__LOG_UNIT,*) " Beta (atomic units):      ", Bath%NonLinCoup_Beta
      __CLOSE_LOG_FILE
#endif

   END SUBROUTINE SetNonLinearCoupling


!===============================================================================================================



!*******************************************************************************
!                     BathPotentialAndForces
!*******************************************************************************
!> Compute the potential and derivatives of a generic bath, considering  
!> the potential terms of the bath and the coupling energy term.
!> The output variables are incremented (and not initialized to zero) to 
!> allow the use of many baths coupled to the same system.
!> The opposite of the derivatives are computed, for convenience (the forces are
!> then computed simply by dividing for the masses)
!>
!> @param  Bath        Bath data type 
!> @param  QCoupl      Value of the linearly coupled coordinate 
!> @param  QBath       Array with the cartesian coordinates for the system and the bath 
!> @param  V           In output, V is incremented Output potential in atomic units
!> @param  CouplForce  In output, CouplForces is incremented with the derivatives of the coupling potential
!> @param  QForces     In output, array is incremented with the derivatives of the potential (in au)
!*******************************************************************************     
   SUBROUTINE BathPotentialAndForces( Bath, QCoupl, QBath, V, CouplForce, QForces, CouplingV, DTimesEffMode ) 
      IMPLICIT NONE
      TYPE(BathData), INTENT(IN)                :: Bath
      REAL, INTENT(IN)                          :: QCoupl
      REAL, DIMENSION(:), INTENT(IN)            :: QBath
      REAL, INTENT(INOUT)                       :: V, CouplForce
      REAL, DIMENSION(:), INTENT(INOUT)         :: QForces 
      REAL, INTENT(OUT), OPTIONAL               :: CouplingV
      REAL, INTENT(OUT), OPTIONAL               :: DTimesEffMode

      INTEGER :: iBath
      REAL, DIMENSION(:), POINTER :: Dn
      REAL  :: TraslCoord, Expon
      REAL    :: EffMode

      ! Check if the bath is setup
      CALL ERROR( .NOT. Bath%BathIsSetup, "IndependentOscillatorsModel.BathPotentialAndForces: bath is not setup" )

      ! Check the size of the position and forces arrays
      CALL ERROR( size(QBath) /= size(QForces), &
                     "IndependentOscillatorsModel.BathPotentialAndForces: array dimension mismatch" )
      CALL ERROR( size(QBath) /= Bath%BathSize, &
                     "IndependentOscillatorsModel.BathPotentialAndForces: wrong bath size" )

      IF ( Bath%BathType == CHAIN_BATH ) THEN

         Dn => Bath%Couplings(2:Bath%BathSize)

         ! COUPLING POTENTIAL AND DISTORSION
         IF ( PRESENT( CouplingV ) ) THEN
            V = V + 0.5 * Bath%DistorsionForce * QCoupl**2
            CouplingV = CouplingV - Bath%Couplings(1) * QCoupl * QBath(1) 
         ELSE
            V = V - Bath%Couplings(1) * QCoupl * QBath(1) + 0.5 * Bath%DistorsionForce * QCoupl**2
         END IF

         ! POTENTIAL OF THE BATH OSCILLATORS
         DO iBath = 1, Bath%BathSize - 1
            V = V + 0.5 * Bath%OscillatorsMass * ( Bath%Frequencies(iBath) * QBath(iBath) )**2 - &
                     Dn(iBath) * QBath(iBath) * QBath(iBath+1)
         END DO
         V = V + 0.5 * Bath%OscillatorsMass *( Bath%Frequencies(Bath%BathSize) * QBath(Bath%BathSize) )**2

         ! DERIVATIVES OF THE BATH POTENTIAL WITH RESPECT TO BATH COORDINATES
         QForces(1) = QForces(1) - Bath%OscillatorsMass * Bath%Frequencies(1)**2 * QBath(1) + Dn(1) * QBath(2)
         DO iBath = 2, Bath%BathSize-1
            QForces(iBath) = QForces(iBath) - Bath%OscillatorsMass * Bath%Frequencies(iBath)**2 * QBath(iBath) + &
                              Dn(iBath-1) * QBath(iBath-1) + Dn(iBath) * QBath(iBath+1)
         END DO
         QForces(Bath%BathSize) = QForces(Bath%BathSize) - Bath%OscillatorsMass * Bath%Frequencies(Bath%BathSize)**2 * &
                  QBath(Bath%BathSize)  + Dn(Bath%BathSize-1) * QBath(Bath%BathSize-1)

         ! DERIVATIVE OF THE COUPLING WITH RESPECT TO THE BATH POTENTIAL
         QForces(1) = QForces(1) + Bath%Couplings(1) * QCoupl 

         ! DERIVATIVE OF THE COUPLING+DISTORSION CORRECTION W.R.T. THE COUPLING COORDINATE
         CouplForce = CouplForce + Bath%Couplings(1) * QBath(1) - Bath%DistorsionForce * QCoupl

         IF ( PRESENT(DTimesEffMode)) DTimesEffMode = Bath%Couplings(1) * QBath(1)

      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN

         ! POTENTIAL OF THE BATH OSCILLATORS AND DERIVATIVES WITH RESPECT TO BATH COORDS
         EffMode = 0.0
         DO iBath = 1, Bath%BathSize
            V = V + 0.5 * Bath%OscillatorsMass * ( Bath%Frequencies(iBath) * QBath(iBath) )**2
            QForces(iBath) = QForces(iBath) - Bath%OscillatorsMass * ( Bath%Frequencies(iBath) )**2 * QBath(iBath) &
                                   + QCoupl * Bath%Couplings(iBath)
            EffMode = EffMode + Bath%Couplings(iBath) * QBath(iBath)
         END DO

         IF ( PRESENT( CouplingV ) ) THEN
            V = V + 0.5 * Bath%DistorsionForce * QCoupl**2
            CouplingV = CouplingV - EffMode * QCoupl
         ELSE
            V = V - EffMode * QCoupl + 0.5 * Bath%DistorsionForce * QCoupl**2
         END IF
         CouplForce = CouplForce + EffMode - Bath%DistorsionForce * QCoupl

         IF ( PRESENT(DTimesEffMode)) DTimesEffMode = EffMode

         IF ( Bath%NonLinearCoup ) THEN

            ! Compute some useful quantities
            EffMode = EffMode / Bath%NonLinCoup_D0
            TraslCoord = EffMode - Bath%NonLinCoup_Beta * QCoupl
            Expon = EXP( Bath%NonLinCoup_Alpha * TraslCoord )
            ! Add potential correction
            IF ( PRESENT( CouplingV ) ) THEN 
               CouplingV = CouplingV + 0.5 * Bath%OscillatorsMass * Bath%NonLinCoup_OmegaBar**2 * &
                                        ( (Expon-1.0)**2/Bath%NonLinCoup_Alpha**2 - TraslCoord**2 )
            ELSE
               V = V  + 0.5 * Bath%OscillatorsMass * Bath%NonLinCoup_OmegaBar**2 * &
                                        ( (Expon-1.0)**2/Bath%NonLinCoup_Alpha**2 - TraslCoord**2 )
            END IF
            ! Add force correction to the system force
            CouplForce = CouplForce + Bath%OscillatorsMass * Bath%NonLinCoup_OmegaBar**2 * &
                       Bath%NonLinCoup_Beta * ( Expon * (Expon-1.0)/Bath%NonLinCoup_Alpha - TraslCoord )
            ! Add force correction to the bath oscillators force
            DO iBath = 1, Bath%BathSize
               QForces(iBath) = QForces(iBath) - Bath%OscillatorsMass * Bath%NonLinCoup_OmegaBar**2 * &
                   Bath%Couplings(iBath) * ( Expon * (Expon-1.0)/Bath%NonLinCoup_Alpha - TraslCoord )
            END DO
         END IF

      END IF

   END SUBROUTINE BathPotentialAndForces


!===============================================================================================================

!*******************************************************************************
!               BathOfRingsThermalConditions
!*******************************************************************************
!> Initial conditions for a bath of ring polymers, corresponding to a bath in 
!> normal form. The coordinates are sampled from the classical thermal distribution
!> at given T.
!>
!> @param Bath            Bath data type 
!> @param NBeads          Nr of beads of each ring polymer
!> @param BeadsFrequency  Frequency of the interbeads harmonic potential
!> @param Temperature     Temperature of the thermal p and q distribution 
!> @param Q               In output, random initial coordinates of the bath
!> @param V               In output, random initial velocities of the bath
!*******************************************************************************  
   SUBROUTINE BathOfRingsThermalConditions( Bath, NBeads, BeadsFrequency, Temperature, Q, V,  &
                                                                           PotEnergy, KinEnergy, RandomNr, RingNormalModes )
      IMPLICIT NONE
      TYPE(BathData), INTENT(IN)                         :: Bath
      INTEGER, INTENT(IN)                                :: NBeads
      REAL, INTENT(IN)                                   :: BeadsFrequency, Temperature
      REAL, DIMENSION(Bath%BathSize,NBeads), INTENT(OUT) :: Q, V 
      REAL, INTENT(OUT)                                  :: PotEnergy, KinEnergy
      TYPE(RNGInternalState), INTENT(INOUT)              :: RandomNr
      TYPE(FFTHalfComplexType), INTENT(INOUT)            :: RingNormalModes
      REAL, DIMENSION(NBeads)                            :: RingEigenvalues
      REAL, DIMENSION(Bath%BathSize,Bath%BathSize)       :: BathPotentialMatrix, BathNormalModes
      REAL, DIMENSION(Bath%BathSize)                     :: BathEigenvalues
      REAL, DIMENSION(Bath%BathSize,NBeads)              :: NormalQ, NormalV 
      INTEGER :: iBath, iBead
      REAL    :: SigmaV

      ! Setup bath potential matrix
      BathPotentialMatrix(:,:) = 0.0
      IF ( Bath%BathType == CHAIN_BATH ) THEN
         ! Set potential matrix
         DO iBath = 1, Bath%BathSize-1
            BathPotentialMatrix(iBath,iBath+1) = - Bath%Couplings(iBath+1) / Bath%OscillatorsMass
            BathPotentialMatrix(iBath+1,iBath) = - Bath%Couplings(iBath+1) / Bath%OscillatorsMass
            BathPotentialMatrix(iBath,iBath)   =   Bath%Frequencies(iBath)**2
         END DO
         BathPotentialMatrix(Bath%BathSize,Bath%BathSize) = Bath%Frequencies(Bath%BathSize)**2
         ! Diagonalize
         CALL TheOneWithDiagonalization(BathPotentialMatrix, BathNormalModes, BathEigenvalues)

      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN
         ! Set potential matrix and its trivial diagonalization
         BathNormalModes(:,:) = 0.0
         DO iBath = 1, Bath%BathSize
            BathEigenvalues(iBath) = Bath%Frequencies(iBath)**2
            BathPotentialMatrix(iBath,iBath) = Bath%Frequencies(iBath)**2
            BathNormalModes(iBath,iBath) = 1.0
         END DO
      END IF

      ! Compute normal mode frequencies
      DO iBead = 1, NBeads
         RingEigenvalues(iBead) = ( 2.0 * BeadsFrequency * SIN( MyConsts_PI * real(iBead-1) / real(NBeads) ) )**2
      END DO

      ! Set sigma of the maxwell boltzmann distribution
      SigmaV = sqrt( Temperature / Bath%OscillatorsMass )
      ! Initialize kinetic energy and potential energy
      PotEnergy = 0.0
      KinEnergy = 0.0

      DO iBead = 1, NBeads
         DO iBath = 1, Bath%BathSize
            ! Define random coordinates and velocities in the normal modes representation
            NormalQ(iBath,iBead) = GaussianRandomNr(RandomNr) * SigmaV / SQRT( BathEigenvalues(iBath) + RingEigenvalues(iBead) )
            NormalV(iBath,iBead) = GaussianRandomNr(RandomNr) * SigmaV            
            ! Compute kinetic and potential energies in normal mode representation
            PotEnergy = PotEnergy + 0.5 * Bath%OscillatorsMass * &
                                           (BathEigenvalues(iBath) + RingEigenvalues(iBead) ) * NormalQ(iBath,iBead)**2
            KinEnergy = KinEnergy + 0.5 * Bath%OscillatorsMass * NormalV(iBath,iBead)**2
         END DO
      END DO

      ! Transform normal modes of the ring polymer to ring polymer coordinates
      DO iBath = 1, Bath%BathSize
         CALL ExecuteFFT( RingNormalModes, NormalQ(iBath,:), INVERSE_FFT ) 
         CALL ExecuteFFT( RingNormalModes, NormalV(iBath,:), INVERSE_FFT ) 
      END DO

      ! Transform normal modes of the bath to original representation
      IF ( Bath%BathType == CHAIN_BATH ) THEN
         Q = TheOneWithMatrixMultiplication( BathNormalModes, NormalQ )
         V = TheOneWithMatrixMultiplication( BathNormalModes, NormalV )
      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN
         Q = NormalQ
         V = NormalV
      END IF

   END SUBROUTINE BathOfRingsThermalConditions

!===============================================================================================================

!*******************************************************************************
!               ThermalEquilibriumBathConditions
!*******************************************************************************
!> Initial conditions for an independent oscillators model in normal or chain form
!> The coordinates are sampled from the classical thermal distribution at given T
!> using an exact canonical distribution for zero coupling.
!> Note that for a linear chain, also the couplings within the chain are assumed
!> to be zero.
!>
!> @param Bath         Bath data type 
!> @param Positions    In output, random initial coordinates of the bath
!> @param Velocities   In output, random initial velocities of the bath
!> @param Temperature  Temperature of the classical distribution 
!*******************************************************************************     
   SUBROUTINE ThermalEquilibriumBathConditions( Bath, Positions, Velocities, Temperature, RandomNr )
      IMPLICIT NONE

      TYPE(BathData), INTENT(IN)            :: Bath
      REAL, DIMENSION(:), INTENT(OUT)       :: Positions, Velocities
      REAL, INTENT(IN)                      :: Temperature
      TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr

      REAL :: SigmaQ, SigmaV
      INTEGER :: iBath
      
      ! Check if the bath is setup
      CALL ERROR( .NOT. Bath%BathIsSetup, "IndependentOscillatorsModel.ThermalEquilibriumBathConditions: bath not setup" )

      ! Check the size of the position and forces arrays
      CALL ERROR( size(Positions) /= size(Velocities), &
                     "IndependentOscillatorsModel.ThermalEquilibriumBathConditions: array dimension mismatch" )
      CALL ERROR( size(Positions) /= Bath%BathSize, &
                     "IndependentOscillatorsModel.ThermalEquilibriumBathConditions: wrong bath size" )

      ! THE OSCILLATORS IN A CORRECT CANONICAL DISTRIBUTION FOR ZERO COUPLING
      SigmaV = sqrt( Temperature / Bath%OscillatorsMass )
      DO iBath = 1, Bath%BathSize
         SigmaQ = SigmaV / Bath%Frequencies(iBath)
         Positions(iBath) = GaussianRandomNr(RandomNr) * SigmaQ
         Velocities(iBath) = GaussianRandomNr(RandomNr) * SigmaV
      END DO
   
   END SUBROUTINE ThermalEquilibriumBathConditions


!===============================================================================================================


!*******************************************************************************
!                     ZeroKelvinBathConditions
!*******************************************************************************
!> Setup initial conditions for the system plus bath for a simulation of 
!> vibrational relaxation. The bath can be fixed in the equilibrium position
!> with no momentum ( classical 0K ) or in a quasiclassical state with the 
!> quantum initial zero point energy. 
!> Data are initialized in ATOMIC UNITS.
!>
!> @param Bath         Bath data type 
!> @param Positions    In output, random initial coordinates of the bath
!> @param Velocities   In output, random initial velocities of the bath
!*******************************************************************************     
   SUBROUTINE ZeroKelvinBathConditions( Bath, Positions, Velocities, ZeroPointEnergy, RandomNr )
      IMPLICIT NONE

      TYPE(BathData), INTENT(IN)            :: Bath
      REAL, DIMENSION(:), INTENT(OUT)       :: Positions, Velocities
      LOGICAL, INTENT(IN)                   :: ZeroPointEnergy
      TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
      INTEGER :: iCoord
      REAL :: Value

      ! Check if the bath is setup
      CALL ERROR( .NOT. Bath%BathIsSetup, "IndependentOscillatorsModel.ZeroKelvinBathConditions: bath not setup" )

      ! Check the size of the position and forces arrays
      CALL ERROR( size(Positions) /= size(Velocities), &
                     "IndependentOscillatorsModel.ZeroKelvinBathConditions: array dimension mismatch" )
      CALL ERROR( size(Positions) /= Bath%BathSize, &
                     "IndependentOscillatorsModel.ZeroKelvinBathConditions: wrong bath size" )

      IF ( .NOT. ZeroPointEnergy ) THEN
         Positions(:) = 0.0
         Velocities(:) = 0.0

      ELSE IF ( ZeroPointEnergy ) THEN

         IF ( Bath%BathType == CHAIN_BATH ) THEN
            CALL ShowWarning( "ZeroKelvinBathConditions: quasi-classical conditions not implemented yet" )
            Positions(:) = 0.0
            Velocities(:) = 0.0
         ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN
            DO iCoord = 1, Bath%BathSize
               Value = UniformRandomNr( RandomNr )*2.0*MyConsts_PI
               Positions(iCoord) = COS(Value)/SQRT(Bath%OscillatorsMass*Bath%Frequencies(iCoord))
               Velocities(iCoord) = SIN(Value)*SQRT(Bath%Frequencies(iCoord)/Bath%OscillatorsMass)
            END DO
         END IF

      END IF

   END SUBROUTINE ZeroKelvinBathConditions


!===============================================================================================================


!*******************************************************************************
!                     EnergyOfTheBath
!*******************************************************************************
!> Compute the potential energy of the bath, the coupling energy and optionally
!> also the first effective mode (not normalized).
!>
!> @param Bath                   Bath data type 
!> @param QCoupl                 Coupling coordinates (shifted if necessary)
!> @param QBath                  Bath coordinates
!> @param VCoupling              After execution, stores the coupling energy
!> @param VCoupling              After execution, stores the energy of the bath
!> @param FirstEffectiveMode     (OPTIONAL) After execution, 
!>                               stores the 1st effective mode, not normalized
!*******************************************************************************     
   SUBROUTINE EnergyOfTheBath( Bath, QCoupl, QBath, VCoupling, VBath, FirstEffectiveMode )
      IMPLICIT NONE
      TYPE(BathData), INTENT(IN)     :: Bath
      REAL, INTENT(IN)               :: QCoupl 
      REAL, INTENT(IN), DIMENSION(:) :: QBath 
      REAL, INTENT(OUT)              :: VCoupling, VBath
      REAL, OPTIONAL, INTENT(OUT)    :: FirstEffectiveMode
      REAL :: Coupl
      INTEGER :: iBath
      REAL    :: EffMode
      REAL  :: TraslCoord, Expon

      ! Check if the bath is setup
      CALL ERROR( .NOT. Bath%BathIsSetup, "IndependentOscillatorsModel.EnergyOfTheBath: bath is not setup" )
      ! Check the size of the position array
      CALL ERROR( size(QBath) /= Bath%BathSize, "IndependentOscillatorsModel.EnergyOfTheBath: wrong bath size" )

      VCoupling = 0.0
      VBath = 0.0

      IF ( Bath%BathType == CHAIN_BATH ) THEN

         ! Bath-system coupling
         VCoupling = - Bath%Couplings(1) * QCoupl * QBath(1)
         ! Potential of the bath
         DO iBath = 1, Bath%BathSize - 1
            VBath = VBath + 0.5 * Bath%OscillatorsMass * ( Bath%Frequencies(iBath) * QBath(iBath) )**2 - &
                     Bath%Couplings(iBath+1) * QBath(iBath) * QBath(iBath+1)
         END DO
         VBath = VBath + 0.5 * Bath%OscillatorsMass * ( Bath%Frequencies(Bath%BathSize) * QBath(Bath%BathSize) )**2
         VBath = VBath + 0.5 * Bath%DistorsionForce * QCoupl**2
         IF( PRESENT( FirstEffectiveMode ) ) FirstEffectiveMode = Bath%Couplings(1)*QBath(1)

      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN

         Coupl = 0.0
         DO iBath = 1, Bath%BathSize
            Coupl = Coupl +  Bath%Couplings(iBath) * QBath(iBath)
            VBath = VBath + 0.5 * Bath%OscillatorsMass * ( Bath%Frequencies(iBath) * QBath(iBath) )**2
         END DO
         VCoupling = - Coupl * QCoupl 
         VBath = VBath + 0.5 * Bath%DistorsionForce * QCoupl**2
         IF( PRESENT( FirstEffectiveMode ) ) FirstEffectiveMode = Coupl

      END IF

      IF ( Bath%NonLinearCoup ) THEN
         ! Compute some useful quantities
         EffMode = Coupl / Bath%NonLinCoup_D0
         TraslCoord = EffMode - Bath%NonLinCoup_Beta * QCoupl
         Expon = EXP( Bath%NonLinCoup_Alpha * TraslCoord )
         ! Add potential correction
         VCoupling = VCoupling + 0.5 * Bath%OscillatorsMass * Bath%NonLinCoup_OmegaBar**2 * &
                                  ( (Expon-1.0)**2/Bath%NonLinCoup_Alpha - TraslCoord**2 )
      END IF

   END SUBROUTINE EnergyOfTheBath


!===============================================================================================================

!*******************************************************************************
!                     GetDistorsionForce
!*******************************************************************************
!> Give the distorsion force, stored in the bath data type.
!>
!> @param   Bath     Bath data type 
!> @result  Dist     Distorsion force of the bath     
!*******************************************************************************     
   REAL FUNCTION GetDistorsionForce( Bath ) RESULT( Dist)
      IMPLICIT NONE
      TYPE(BathData), INTENT(IN)     :: Bath

      Dist = Bath%DistorsionForce
   END FUNCTION GetDistorsionForce

!===============================================================================================================

!*******************************************************************************
!                     GetFirstEffectiveMode
!*******************************************************************************
!> Give the distorsion force, stored in the bath data type.
!>
!> @param   Bath     Bath data type 
!> @result  Coords   Coordinates of the bath
!*******************************************************************************     
   REAL FUNCTION GetFirstEffectiveMode( Bath, QBath )
      IMPLICIT NONE
      TYPE(BathData), INTENT(IN)     :: Bath
      REAL, INTENT(IN), DIMENSION(:) :: QBath 
      REAL :: Coupl, DSq
      INTEGER  :: iBath

      ! Check if the bath is setup
      CALL ERROR( .NOT. Bath%BathIsSetup, "IndependentOscillatorsModel.GetFirstEffectiveMode: bath is not setup" )
      ! Check the size of the position array
      CALL ERROR( size(QBath) /= Bath%BathSize, "IndependentOscillatorsModel.GetFirstEffectiveMode: wrong bath size" )

      IF ( Bath%BathType == CHAIN_BATH ) THEN

         GetFirstEffectiveMode = Bath%Couplings(1)*QBath(1)

      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN

         Coupl = 0.0
         DSq = 0.0
         DO iBath = 1, Bath%BathSize
            Coupl = Coupl +  Bath%Couplings(iBath) * QBath(iBath)
            DSq   = DSq   +  Bath%Couplings(iBath)**2
         END DO
         GetFirstEffectiveMode = Coupl / SQRT(DSq)

      END IF

   END FUNCTION GetFirstEffectiveMode

!===============================================================================================================

!*******************************************************************************
!             CouplingAndDistortionHessian
!*******************************************************************************
!> Give the contribution to the hessian which is due to the coupling
!> potential and the distortion force.
!>  NOTE THAT RESULTS ARE NOT MASS SCALED
!>
!> @param   Bath     Bath data type 
!> @result  Dist     Distortion force of the bath     
!*******************************************************************************     
   SUBROUTINE CouplingAndDistortionHessian( Bath, CouplingHessian, DistortionHessian ) 
      IMPLICIT NONE
      TYPE(BathData), INTENT(IN)                    :: Bath
      REAL, DIMENSION( Bath%BathSize ), INTENT(OUT) :: CouplingHessian
      REAL, INTENT(OUT)                             :: DistortionHessian

      DistortionHessian = Bath%DistorsionForce
      CouplingHessian(:) = 0.0

      IF ( Bath%BathType == CHAIN_BATH ) THEN
         CouplingHessian(1) = - Bath%Couplings(1)
      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN
         CouplingHessian(:) = - Bath%Couplings(:)
      END IF

   END SUBROUTINE CouplingAndDistortionHessian


!===============================================================================================================

!*******************************************************************************
!                     HessianOfTheBath
!*******************************************************************************
!> Compute the Hessian of the potential of the bath only
!> divided by the masses to work with mass weighted coordianates
!>
!> @param Bath       Bath data type 
!> @param Hessian    output matrix of size NxN, where N is the size of the bath
!*******************************************************************************     
   SUBROUTINE HessianOfTheBath( Bath, Hessian ) 
      IMPLICIT NONE
      TYPE(BathData), INTENT(IN)  :: Bath
      REAL, DIMENSION(Bath%BathSize, Bath%BathSize), INTENT(INOUT) :: Hessian
      INTEGER :: iBath

      ! Check if the bath is setup
      CALL ERROR( .NOT. Bath%BathIsSetup, " IndependentOscillatorsModel.HessianOfTheBath: bath is not setup" )

      ! Initialize
      Hessian(:,:) = 0.0

      IF ( Bath%BathType == CHAIN_BATH ) THEN
         ! Diagonal elements: quadratic terms of the potential
         DO iBath = 1, Bath%BathSize
            Hessian(iBath,iBath) = Bath%Frequencies(iBath)**2
         END DO
         ! off-diagonal elements: couplings
         DO iBath = 1, Bath%BathSize-1
            Hessian(iBath,iBath+1) = - Bath%Couplings(iBath+1) / Bath%OscillatorsMass
            Hessian(iBath+1,iBath) = - Bath%Couplings(iBath+1) / Bath%OscillatorsMass
         END DO

      ELSE IF ( Bath%BathType == STANDARD_BATH ) THEN
         ! Diagonal elements: quadratic terms of the potential
         DO iBath = 1, Bath%BathSize
            Hessian(iBath,iBath) = Bath%Frequencies(iBath)**2
         END DO
      END IF

   END SUBROUTINE HessianOfTheBath


!===============================================================================================================


!*******************************************************************************
! DisposeIndepOscillatorsModel
!*******************************************************************************
!>  Deallocate memory used by the Bath data type

!> @param Bath       Bath data type 
!*******************************************************************************
   SUBROUTINE DisposeIndepOscillatorsModel( Bath )
      IMPLICIT NONE
      TYPE(BathData), INTENT(INOUT)    :: Bath

      ! Do nothing if bath is not yet setup
      IF (.NOT. Bath%BathIsSetup)  RETURN

      ! deallocate memory
      DEALLOCATE( Bath%Frequencies, Bath%Couplings )

      ! Set linear coupling
      Bath%NonLinearCoup = .FALSE.

      ! Module is no longer setup
      Bath%BathIsSetup = .FALSE.

#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " Independent oscillator model data has been disposed"
      __CLOSE_LOG_FILE
#endif

   END SUBROUTINE DisposeIndepOscillatorsModel


!===============================================================================================================

END MODULE IndependentOscillatorsModel
