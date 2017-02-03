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
!>  \version    1.1
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
!>  \arg 10 Feb 2017: Hessian computation and potential optimization have
!>                    been moved outside the module in other general purpose
!>                    modules, which are now used (FiniteDifference, Optimize)
!
!>  \todo N.A.
!
!***************************************************************************************

MODULE PotentialModule
#include "preprocessoptions.cpp"
   USE RandomNumberGenerator
   USE FiniteDifference
   USE Optimize
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetupPotential                                        !< setup subroutine
   PUBLIC :: GetXLabel, GetSystemDimension, PESIsCollinear         !< info subroutines
   PUBLIC :: GetPotential, GetPotAndForces                         !< get potential and forces
   PUBLIC :: GetVPartitions                                        !< get pot energy partitioned according to some relevant scheme
   PUBLIC :: StartSystemForScattering                              !< system initial conditions subroutines
   PUBLIC :: GetNrChannels, GetChannelLabel, GetCurrentChannel     !< trajectory analysis
   
   !> Setup variable for the potential
   LOGICAL, SAVE :: PotentialModuleIsSetup = .FALSE.

   !> Number of dimensions of the potential
   INTEGER, SAVE :: NDim

   !> Labels of the potential coordinates
   CHARACTER(5), DIMENSION(:), ALLOCATABLE, SAVE :: CoordLabels

   !> Step for computing the numerical derivatives with finite differences
   REAL, PARAMETER :: SmallDelta = 0.00001

   !> Logical control variable for the log of subroutine StartSystemForScattering
   LOGICAL, SAVE  :: LogScattInit = .TRUE.
   
#if defined(LOG_FILE)
   CHARACTER(17), SAVE :: LogStr = " SetupPotential |"
#endif

   !> Variable to define the type of reduced dimensionality calculation
   INTEGER, SAVE              :: VReducedDim
   INTEGER, PARAMETER, PUBLIC :: FULLPOT    = 0,  & ! potential used in its full dimensionality
                                 SUDDEN     = 1,  & ! vibrational dof is frozen in initial position
                                 ADIABATIC  = 2     ! vibrational dof adiabatically adjusts to other degrees of freedom

   !> Variable to temprarary store optimize value of zc for adiabatic potential
   REAL, SAVE  :: OptZc = 0.0

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
!> @param ReducedDim    Integer variable to set adiabatic/sudden approx
!**************************************************************************************
      SUBROUTINE SetupPotential( ReducedDim, DerivTest )
         IMPLICIT NONE
         INTEGER, INTENT(IN)           :: ReducedDim
         LOGICAL, INTENT(IN), OPTIONAL :: DerivTest
         REAL, DIMENSION(:), ALLOCATABLE :: CoordMin, CoordMax         !< Coordinate intervals where to check the derivatives
         INTEGER :: i
         
         ! exit if module is setup
         IF ( PotentialModuleIsSetup ) RETURN
 
         ! Potential module is set up
         PotentialModuleIsSetup = .TRUE.
 
         ! Store the info on static approximation of the potential
         VReducedDim = ReducedDim

         ! Set the number of dimensions
         IF ( VReducedDim ==  FULLPOT ) THEN
            NDim = 3
         ELSE IF ( VReducedDim == SUDDEN .OR. VReducedDim == ADIABATIC ) THEN
            NDim = 2
         ELSE
            CALL AbortWithError( "PotentialModule.SetupPotential: wrong ReducedDim value " )
         END IF

         ! Set the labels of the coordinates
         ALLOCATE( CoordLabels(NDim) )
         CoordLabels(1) = "Hi_z"
         CoordLabels(2) = "Ht_z"
         IF ( VReducedDim ==  FULLPOT )  CoordLabels(3) = "C_z"

#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,"(/,A,A,/)")  LogStr," System potential has been setup"
      IF ( VReducedDim ==  FULLPOT ) THEN
         WRITE(__LOG_UNIT,"(A,A,/)") LogStr," Setup PES is a 3D pot for Eley-Rideal abstraction of H2 on a graphene surface."
      ELSE IF ( VReducedDim == SUDDEN ) THEN
         WRITE(__LOG_UNIT,"(A,A,/)") LogStr," Setup PES is a 2D sudden pot for Eley-Rideal abstraction of H2 on a graphene surface."
      ELSE IF ( VReducedDim == ADIABATIC ) THEN
         WRITE(__LOG_UNIT,"(A,A,/)") LogStr," Setup PES is a 2D adiab pot for Eley-Rideal abstraction of H2 on a graphene surface."
      END IF
      WRITE(__LOG_UNIT,"(A,A,I2,A)") LogStr," PES is a function  of ",NDim," coordinates:"
      DO i = 1, NDim 
         WRITE(__LOG_UNIT,"(A,A,I2,A,A)") LogStr," * Coord. ",i," - ",CoordLabels(i)
      END DO
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
               CALL TestAnalyticForces( GetPotAndForces, CoordMin, CoordMax, SmallDelta ) 

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
         CALL ERROR( (iCoord > NDim) .OR. (iCoord < 1), "PotentialModule.GetXLabel : wrong coordinate number" )
         Label = CoordLabels(iCoord)
 
      END FUNCTION GetXLabel

!===============================================================================================================================

!**************************************************************************************
!> Function that returns the number of dimensions of the system.
!> @returns   Nr of degrees of freedom on which the system potential depends
!**************************************************************************************
      INTEGER FUNCTION GetSystemDimension( )
         IMPLICIT NONE

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetSystemDimension : Module not Setup" )
         GetSystemDimension = NDim
 
      END FUNCTION GetSystemDimension

!===============================================================================================================================

!**************************************************************************************
!> Inquire whether the potential is collinear or not
!> @returns   Logical value: .TRUE. for a collinear PES, .FALSE. otherwise
!**************************************************************************************
      LOGICAL FUNCTION PESIsCollinear( )
         IMPLICIT NONE

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.PESIsCollinear : Module not Setup" )
         PESIsCollinear = .TRUE.
 
      END FUNCTION PESIsCollinear

!===============================================================================================================================

!**************************************************************************************
!> Function that returns the number of final channels in the trajectory analysis
!> @returns   Nr of channels
!**************************************************************************************
      INTEGER FUNCTION GetNrChannels( )
         IMPLICIT NONE

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetNrChannels : Module not Setup" )
         GetNrChannels = 4
 
      END FUNCTION GetNrChannels

!===============================================================================================================================

!**************************************************************************************
!> Function that returns a string describing the i-th channel, where i is the integer input
!> @param ChannelIdNumber  ordinal number of the channel
!> @returns                string describing the ChannelIdNumber-th channel
!**************************************************************************************
      CHARACTER(20) FUNCTION GetChannelLabel( ChannelIdNumber )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ChannelIdNumber

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetChannelLabel : Module not Setup" )

         SELECT CASE ( ChannelIdNumber )
            CASE ( 0 )
               GetChannelLabel = "interaction"
            CASE ( 1 )
               GetChannelLabel = "reflection"
            CASE ( 2 )
               GetChannelLabel = "reaction"
            CASE ( 3 )
               GetChannelLabel = "cid"
            CASE DEFAULT
               CALL AbortWithError( "PotentialModule.GetChannelLabel : channel ID does not exist" )
         END SELECT
 
      END FUNCTION GetChannelLabel

!===============================================================================================================================

!**************************************************************************************
!> Function that returns the current channel for the given coordinate values
!> @param Position  given position
!> @returns         current channel corresponding to Position
!**************************************************************************************
      INTEGER FUNCTION GetCurrentChannel( Position )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(IN) :: Position

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetCurrentChannel : Module not Setup" )
         CALL ERROR( size(Position) > NDim .OR. size(Position) < 1, "PotentialModule.GetCurrentChannel : wrong coordinate number" )

         IF ( Position(1) > 20. .AND. Position(2) < 5. ) THEN
            GetCurrentChannel = 1         ! projectile is reflected to the gas phase, target still bound
         ELSE IF ( Position(1) > 20. .AND. Position(2) > 20. .AND. abs(Position(1)-Position(2)) < 4.0 ) THEN
            GetCurrentChannel = 2         ! projectile and target are reflected to the gas phase in a bound state
         ELSE IF ( Position(1) > 20. .AND. Position(2) > 20. .AND. abs(Position(1)-Position(2)) > 4.0 ) THEN
            GetCurrentChannel = 3         ! projectile and target are reflected to the gas phase, without being in a bound state
         ELSE
            GetCurrentChannel = 0
         END IF
 
      END FUNCTION GetCurrentChannel

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
         REAL, DIMENSION(:), INTENT(IN)  :: Positions
         REAL, DIMENSION(3) :: Dummy
         REAL               :: FirstDer, SecDer, DeltaF
         INTEGER            :: NIter, k

         REAL, DIMENSION(4), PARAMETER :: Deltas = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
         REAL, DIMENSION(4), PARAMETER :: Coeffs = (/ +1./12., -8./12., +8./12., -1./12. /) 
         INTEGER :: NMaxIter = 10000
         REAL    :: GradThresh = 0.00001

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
         IF ( VReducedDim ==  FULLPOT ) THEN
           CALL ER_3D( Positions(1), Positions(2), Positions(3), V, Dummy(1), Dummy(2), Dummy(3) )

         ELSE IF ( VReducedDim == SUDDEN ) THEN
            CALL ER_3D( Positions(1), Positions(2), 0.691576, V, Dummy(1), Dummy(2), Dummy(3) )

         ELSE IF ( VReducedDim == ADIABATIC ) THEN
            CALL AbortWithError(" Adiabatic not yet implemented" )

            ! Optimize carbon position with newton method (numerical sec derivs)
            Iterations: DO NIter = 1, NMaxIter
               CALL ER_3D( Positions(1), Positions(2), OptZc, V, Dummy(1), Dummy(2), FirstDer )
               IF ( ABS(FirstDer) > 0.01 ) THEN
                 OptZc = OptZc - FirstDer
               ELSE 
                 SecDer = 0.0
                 DO k = 1, size(Deltas)
                    CALL ER_3D( Positions(1), Positions(2), OptZc+Deltas(k)*SmallDelta, V, Dummy(1), Dummy(2), DeltaF )
                    SecDer = SecDer + Coeffs(k)*DeltaF
                 END DO
                 SecDer = SecDer / SmallDelta
                 IF (ABS(SecDer) < 0.001) THEN
                    OptZc = OptZc - FirstDer
                 ELSE
                    OptZc = OptZc - FirstDer/SecDer
                 ENDIF
               END IF
               IF ( ABS(FirstDer) < GradThresh ) EXIT Iterations
            END DO Iterations
            CALL ER_3D( Positions(1), Positions(2), OptZc, V, Dummy(1), Dummy(2), FirstDer )
         END IF

      END FUNCTION GetPotential

!===============================================================================================================================

!******************************************************************************
!>  Subroutine to compute special potential expectation values, 
!>  corresponding to a specific partitioning of the potential energy.
!>
!>  @param Positions   array with current coordinates of the system
!>  @returns Expect    array with the potential expectations values
!******************************************************************************
      FUNCTION GetVPartitions( Positions ) RESULT( VPart ) 
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(IN)           :: Positions    
         REAL, DIMENSION(3)                       :: VPart 

         REAL, DIMENSION(3) :: Dummy
         REAL, PARAMETER :: H2EqD = 1.51178089965204958905d0  !< equilibrium distance of H2 / Bohr
         REAL :: zC, E1, E2
         
         INTERFACE
            SUBROUTINE ER_3D (zi, zt, zc, vv, dv_zi, dv_zt, dv_zc)
               REAL*8, INTENT (IN) :: zi, zt, zc
               REAL*8, INTENT(OUT) :: vv, dv_zi, dv_zt, dv_zc
            END SUBROUTINE
         END INTERFACE

         ! Error if module not have been setup yet
         CALL ERROR(.NOT. PotentialModuleIsSetup,"PotentialModule.GetVPartitions: module not set")
         ! Check the number of degree of freedom
         CALL ERROR( size(Positions) /= NDim, "PotentialModule.GetVPartitions: input array dimension mismatch" )

         IF ( VReducedDim ==  FULLPOT ) THEN
            zC = Positions(3)
         ELSE IF ( VReducedDim == SUDDEN ) THEN
            zC = 0.691576
         ELSE IF ( VReducedDim == ADIABATIC ) THEN
            CALL AbortWithError(" Adiabatic not yet implemented" )
         END IF

         ! 3 expectation values are computed:
         ! 1) Compute energy of the carbon for H1 and H2 far from surface at eq position
         CALL ER_3D (100., 100.+H2EqD, zC,  E1, Dummy(1), Dummy(2), Dummy(3)) 
         CALL ER_3D (100., 100.+H2EqD, 0.0, E2, Dummy(1), Dummy(2), Dummy(3)) 
         VPart(1) = E1-E2
         ! 2) Compute energy of H-H far from the surface, for the carbon planar
         CALL ER_3D (Positions(1)+100., Positions(2)+100., 0.0, VPart(2), Dummy(1), Dummy(2), Dummy(3))
         ! 3) Compute energy of the C-H for the other H far from the surface
         CALL ER_3D (100., Positions(2), zC, VPart(3), Dummy(1), Dummy(2), Dummy(3))
         
      END FUNCTION GetVPartitions

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
         REAL, DIMENSION(:), INTENT(IN)    :: Positions     
         REAL, DIMENSION(:), INTENT(OUT)   :: Forces 
         REAL               :: FirstDer, SecDer, DeltaF
         INTEGER            :: NIter, k

         REAL, DIMENSION(4), PARAMETER :: Deltas = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
         REAL, DIMENSION(4), PARAMETER :: Coeffs = (/ +1./12., -8./12., +8./12., -1./12. /) 
         INTEGER :: NMaxIter = 10000
         REAL    :: GradThresh = 0.00001

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
         CALL ERROR( size(Forces) /= NDim, "PotentialModule.GetPotAndForces: output array dimension mismatch" )

         ! Compute energy
         IF ( VReducedDim ==  FULLPOT ) THEN
           CALL ER_3D( Positions(1), Positions(2), Positions(3), V, Forces(1), Forces(2), Forces(3) )

         ELSE IF ( VReducedDim == SUDDEN ) THEN
            CALL ER_3D( Positions(1), Positions(2), 0.691576, V, Forces(1), Forces(2), FirstDer )

         ELSE IF ( VReducedDim == ADIABATIC ) THEN
            CALL AbortWithError(" Adiabatic not yet implemented" )

            ! Optimize carbon position with newton method (numerical sec derivs)
            Iterations: DO NIter = 1, NMaxIter
               CALL ER_3D( Positions(1), Positions(2), OptZc, V, Forces(1), Forces(2), FirstDer )
               IF ( ABS(FirstDer) > 0.01 ) THEN
                 OptZc = OptZc - FirstDer
               ELSE 
                 SecDer = 0.0
                 DO k = 1, size(Deltas)
                    CALL ER_3D( Positions(1), Positions(2), OptZc+Deltas(k)*SmallDelta, V, Forces(1), Forces(2), DeltaF )
                    SecDer = SecDer + Coeffs(k)*DeltaF
                 END DO
                 SecDer = SecDer / SmallDelta
                 IF (ABS(SecDer) < 0.001) THEN
                    OptZc = OptZc - FirstDer
                 ELSE
                    OptZc = OptZc - FirstDer/SecDer
                 ENDIF
               END IF
               IF ( ABS(FirstDer) < GradThresh ) EXIT Iterations
            END DO Iterations
            CALL ER_3D( Positions(1), Positions(2), OptZc, V, Forces(1), Forces(2), FirstDer )
         END IF

         ! Compute forces from derivatives
         Forces(:) = - Forces(:)

      END FUNCTION GetPotAndForces

 
!===============================================================================================================================


!**************************************************************************************
!> Setup initial condition for the system, in the case of a scattering simulation.
!> The initial state can be chosen to be either a finite temperature classical 
!> state or a quasi-classical T=0 K state.
!>
!> @param X              Output vector with the initial positions of the system
!> @param V              Output vector with the initial Velocities of the system
!> @param M              Input vector with the masses of the system coordinates
!> @param InitDist       Input real value of the scatterer distace
!> @param InitEKin       Input real value of the scatterer kinetic energy
!> @param ImpactPar      Input real value of the impact parameter
!> @param Temperature    Input real value of the temperature
!> @param RandomNr       Internal state of the random number generator
!> @param Task           Integer flag to set the task to perform
!**************************************************************************************
      SUBROUTINE StartSystemForScattering( X, V, M, InitDist, InitEKin, ImpactPar, Temperature, RandomNr, Task )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(OUT)       :: X
         REAL, DIMENSION(size(X)), INTENT(OUT) :: V
         REAL, DIMENSION(size(X)), INTENT(IN)  :: M
         REAL                                  :: InitDist, InitEKin, ImpactPar
         REAL                                  :: Temperature
         TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
         INTEGER                               :: Task

         INTEGER, PARAMETER   ::  EQUILIBRATION_CLASSICAL = 1
         INTEGER, PARAMETER   ::  SCATTERING_CLASSICAL    = 2
         INTEGER, PARAMETER   ::  QUASICLASSICAL          = 3
         
         REAL, DIMENSION(size(X),size(X)) ::  Hessian, NormalModesVec
         REAL, DIMENSION(size(X))         ::  NormalModesVal
         INTEGER :: i,j
         REAL    :: Value, CarbonFreq, SigmaV

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.StartSystemForScattering : module not set" )
         ! Check the number of degree of freedom
         CALL ERROR( size(X) /= NDim, "PotentialModule.StartSystemForScattering: input array dimension mismatch" )

         ! In case of a SCATTERING_CLASSICAL task, coordinates of the substrate are left untouched
         IF ( Task == EQUILIBRATION_CLASSICAL .OR. Task ==  QUASICLASSICAL ) THEN
         
            ! SET INITIAL GEOMETRY, OPTIMIZING THE NON-SCATTERING COORDINATES
            X(1) = 10.0000 / MyConsts_Bohr2Ang
            X(2) = 1.45000 / MyConsts_Bohr2Ang
            IF (NDim == 3) X(3) = 0.35000 / MyConsts_Bohr2Ang
            X = NewtonLocator( GetPotAndForces, X, 10**3, 1.E-6, 1.E-6, SmallDelta )
            
            ! COMPUTE HESSIAN, DIVIDE IT BY THE MASSES AND FIND NORMAL MODES
            ! Use subroutine GetHessianFromForces to compute the matrix of the 2nd derivatives
            Hessian(:,:) = GetHessianFromForces( X, GetPotAndForces, SmallDelta ) 
            ! Numerical hessian of the potential in mass weighted coordinates
            DO j = 1, NDim
               DO i = 1, NDim
                  Hessian(i,j) = Hessian(i,j) / SQRT( M(i)*M(j) )
               END DO
            END DO
            ! Diagonalize
            CALL TheOneWithDiagonalization(Hessian, NormalModesVec, NormalModesVal)

            ! WEIGHT STARTING COORDS BY MASS AND TRANSFORM TO THE NORMAL MODES REFERENCE
            DO i = 1, NDim
               X(i) = X(i) * SQRT( M(i) )
            END DO
            X = TheOneWithMatrixVectorProduct( TheOneWithTransposeMatrix(NormalModesVec), X )

            ! SET VELOCITIES AND MOMENTA OF BOUND COORDINATES

            ! cycle over the coordinates
            DO i = 1, NDim

               ! Set sigma of the maxwell boltzmann distribution
               SigmaV = sqrt( Temperature )

               ! the coordinate is bound
               IF ( NormalModesVal(i) > 0. ) THEN
                  ! set frequency of the normal mode
                  CarbonFreq = SQRT( NormalModesVal(i) )
                  ! Classical distribution
                  IF ( Task == EQUILIBRATION_CLASSICAL ) THEN
                     X(i) = X(i) + GaussianRandomNr(RandomNr) * SigmaV / CarbonFreq
                     V(i) = GaussianRandomNr(RandomNr) * SigmaV
                  ! Quasiclassical distribution
                  ELSE IF (  Task ==  QUASICLASSICAL ) THEN
                     Value = UniformRandomNr( RandomNr )*2.0*MyConsts_PI
                     X(i) = X(i) + COS(Value) / SQRT(CarbonFreq)
                     V(i) = SIN(Value) * SQRT(CarbonFreq)
                  END IF

               ! unbound coordinate (set velocity according to maxwell-boltzman) 
               ELSE
                  X(i) = X(i)
                  V(i) = GaussianRandomNr(RandomNr) * SigmaV

               ENDIF
            END DO
               
            ! TRASFORM BACK TO ORIGINAL FRAME AND TO NOT-MASS-WEIGHTED COORDINATES
            X = TheOneWithMatrixVectorProduct( NormalModesVec, X )
            V = TheOneWithMatrixVectorProduct( NormalModesVec, V )
            DO i = 1, NDim
               X(i) = X(i) / SQRT( M(i) )
               V(i) = V(i) / SQRT( M(i) )
            END DO

         END IF 
         
         ! Set initial coord and vel of the scatterer
         SELECT CASE ( Task )
            CASE ( EQUILIBRATION_CLASSICAL )
               X(1) = X(1)
               V(1) = 0.0
            CASE ( QUASICLASSICAL, SCATTERING_CLASSICAL )
               X(1) = InitDist
               V(1) = -SQRT( 2.0 * InitEKin / M(1) )
            CASE DEFAULT
               CALL AbortWithError( "PotentialModule.StartSystemForScattering : Task is not defined" )
         END SELECT

         IF ( LogScattInit ) THEN

              ! PRINT INFO
              ! system starts with random displacement around this minimum geometry and random momenta
              ! computed in the frame of the normal modes
              ! frequencies for the bound states are computed with a normal modes analysis at the minimum
              ! here are the frequencies: 
              ! i-th normal mode is a bound coordinate of frequency tot
              ! j-th degree is a free coordinate of imaginary frequency tot
              ! displacement and momenta along the normal modes follow a classical boltzmann distribution for a harmonic Oscillators
              ! or a quasi-classical distribution for a harmonic oscillator at a classical orbit of energy
              ! given by the ZPE hbar*omega
              
              ! -----------------------
              
!             DO i = 1, NDim
!                IF ( NormalModesVal(i) > 0. ) THEN
!                   PRINT*, SQRT(NormalModesVal(i))*FreqConversion(InternalUnits, InputUnits)
!                ELSE
!                   PRINT*, SQRT(-NormalModesVal(i))*FreqConversion(InternalUnits, InputUnits), " *i"
!                ENDIF
!             END DO

              LogScattInit = .FALSE.
              
          END IF 

      END SUBROUTINE StartSystemForScattering

!===============================================================================================================================

END MODULE PotentialModule

