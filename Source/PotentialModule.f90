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
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetupPotential                                        !< setup subroutine
   PUBLIC :: GetXLabel, GetSystemDimension, PESIsCollinear         !< info subroutines
   PUBLIC :: GetPotential, GetPotAndForces                         !< get potential and forces
   PUBLIC :: GetVPartitions                                        !< get pot energy partitioned according to some relevant scheme
   PUBLIC :: StartSystemForScattering, GetInitialAsymptoteMask     !< system initial conditions subroutines
   PUBLIC :: GetNrChannels, GetChannelLabel, GetCurrentChannel     !< trajectory analysis
   
   !> Setup variable for the potential
   LOGICAL, SAVE :: PotentialModuleIsSetup = .FALSE.

   !> Number of dimensions of the potential
   INTEGER, SAVE :: NDim

   !> Labels of the potential coordinates
   CHARACTER(5), DIMENSION(:), ALLOCATABLE, SAVE :: CoordLabels

   !> Step for computing the numerical derivatives with finite differences
   REAL, PARAMETER :: SmallDelta = 0.00001
   
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

   ! Reference geometry values of the potential
   REAL, PARAMETER  ::  Minimum_ZHTar = 2.7955073359492473  !< equilibrium cartesian Z of HTargon in initial asymptote / bohr
   REAL, PARAMETER  ::  Minimum_ZCarb = 0.69156392247796172 !< equilibrium cartesian Z of Carbon in initial asymptote / bohr
   REAL, PARAMETER  ::  H2EqD = 1.51178089965204958905       !< equilibrium distance of H2 in final asymptote / bohr

   
   
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
!> Setup appropriate initial condition for the subsystem.
!>
!> @param X              Output vector with the initial positions of the system
!> @param V              Output vector with the initial Velocities of the system
!> @param M              Input vector with the masses of the system coordinates
!> @param InitDist       Input real value of the scatterer distace
!> @param InitEKin       Input real value of the scatterer kinetic energy
!> @param ImpactPar      Input real value of the impact parameter
!> @param Task           Input integer to define the task of the initialization
!**************************************************************************************
      SUBROUTINE StartSystemForScattering( X, V, M, InitDist, InitEKin, ImpactPar, Task )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(OUT)       :: X
         REAL, DIMENSION(size(X)), INTENT(OUT) :: V
         REAL, DIMENSION(size(X)), INTENT(IN)  :: M
         REAL                                  :: InitDist, InitEKin, ImpactPar
         INTEGER                               :: Task

         INTEGER, PARAMETER   ::  SCATTERING       = 1
         INTEGER, PARAMETER   ::  ASYMPTOTIC_START = 2
        
         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.StartSystemForScattering : module not set" )
         ! Check the number of degree of freedom
         CALL ERROR( size(X) /= NDim, "PotentialModule.StartSystemForScattering: input array dimension mismatch" )

         SELECT CASE ( Task )
         
            ! In case of a ASYMPTOTIC_START task, the scatterer is placed far from the target, with no incoming velocity
            ! the substrate is placed in equilibrium geometry
            CASE ( ASYMPTOTIC_START )       
                X(1) = InitDist
                X(2) = Minimum_ZHTar
                IF (NDim == 3) X(3) = Minimum_ZCarb
                
            ! In case of a SCATTERING task, the scatterer is placed at the input initial conditions 
            ! the substrate is left untouched
            CASE ( SCATTERING )
               X(1) = InitDist
               V(1) = -SQRT( 2.0 * InitEKin / M(1) )
               
            CASE DEFAULT
               CALL AbortWithError( "PotentialModule.StartSystemForScattering : Task is not defined" )
               
         END SELECT

      END SUBROUTINE StartSystemForScattering

      
!===============================================================================================================================

!**************************************************************************************
!> Function that returns a logical mask for initial asymptote optimization 
!> @returns   Logical vector with .FALSE. corresponding to the scattering coordinate
!**************************************************************************************
      FUNCTION GetInitialAsymptoteMask( ) RESULT(Mask)
         IMPLICIT NONE
         LOGICAL, DIMENSION(NDim) :: Mask
         
         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GetInitialAsymptoteMask : Module not Setup" )

         Mask(1) = .FALSE.
         Mask(2) = .TRUE.
         IF ( VReducedDim ==  FULLPOT )  Mask(3) = .TRUE.
 
      END FUNCTION GetInitialAsymptoteMask
      
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
            CALL ER_3D( Positions(1), Positions(2), Minimum_ZCarb, V, Dummy(1), Dummy(2), Dummy(3) )

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
            zC = Minimum_ZCarb
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
            CALL ER_3D( Positions(1), Positions(2), Minimum_ZCarb, V, Forces(1), Forces(2), FirstDer )

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


END MODULE PotentialModule

