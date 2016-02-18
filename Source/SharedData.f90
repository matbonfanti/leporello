!***************************************************************************************
!*                              MODULE SharedData
!***************************************************************************************
!
!>  \brief     Common data
!>  \details   This module include the data to be shared by the
!>             all the other modules of the code
!
!***************************************************************************************
MODULE SharedData
#include "preprocessoptions.cpp"
   USE IndependentOscillatorsModel

   IMPLICIT NONE

   PUBLIC          ! This module contain data that is supposed to be shared

!=============================================================================================================

   ! PARAMETERS 
   
   ! Variable to define adsorption condition
   REAL, PARAMETER :: AdsorpLimit = 5.0

!=============================================================================================================

   ! VARIABLES SET FROM INPUT
   
   !> Variable to define which kind of calculation is required 
   INTEGER :: RunType
   INTEGER, PARAMETER :: EQUILIBRIUM      = 1,  & ! Equilibrium calculation with H already adsorbed
                         RELAXATION       = 2,  & ! Relaxation dynamics of a CH bound state, with the bath at 0K
                         SCATTERING       = 4,  & ! Scattering calculation with H coming from gas-phase
                         POTENTIAL        = 10    ! Static analysis of the potential 

   !> Variable to set the print level of the calculation
   INTEGER :: PrintType
   INTEGER, PARAMETER :: DEBUG       = 3,  &   ! fully detailed information about the trajs
                         FULL        = 2,  &   ! files to plot the make animations, averages for each traj
                         MINIMAL     = 1       ! minimal level of output, only final averages

   !> Variable to set the kind of bath included in the dynamics
   INTEGER :: BathType
   INTEGER, PARAMETER :: NORMAL_BATH    = 3, &   ! bath of HO on a regular freq grid, all coupled to the system
                         CHAIN_BATH     = 2, &   ! bath of HO in a linear chain form
                         LANGEVIN_DYN   = 1      ! effective relaxation dynamics, with langevin eq of motion

   !> Gamma of the relaxation during dynamics (its meaning depends on the bath representation)
   REAL    :: DynamicsGamma             

   !> Logical variable to perform the tests on the derivatives of the potential !!! the code is terminated after the tests !!!
   LOGICAL :: DerivTesting

   ! INFORMATION ON THE SYSTEM

   REAL    :: MassHTar                  !< Mass of the target hydrogen atom
   REAL    :: MassHInc                  !< Mass of the incident hydrogen atom
   REAL    :: MassC                     !< Mass of the carbon atoms

   ! THE FOLLOWING INPUT DATA ARE RELEVANT FOR A NORMAL/CHAIN BATH DEFINITION

   INTEGER        :: NBath                     !< Nr of bath degrees of freedom
   REAL           :: MassBath                  !< Mass of the HO in the bath
   REAL           :: BathCutOffFreq            !< cutoff frequency of the bath
   REAL           :: BathLowerCutOffFreq       !< lower cutoff frequency of the bath
   CHARACTER(100) :: SpectralDensityFile       !< spectral density file name
   TYPE(BathData), SAVE :: Bath                      !< derived datatype to define a single bath
   LOGICAL        :: ZPECorrection             !< ZeroPointEnergy correction in the initial conditions of the bath (at 0 K)
   REAL           :: OhmicGammaTimesMass       !< Gamma of an ohmic spectral density of the bath
   LOGICAL        :: NonLinearCoupling         !< system - bath coupling is non linear
   REAL           :: AlphaCoupling             !< paramter in non linear system-bath coupling
 
   ! POSITION, VELOCITY, ACCELERATION 

   REAL, DIMENSION(:), ALLOCATABLE, TARGET :: X    !< Position at given timestep
   REAL, DIMENSION(:), ALLOCATABLE, TARGET :: V    !< Velocity at given timestep
   REAL, DIMENSION(:), ALLOCATABLE, TARGET :: A    !< Acceleration at ginve timestep

   ! VECTOR WITH THE MASSES

   REAL, DIMENSION(:), ALLOCATABLE :: MassVector         !< Vector with the masses of the system
   
CONTAINS

   !> Subroutine to check the availability of a given runtype option
   SUBROUTINE CheckRunType( IntNr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IntNr
      LOGICAL :: Check 

      Check = ( IntNr /= EQUILIBRIUM .AND. &
                IntNr /= RELAXATION .AND. &
                IntNr /= SCATTERING .AND. &
                IntNr /= POTENTIAL )
      CALL ERROR( Check, " SharedData.CheckRunType: Invalid RunType option " )
   END SUBROUTINE CheckRunType

   !> Subroutine to check the availability of a given printtype option
   SUBROUTINE CheckPrintType( IntNr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IntNr
      LOGICAL :: Check 

      Check = ( IntNr /= DEBUG .AND. &
                IntNr /= FULL .AND. &
                IntNr /= MINIMAL  )
      CALL ERROR( Check, " SharedData.CheckPrintType: Invalid PrintType option " )
   END SUBROUTINE CheckPrintType

   !> Subroutine to check the availability of a given printtype option
   SUBROUTINE CheckBathType( IntNr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IntNr
      LOGICAL :: Check 

      Check = ( IntNr /= NORMAL_BATH .AND. &
                IntNr /= CHAIN_BATH .AND. &
                IntNr /= LANGEVIN_DYN  )
      CALL ERROR( Check, " SharedData.CheckBathType: Invalid BathType option " )
   END SUBROUTINE CheckBathType

END MODULE SharedData