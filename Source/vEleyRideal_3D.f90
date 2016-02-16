
  !-------------------ER_3D------------------------------ 
  !    Subroutine to compute the 3D potential 
  !        for the collinear ER reaction.
  !  DOFs = incidon height, targon height, C height
  !------------------------------------------------------

  SUBROUTINE ER_3D (zi, zt, zc, vv)
  IMPLICIT NONE

  REAL*8, INTENT (IN) :: zi, zt, zc
  REAL*8 :: ztc, ri, rt
  REAL*8 :: vv, vdiab, vstick, vmorse, vc, vh, vcoup

  !LEPS parameters for the incident atom
  REAL*8, PARAMETER :: di = 0.00775d0	!incidon potential well, in eV
  REAL*8, PARAMETER :: alphai = 0.954d0	!incidon curvature, in Ang(-1)
  REAL*8, PARAMETER :: ai = 4.01d0	!incidon eq. position, in Ang
  REAL*8, PARAMETER :: deltai = -0.80d0	

  !LEPS parameters for the target atom
  REAL*8, PARAMETER :: dt = 1.55d0	!targon potential well, in eV
  REAL*8, PARAMETER :: alphat = 2.34d0	!targon curvature, in Ang(-1)
  REAL*8, PARAMETER :: at = 1.52d0 	!targon eq. position, in Ang
  REAL*8, PARAMETER :: deltat = 0.035d0
  REAL*8, PARAMETER :: alpha2t = 1.73d0	!in Ang(-1)
  REAL*8, PARAMETER :: betat = 5.5d0	!in Ang(-1)
  REAL*8, PARAMETER :: bt = 2.10d0	!in Ang

  !LEPS parameters for the molecule
  REAL*8, PARAMETER :: dm = 4.58d0	!H2 potential well, in eV
  REAL*8, PARAMETER :: alpham = 1.81d0	!H2 curvature, in Ang(-1)
  REAL*8, PARAMETER :: am = 0.80d0	!H2 eq position, in Ang
  REAL*8, PARAMETER :: deltam = 0.18d0
  REAL*8, PARAMETER :: alpha2m = 0.96d0	!in Ang(-1)
  REAL*8, PARAMETER :: betam = 5.5d0	!in Ang(-1)
  REAL*8, PARAMETER :: bm = 1.89d0	!in Ang

  !temporary variables to store the different terms of the potential
  REAL*8 :: ft, fm, ui, ut, um, qi, qt, qm, r

  !Useful conversion factors
  REAL*8, PARAMETER :: MyConsts_Bohr2Ang = 0.52917721092d0	!conversion factor from Bohr to Ang
  REAL*8, PARAMETER :: MyConsts_Hartree2eV = 27.21138505d0	!conversion factor from Hartree to eV

  INTERFACE
   SUBROUTINE hstick_carbon( rcz, vv )
      REAL*8  rcz, vv
   END SUBROUTINE
  END INTERFACE

  INTERFACE   
    SUBROUTINE hstick_hydro( rhx, rhy, rhz, vv )
      REAL*8  rhx, rhy, rhz, vv
    END SUBROUTINE
  END INTERFACE

  INTERFACE
    SUBROUTINE hstick_coupling( rhx, rhy, rhz, rcz ,vv )
      REAL*8  rhx, rhy, rhz, rcz, vv
    END SUBROUTINE
  END INTERFACE

!LEPS terms:
 !1) full diabatic potential
  !convert the positions from a.u. in Ang
  ri = (zi - zc) * MyConsts_Bohr2Ang
  rt = (zt - zc) * MyConsts_Bohr2Ang

  !compute the incident atom-surface terms
  ui = di/(4.0d0*(1+deltai))*((3.0d0+deltai)*exp(-2.0d0*alphai*(ri-ai))&
       -(2.0d0+6.0d0*deltai)*exp(-alphai*(ri-ai)))

  qi = di/(4.0d0*(1.0d0+deltai))*((1.0d0+3.0d0*deltai)  &
       *exp(-2.0d0*alphai*(ri-ai))-(6.0d0+2.0d0*deltai)  &
       *exp(-alphai*(ri-ai)))

  !compute the target atom-surface terms
  ft = exp(-alpha2t*(rt-at)/(1.0d0+exp(-betat*(rt-bt))))

  ut = dt/(4.0d0*(1.0d0+deltat))*((3.0d0+deltat)  &
       *exp(-2.d0*alphat*(rt-at))-(2.0d0+6.0d0*deltat)  &
       *ft*exp(-alphat*(rt-at)))

  qt = dt/(4.0d0*(1.0d0+deltat))*((1.0d0+3.0d0*deltat)  &
       *exp(-2.0d0*alphat*(rt-at))-(6.0d0+2.0d0*deltat)  &
       *ft*exp(-alphat*(rt-at)))

  !compute the molecular terms
  r = abs(ri-rt)

  fm = exp(-alpha2m*(r-am)/(1+exp(-betam*(r-bm))))

  um = dm/(4.0d0*(1.0d0+deltam))*((3.0d0+deltam)  &
       *exp(-2.0d0*alpham*(r-am))-(2.0d0+6.0d0*deltam)*fm  &
       *exp(-alpham*(r-am)))

  qm = dm/(4.0d0*(1.0d0+deltam))*((1.0d0+3.0d0*deltam)  &
       *exp(-2*alpham*(r-am))-(6.0d0+2.0d0*deltam)*fm  &
       *exp(-alpham*(r-am)))

  !compute the full diabatic potential
  vdiab = ui+ut+um-sqrt((qm**2+(qi+qt)**2-qm*(qi+qt)))

!   ! Upper and lower energy cutoff
!   if ( vdiab < -5.0 ) vdiab = -5.0
!   if ( vdiab > 20.0 ) vdiab = 20.0

  !convert the full potential in a.u.
  vdiab = vdiab / MyConsts_Hartree2eV

 !2) 1D targon term

  vmorse = ut - abs(qt)

!   ! Upper and lower energy cutoff
!   if ( vmorse < -5.0 ) vmorse = -5.0
!   if ( vmorse > 20.0 ) vmorse = 20.0

  !convert the potential in a.u.
  vmorse = vmorse / MyConsts_Hartree2eV

!N.B. everything is in a.u.

  ztc = zt - zc


!sticking potential for targon and C: it has three parts
  !carbon term
  CALL hstick_carbon( zc, vc )

  !hydrogen (targon) term
  CALL hstick_hydro( 0.0d0, 0.0d0, zt, vh )

  !coupling term
  CALL hstick_coupling( 0.0d0, 0.0d0, zt, zc ,vcoup )

  !full sticking potential
  vstick = vc + vh + vcoup

!full potential
  vv = vdiab + vstick - vmorse

!Upper and lower energy cutoff
  if ( vv < -0.18d0 ) vv = -0.18d0
  if ( vv > 0.735d0 ) vv = 0.735d0

!convert the full potential in eV
!  vv = vv * MyConsts_Hartree2eV
  

  END SUBROUTINE
