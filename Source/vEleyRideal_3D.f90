
  !-------------------ER_3D------------------------------ 
  !    Subroutine to compute the 3D potential 
  !        for the collinear ER reaction.
  !  DOFs = incidon height, targon height, C height
  !------------------------------------------------------

  SUBROUTINE ER_3D (zi, zt, zc, vv, dv_zi, dv_zt, dv_zc)
  IMPLICIT NONE

  REAL*8, INTENT (IN) :: zi, zt, zc
  REAL*8, INTENT(OUT) :: vv, dv_zi, dv_zt, dv_zc

  REAL*8 :: ri, rt
  REAL*8 :: vdiab, vstick, vmorse

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

  !potential derivatives
  REAL*8 :: d_ui, d_qi, d_ut, d_qt, d_um, d_qm
  REAL*8 :: dv_ri, dv_rt, dv_r, d_ft, d_fm, d_vmorse
  REAL*8 :: dvstick_dxt, dvstick_dyt, dvstick_dzt, dvstick_dzc

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

  !compute the incident atom-surface terms and their derivatives respect to ri

  ui = di/(4.0d0*(1.0d0+deltai))*((3.0d0+deltai)*exp(-2.0d0*alphai* &
       (ri-ai))-(2.0d0+6.0d0*deltai)*exp(-alphai*(ri-ai)))

  d_ui = di*alphai/4.0d0*(1.0d0+deltai)*(-2.0d0*(3.0d0+deltai)* &
	 exp(-2.0d0*alphai*(ri-ai))+(2.0d0+6.0d0*deltai)* &
	 exp(-alphai*(ri-ai)))

  qi = di/(4.0d0*(1.0d0+deltai))*((1.0d0+3.0d0*deltai)  &
       *exp(-2.0d0*alphai*(ri-ai))-(6.0d0+2.0d0*deltai)  &
       *exp(-alphai*(ri-ai)))

  d_qi = di*alphai/(4.0d0*(1.0d0+deltai))*(-2.0d0*(1.0d0+3.0d0*deltai) &
	 *exp(-2.0d0*alphai*(ri-ai))+(6.0d0+2.0d0*deltai)  &
         *exp(-alphai*(ri-ai)))


  !compute the target atom-surface terms and their derivatives respect to rt

  ft = exp(-alpha2t*(rt-at)/(1.0d0+exp(-betat*(rt-bt))))

  d_ft = ft*((-alpha2t*(1.0d0+exp(-betat*(rt-bt)))+alpha2t*(rt-at)  &
	 *exp(-betat*(rt-bt))*(-betat))/(1+exp(-betat*(rt-bt)))**2)

  ut = dt/(4.0d0*(1.0d0+deltat))*((3.0d0+deltat)  &
       *exp(-2.d0*alphat*(rt-at))-(2.0d0+6.0d0*deltat)  &
       *ft*exp(-alphat*(rt-at)))

  d_ut = dt/(4.0d0*(1.0d0+deltat))*(-2.0d0*alphat*(3.0d0+deltat)  &
	 *exp(-2.d0*alphat*(rt-at))-(2.0d0+6.0d0*deltat)*(-alphat  &
	 *exp(-alphat*(rt-at))*ft+exp(-alphat*(rt-at))*d_ft))

  qt = dt/(4.0d0*(1.0d0+deltat))*((1.0d0+3.0d0*deltat)  &
       *exp(-2.0d0*alphat*(rt-at))-(6.0d0+2.0d0*deltat)  &
       *ft*exp(-alphat*(rt-at)))

  d_qt = dt/(4.0d0*(1.0d0+deltat))*(-2.0d0*alphat*(1.0d0+3.0d0*deltat)&
	 *exp(-2.0d0*alphat*(rt-at))-(6.0d0+2.0d0*deltat)*(-alphat  &
	 *exp(-alphat*(rt-at))*ft+exp(-alphat*(rt-at))*d_ft))

  !compute the molecular terms and their derivatives respect to r

  r = abs(ri-rt)

  fm = exp(-alpha2m*(r-am)/(1+exp(-betam*(r-bm))))

  d_fm = fm*((-alpha2m*(1+exp(-betam*(r-bm)))+alpha2m*(r-am)  &
	 *exp(-betam*(r-bm))*(-betam))/(1+exp(-betam*(r-bm)))**2)

  um = dm/(4.0d0*(1.0d0+deltam))*((3.0d0+deltam)  &
       *exp(-2.0d0*alpham*(r-am))-(2.0d0+6.0d0*deltam)*fm  &
       *exp(-alpham*(r-am)))

  d_um = dm/(4.0d0*(1.0d0+deltam))*(-2.0d0*alpham*(3.0d0+deltam)  &
	 *exp(-2.0d0*alpham*(r-am))-(2.0d0+6.0d0*deltam)*(-alpham  &
	 *exp(-alpham*(r-am))*fm+exp(-alpham*(r-am))*d_fm))

  qm = dm/(4.0d0*(1.0d0+deltam))*((1.0d0+3.0d0*deltam)  &
       *exp(-2*alpham*(r-am))-(6.0d0+2.0d0*deltam)*fm  &
       *exp(-alpham*(r-am)))

  d_qm = dm/(4.0d0*(1.0d0+deltam))*(-2.0d0*alpham*(1.0d0+3.0d0*deltam)&
	 *exp(-2*alpham*(r-am))-(6.0d0+2.0d0*deltam)*(-alpham  &
	 *exp(-alpham*(r-am))*fm+exp(-alpham*(r-am))*d_fm))

  !compute the full diabatic potential and its partial derivatives
  vdiab = ui+ut+um-sqrt((qm*qt*2+(qi+qt)**2-qm*(qi+qt)))

  dv_ri = d_ui-0.5d0*d_qi*((2.0d0*(qi+qt)-qm)/(sqrt(qm**2-(qi+qt)**2 &
	  -qm*(qi+qt))))

  dv_rt = d_ut-0.5d0*d_qt*((2.0d0*(qi+qt)-qm)/(sqrt(qm**2-(qi+qt)**2 &
	  -qm*(qi+qt))))

  dv_r = d_um-0.5d0*d_qm*((2.0d0*qm-qi-qt)/(sqrt(qm**2-(qi+qt)**2 &
	  -qm*(qi+qt))))

  !convert the full potential in a.u.
  vdiab = vdiab / MyConsts_Hartree2eV
  !convert the derivatives to a.u.
  dv_ri = dv_ri / MyConsts_Hartree2eV * MyConsts_Bohr2Ang
  dv_rt = dv_rt / MyConsts_Hartree2eV * MyConsts_Bohr2Ang
  dv_r = dv_r / MyConsts_Hartree2eV * MyConsts_Bohr2Ang

 !2) 1D targon term

  vmorse = ut - abs(qt)
  d_vmorse = d_ut - d_qt * SIGN(1.0,qt)

  !convert the potential in a.u.
  vmorse = vmorse / MyConsts_Hartree2eV
  !convert the derivative to a.u.
  d_vmorse = d_vmorse / MyConsts_Hartree2eV * MyConsts_Bohr2Ang
  ! d_vmorse is derivative with respect to rt, sum to the corresponing term (with minus sign V = Vdiab - Vmorse)
  dv_rt = dv_rt - d_vmorse

!N.B. everything is in a.u.

  !sticking potential for targon and C: it has three parts
  CALL hstick_coupling( 0.0d0, 0.0d0, zt, zc, vstick, dvstick_dxt, dvstick_dyt, dvstick_dzt, dvstick_dzc )

  !full potential
  vv = vdiab - vmorse + vstick 

  ! sum derivatives
  dv_zi = +dv_ri +dv_r
  dv_zt = +dv_rt -dv_r  + dvstick_dzt
  dv_zc = -dv_ri -dv_rt + dvstick_dzc

!Upper and lower energy cutoff
  if ( vv < -0.18d0 ) vv = -0.18d0
  if ( vv > 0.735d0 ) vv = 0.735d0

!convert the full potential in eV
!  vv = vv * MyConsts_Hartree2eV
  

  END SUBROUTINE
