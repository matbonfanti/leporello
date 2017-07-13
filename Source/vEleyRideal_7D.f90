
  !-------------------ER_3D------------------------------
  !    Subroutine to compute the 7D potential
  !        for the ER reaction.
  !        DOFs = incidon x, y, z;
  !               targon x, y, z;
  !               carbon z.
  !------------------------------------------------------

  SUBROUTINE ER_7D (dof, vv, dv)
  IMPLICIT NONE

  REAL*8, DIMENSION(7), INTENT(IN) :: dof
  REAL*8, INTENT(OUT) :: vv
  REAL*8, DIMENSION(7), INTENT(OUT) :: dv

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
  REAL*8 :: ft, fm, ui, ut, um, qi, qt, qm, r, r_bohr

  !potential derivatives
  REAL*8 :: d_ui, d_qi, d_ut, d_qt, d_um, d_qm
  REAL*8 :: dv_ri, dv_rt, dv_r, d_ft, d_fm, d_vmorse
  REAL*8 :: dvstick_dxt, dvstick_dyt, dvstick_dzt, dvstick_dzc
  REAL*8 :: dr_dxi, dr_dyi, dr_dzi, dr_dxt, dr_dyt, dr_dzt

  !Useful conversion factors
  REAL*8, PARAMETER :: MyConsts_Bohr2Ang = 0.52917721092d0	!conversion factor from Bohr to Ang
  REAL*8, PARAMETER :: MyConsts_Hartree2eV = 27.21138505d0	!conversion factor from Hartree to eV

!LEPS terms:
 !1) full diabatic potential
  !convert the positions from a.u. in Ang

  !ri is the heigth of the incidon above the surface (zi - zc)
  ri = (dof(3) - dof(7)) * MyConsts_Bohr2Ang

  !rt is the heigth of the targon above the surface (zt - zc)
  rt = (dof(6) - dof(7)) * MyConsts_Bohr2Ang

  !compute the incident atom-surface terms and their derivatives respect to ri

  ui = di/(4.0d0*(1.0d0+deltai))*((3.0d0+deltai)*exp(-2.0d0*alphai* &
       (ri-ai))-(2.0d0+6.0d0*deltai)*exp(-alphai*(ri-ai)))

  d_ui = di*alphai/4.0d0/(1.0d0+deltai)*(-2.0d0*(3.0d0+deltai)* &
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

  !definition of r = sqrt( (xi-xt)**2 + (yi-yt)**2 + (zi-zt)**2 )
  r_bohr = sqrt( (dof(1)-dof(4))**2+(dof(2)-dof(5))**2+(dof(3)-dof(6))**2)
  r = r_bohr * MyConsts_Bohr2Ang

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
  vdiab = ui+ut+um-sqrt((qm**2+(qi+qt)**2-qm*(qi+qt)))

  dv_ri = d_ui-0.5d0*d_qi*((2.0d0*(qi+qt)-qm)/(sqrt(qm**2+(qi+qt)**2 &
	  -qm*(qi+qt))))

  dv_rt = d_ut-0.5d0*d_qt*((2.0d0*(qi+qt)-qm)/(sqrt(qm**2+(qi+qt)**2 &
	  -qm*(qi+qt))))

  dv_r = d_um-0.5d0*d_qm*((2.0d0*qm-qi-qt)/(sqrt(qm**2+(qi+qt)**2 &
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

  !sticking potential for targon and C
  CALL hstick( dof(4),dof(5),dof(6),dof(7), vstick, dvstick_dxt, dvstick_dyt, dvstick_dzt, dvstick_dzc )

  !full potential
  vv = vdiab - vmorse + vstick


  ! derivatives of r wrt the subsystem DOFs

  dr_dxi = (dof(1) - dof(4)) / MAX(r_bohr,0.0001)
  dr_dyi = (dof(2) - dof(5)) / MAX(r_bohr,0.0001)
  dr_dzi = (dof(3) - dof(6)) / MAX(r_bohr,0.0001)
  dr_dxt = (dof(4) - dof(1)) / MAX(r_bohr,0.0001)
  dr_dyt = (dof(5) - dof(2)) / MAX(r_bohr,0.0001)
  dr_dzt = (dof(6) - dof(3)) / MAX(r_bohr,0.0001)

  ! sum derivatives
  !derivatives of the potential wrt subsytem dofs
  dv(1) =         + dv_r * dr_dxi                 !dV/dxi = dV/dr * dr/dxi

  dv(2) =         + dv_r * dr_dyi                 !dV/dyi = dV/dr * dr/dyi

  dv(3) = + dv_ri + dv_r * dr_dzi                 !dV/dzi = dV/dri*dri/dzi + dV/dr*dr/dzi

  dv(4) =         + dv_r * dr_dxt + dvstick_dxt   !dV/dxt = dV/dr*dr/dxt + dVstick/dxt

  dv(5) =         + dv_r * dr_dyt + dvstick_dyt   !dV/dyt = dV/dr * dr/dyt + dVstick/dyt

  dv(6) = + dv_rt + dv_r * dr_dzt + dvstick_dzt   !dV/dzt

  dv(7) = - dv_ri - dv_rt         + dvstick_dzc   !dV/dzc

!Upper and lower energy cutoff
  if ( vv < -0.18d0 ) then
       vv = -0.18d0
       dv = 0.
  else if ( vv > 0.367d0 ) then
       vv = 0.367d0
       dv = 0.
  endif

!convert the full potential in eV
!  vv = vv * MyConsts_Hartree2eV


  CONTAINS

   SUBROUTINE hstick( rhx, rhy, rhz, rcz ,vv, dv_rhx, dv_rhy, dv_rhz, dv_rcz )
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: rhx, rhy, rhz, rcz
      REAL*8, INTENT(OUT) :: vv, dv_rhx, dv_rhy, dv_rhz, dv_rcz

      REAL*8, PARAMETER    :: MyConsts_Bohr2Ang       = 0.52917721092d0  !< Conversion factor from Bohr to Ang (from NIST )
      REAL*8, PARAMETER    :: MyConsts_Hartree2eV     = 27.21138505d0    !< Conversion factor from Hartree to eV (from NIST )
      REAL*8, PARAMETER :: di = 0.00775d0           !< potential well (eV)
      REAL*8, PARAMETER :: alphi = 0.954d0          !< curvature (Ang-1)
      REAL*8, PARAMETER :: ai = 4.01d0              !< equilibrium position (Ang)
      ! parameters to describe the behaviour along rho at large rho
      REAL*8, PARAMETER :: alp = 1.4d0              !< a_f in the text (Ang-1)
      REAL*8, PARAMETER :: alp2 = 1.65d0            !< b_f in the text (Ang-1)
      REAL*8, PARAMETER :: ba = 2.9d0               !< c_f in the text (Ang-1)
      REAL*8, PARAMETER :: rhoa = 1.0d0             !< rho_f in the text (Ang)
      ! switching function from small rho to large rho
      REAL*8, PARAMETER :: bs = 6.0d0               !< b_s in the text (Ang-1)
      REAL*8, PARAMETER :: rhos = 0.53d0            !< rho_s in the text (Ang)

      REAL*8, PARAMETER :: delt = 2.651d0           !< twisting potential parameter delta (eV)
      REAL*8, PARAMETER :: gam2 = 1.958d0           !< puckering paramter gamma_2 (eV)
      REAL*8, PARAMETER :: bndprm = 2.4612d0        !< lattice constant (in Ang)

      REAL*8 :: xh, yh, rkc
      REAL*8 :: a, b, c1, c2, d0

      REAL*8 :: dbdzc, dc1dzc, dc2dzc, dd0dzc, df2dr, dkrdzc, dkrdzh
      REAL*8 :: dswdr, dvdzc, dvdzh, dvidzc, dvidzh, dvqdr, dvqdzc
      REAL*8 :: dvqdzh, dvtdrho, dvtds1, dvtds2
      REAL*8 :: dzgdzc, dzmdzc

      REAL*8 :: fexp, ff1, ff2, ff3
      REAL*8 :: rho, rkrho, sub1, sub2, sw
      REAL*8 :: v, vi, vq, vt
      REAL*8 :: zg, zm

      INTEGER :: nn

      ! Coordinates in atomic units
      xh = rhx*MyConsts_Bohr2Ang
      yh = rhy*MyConsts_Bohr2Ang
      rho=sqrt(rhx**2+rhy**2)* MyConsts_Bohr2Ang
      sub1 = rhz * MyConsts_Bohr2Ang
      sub2 = rcz * MyConsts_Bohr2Ang

      rkc = (36.0*gam2+6.0*delt)/(bndprm**2)

! **************************************************************************************************
!                           POTENTIAL FOR THE C-H SYSTEM
! **************************************************************************************************

      ! Compute the parameters for the morse + gaussian
      ! functional form for the collinear potential

      ! D_0 is the morse depth (eV)
      d0=0.474801+0.9878257*sub2-1.3921499*sub2**2              &
      +0.028278*sub2**3-1.8879928*sub2**4                       &
      +0.11*exp(-8.0*(sub2-0.28)**2)
      dd0dzc=0.9878257-2.7842998*sub2+0.084834*sub2**2-         &
      7.5519712*sub2**3+(                                       &
      -1.76*(sub2-0.28)*exp(-8.0*(sub2-0.28)**2))

      ! A is the morse curvature (Ang)
      a=2.276211

      ! Z_M is the morse equilibrium distance (Ang)
      zm=0.87447*(sub2)+1.17425
      dzmdzc=0.87447

      ! C_1 is the asympotic harmonic potential for C1 vibrations (eV)
      c1=0.5*rkc*(sub2)**2-0.00326
      dc1dzc=rkc*(sub2)

      ! C_2 is the gaussian height (eV)
      c2= 0.3090344*exp(-2.741813*(sub2-0.2619756))+            &
      0.03113325*exp(3.1844857*(sub2-0.186741))+                &
      0.02*exp(-20.0*(sub2-0.1)**2)
      dc2dzc=-0.8473145*exp(-2.741813*(sub2-0.2619756))+        &
      0.0991434*exp(3.1844857*(sub2-0.186741))+(                &
      -0.8*(sub2-0.1)*exp(-20.0*(sub2-0.1)**2))

      ! B is the gaussian curvature (Ang-2)
      b=4.00181*exp(1.25965*(sub2-0.58729))
      dbdzc=5.0408799*exp(1.25965*(sub2-0.58729))

      ! Z_G is the center of the gaussian (Ang)
      zg=1.99155*(sub2)+1.46095
      dzgdzc=1.99155

      ! Compute the potential and the derivatives of the
      ! collinear potential V_0

      v=c1+(c1+d0)*(exp(-2.0*a*(sub1-zm))-2.0*exp(-a*(sub1-zm)))+       &
      c2*exp(-b*(sub1-zg)**2)

      dvdzh=-b*c2*(sub1-zg)*2.0*exp(-b*(sub1-zg)**2) +                  &
      (c1+d0)*2.0*(-a*exp(-2.0*a*(sub1-zm))+a*exp(-a*(sub1-zm)))

      dvdzc=dc1dzc+(dc1dzc+dd0dzc)*                             &
      (exp(-2.0*a*(sub1-zm))-2.0*exp(-a*(sub1-zm)))+            &
      (c1+d0)*(a*dzmdzc*2.0*exp(-2.0*a*(sub1-zm))+              &
      (-a*dzmdzc)*2.0*exp(-a*(sub1-zm)))+                       &
      dc2dzc*exp(-b*(sub1-zg)**2)+                              &
      (-c2*dbdzc*(sub1-zg)**2+c2*b*2.0*(sub1-zg)*dzgdzc)*       &
      exp(-b*(sub1-zg)**2)

      ! Compute the force constant (and derivatives) for small rho
      ! potential rkrho(zh-q,zc-q)

      rkrho=3.866259*exp(-17.038588*(sub2-0.537621)**2+         &
      0.312355*(sub2-0.537621)*(sub1-2.003753)-                 &
      4.479864*(sub1-2.003753)**2)+                             &
      4.317415*exp(-11.931770*(sub2-0.286858)**2+               &
      18.540974*(sub2-0.286858)*(sub1-1.540947)-                &
      14.537321*(sub1-1.540947)**2)

      dkrdzc=(-34.077176*(sub2-0.537621)+0.312355*              &
      (sub1-2.003753))                                          &
      *3.866259*exp(-17.038588*(sub2-0.537621)**2+              &
      0.312355*(sub2-0.537621)*(sub1-2.003753)-                 &
      4.479864*(sub1-2.003753)**2)+                             &
      (-23.86354*(sub2-0.286858)+18.540974*(sub1-1.540947))     &
      *4.317415*exp(-11.931770*(sub2-0.286858)**2+              &
      18.540974*(sub2-0.286858)*(sub1-1.540947)-                &
      14.537321*(sub1-1.540947)**2)

      dkrdzh=(0.312355*(sub2-0.537621)-8.959728*(sub1-2.003753))        &
      *3.866259*exp(-17.038588*(sub2-0.537621)**2+                      &
      0.312355*(sub2-0.537621)*(sub1-2.003753)-                         &
      4.479864*(sub1-2.003753)**2)+                                     &
      (18.540974*(sub2-0.286858)-29.074642*(sub1-1.540947))             &
      *4.317415*exp(-11.931770*(sub2-0.286858)**2+                      &
      18.540974*(sub2-0.286858)*(sub1-1.540947)-                        &
      14.537321*(sub1-1.540947)**2)

      ! Compute the small rho potential/derivatives

      vq=v+0.5*rkrho*rho**2
      dvqdzh=dvdzh+0.5*dkrdzh*rho**2
      dvqdzc=dvdzc+0.5*dkrdzc*rho**2
      dvqdr=rkrho*rho

      ! Compute the  "infinite" rho potential/derivatives

      vi=0.5*rkc*sub2**2-0.00326+                                       &
      di*(exp(-2.0*alphi*(sub1-ai))-2.0*exp(-alphi*(sub1-ai)))
      dvidzh=di*(-2.0*alphi*exp(-2.0*alphi*(sub1-ai))+                  &
            2.0*alphi*exp(-alphi*(sub1-ai)))
      dvidzc=rkc*sub2

      ! Switching function and associated functions

      fexp=exp(-ba*(rho-rhoa))
      ff1=1.0+fexp
      ff2=exp(-2.0*alp*rho)-2.0*exp(-alp*rho)*exp(-alp2*rho/ff1)
      ff3=(vi-v)*ff2+vi
      sw=1.0/(1.0+exp(-bs*(rho-rhos)))

      ! Total H,C1 potential/derivatives

      vt=vq*(1.0-sw)+ff3*sw

      df2dr=-2.0*alp*exp(-2.0*alp*rho)-                &
         2.0*exp(-alp*rho)*exp(-alp2*rho/ff1)*       &
         (-alp-(alp2/ff1)-alp2*rho*ba*fexp/(ff1**2))
      dswdr=(bs*exp(-bs*(rho-rhos)))/((1.0+exp(-bs*(rho-rhos)))**2)

      dvtds1=dvqdzh*(1.0-sw)+sw*((dvidzh-dvdzh)*ff2+dvidzh)
      dvtds2=dvqdzc*(1.0-sw)+sw*((dvidzc-dvdzc)*ff2+dvidzc)
      dvtdrho=dvqdr*(1.0-sw)+vq*(-dswdr)+sw*(vi-v)*df2dr+ff3*dswdr

! **************************************************************************************************
!            OUTPUT VALUES
! **************************************************************************************************

      ! Convert total potential to AU
      vv=vt / MyConsts_Hartree2eV

      ! potential derivatives
      IF ( rho  < 0.001 ) THEN
         dv_rhx = 0.0
         dv_rhy = 0.0
      ELSE
         dv_rhx = dvtdrho*xh/rho
         dv_rhy = dvtdrho*yh/rho
      ENDIF
      dv_rhz=dvtds1
      dv_rcz=dvtds2

      ! Transform forces in atomic units (from eV Ang^-1 to Hartree Bohr^-1)
      dv_rhx = dv_rhx * MyConsts_Bohr2Ang / MyConsts_Hartree2eV
      dv_rhy = dv_rhy * MyConsts_Bohr2Ang / MyConsts_Hartree2eV
      dv_rhz = dv_rhz * MyConsts_Bohr2Ang / MyConsts_Hartree2eV
      dv_rcz = dv_rcz * MyConsts_Bohr2Ang / MyConsts_Hartree2eV

   END SUBROUTINE hstick
END SUBROUTINE

