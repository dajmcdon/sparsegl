! --------------------------------------------------
SUBROUTINE ls_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
  ! --------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: big=9.9E30
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
  INTEGER, PARAMETER :: mnlam = 6
  INTEGER:: mnl
  INTEGER:: bn
  INTEGER::bs(bn)
  INTEGER::ix(bn)
  INTEGER::iy(bn)
  INTEGER:: nobs
  INTEGER::nvars
  INTEGER::dfmax
  INTEGER::pmax
  INTEGER::nlam
  INTEGER::nalam
  INTEGER::npass
  INTEGER::jerr
  INTEGER::maxit
  INTEGER::intr
  INTEGER:: idx(pmax)
  INTEGER::nbeta(nlam)
  DOUBLE PRECISION:: flmin
  DOUBLE PRECISION::eps
  DOUBLE PRECISION:: x(nobs,nvars)
  DOUBLE PRECISION::y(nobs)
  DOUBLE PRECISION::pf(bn)
  DOUBLE PRECISION::ulam(nlam)
  DOUBLE PRECISION::gam(bn)
  DOUBLE PRECISION:: b0(nlam)
  DOUBLE PRECISION::beta(nvars,nlam)
  DOUBLE PRECISION::alam(nlam)
  ! - - - local declarations - - -
  DOUBLE PRECISION:: max_gam
  DOUBLE PRECISION::d
  DOUBLE PRECISION::t
  DOUBLE PRECISION::dif
  DOUBLE PRECISION::unorm
  DOUBLE PRECISION::al
  DOUBLE PRECISION::alf
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
  INTEGER:: g
  INTEGER::j
  INTEGER::l
  INTEGER::ni
  INTEGER::me
  INTEGER::startix
  INTEGER::endix
  ! - - - begin - - -
  ! - - - local declarations - - -
  DOUBLE PRECISION:: tlam
  INTEGER:: jx
  INTEGER:: jxx(bn)
  DOUBLE PRECISION:: ga(bn)
  DOUBLE PRECISION:: vl(nvars)
  DOUBLE PRECISION:: al0
  ! - - - allocate variables - - -
  ALLOCATE(b(0:nvars))
  ALLOCATE(oldbeta(0:nvars))
  ALLOCATE(r(1:nobs))
  ALLOCATE(oidx(1:bn))
  ! - - - checking pf - - -
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  jxx = 0
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  idx = 0
  oidx = 0
  npass = 0
  ni = npass
  alf = 0.0D0
  ! --------- lambda loop ----------------------------
  IF(flmin < 1.0D0) THEN ! flmin <1 means NOT user-supplied, so this is default
     flmin = Max (mfl, flmin)
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  vl = matmul(r, x)/nobs
  DO g = 1,bn
     ALLOCATE(u(bs(g)))
     u = vl(ix(g):iy(g))
     ga(g) = sqrt(dot_product(u,u))
     DEALLOCATE(u)
  ENDDO
  DO l=1,nlam
     al0 = al
     IF(flmin>=1.0D0) THEN
        al=ulam(l)
     ELSE
        IF(l > 2) THEN
           al=al*alf
        ELSE IF(l==1) THEN
           al=big
        ELSE IF(l==2) THEN
           al0 = 0.0D0
           DO g = 1,bn
              IF(pf(g)>0.0D0) THEN
                 al0 = max(al0, ga(g) / pf(g))
              ENDIF
           ENDDO
           al = al0 * alf
        ENDIF
     ENDIF
     tlam = (2.0*al-al0)
     DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
     ENDDO
     ! --------- outer loop ----------------------------
     DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
           DO j=1,ni
              g=idx(j)
              oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
           npass=npass+1
           dif=0.0D0
           DO g=1,bn
              IF(jxx(g) == 0) CYCLE
              startix=ix(g)
              endix=iy(g)
              ALLOCATE(u(bs(g)))
              ALLOCATE(dd(bs(g)))
              ALLOCATE(oldb(bs(g)))
              oldb=b(startix:endix)
              u=matmul(r,x(:,startix:endix))/nobs
              u=gam(g)*b(startix:endix)+u
              unorm=sqrt(dot_product(u,u))
              t=unorm-pf(g)*al
              IF(t>0.0D0) THEN
                 b(startix:endix)=u*t/(gam(g)*unorm)
              ELSE
                 b(startix:endix)=0.0D0
              ENDIF
              dd=b(startix:endix)-oldb
              IF(any(dd/=0.0D0)) THEN
                 dif=max(dif,gam(g)**2*dot_product(dd,dd))
                 r=r-matmul(x(:,startix:endix),dd)
                 IF(oidx(g)==0) THEN
                    ni=ni+1
                    IF(ni>pmax) EXIT
                    oidx(g)=ni
                    idx(ni)=g
                 ENDIF
              ENDIF
              DEALLOCATE(u,dd,oldb)
           ENDDO
           IF(intr /= 0) THEN
              d=sum(r)/nobs
              IF(d/=0.0D0) THEN
                 b(0)=b(0)+d
                 r=r-d
                 dif=max(dif,d**2)
              ENDIF
           ENDIF
           IF (ni > pmax) EXIT
           IF (dif < eps) EXIT
           IF(npass > maxit) THEN
              jerr=-l
              RETURN
           ENDIF
           ! --inner loop----------------------
           DO
              npass=npass+1
              dif=0.0D0
              DO j=1,ni
                 g=idx(j)
                 startix=ix(g)
                 endix=iy(g)
                 ALLOCATE(u(bs(g)))
                 ALLOCATE(dd(bs(g)))
                 ALLOCATE(oldb(bs(g)))
                 oldb=b(startix:endix)
                 u=matmul(r,x(:,startix:endix))/nobs
                 u=gam(g)*b(startix:endix)+u
                 unorm=sqrt(dot_product(u,u))
                 t=unorm-pf(g)*al
                 IF(t>0.0D0) THEN
                    b(startix:endix)=u*t/(gam(g)*unorm)
                 ELSE
                    b(startix:endix)=0.0D0
                 ENDIF
                 dd=b(startix:endix)-oldb
                 IF(any(dd/=0.0D0)) THEN
                    dif=max(dif,gam(g)**2*dot_product(dd,dd))
                    r=r-matmul(x(:,startix:endix),dd)
                 ENDIF
                 DEALLOCATE(u,dd,oldb)
              ENDDO
              IF(intr /= 0) THEN
                 d=sum(r)/nobs
                 IF(d/=0.0D0) THEN
                    b(0)=b(0)+d
                    r=r-d
                    dif=max(dif,d**2)
                 ENDIF
              ENDIF
              IF(dif<eps) EXIT
              IF(npass > maxit) THEN
                 jerr=-l
                 RETURN
              ENDIF
           ENDDO ! end inner loop
        ENDDO ! end middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
        IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        vl = matmul(r, x)/nobs
        DO g = 1, bn
           IF(jxx(g) == 1) CYCLE
           ALLOCATE(u(bs(g)))
           u = vl(ix(g):iy(g))
           ga(g) = sqrt(dot_product(u,u))
           IF(ga(g) > al*pf(g))THEN
              jxx(g) = 1
              jx = 1
           ENDIF
           DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
     ENDDO ! End outer loop
     !---------- final update variable and save results------------
     IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
     ENDIF
     IF(ni>0) THEN
        DO j=1,ni
           g=idx(j)
           beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
     ENDIF
     nbeta(l)=ni
     b0(l)=b(0)
     alam(l)=al
     nalam=l
     IF (l < mnl) CYCLE
     me=0
     DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,oidx)
  RETURN
END SUBROUTINE ls_f

! --------------------------------------------------
SUBROUTINE ls_f_sparse (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr,alsparse)
  ! --------------------------------------------------
  ! to do: remove intr, b0,
  
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: big=9.9E30
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
  INTEGER, PARAMETER :: mnlam = 6
  INTEGER:: mnl
  INTEGER:: bn
  INTEGER::bs(bn)
  INTEGER::ix(bn)
  INTEGER::iy(bn)
  INTEGER:: nobs
  INTEGER::nvars
  INTEGER::dfmax
  INTEGER::pmax
  INTEGER::nlam
  INTEGER::nalam
  INTEGER::npass
  INTEGER::jerr
  INTEGER::maxit
  INTEGER::intr ! dont need you
  INTEGER:: idx(pmax)
  INTEGER::nbeta(nlam)
  DOUBLE PRECISION:: flmin
  DOUBLE PRECISION::eps
  DOUBLE PRECISION:: x(nobs,nvars)
  DOUBLE PRECISION::y(nobs)
  DOUBLE PRECISION::pf(bn)
  DOUBLE PRECISION::ulam(nlam)
  DOUBLE PRECISION::gam(bn)
  DOUBLE PRECISION:: b0(nlam) ! don't need you
  DOUBLE PRECISION::beta(nvars,nlam) ! this will be sparsified
  DOUBLE PRECISION::alam(nlam)
  DOUBLE PRECISION::alsparse ! Should be alpha for the sparsity weight
  ! - - - local declarations - - -
  DOUBLE PRECISION:: max_gam
  DOUBLE PRECISION::d
  ! DOUBLE PRECISION::t ! No longer using this
  DOUBLE PRECISION::dif
  ! DOUBLE PRECISION::unorm ! No longer using this for ls_new
  DOUBLE PRECISION::al
  DOUBLE PRECISION::alf
  DOUBLE PRECISION::sg
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r ! Residue y-beta_k*x etc
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u ! No longer using this for update step, but still need for other parts
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
  INTEGER:: g
  INTEGER::j
  INTEGER::l
  INTEGER::ni
  INTEGER::me
  INTEGER::startix
  INTEGER::endix
  ! - - - Aaron's declarations
  DOUBLE PRECISION::snorm
  DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma
  DOUBLE PRECISION::tea ! this takes the place of 't' in the update step for ls
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s ! takes the place of 'u' in update for ls
  INTEGER::soft_g ! this is an iterating variable 'vectorizing' the soft thresholding operator
  INTEGER::vl_iter ! for iterating over columns(?) of x*r
  INTEGER::kill_count
  ! - - - begin - - -
  ! - - - local declarations - - -
  DOUBLE PRECISION:: tlam
  DOUBLE PRECISION:: lama
  DOUBLE PRECISION:: lam1ma
  INTEGER:: jx
  INTEGER:: jxx(bn)
  DOUBLE PRECISION:: ga(bn) ! What is this for??
  DOUBLE PRECISION:: vl(nvars) ! What is this for?
  DOUBLE PRECISION:: al0
  ! - - - allocate variables - - -
  ALLOCATE(b(0:nvars)) ! make them 1 index (ditch intercept), to pmax?
  ALLOCATE(oldbeta(0:nvars)) ! make them 1 index (ditch intercept)
  ALLOCATE(r(1:nobs))
  ALLOCATE(oidx(1:bn))
  !    ALLOCATE(al_sparse)
  ! - - - checking pf - ! pf is the relative penalties for each group
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  jxx = 0
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  idx = 0
  kill_count = 0
  oidx = 0
  npass = 0 ! This is a count, correct?
  ni = npass ! This controls so-called "outer loop"
  alf = 0.0D0
  t_for_s = 1/gam ! might need to use a loop if no vectorization.........
  ! --------- lambda loop ----------------------------
  IF(flmin < 1.0D0) THEN ! THIS is the default...
     flmin = Max (mfl, flmin) ! just sets a threshold above zero
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  ! PRINT *, alf
  vl = matmul(r, x)/nobs ! Note r gets updated in middle and inner loop
  al0 = 0.0D0
  DO g = 1,bn ! For each group... !
     ALLOCATE(u(bs(g)))
     u = vl(ix(g):iy(g))
     ga(g) = sqrt(dot_product(u,u))
     DEALLOCATE(u)
  ENDDO
  DO vl_iter = 1,nvars
     al0 = max(al0, abs(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
  ENDDO
  ! PRINT *, alsparse
  al = al0 ! / (1-alsparse) ! this value ensures all betas are 0 , divide by 1-a?
  l = 0
  tlam = 0.0D0
  DO WHILE (l < nlam) !! This is the start of the loop over all lambda values...
     ! print *, "l = ", l
     ! print *, "al = ", al
     ! IF(kill_count > 1000) RETURN
     ! kill_count = kill_count + 1
     al0 = al ! store old al value on subsequent loops, first set to al
     IF(flmin>=1.0D0) THEN ! user supplied lambda value, break out of everything
        l = l+1
        al=ulam(l)
        ! print *, "This is at the flmin step of the while loop"
     ELSE
        IF(l > 1) THEN ! have some active groups
          al=al*alf
          tlam = max((2.0*al-al0), 0.0) ! Here is the strong rule...
          l = l+1
          ! print *, "This is the l>1 step of while loop"
        ELSE IF(l==0) THEN
          al=al*.99
          tlam=al
          ! Trying to find an active group
        ENDIF
     ENDIF
     lama = al*alsparse
     lam1ma = al*(1-alsparse)
     ! This is the start of the algorithm, for a given lambda...
     ! print *, "Here is tlam = ", tlam
     DO g = 1, bn
        IF(jxx(g) == 1)  CYCLE
        IF(ga(g) > pf(g)*tlam*(1-alsparse)) jxx(g) = 1 ! Implementing the strong rule
     ENDDO
     ! --------- outer loop ---------------------------- !
     DO
        oldbeta(0)=b(0) ! intercept, to be removed
        ! print *, "Here is the outer loop, and here's oldbeta:", oldbeta
        IF(ni>0) THEN
           DO j=1,ni ! sparsity alters the copy here
              g=idx(j)
              oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
           ! print *, "This is where we enter the middle loop"
           npass=npass+1
           dif=0.0D0
           DO g=1,bn
              IF(jxx(g) == 0) CYCLE
              startix=ix(g)
              endix=iy(g)
              ALLOCATE(dd(bs(g)))
              ALLOCATE(oldb(bs(g)))
              oldb=b(startix:endix) ! how do we index into this? we do it repeatedly
                ! maybe replace with a temp?
              ALLOCATE(s(bs(g)))
              s = matmul(r,x(:,startix:endix))/nobs
              s = s*t_for_s(g) + b(startix:endix)
              DO soft_g = 1, bs(g)
                 sg = s(soft_g)
                 s(soft_g) = sign(max(abs(s(soft_g))-lama*t_for_s(g), 0.0D0), s(soft_g))
              ENDDO
              snorm = sqrt(dot_product(s,s))
              tea = snorm - t_for_s(g)*lam1ma*pf(g)
              IF(tea>0.0D0) THEN
                 b(startix:endix) = s*tea/snorm
              ELSE
                 b(startix:endix) = 0.0D0
              ENDIF
              dd=b(startix:endix)-oldb
              IF(any(dd/=0.0D0)) THEN
                 dif=max(dif,gam(g)**2*dot_product(dd,dd))
                 r=r-matmul(x(:,startix:endix),dd)
                 IF(oidx(g)==0) THEN ! Here is where middle loop is different; if group g was not in oidx (active), and the
                    ! difference was nonzero, put it in active (ni)
                    ni=ni+1
                    IF(ni>pmax) EXIT
                    oidx(g)=ni
                    idx(ni)=g
                 ENDIF
              ENDIF
              DEALLOCATE(s,dd,oldb)
              ! DEALLOCATE(u,dd,oldb)
           ENDDO ! End middle loop
           IF(intr /= 0) THEN
              d=sum(r)/nobs
              IF(d/=0.0D0) THEN
                 b(0)=b(0)+d
                 r=r-d
                 dif=max(dif,d**2)
              ENDIF
           ENDIF
           IF (ni > pmax) EXIT
           IF (dif < eps) EXIT
           IF(npass > maxit) THEN
              jerr=-l
              RETURN
           ENDIF
           ! --inner loop----------------------
           DO
              ! PRINT *, "Here is where the inner loop starts"
              npass=npass+1
              dif=0.0D0
              DO j=1,ni
                 g=idx(j)
                 startix=ix(g)
                 endix=iy(g)
                 ALLOCATE(s(bs(g)))
                 ALLOCATE(dd(bs(g)))
                 ALLOCATE(oldb(bs(g)))
                 oldb=b(startix:endix) ! same shitck, still indexing in
                 s = matmul(r,x(:,startix:endix))/nobs
                 s = s*t_for_s(g) + b(startix:endix)
                 DO soft_g = 1, bs(g)
                    s(soft_g) = sign(max(abs(s(soft_g))-lama*t_for_s(g), 0.0D0), s(soft_g))
                 ENDDO
                 snorm = sqrt(dot_product(s,s))
                 tea = snorm - t_for_s(g)*lam1ma*pf(g)
                 IF(tea>0.0D0) THEN
                    b(startix:endix) = s*tea/snorm
                 ELSE
                    b(startix:endix) = 0.0D0
                 ENDIF
                 dd=b(startix:endix)-oldb
                 IF(any(dd/=0.0D0)) THEN
                    dif=max(dif,gam(g)**2*dot_product(dd,dd))
                    r=r-matmul(x(:,startix:endix),dd)
                 ENDIF
                 DEALLOCATE(s,dd,oldb)
                 ! DEALLOCATE(u,dd,oldb)
              ENDDO ! END INNER LOOP
              IF(intr /= 0) THEN ! intr is whether to include intercept
                 d=sum(r)/nobs
                 IF(d/=0.0D0) THEN
                    b(0)=b(0)+d
                    r=r-d
                    dif=max(dif,d**2)
                 ENDIF
              ENDIF
              IF(dif<eps) EXIT ! Exit nearest loop. This is till convergence.
              IF(npass > maxit) THEN
                 jerr=-l
                 RETURN
              ENDIF
           ENDDO ! End Inner loop
        ENDDO ! End middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------ ! This checks which violate KKT condition
        ! PRINT *, "Here is where the final check starts"
        jx = 0
        max_gam = maxval(gam)
        ! need to do sp difference times scalars, sp division
        IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        vl = matmul(r, x)/nobs
        DO g = 1, bn
           IF(jxx(g) == 1) CYCLE
           startix=ix(g)
           endix=iy(g)
           ALLOCATE(s(bs(g)))
           s = matmul(r,x(:,startix:endix))/nobs
           DO soft_g = 1, bs(g)
              s(soft_g) = sign(max(abs(s(soft_g))-lama, 0.0D0), s(soft_g))
           ENDDO
           snorm = sqrt(dot_product(s,s))
           ga(g) = snorm
           DEALLOCATE(s)
           IF(ga(g) > pf(g)*lam1ma)THEN
              jxx(g) = 1
              jx = 1
           ENDIF
        ENDDO
        IF(jx == 1) CYCLE ! This goes all the way back to outer loop
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     IF(l==0) THEN
        IF(maxval(jxx)==0) THEN
           CYCLE ! don't save anything, we're still decrementing lambda
        ELSE
           l=2
           alam(1) = al / max(alf,.99) ! store previous, larger value
        ENDIF
     ENDIF
     ! PRINT *, "Here is where the final update starts"
     IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
     ENDIF
     IF(ni>0) THEN
        DO j=1,ni
           g=idx(j)
           beta(ix(g):iy(g),l)=b(ix(g):iy(g)) ! here's the copy
        ENDDO
     ENDIF
     nbeta(l)=ni
     b0(l)=b(0) ! delete me
     alam(l)=al
     nalam=l
     IF (l < mnl) CYCLE
     me=0
     DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1 ! just checks active groups
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,oidx)
  RETURN
END SUBROUTINE ls_f_sparse

SUBROUTINE strong_rule (njxx, pmax, jxx, bn, ga, pf, tlam, alsparse)
  IMPLICIT NONE
  
  ! This is the strong rule check, unfortunately, jxx should be sorted...
  INTEGER ::  njxx, pmax, bn, posxx, g
  INTEGER :: jxx(pmax)
  DOUBLE PRECISION:: ga(bn) ! this is possibly big, and redundant. Use Dmat structure
  DOUBLE PRECISION::pf(bn)
  DOUBLE PRECISION:: tlam, alsparse
  IF(njxx < 1) RETURN
  posxx = 1
  DO g = 1, bn
    IF(jxx(posxx) == g) THEN
      posxx = posxx + 1
      CYCLE
    ENDIF
    IF(ga(g) > pf(g)*tlam*(1-alsparse)) THEN
      njxx = njxx + 1
      IF(njxx > pmax) RETURN
      IF(posxx > 1) THEN
        jxx(1:njxx) = (/jxx(1:posxx-1),g,jxx(posxx:njxx-1)/) ! shift and insert
      ELSE
        jxx(1:njxx) = (/ g, jxx(1:njxx-1) /)
      ENDIF
      posxx = posxx + 1
    ENDIF
  ENDDO
END SUBROUTINE strong_rule

SUBROUTINE softthresh(vec, thresh, n)
  INTEGER :: n, it
  DOUBLE PRECISION :: vec(n)
  DOUBLE PRECISION :: sg
  DO it=1,n
    sg = vec(it)
    vec(it) = sign(max(abs(sg) - thresh, 0.0D0), sg)
  ENDDO
END SUBROUTINE softthresh



SUBROUTINE kkt_check (strong, nstrong, active, nactive, ix, iy, bs, bn,&
     lama, r, x, nobs, nvars, pf, ga, lam1ma, addactive, dfmax, checkall)
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::nstrong, nobs, nvars, dfmax, bn, checkall
  INTEGER, INTENT(IN) ::ix(bn)
  INTEGER, INTENT(IN) ::iy(bn)
  INTEGER, INTENT(IN) ::bs(bn)
  INTEGER, INTENT(IN) ::strong(dfmax)
  DOUBLE PRECISION, INTENT(IN) :: lama, lam1ma
  DOUBLE PRECISION, INTENT(IN) :: pf(bn)
  DOUBLE PRECISION, INTENT(IN) :: ga(bn) ! this is possibly big, and redundant. Use Dmat structure
  DOUBLE PRECISION, INTENT(IN) :: r(nobs)
  DOUBLE PRECISION, INTENT(IN) :: x(nobs,nvars) ! supposedly means x constant
  INTEGER :: addactive, nactive
  INTEGER :: active(dfmax)
  INTEGER :: posxx, startix, endix, gidx, g, stopper
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s 
  
  posxx = 1
  IF(checkall == 1) THEN
     stopper = bn
  ELSE
     stopper = nstrong
  ENDIF
  DO g = 1, stopper
     IF(checkall==1) THEN
        gidx = g
     ELSE
        gidx = strong(g)
     ENDIF
     IF(active(posxx) == gidx) THEN
        posxx = posxx + 1
        CYCLE
     ENDIF
     startix=ix(gidx)
     endix=iy(gidx)
     ALLOCATE(s(bs(gidx)))
     s = matmul(r,x(:,startix:endix))/nobs
     softthresh(s, lama)
     snorm = sqrt(dot_product(s,s))
     ga(gidx) = snorm
     DEALLOCATE(s)
     IF(ga(gidx) > pf(gidx)*lam1ma) THEN
        nactive = nactive + 1
        addactive = 1
        IF(nactive > dfmax) RETURN
        active(posxx:nactive) = (/gidx, active(posxx:nactive-1)/)
        posxx = posxx + 1
     ENDIF
  ENDDO
END SUBROUTINE kkt_check



SUBROUTINE LOOP_ACTIVE (nactive, active_set, ix, iy, bs, bn, b,&
      bspot, lb, r, x, nobs, ts, lama, lam1ma, pf, gam, dif)

  IMPLICIT NONE
  INTEGER :: nactive, gidx, g, bn, startix, endix, nobs, lb
  INTEGER :: active_set(nactive)
  INTEGER, INTENT(IN) :: bs(bn), ix(bn), iy(bn)
  INTEGER, INTENT(IN) :: bspot(nactive)
  DOUBLE PRECISION :: b(lb)
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: tempb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
  DOUBLE PRECISION :: tea, snorm, ts, lama, lam1ma, dif
  DOUBLE PRECISION :: r(nobs)
  DOUBLE PRECISION, INTENT(IN) :: x(nobs, nvars)
  DOUBLE PRECISION, INTENT(IN)::pf(bn)
  DOUBLE PRECISION, INTENT(IN)::gam(bn)

  DO g=1,nactive
    gidx = active_set(g)
    startix = ix(gidx)
    endix = iy(gidx)
    ALLOCATE(dd(bs(gidx)))
    ALLOCATE(oldb(bs(gidx)))
    ALLOCATE(tempb(bs(gidx)))
    ALLOCATE(s(bs(gidx)))
    tempb = b(bspot(g):bspot(g) + bs(gidx))
    oldb = tempb
    s = matmul(r, x(:,startix:endix))/nobs
    s = s*ts(gidx) + tempb
    softthresh(s, lama*ts(gidx), bs(gidx))
    snorm = sqrt(dot_product(s,s))
    tea = snorm - ts(gidx)*lam1ma*pf(gidx)
    IF(tea > 0.0D0) THEN
      tempb = s*tea/snorm
    ELSE
      tempb = 0.0D0
    ENDIF
    b(bspot(g):bspot(g)+bs(gidx)) = tempb
    dd = tempb - oldb
    IF(any(dd .ne. 0.0D0)) THEN
      dif = max(dif, gam(gidx)**2*dot_(dd,dd))
      r=r-matmul(x(:,startix:endix),dd)
    ENDIF
    DEALLOCATE(dd,oldb,tempb,s)
  ENDDO
END SUBROUTINE LOOP_ACTIVE

SUBROUTINE new_loss ()
  IMPLICIT NONE

  ! Declarations

  ! Initialization

  ! Lambda loop
   IF(flmin < 1.0D0) THEN ! THIS is the default...
     flmin = Max (mfl, flmin) ! just sets a threshold above zero
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  vl = matmul(r, x)/nobs ! Note r gets updated in middle and inner loop, big, but we need it
  al0 = 0.0D0
  DO g = 1,bn ! For each group... !
    ALLOCATE(u(bs(g)))
    u = vl(ix(g):iy(g))
    ga(g) = sqrt(dot_product(u,u))
    DEALLOCATE(u)
  ENDDO
  DO vl_iter = 1,nvars
    al0 = max(al0, abs(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
  ENDDO
  al = al0 ! / (1-alsparse) ! this value ensures all betas are 0
  l = 0
  tlam = 0.0D0
  DO WHILE (l < nlam) !! This is the start of the loop over all lambda values...
    al0 = al ! store old al value on subsequent loops, first set to al
    IF(flmin>=1.0D0) THEN ! user supplied lambda value, break out of everything
      l = l+1
      al=ulam(l)
    ELSE
      IF(l > 1) THEN ! have some active groups
        al=al*alf
        tlam = max((2.0*al-al0), 0.0) ! Here is the strong rule...
        l = l+1
      ELSE IF(l==0) THEN ! Trying to find an active group
        al=al*.99
        tlam=al
      ENDIF
    ENDIF
    lama = al*alsparse
    lam1ma = al*(1-alsparse)
    ! This is the start of the algorithm, for a given lambda...
    ! We no longer check the strong rule first, instead, we iterate over A(lam) (ever-active set)
    DO
       IF(nactive>0) THEN
          oldbeta(:,1:ni) = b(:,1:ni) ! here's a copy to sparsify
       ENDIF
       DO ! print *, "This is where we enter the middle loop"
          npass=npass+1
          dif=0.0D0
          LOOP_ACTIVE (nactive, active, ix, iy, bs, bn, b,&
               bspot, lb, r, x, nobs, ts, lama, lam1ma, pf, gam, dif)
          IF (nactive > pmax) EXIT
          IF (dif < eps) EXIT
          IF(npass > maxit) THEN
             jerr=-l
             RETURN
          ENDIF
  

END SUBROUTINE NEW_LOSS


!FUNCTION mulvecD(b, D, drows, dcols, belem)
!  IMPLICIT NONE
!
!  INTEGER, INTENT(IN) :: drows, dcols, belem
!  DOUBLE PRECISION, INTENT(IN) :: D(drows, dcols)
!  DOUBLE PRECISION, INTENT(IN) :: b(belem)
!  DOUBLE PRECISION, INTENT(OUT) :: z(drows)
!
!  INTEGER :: i,j,k
!
!  z = 0
!
!  DO j=0,belem-1
!    k = MOD(j, dcols) + 1
!    DO i=1,drows
!      z(i) += D(i,k) * b(j)
!    ENDDO
!  ENDDO
!
!  RETURN
!END FUNCTION






! --------------------------------------------------
SUBROUTINE ls_f_sparse_beta (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,nalam,beta,idx,nbeta,alam,npass,jerr,alsparse)
  ! --------------------------------------------------
  
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: big=9.9E30
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
  INTEGER, PARAMETER :: mnlam = 6
  INTEGER:: mnl
  INTEGER:: bn
  INTEGER::bs(bn)
  INTEGER, INTENT(IN) ::ix(bn)
  INTEGER, INTENT(IN) ::iy(bn)
  INTEGER, INTENT(IN) :: nobs
  INTEGER, INTENT(IN) ::nvars
  INTEGER, INTENT(IN) ::dfmax
  INTEGER, INTENT(IN) ::pmax
  INTEGER, INTENT(IN) ::nlam
  INTEGER::nalam
  INTEGER::npass
  INTEGER::jerr
  INTEGER::maxit
  INTEGER:: idx(pmax)
  INTEGER::nbeta(nlam)
  DOUBLE PRECISION:: flmin
  DOUBLE PRECISION::eps
  DOUBLE PRECISION, INTENT(IN) :: x(nobs,nvars) ! supposedly means x constant
  DOUBLE PRECISION, INTENT(IN) :: y(nobs)
  DOUBLE PRECISION, INTENT(IN) ::pf(bn)
  DOUBLE PRECISION::ulam(nlam)
  DOUBLE PRECISION, INTENT(IN) ::gam(bn)
  DOUBLE PRECISION::beta(nvars,nlam) ! this will be sparsified
  DOUBLE PRECISION::alam(nlam)
  DOUBLE PRECISION::alsparse ! Should be alpha for the sparsity weight
  ! - - - local declarations - - -
  DOUBLE PRECISION:: max_gam
  DOUBLE PRECISION::d
  ! DOUBLE PRECISION::t ! No longer using this
  DOUBLE PRECISION::dif
  ! DOUBLE PRECISION::unorm ! No longer using this for ls_new
  DOUBLE PRECISION::al
  DOUBLE PRECISION::alf
  DOUBLE PRECISION::sg
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r ! Residue y-beta_k*x etc
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: tempb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u ! No longer using this for update step, but still need for other parts
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
  INTEGER:: g
  INTEGER::j
  INTEGER::l
  INTEGER::ni
  INTEGER::me
  INTEGER::startix
  INTEGER::endix
  ! - - - Aaron's declarations
  DOUBLE PRECISION::snorm
  DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma
  DOUBLE PRECISION::tea ! this takes the place of 't' in the update step for ls
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s ! takes the place of 'u' in update for ls
  INTEGER::soft_g ! this is an iterating variable 'vectorizing' the soft thresholding operator
  INTEGER::vl_iter ! for iterating over columns(?) of x*r
  INTEGER::kill_count
  ! - - - begin - - -
  ! - - - local declarations - - -
  DOUBLE PRECISION:: tlam
  DOUBLE PRECISION:: lama
  DOUBLE PRECISION:: lam1ma
  INTEGER:: jx, bsmax, grind
  INTEGER:: jxx(dfmax) ! was length bn, 0/1 if active, making list of actives
  INTEGER:: njxx ! place to store how many are nonzero
  DOUBLE PRECISION:: ga(bn)
  DOUBLE PRECISION:: vl(nvars) ! this is big, what can we do?
  DOUBLE PRECISION:: al0
  bsmax = maxval(bs)
  ! - - - allocate variables - - -
  ALLOCATE(b(1:dfmax*bsmax)) ! make them 1 index (ditch intercept), to dfmax
  ALLOCATE(oldbeta(1:dfmax*bsmax)) ! make them 1 index (ditch intercept)
  ALLOCATE(r(1:nobs))
  ALLOCATE(oidx(1:bn)) ! what's the deal here?
  ! - - - checking pf - ! pf is the relative penalties for each group
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  max_gam = maxval(gam)
  jxx = 0
  njxx = 0
  onjxx = 0
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  idx = 0
  kill_count = 0
  oidx = 0
  npass = 0 ! This is a count, correct?
  ni = npass ! This controls so-called "outer loop"
  alf = 0.0D0
  t_for_s = 1/gam
  ! --------- lambda loop ----------------------------
  IF(flmin < 1.0D0) THEN ! THIS is the default...
     flmin = Max (mfl, flmin) ! just sets a threshold above zero
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  ! PRINT *, alf
  vl = matmul(r, x)/nobs ! Note r gets updated in middle and inner loop, big, but we need it
  al0 = 0.0D0
  DO g = 1,bn ! For each group... !
    ALLOCATE(u(bs(g)))
    u = vl(ix(g):iy(g))
    ga(g) = sqrt(dot_product(u,u))
    DEALLOCATE(u)
  ENDDO
  DO vl_iter = 1,nvars
    al0 = max(al0, abs(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
  ENDDO
  al = al0 ! / (1-alsparse) ! this value ensures all betas are 0
  l = 0
  tlam = 0.0D0
  DO WHILE (l < nlam) !! This is the start of the loop over all lambda values...
    al0 = al ! store old al value on subsequent loops, first set to al
    IF(flmin>=1.0D0) THEN ! user supplied lambda value, break out of everything
      l = l+1
      al=ulam(l)
    ELSE
      IF(l > 1) THEN ! have some active groups
        al=al*alf
        tlam = max((2.0*al-al0), 0.0) ! Here is the strong rule...
        l = l+1
      ELSE IF(l==0) THEN ! Trying to find an active group
        al=al*.99
        tlam=al
      ENDIF
    ENDIF
    lama = al*alsparse
    lam1ma = al*(1-alsparse)
    ! This is the start of the algorithm, for a given lambda...
    ! We no longer check the strong rule first, instead, we iterate over A(lam) (ever-active set)
    ! OLD: This is the strong rule check, unfortunately, jxx should be sorted...
    ! OLD: strong_rule(njxx, dfmax, jxx, bn, ga, pf, tlam, alsparse)
    ! --------- outer loop ---------------------------- !
    DO
      IF(ni>0) THEN
        oldbeta(:,1:ni) = b(:,1:ni)
      ENDIF
      ! --middle loop-------------------------------------
      DO ! print *, "This is where we enter the middle loop"
        npass=npass+1
        dif=0.0D0
        DO g=1,njxx
          grind = jxx(g)
          startix=ix(grind)
          endix=iy(grind)
          ALLOCATE(dd(bs(grind)))
          ALLOCATE(oldb(bs(grind)))
          ALLOCATE(tempb(bs(grind)))
          tempb=b(1:bs(grind),g)
          oldb=tempb
          ALLOCATE(s(bs(grind)))
          s = matmul(r,x(:,startix:endix))/nobs
          s = s*t_for_s(grind) + tempb
          DO soft_g = 1, bs(grind)
            sg = s(soft_g)
            s(soft_g) = sign(max(abs(sg)-lama*t_for_s(grind), 0.0D0), sg)
          ENDDO
          snorm = sqrt(dot_product(s,s))
          tea = snorm - t_for_s(grind)*lam1ma*pf(grind)
          IF(tea>0.0D0) THEN
            tempb = s*tea/snorm
          ELSE
            tempb = 0.0D0
          ENDIF
          dd=tempb-oldb
          IF(any(dd/=0.0D0)) THEN
            dif=max(dif,gam(grind)**2*dot_product(dd,dd))
            r=r-matmul(x(:,startix:endix),dd)
            IF(oidx(g)==0) THEN ! Here is where middle loop is different;
              ! if group g was not in oidx (active), and the
              ! difference was nonzero, put it in active (ni)
              ni=ni+1
              IF(ni>pmax) EXIT
              oidx(g)=ni
              idx(ni)=g
            ENDIF
          ENDIF
          DEALLOCATE(s,dd,oldb)
        ENDDO ! End middle loop
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
          jerr=-l
          RETURN
        ENDIF
        ! --inner loop----------------------
        DO ! PRINT *, "Here is where the inner loop starts"
          npass=npass+1
          dif=0.0D0
          DO j=1,ni
            g=idx(j)
            startix=ix(g)
            endix=iy(g)
            ALLOCATE(s(bs(g)))
            ALLOCATE(dd(bs(g)))
            ALLOCATE(oldb(bs(g)))
            oldb=b(startix:endix) ! same shitck, still indexing in
            s = matmul(r,x(:,startix:endix))/nobs
            s = s*t_for_s(g) + b(startix:endix)
            DO soft_g = 1, bs(g)
              s(soft_g) = sign(max(abs(s(soft_g))-lama*t_for_s(g), 0.0D0), s(soft_g))
            ENDDO
            snorm = sqrt(dot_product(s,s))
            tea = snorm - t_for_s(g)*lam1ma*pf(g)
            IF(tea>0.0D0) THEN
              b(startix:endix) = s*tea/snorm
            ELSE
              b(startix:endix) = 0.0D0
            ENDIF
            dd=b(startix:endix)-oldb
            IF(any(dd/=0.0D0)) THEN
              dif=max(dif,gam(g)**2*dot_product(dd,dd))
              r=r-matmul(x(:,startix:endix),dd)
            ENDIF
            DEALLOCATE(s,dd,oldb)
          ENDDO
          IF(dif<eps) EXIT ! Exit nearest loop. This is till convergence.
          IF(npass > maxit) THEN
            jerr=-l
            RETURN
          ENDIF
        ENDDO ! End Inner loop
      ENDDO ! End middle loop
      IF(ni>pmax) EXIT
      !--- final check ------------------------ ! This checks which violate KKT condition
      jx = 0
      
      ! need to do sp difference times scalars, sp division
      IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
      IF (jx /= 0) CYCLE
      vl = matmul(r, x)/nobs
      kkt_check (bn, ix, iy, bs, lama, r, x, nobs, nvars, pf, ga, lam1ma, jx, jxx, njxx, dfmax)
      IF(njxx > dfmax) EXIT
      IF(jx == 1) CYCLE ! This goes all the way back to outer loop
      EXIT
    ENDDO ! Ends outer loop
    !---------- final update variable and save results------------
    IF(l==0) THEN
      IF(njxx==0) THEN
        CYCLE ! don't save anything, we're still decrementing lambda
      ELSE
        l=2
        alam(1) = al / max(alf,.99) ! store previous, larger value
      ENDIF
    ENDIF
    ! PRINT *, "Here is where the final update starts"
    IF(ni>pmax .or. njxx > dfmax) THEN
      jerr=-10000-l
      EXIT
    ENDIF
    IF(ni>0) THEN
      DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g)) ! here's the copy
      ENDDO
    ENDIF
    nbeta(l)=ni
    alam(l)=al
    nalam=l
    IF (l < mnl) CYCLE
    me=0
    DO j=1,ni
      g=idx(j)
      IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1 ! just counts active groups
    ENDDO
    IF(me > dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,oidx)
  RETURN
END SUBROUTINE ls_f_sparse
