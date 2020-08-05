module mySubby

contains
!
!----------------------------------------------------
subroutine strong_rule (is_in_E_set, ga, pf, tlam, alsparse)
        !-------------------------------------------
        implicit none
        integer :: g
        integer, dimension (:), intent(inout) :: is_in_E_set
        double precision, dimension (:), intent(in) :: ga
        double precision, dimension (:), intent(in) :: pf
        double precision, intent(in) :: tlam, alsparse
        !------------------------------------------
        do g = 1, size(is_in_E_set)
                if(is_in_E_set(g) == 1) cycle
                if(ga(g) > pf(g)*tlam*(1-alsparse)) is_in_E_set(g) = 1
        enddo
        RETURN
end subroutine strong_rule


!--------------------------------------------
SUBROUTINE softthresh(vec, thresh, n)
        !------------------------------------
  implicit none
  INTEGER :: n, it
  double precision :: sg
  DOUBLE PRECISION :: vec(n)
  double precision :: thresh
  DO it=1,n
    sg = vec(it)
    vec(it) = sign(max(abs(sg) - thresh, 0.0D0), sg)
  ENDDO
  RETURN
END SUBROUTINE softthresh



!----------------------------------------------
subroutine kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama)
        implicit none
        integer :: g, startix, endix
        integer, intent(in) :: bn
        INTEGER, intent(in) ::bs(bn)
        integer, intent(in) :: ix(bn), iy(bn)
        integer, dimension(:), intent(inout) :: is_in_E_set
        double precision, dimension (:), intent(in) :: vl
        double precision, dimension (:), allocatable :: s
        double precision :: snorm
        double precision, intent(in) :: pf(bn)
        integer, intent(inout) :: violation
        double precision, intent(in) :: lam1ma, lama
        !------------------------
        do g = 1, bn
                if(is_in_E_set(g) ==1) cycle
                startix = ix(g)
                endix = iy(g)
                allocate(s(bs(g)))
                s = vl(startix:endix)
                call softthresh(s, lama, bs(g))
                snorm = sqrt(dot_product(s,s))
                if(snorm > pf(g)*lam1ma) then
                        is_in_E_set(g) = 1
                        violation = 1
                endif
                deallocate(s)
        enddo
        RETURN
end subroutine kkt_check

!---------------------------------------
subroutine update_step(bsg, startix, endix, b, lama, t_for_sg, pfg, lam1ma, x,&
    isDifZero, nobs, r, gamg, maxDif,nvars)
        implicit none
        integer, intent(in) :: bsg, nobs, nvars
        integer, intent(in) :: startix, endix
        double precision :: gamg
        double precision, intent(inout) :: maxDif
        double precision, dimension (:), allocatable :: oldb, s, dd
        double precision, dimension (:), intent(inout) :: b, r
        double precision :: snorm, tea
        double precision, intent(in) :: lama, t_for_sg, pfg, lam1ma
        double precision, dimension (:), intent(in) :: x(nobs,nvars)
        integer, intent(inout) :: isDifZero
        !------------------------------
        allocate(s(bsg))
        allocate(oldb(bsg))
        isDifZero = 0
        oldb = b(startix:endix)
        s = matmul(r, x(:, startix:endix))/nobs
        s = s*t_for_sg + b(startix:endix)
        call softthresh(s, lama*t_for_sg, bsg)
        snorm = sqrt(dot_product(s,s))
        tea = snorm - t_for_sg*lam1ma*pfg
        if(tea > 0.0D0) then
                b(startix:endix) = s*tea/snorm
        else
                b(startix:endix) = 0.0D0
        endif
        allocate(dd(bsg))
        dd = b(startix:endix) - oldb
        if(any(dd/=0.0D0)) then
                maxDif = max(maxDif,gamg**2*dot_product(dd,dd))
                r=r-matmul(x(:,startix:endix),dd)
                isDifZero = 1
        endif
        deallocate(s, oldb, dd)
        RETURN
end subroutine update_step

end module mySubby

! --------------------------------------------------
SUBROUTINE sparse_three (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,nalam,beta,activeGroup,nbeta,alam,npass,jerr,alsparse)
  ! --------------------------------------------------
  USE mySubby
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: big=9.9E30
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
  INTEGER, PARAMETER :: mnlam = 6
  INTEGER:: isDifZero
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
  INTEGER:: activeGroup(pmax)
  INTEGER::nbeta(nlam)
  DOUBLE PRECISION:: flmin
  DOUBLE PRECISION::eps
  DOUBLE PRECISION:: x(nobs,nvars)
  DOUBLE PRECISION::y(nobs)
  DOUBLE PRECISION::pf(bn)
  DOUBLE PRECISION::ulam(nlam)
  DOUBLE PRECISION::gam(bn)
  DOUBLE PRECISION::beta(nvars,nlam)
  DOUBLE PRECISION::alam(nlam)
  DOUBLE PRECISION::alsparse ! Should be alpha for the sparsity weight
  ! - - - local declarations - - -
  DOUBLE PRECISION:: max_gam
  ! DOUBLE PRECISION::d
  ! DOUBLE PRECISION::t ! No longer using this
  DOUBLE PRECISION::maxDif
  ! DOUBLE PRECISION::unorm ! No longer using this for ls_new
  DOUBLE PRECISION::al
  DOUBLE PRECISION::alf
  ! DOUBLE PRECISION::sg
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r ! Residue y-beta_k*x etc
  ! DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u ! No longer using this for update step, but still need for other parts
  ! DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  INTEGER, DIMENSION (:), ALLOCATABLE :: activeGroupIndex
  INTEGER:: g
  INTEGER::j
  INTEGER::l
  INTEGER::ni
  INTEGER::me
  INTEGER::startix
  INTEGER::endix
  ! - - - Aaron's declarations
  ! DOUBLE PRECISION::snorm
  DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma
  ! DOUBLE PRECISION::tea ! this takes the place of 't' in the update step for ls
  ! DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s ! takes the place of 'u' in update for ls
  ! INTEGER::soft_g ! this is an iterating variable 'vectorizing' the soft thresholding operator
  INTEGER::vl_iter ! for iterating over columns(?) of x*r
  INTEGER::kill_count
  ! - - - begin - - -
  ! - - - local declarations - - -
  DOUBLE PRECISION:: tlam
  DOUBLE PRECISION:: lama
  DOUBLE PRECISION:: lam1ma
  INTEGER:: violation
  INTEGER:: is_in_E_set(bn)
  DOUBLE PRECISION:: ga(bn) ! What is this for??
  DOUBLE PRECISION:: vl(nvars) ! What is this for?
  DOUBLE PRECISION:: al0
  ! - - - allocate variables - - -
  ALLOCATE(b(1:nvars))
  ALLOCATE(oldbeta(1:nvars))
  ALLOCATE(r(1:nobs))
  ALLOCATE(activeGroupIndex(1:bn))
  !    ALLOCATE(al_sparse)
  ! - - - checking pf - ! pf is the relative penalties for each group
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  is_in_E_set = 0
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  activeGroup = 0
  kill_count = 0
  activeGroupIndex = 0
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
  DO g = 1,bn ! For each group...
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
     call strong_rule (is_in_E_set, ga, pf, tlam, alsparse) !implementing strong rule, updates is_in_E_set
     ! --------- outer loop ---------------------------- !
     DO
        IF(ni>0) THEN
           DO j=1,ni
              g=activeGroup(j)
              oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
           ! print *, "This is where we enter the middle loop"
           npass=npass+1
           maxDif=0.0D0
           ! isDifZero = 0 !Boolean to check if b-oldb nonzero. Unnec, in fn.
           DO g=1,bn
              IF(is_in_E_set(g) == 0) CYCLE
              startix=ix(g)
              endix=iy(g)
              call update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g), lam1ma, x,&
                  isDifZero, nobs, r, gam(g), maxDif, nvars)
              IF(activeGroupIndex(g)==0 .and. isDifZero == 1) THEN
                    ni=ni+1
                    IF(ni>pmax) EXIT
                    activeGroupIndex(g)=ni
                    activeGroup(ni)=g
              ENDIF
           ENDDO ! End middle loop
           IF (ni > pmax) EXIT
           IF (maxDif < eps) EXIT
           IF(npass > maxit) THEN !Is this needed?
              jerr=-l
              RETURN
           ENDIF
           ! --inner loop----------------------
          DO
             ! PRINT *, "Here is where the inner loop starts"
             npass=npass+1
             maxDif=0.0D0
             isDifZero = 0
             DO j=1,ni
                g=activeGroup(j)
                startix=ix(g)
                endix=iy(g)
                call update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g),&
                   lam1ma, x, isDifZero, nobs, r, gam(g), maxDif, nvars)
             ENDDO ! END INNER LOOP
             IF(maxDif<eps) EXIT ! Exit nearest loop. This is till convergence.
             IF(npass > maxit) THEN
                jerr=-l
                RETURN
             ENDIF
          ENDDO ! End Inner loop
        ENDDO ! End middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------ ! This checks which violate KKT condition
        ! PRINT *, "Here is where the final check starts"
        violation = 0
        max_gam = maxval(gam)
        IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) violation = 1
        IF (violation == 1) CYCLE
        vl = matmul(r, x)/nobs
        call kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama) ! kkt subroutine
        IF(violation == 1) CYCLE ! This goes all the way back to outer loop
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     IF(l==0) THEN
        IF(maxval(is_in_E_set)==0) THEN
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
           g=activeGroup(j)
           beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
     ENDIF
     nbeta(l)=ni
     alam(l)=al
     nalam=l
     IF (l < mnl) CYCLE
     me=0
     DO j=1,ni
        g=activeGroup(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,activeGroupIndex)
  RETURN
END SUBROUTINE sparse_three

! --------------------------------------------------
SUBROUTINE sparse_four (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,nalam,beta,activeGroup,nbeta,alam,npass,jerr,alsparse)
  ! --------------------------------------------------
  USE mySubby
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: big=9.9E30
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
  INTEGER, PARAMETER :: mnlam = 6
  INTEGER:: isDifZero
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
  INTEGER:: activeGroup(pmax)
  INTEGER::nbeta(nlam)
  DOUBLE PRECISION:: flmin
  DOUBLE PRECISION::eps
  DOUBLE PRECISION:: x(nobs,nvars)
  DOUBLE PRECISION::y(nobs)
  DOUBLE PRECISION::pf(bn)
  DOUBLE PRECISION::ulam(nlam)
  DOUBLE PRECISION::gam(bn)
  DOUBLE PRECISION::beta(nvars,nlam)
  DOUBLE PRECISION::alam(nlam)
  DOUBLE PRECISION::alsparse ! Should be alpha for the sparsity weight
  ! - - - local declarations - - -
  DOUBLE PRECISION:: max_gam
  ! DOUBLE PRECISION::d
  ! DOUBLE PRECISION::t ! No longer using this
  DOUBLE PRECISION::maxDif
  ! DOUBLE PRECISION::unorm ! No longer using this for ls_new
  DOUBLE PRECISION::al
  DOUBLE PRECISION::alf
  ! DOUBLE PRECISION::sg
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r ! Residue y-beta_k*x etc
  ! DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u ! No longer using this for update step, but still need for other parts
  ! DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  INTEGER, DIMENSION (:), ALLOCATABLE :: activeGroupIndex
  INTEGER:: g
  INTEGER::j
  INTEGER::l
  INTEGER::ni
  INTEGER::me
  INTEGER::startix
  INTEGER::endix
  ! - - - Aaron's declarations
  ! DOUBLE PRECISION::snorm
  DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma
  ! DOUBLE PRECISION::tea ! this takes the place of 't' in the update step for ls
  ! DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s ! takes the place of 'u' in update for ls
  ! INTEGER::soft_g ! this is an iterating variable 'vectorizing' the soft thresholding operator
  INTEGER::vl_iter ! for iterating over columns(?) of x*r
  INTEGER::kill_count
  ! - - - begin - - -
  ! - - - local declarations - - -
  DOUBLE PRECISION:: tlam
  DOUBLE PRECISION:: lama
  DOUBLE PRECISION:: lam1ma
  INTEGER:: violation
  INTEGER:: is_in_E_set(bn)
  DOUBLE PRECISION:: ga(bn) ! What is this for??
  DOUBLE PRECISION:: vl(nvars) ! What is this for?
  DOUBLE PRECISION:: al0
  ! - - - allocate variables - - -
  ALLOCATE(b(1:nvars))
  ALLOCATE(oldbeta(1:nvars))
  ALLOCATE(r(1:nobs))
  ALLOCATE(activeGroupIndex(1:bn))
  !    ALLOCATE(al_sparse)
  ! - - - checking pf - ! pf is the relative penalties for each group
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  is_in_E_set = 0
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  activeGroup = 0
  kill_count = 0
  activeGroupIndex = 0
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
  DO g = 1,bn ! For each group...
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
     call strong_rule (is_in_E_set, ga, pf, tlam, alsparse) !implementing strong rule, updates is_in_E_set
     ! --------- outer loop ---------------------------- !
     DO
        IF(ni>0) THEN
           DO j=1,ni
              g=activeGroup(j)
              oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
           ! print *, "This is where we enter the middle loop"
           npass=npass+1
           maxDif=0.0D0
           isDifZero = 0 !Boolean to check if b-oldb nonzero
           DO g=1,bn
              IF(is_in_E_set(g) == 0) CYCLE
              startix=ix(g)
              endix=iy(g)
              call update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g), lam1ma, x,&
                  isDifZero, nobs, r, gam(g), maxDif, nvars)
              IF(activeGroupIndex(g)==0 .and. isDifZero == 1) THEN
                    ni=ni+1
                    IF(ni>pmax) EXIT
                    activeGroupIndex(g)=ni
                    activeGroup(ni)=g
              ENDIF
           ENDDO ! End middle loop
           IF (ni > pmax) EXIT
           IF (maxDif < eps) EXIT
           IF(npass > maxit) THEN !Is this needed?
              jerr=-l
              RETURN
           ENDIF
           ! --inner loop----------------------
          DO
             ! PRINT *, "Here is where the inner loop starts"
             npass=npass+1
             maxDif=0.0D0
             isDifZero = 0
             DO j=1,ni
                g=activeGroup(j)
                startix=ix(g)
                endix=iy(g)
                call update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g),&
                   lam1ma, x, isDifZero, nobs, r, gam(g), maxDif, nvars)
             ENDDO ! END INNER LOOP
             IF(maxDif<eps) EXIT ! Exit nearest loop. This is till convergence.
             IF(npass > maxit) THEN
                jerr=-l
                RETURN
             ENDIF
          ENDDO ! End Inner loop
        ENDDO ! End middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------ ! This checks which violate KKT condition
        ! PRINT *, "Here is where the final check starts"
        violation = 0
        max_gam = maxval(gam)
        IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) violation = 1
        IF (violation == 1) CYCLE
        vl = matmul(r, x)/nobs
        call kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama) ! kkt subroutine
        IF(violation == 1) CYCLE ! This goes all the way back to outer loop
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     IF(l==0) THEN
        IF(maxval(is_in_E_set)==0) THEN
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
           g=activeGroup(j)
           beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
     ENDIF
     nbeta(l)=ni
     alam(l)=al
     nalam=l
     IF (l < mnl) CYCLE
     me=0
     DO j=1,ni
        g=activeGroup(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,activeGroupIndex)
  RETURN
END SUBROUTINE sparse_four

! --------------------------------------------------
SUBROUTINE gglasso (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,intr,nalam,b0,beta,activeGroup,nbeta,alam,npass,jerr)
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
  INTEGER:: activeGroup(pmax)
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
  DOUBLE PRECISION::maxDif
  DOUBLE PRECISION::unorm
  DOUBLE PRECISION::al
  DOUBLE PRECISION::alf
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  INTEGER, DIMENSION (:), ALLOCATABLE :: activeGroupIndex
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
  INTEGER:: violation
  INTEGER:: is_in_E_set(bn)
  DOUBLE PRECISION:: ga(bn)
  DOUBLE PRECISION:: vl(nvars)
  DOUBLE PRECISION:: al0
  ! - - - allocate variables - - -
  ALLOCATE(b(0:nvars))
  ALLOCATE(oldbeta(0:nvars))
  ALLOCATE(r(1:nobs))
  ALLOCATE(activeGroupIndex(1:bn))
  ! - - - checking pf - - -
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  is_in_E_set = 0
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  activeGroup = 0
  activeGroupIndex = 0
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
        IF(is_in_E_set(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) is_in_E_set(g) = 1
     ENDDO
     ! --------- outer loop ----------------------------
     DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
           DO j=1,ni
              g=activeGroup(j)
              oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
           npass=npass+1
           maxDif=0.0D0
           DO g=1,bn
              IF(is_in_E_set(g) == 0) CYCLE
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
                 maxDif=max(maxDif,gam(g)**2*dot_product(dd,dd))
                 r=r-matmul(x(:,startix:endix),dd)
                 IF(activeGroupIndex(g)==0) THEN
                    ni=ni+1
                    IF(ni>pmax) EXIT
                    activeGroupIndex(g)=ni
                    activeGroup(ni)=g
                 ENDIF
              ENDIF
              DEALLOCATE(u,dd,oldb)
           ENDDO
           IF(intr /= 0) THEN
              d=sum(r)/nobs
              IF(d/=0.0D0) THEN
                 b(0)=b(0)+d
                 r=r-d
                 maxDif=max(maxDif,d**2)
              ENDIF
           ENDIF
           IF (ni > pmax) EXIT
           IF (maxDif < eps) EXIT
           IF(npass > maxit) THEN
              jerr=-l
              RETURN
           ENDIF
           ! --inner loop----------------------
           DO
              npass=npass+1
              maxDif=0.0D0
              DO j=1,ni
                 g=activeGroup(j)
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
                    maxDif=max(maxDif,gam(g)**2*dot_product(dd,dd))
                    r=r-matmul(x(:,startix:endix),dd)
                 ENDIF
                 DEALLOCATE(u,dd,oldb)
              ENDDO
              IF(intr /= 0) THEN
                 d=sum(r)/nobs
                 IF(d/=0.0D0) THEN
                    b(0)=b(0)+d
                    r=r-d
                    maxDif=max(maxDif,d**2)
                 ENDIF
              ENDIF
              IF(maxDif<eps) EXIT
              IF(npass > maxit) THEN
                 jerr=-l
                 RETURN
              ENDIF
           ENDDO ! end inner loop
        ENDDO ! end middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        violation = 0
        max_gam = maxval(gam)
        IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) violation = 1
        IF (violation /= 0) CYCLE
        vl = matmul(r, x)/nobs
        DO g = 1, bn
           IF(is_in_E_set(g) == 1) CYCLE
           ALLOCATE(u(bs(g)))
           u = vl(ix(g):iy(g))
           ga(g) = sqrt(dot_product(u,u))
           IF(ga(g) > al*pf(g))THEN
              is_in_E_set(g) = 1
              violation = 1
           ENDIF
           DEALLOCATE(u)
        ENDDO
        IF(violation == 1) CYCLE
        EXIT
     ENDDO ! End outer loop
     !---------- final update variable and save results------------
     IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
     ENDIF
     IF(ni>0) THEN
        DO j=1,ni
           g=activeGroup(j)
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
        g=activeGroup(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,activeGroupIndex)
  RETURN
END SUBROUTINE gglasso

! --------------------------------------------------
SUBROUTINE sparse_orig (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,nalam,beta,activeGroup,nbeta,alam,npass,jerr,alsparse)
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
  INTEGER:: activeGroup(pmax)
  INTEGER::nbeta(nlam)
  DOUBLE PRECISION:: flmin
  DOUBLE PRECISION::eps
  DOUBLE PRECISION:: x(nobs,nvars)
  DOUBLE PRECISION::y(nobs)
  DOUBLE PRECISION::pf(bn)
  DOUBLE PRECISION::ulam(nlam)
  DOUBLE PRECISION::gam(bn)
  DOUBLE PRECISION::beta(nvars,nlam)
  DOUBLE PRECISION::alam(nlam)
  DOUBLE PRECISION::alsparse ! Should be alpha for the sparsity weight
  ! - - - local declarations - - -
  DOUBLE PRECISION:: max_gam
  ! DOUBLE PRECISION::d
  ! DOUBLE PRECISION::t ! No longer using this
  DOUBLE PRECISION::maxDif
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
  INTEGER, DIMENSION (:), ALLOCATABLE :: activeGroupIndex
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
  INTEGER:: violation
  INTEGER:: is_in_E_set(bn)
  DOUBLE PRECISION:: ga(bn) ! What is this for??
  DOUBLE PRECISION:: vl(nvars) ! What is this for?
  DOUBLE PRECISION:: al0
  ! - - - allocate variables - - -
  ALLOCATE(b(1:nvars))
  ALLOCATE(oldbeta(1:nvars))
  ALLOCATE(r(1:nobs))
  ALLOCATE(activeGroupIndex(1:bn))
  !    ALLOCATE(al_sparse)
  ! - - - checking pf - ! pf is the relative penalties for each group
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  is_in_E_set = 0
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  activeGroup = 0
  kill_count = 0
  activeGroupIndex = 0
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
  DO g = 1,bn ! For each group...
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
        IF(is_in_E_set(g) == 1)  CYCLE
        IF(ga(g) > pf(g)*tlam*(1-alsparse)) is_in_E_set(g) = 1 ! Implementing the strong rule
     ENDDO
     ! --------- outer loop ---------------------------- !
     DO
        ! print *, "Here is the outer loop, and here's oldbeta:", oldbeta
        IF(ni>0) THEN
           DO j=1,ni
              g=activeGroup(j)
              oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
           ! print *, "This is where we enter the middle loop"
           npass=npass+1
           maxDif=0.0D0
           DO g=1,bn
              IF(is_in_E_set(g) == 0) CYCLE
              startix=ix(g)
              endix=iy(g)
              ALLOCATE(dd(bs(g)))
              ALLOCATE(oldb(bs(g)))
              oldb=b(startix:endix)
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
                 maxDif=max(maxDif,gam(g)**2*dot_product(dd,dd))
                 r=r-matmul(x(:,startix:endix),dd)
                 IF(activeGroupIndex(g)==0) THEN ! Here is where middle loop is different; if group g was not in activeGroupIndex (active), and the
                    ! difference was nonzero, put it in active (ni)
                    ni=ni+1
                    IF(ni>pmax) EXIT
                    activeGroupIndex(g)=ni
                    activeGroup(ni)=g
                 ENDIF
              ENDIF
              DEALLOCATE(s,dd,oldb)
              ! DEALLOCATE(u,dd,oldb)
           ENDDO ! End middle loop
           IF (ni > pmax) EXIT
           IF (maxDif < eps) EXIT
           IF(npass > maxit) THEN
              jerr=-l
              RETURN
           ENDIF
           ! --inner loop----------------------
           DO
              ! PRINT *, "Here is where the inner loop starts"
              npass=npass+1
              maxDif=0.0D0
              DO j=1,ni
                 g=activeGroup(j)
                 startix=ix(g)
                 endix=iy(g)
                 ALLOCATE(s(bs(g)))
                 ALLOCATE(dd(bs(g)))
                 ALLOCATE(oldb(bs(g)))
                 oldb=b(startix:endix)
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
                    maxDif=max(maxDif,gam(g)**2*dot_product(dd,dd))
                    r=r-matmul(x(:,startix:endix),dd)
                 ENDIF
                 DEALLOCATE(s,dd,oldb)
                 ! DEALLOCATE(u,dd,oldb)
              ENDDO ! END INNER LOOP
              IF(maxDif<eps) EXIT ! Exit nearest loop. This is till convergence.
              IF(npass > maxit) THEN
                 jerr=-l
                 RETURN
              ENDIF
           ENDDO ! End Inner loop
        ENDDO ! End middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------ ! This checks which violate KKT condition
        ! PRINT *, "Here is where the final check starts"
        violation = 0
        max_gam = maxval(gam)
        IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) violation = 1
        IF (violation /= 0) CYCLE
        vl = matmul(r, x)/nobs
        DO g = 1, bn
           IF(is_in_E_set(g) == 1) CYCLE
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
              is_in_E_set(g) = 1
              violation = 1
           ENDIF
        ENDDO
        IF(violation == 1) CYCLE ! This goes all the way back to outer loop
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     IF(l==0) THEN
        IF(maxval(is_in_E_set)==0) THEN
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
           g=activeGroup(j)
           beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
     ENDIF
     nbeta(l)=ni
     alam(l)=al
     nalam=l
     IF (l < mnl) CYCLE
     me=0
     DO j=1,ni
        g=activeGroup(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,activeGroupIndex)
  RETURN
END SUBROUTINE sparse_orig





