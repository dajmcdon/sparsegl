MODULE mySubby

CONTAINS
  SUBROUTINE strong_rule (is_in_E_set, ga, pf, tlam, alsparse)
    IMPLICIT NONE
    INTEGER :: g
    INTEGER, DIMENSION (:), INTENT(inout) :: is_in_E_set
    DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: ga
    DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: pf
    DOUBLE PRECISION, INTENT(in) :: tlam, alsparse
    !------------------------------------------
    DO g = 1, SIZE(is_in_E_set)
       IF(is_in_E_set(g) == 1) CYCLE
       IF(ga(g) > pf(g)*tlam*(1-alsparse)) is_in_E_set(g) = 1
    ENDDO
    RETURN
  END SUBROUTINE strong_rule


  SUBROUTINE softthresh(vec, thresh, n)
    IMPLICIT NONE
    INTEGER :: n, it
    DOUBLE PRECISION :: sg
    DOUBLE PRECISION :: vec(n)
    DOUBLE PRECISION :: thresh
    DO it=1,n
       sg = vec(it)
       vec(it) = SIGN(MAX(ABS(sg) - thresh, 0.0D0), sg)
    ENDDO
    RETURN
  END SUBROUTINE softthresh



  SUBROUTINE kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama, ga)
    IMPLICIT NONE
    INTEGER :: g, startix, endix
    INTEGER, INTENT(in) :: bn
    INTEGER, INTENT(in) ::bs(bn)
    INTEGER, INTENT(in) :: ix(bn), iy(bn)
    INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
    DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
    DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: vl
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
    DOUBLE PRECISION :: snorm
    DOUBLE PRECISION, INTENT(in) :: pf(bn)
    INTEGER, INTENT(inout) :: violation
    DOUBLE PRECISION, INTENT(in) :: lam1ma, lama

    DO g = 1, bn
       IF(is_in_E_set(g) ==1) CYCLE
       startix = ix(g)
       endix = iy(g)
       ALLOCATE(s(bs(g)))
       s = vl(startix:endix)
       CALL softthresh(s, lama, bs(g))
       snorm = SQRT(dot_PRODUCT(s,s))
       ga(g) = snorm
       IF(ga(g) > pf(g)*lam1ma) THEN
          is_in_E_set(g) = 1
          violation = 1
       ENDIF
       DEALLOCATE(s)
    ENDDO
    RETURN
  END SUBROUTINE kkt_check


  SUBROUTINE update_step(bsg, startix, endix, b, lama, t_for_sg, pfg, lam1ma, x,&
       isDifZero, nobs, r, gamg, maxDif,nvars)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: bsg, nobs, nvars
    INTEGER, INTENT(in) :: startix, endix
    DOUBLE PRECISION :: gamg
    DOUBLE PRECISION, INTENT(inout) :: maxDif
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb, s, dd
    DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: b, r
    DOUBLE PRECISION :: snorm, tea
    DOUBLE PRECISION, INTENT(in) :: lama, t_for_sg, pfg, lam1ma
    DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: x(nobs,nvars)
    INTEGER, INTENT(inout) :: isDifZero

    ALLOCATE(s(bsg))
    ALLOCATE(oldb(bsg))
    isDifZero = 0
    oldb = b(startix:endix)
    s = MATMUL(r, x(:, startix:endix))/nobs
    s = s*t_for_sg + b(startix:endix)
    CALL softthresh(s, lama*t_for_sg, bsg)
    snorm = SQRT(dot_PRODUCT(s,s))
    tea = snorm - t_for_sg*lam1ma*pfg
    IF(tea > 0.0D0) THEN
       b(startix:endix) = s*tea/snorm
    ELSE
       b(startix:endix) = 0.0D0
    ENDIF
    ALLOCATE(dd(bsg))
    dd = b(startix:endix) - oldb
    IF(ANY(dd/=0.0D0)) THEN
       maxDif = MAX(maxDif,gamg**2*dot_PRODUCT(dd,dd))
       r=r-MATMUL(x(:,startix:endix),dd)
       isDifZero = 1
    ENDIF
    DEALLOCATE(s, oldb, dd)
    RETURN
  END SUBROUTINE update_step


  SUBROUTINE strong_kkt_check(is_in_E_set,violation,bn,ix,iy,pf,lam1ma,bs,&
       lama,ga,is_in_S_set,x,r,nobs,nvars,vl)
    IMPLICIT NONE
    INTEGER, INTENT(in)::nobs
    INTEGER, INTENT(in)::nvars
    DOUBLE PRECISION,INTENT(in):: x(nobs, nvars)
    DOUBLE PRECISION, INTENT(in):: r(nobs)
    DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: vl
    INTEGER :: g, startix, endix
    INTEGER, INTENT(in) :: bn
    INTEGER, INTENT(in) ::bs(bn)
    INTEGER, INTENT(in) :: ix(bn), iy(bn)
    INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
    INTEGER, DIMENSION(:), INTENT(in) :: is_in_S_set
    DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
    DOUBLE PRECISION :: snorm
    DOUBLE PRECISION, INTENT(in) :: pf(bn)
    INTEGER, INTENT(inout) :: violation
    DOUBLE PRECISION, INTENT(in) :: lam1ma, lama

    violation = 0
    DO g = 1, bn
       IF(is_in_S_set(g) == 1) THEN
          startix = ix(g)
          endix = iy(g)
          ALLOCATE(s(bs(g)))
          s = MATMUL(r,x(:,startix:endix))/nobs
          vl(startix:endix) = s
          CALL softthresh(s, lama, bs(g))
          snorm = SQRT(dot_PRODUCT(s,s))
          ga(g) = snorm
          DEALLOCATE(s)
          IF(is_in_E_set(g) == 1) CYCLE
          IF(ga(g) > pf(g)*lam1ma) THEN
             is_in_E_set(g) = 1
             violation = 1
          ENDIF
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE strong_kkt_check

END MODULE mySubby


!--------------------------------------------------------------


MODULE spmatmul

CONTAINS

  ! ----
  ! These functions perform Ax=y and Atx = y for A in CSC form
  ! A is given by (a, ridx, cptr)
  ! a is the value, size nnz
  ! ridx is the row index of each nz
  ! cptr is of length ncol(A)+1 where cptr[j] is the first entry of a in column j, last entry == nnz

  ! to slice into columns j:k, you need a[cptr[j]:(cptrr[k]-1)] and
  ! ridx[cptr[j]:(cptr[k]-1)] and cptr[j:k])
  SUBROUTINE spax (a, ridx, cptr, n, p, nnz, x, y)
    IMPLICIT NONE
    INTEGER n, p, nnz
    DOUBLE PRECISION, INTENT(in) :: a(nnz)
    DOUBLE PRECISION, INTENT(in) :: ridx(nnz)
    DOUBLE PRECISION, INTENT(in) :: cptr(p+1)
    DOUBLE PRECISION, INTENT(in) :: x(p)
    DOUBLE PRECISION, INTENT(inout) :: y(n)

    INTEGER i, j
    y = 0;

    DO i = 1, p
       DO j = cptr(i), (cptr(i+1)-1)
          y(ridx(j)) += x(i)*a(j)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE spax

  SUBROUTINE spatx (a, ridx, cptr, n, p, nnz, x, y)
    IMPLICIT NONE
    INTEGER n, p, nnz
    DOUBLE PRECISION, INTENT(in) :: a(nnz)
    DOUBLE PRECISION, INTENT(in) :: ridx(nnz)
    DOUBLE PRECISION, INTENT(in) :: cptr(p+1)
    DOUBLE PRECISION, INTENT(in) :: x(n)
    DOUBLE PRECISION, INTENT(inout) :: y(p)

    INTEGER i, j
    y = 0;

    DO i = 1, p
       DO j = cptr(i), (cptr(i+1)-1)
          y(i) += x(ridx(j))*a(j)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE spatx

END MODULE spmatmul

!-------------------------------------------------------------


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
  IF(MAXVAL(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=MAX(0.0D0,pf)
  ! - - - some initial setup - - -
  is_in_E_set = 0
  al = 0.0D0
  mnl = MIN (mnlam, nlam)
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
     flmin = MAX (mfl, flmin) ! just sets a threshold above zero
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  ! PRINT *, alf
  vl = MATMUL(r, x)/nobs ! Note r gets updated in middle and inner loop
  al0 = 0.0D0
  DO g = 1,bn ! For each group...
     ALLOCATE(u(bs(g)))
     u = vl(ix(g):iy(g))
     ga(g) = SQRT(dot_PRODUCT(u,u))
     DEALLOCATE(u)
  ENDDO
  DO vl_iter = 1,nvars
     al0 = MAX(al0, ABS(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
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
           tlam = MAX((2.0*al-al0), 0.0) ! Here is the strong rule...
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
     CALL strong_rule (is_in_E_set, ga, pf, tlam, alsparse) !implementing strong rule, updates is_in_E_set
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
              CALL update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g), lam1ma, x,&
                   isDifZero, nobs, r, gam(g), maxDif, nvars)
              IF(activeGroupIndex(g)==0 .AND. isDifZero == 1) THEN
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
                 CALL update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g),&
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
        max_gam = MAXVAL(gam)
        IF(ANY((max_gam*(b-oldbeta)/(1+ABS(b)))**2 >= eps)) violation = 1
        IF (violation == 1) CYCLE
        vl = MATMUL(r, x)/nobs
        CALL kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama, ga) ! kkt subroutine
        IF(violation == 1) CYCLE ! This goes all the way back to outer loop
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     IF(l==0) THEN
        IF(MAXVAL(is_in_E_set)==0) THEN
           CYCLE ! don't save anything, we're still decrementing lambda
        ELSE
           l=2
           alam(1) = al / MAX(alf,.99) ! store previous, larger value
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
        IF(ANY(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,activeGroupIndex)
  RETURN
END SUBROUTINE sparse_three

! --------------------------------------------------
SUBROUTINE sparse_three_alt (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
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
  ! integer :: i ! Just to check how many final checks etc.
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
  IF(MAXVAL(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=MAX(0.0D0,pf)
  ! - - - some initial setup - - -
  ! i = 0
  is_in_E_set = 0
  al = 0.0D0
  mnl = MIN (mnlam, nlam)
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
     flmin = MAX (mfl, flmin) ! just sets a threshold above zero
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  ! PRINT *, alf
  vl = MATMUL(r, x)/nobs ! Note r gets updated in middle and inner loop
  al0 = 0.0D0
  DO g = 1,bn ! For each group...
     ALLOCATE(u(bs(g)))
     u = vl(ix(g):iy(g))
     ga(g) = SQRT(dot_PRODUCT(u,u))
     DEALLOCATE(u)
  ENDDO
  DO vl_iter = 1,nvars
     al0 = MAX(al0, ABS(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
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
           tlam = MAX((2.0*al-al0), 0.0) ! Here is the strong rule...
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
     CALL strong_rule (is_in_E_set, ga, pf, tlam, alsparse) !implementing strong rule, updates is_in_E_set
     ! --------- outer loop ---------------------------- !
     DO
        ! print *, is_in_E_set
        IF(ni>0) THEN
           ! print *, "ni > 0"
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
              CALL update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g), lam1ma, x,&
                   isDifZero, nobs, r, gam(g), maxDif, nvars)
              IF(activeGroupIndex(g)==0 .AND. isDifZero == 1) THEN
                 ni=ni+1
                 IF(ni>pmax) EXIT
                 activeGroupIndex(g)=ni
                 activeGroup(ni)=g
              ENDIF
           ENDDO
           IF (ni > pmax) EXIT
           IF (maxDif < eps) EXIT
           IF(npass > maxit) THEN
              jerr=-l
              RETURN
           ENDIF
        ENDDO ! End middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------ ! This checks which violate KKT condition
        ! PRINT *, "Here is where the final check starts"
        ! print *, i
        ! i = i + 1
        violation = 0
        max_gam = MAXVAL(gam)
        IF(ANY((max_gam*(b-oldbeta)/(1+ABS(b)))**2 >= eps)) violation = 1 !has beta moved globally
        IF (violation == 1) CYCLE
        vl = MATMUL(r, x)/nobs
        CALL kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama, ga) ! kkt subroutine
        IF(violation == 1) CYCLE ! This goes all the way back to outer loop
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     IF(l==0) THEN
        IF(MAXVAL(is_in_E_set)==0) THEN
           CYCLE ! don't save anything, we're still decrementing lambda
        ELSE
           l=2
           alam(1) = al / MAX(alf,.99) ! store previous, larger value
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
        IF(ANY(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  ! print *, is_in_E_set
  DEALLOCATE(b,oldbeta,r,activeGroupIndex)
  RETURN
END SUBROUTINE sparse_three_alt

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
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s !need for sparse_four
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
  ! integer:: i !Just to check how many final checks etc.
  DOUBLE PRECISION::snorm
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
  INTEGER:: is_in_S_set(bn) !this is for 4-step alg
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
  IF(MAXVAL(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=MAX(0.0D0,pf)
  ! - - - some initial setup - - -
  ! i = 0
  is_in_E_set = 0
  al = 0.0D0
  mnl = MIN (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  activeGroup = 0
  kill_count = 0
  activeGroupIndex = 0
  npass = 0 ! This is a count, correct?
  ni = npass ! This controls so-called "outer loop"
  alf = 0.0D0
  max_gam = MAXVAL(gam) !should be outside of loop...
  t_for_s = 1/gam ! might need to use a loop if no vectorization.........
  ! --------- lambda loop ----------------------------
  IF(flmin < 1.0D0) THEN ! THIS is the default...
     flmin = MAX (mfl, flmin) ! just sets a threshold above zero
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  ! PRINT *, alf
  vl = MATMUL(r, x)/nobs ! Note r gets updated in middle and inner loop
  al0 = 0.0D0
  DO g = 1,bn ! For each group...
     ALLOCATE(u(bs(g)))
     u = vl(ix(g):iy(g))
     ga(g) = SQRT(dot_PRODUCT(u,u))
     DEALLOCATE(u)
  ENDDO
  DO vl_iter = 1,nvars
     al0 = MAX(al0, ABS(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
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
           tlam = MAX((2.0*al-al0), 0.0) ! Here is the strong rule...
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
     CALL strong_rule (is_in_S_set, ga, pf, tlam, alsparse) !uses s_set instead of e_set...
     ! --------- outer loop ---------------------------- !
     DO
        ! print *, is_in_E_set
        IF(ni>0) THEN
           ! print *, "ni > 0"
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
           isDifZero = 0 !Boolean to check if b-oldb nonzero. Unnec, in fn.
           DO g=1,bn
              IF(is_in_E_set(g)==0) CYCLE
              startix=ix(g)
              endix=iy(g)
              CALL update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g), lam1ma, x,&
                   isDifZero, nobs, r, gam(g), maxDif, nvars)
              IF(activeGroupIndex(g)==0 .AND. isDifZero == 1) THEN
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
        ENDDO ! End middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------ ! This checks which violate KKT condition
        ! PRINT *, "Here is where the final check starts"
        ! print *, i ! Just to check how many final checks...
        ! i = i+1
        violation = 0
        IF(ANY((max_gam*(b-oldbeta)/(1+ABS(b)))**2 >= eps)) violation = 1 !has beta moved globally
        IF (violation == 1) CYCLE
        CALL strong_kkt_check(is_in_E_set, violation, bn, ix, iy, pf, lam1ma,&
             bs, lama, ga, is_in_S_set, x,r, nobs,nvars, vl) ! Step 3
        IF(violation == 1) CYCLE
        ! Need to compute vl/ga for the ones that aren't already updated, before kkt_check
        DO g = 1, bn
           IF(is_in_S_set(g)==0) THEN
              startix = ix(g)
              endix = iy(g)
              ALLOCATE(s(bs(g)))
              s = MATMUL(r,x(:,startix:endix))/nobs
              vl(startix:endix) = s
              CALL softthresh(s, lama, bs(g))
              snorm = SQRT(dot_PRODUCT(s,s))
              ga(g) = snorm
              DEALLOCATE(s)
           ENDIF
        ENDDO
        !IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) violation = 1 !has beta moved globally
        !IF (violation == 1) then
        !        print *, "violation from beta moving globally"
        !        CYCLE
        !endif
        CALL kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama, ga) ! Step 4
        IF(violation == 1) CYCLE
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     IF(l==0) THEN
        IF(MAXVAL(is_in_E_set)==0) THEN
           CYCLE ! don't save anything, we're still decrementing lambda
        ELSE
           l=2
           alam(1) = al / MAX(alf,.99) ! store previous, larger value
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
        IF(ANY(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  ! print *, is_in_E_set
  DEALLOCATE(b,oldbeta,r,activeGroupIndex)
  RETURN
END SUBROUTINE sparse_four
