SUBROUTINE wls_step (bn,bs,ix,iy,nobs,nvars,x,y,pf,pfl1,dfmax,pmax,&
  flmin,lam,eps,maxit,intr,b0,beta,activeGroup,alam,npass,jerr,&
  alsparse,lb,ub,eset,sset)
  USE irlsgl_subfuns
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

  INTEGER:: isDifZero
  INTEGER:: mnl, intr
  INTEGER:: bn
  INTEGER:: bs(bn)
  INTEGER:: ix(bn)
  INTEGER:: iy(bn)
  INTEGER:: nobs, nvars, dfmax, pmax, nlam, nalam, npass, jerr, maxit
  INTEGER:: activeGroup(pmax)
  INTEGER:: nbeta(nlam)
  DOUBLE PRECISION :: flmin, eps, alsparse, max_gam, d, maxDif, al, alf, snorm
  DOUBLE PRECISION, INTENT(in) :: x(nobs,nvars)
  DOUBLE PRECISION, INTENT(in) :: y(nobs)
  DOUBLE PRECISION, INTENT(in) :: pf(bn)
  DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
  DOUBLE PRECISION :: ulam(nlam)
  DOUBLE PRECISION :: gam(bn)
  DOUBLE PRECISION, INTENT(in) :: lb(bn), ub(bn)
  DOUBLE PRECISION :: b0(nlam)
  DOUBLE PRECISION :: beta(nvars,nlam)
  DOUBLE PRECISION :: alam(nlam), mse(nlam)

  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s !need for sparse_four
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r ! Residual
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u ! No longer using this for update step, but still need for other parts

INTEGER, DIMENSION (:), ALLOCATABLE :: activeGroupIndex
INTEGER:: g, j, l, ni, me, startix, endix, vl_iter
DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma

! - - - begin local declarations - - -
DOUBLE PRECISION:: tlam, lama, lam1ma, al0
INTEGER:: violation
INTEGER:: is_in_E_set(bn)
INTEGER:: is_in_S_set(bn) ! this is for 4-step alg
DOUBLE PRECISION:: ga(bn)
DOUBLE PRECISION:: vl(nvars)
! - - - allocate variables - - -
ALLOCATE(b(0:nvars))
ALLOCATE(oldbeta(0:nvars))
ALLOCATE(r(1:nobs))
ALLOCATE(activeGroupIndex(1:bn))
!    ALLOCATE(al_sparse)
! - - - checking pf - ! pf is the relative penalties for each group
IF(MAXVAL(pf) <= 0.0D0) THEN
 jerr = 10000
 RETURN
ENDIF

! - - - some initial setup - - -
is_in_E_set = 0
is_in_S_set = 0
al = 0.0D0
mnl = MIN(mnlam, nlam)
r = y
b = 0.0D0
oldbeta = 0.0D0
activeGroup = 0
activeGroupIndex = 0
npass = 0
ni = 0
alf = 0.0D0
max_gam = MAXVAL(gam)
t_for_s = 1 / gam
! --------- lambda loop ----------------------------
! PRINT *, alf
al0 = 0.0D0
DO g = 1, bn ! For each group...
 ALLOCATE(u(bs(g)))
 u = vl(ix(g):iy(g))
 ga(g) = SQRT(DOT_PRODUCT(u, u))
 DEALLOCATE(u)
ENDDO
CALL rchkusr()
! PRINT *, alsparse
l = 0
tlam = 0.0D0
 
 lama = al * alsparse
 lam1ma = al * (1 - alsparse)
    ! --inner loop-------------------------------------
    DO
       CALL rchkusr()
       ! print *, "This is where we enter the inner loop"
       npass = npass + 1
       maxDif = 0.0D0
       isDifZero = 0 !Boolean to check if b-oldb nonzero. Unnec, in fn.
       DO g = 1, bn
          IF (is_in_E_set(g) == 0) CYCLE
          startix = ix(g)
          endix = iy(g)
          CALL update_step(bs(g), startix, endix, b, lama, t_for_s(g),&
               pf(g), pfl1(startix:endix), lam1ma, x,&
               isDifZero, nobs, r, gam(g), maxDif, nvars, lb(g), ub(g))
          IF (activeGroupIndex(g) == 0 .AND. isDifZero == 1) THEN
             ni = ni + 1
             IF (ni > pmax) EXIT
             activeGroupIndex(g) = ni
             activeGroup(ni) = g
          ENDIF
       ENDDO
       IF (intr .ne. 0) THEN
        d = sum(r) / nobs
        IF (d .ne. 0.0D0) THEN
           b(0) = b(0) + d
           r = r - d
           maxDif = max(maxDif, d**2)
        ENDIF
       ENDIF
       IF (ni > pmax) EXIT
       IF (maxDif < eps) EXIT
       IF (npass > maxit) THEN
          jerr = -l
          RETURN
       ENDIF
    ENDDO ! End inner loop
    IF (ni > pmax) EXIT
    !--- final check ------------------------ ! This checks which violate KKT condition
    ! PRINT *, "Here is where the final check starts"
    ! print *, i ! Just to check how many final checks...
    ! i = i+1
    violation = 0
    IF (ANY((max_gam * (b - oldbeta) / (1 + ABS(b)))**2 >= eps)) violation = 1 !has beta moved globally
    IF (violation == 1) CYCLE
DEALLOCATE(b, oldbeta, r, activeGroupIndex)
RETURN
END SUBROUTINE sparse_four
