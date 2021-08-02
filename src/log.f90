MODULE spmatmul

   IMPLICIT NONE

   CONTAINS
   ! ----
   ! These functions perform Ax=y and Atx = y for A in CSC form
   ! A is given by (a, ridx, cptr)
   ! a is the value, size nnz
   ! ridx is the row index of each nz
   ! cptr is of length ncol(A)+1 where cptr[j] is the first entry of a in column j, last entry == nnz
   ! to slice into columns cj:ck, you need a[cptr[cj]:(cptrr[ck+1]-1)] and
   ! ridx[cptr[cj]:(cptr[ck+1]-1)] and cptr[cj:ck])
   ! We don't actually need Ax->y ever, we need y <- y-A[,cj:ck] x
   SUBROUTINE log_ymspax (a, ridx, cptr, n, p, nnz, x, y, cj, ck, lx)
      IMPLICIT NONE
      INTEGER n, p, nnz, cj, ck, lx
      DOUBLE PRECISION, INTENT(in) :: a(nnz)
      INTEGER, INTENT(in) :: ridx(nnz)
      INTEGER, INTENT(in) :: cptr(p+1)
      DOUBLE PRECISION, INTENT(in) :: x(lx)
      DOUBLE PRECISION, INTENT(inout) :: y(n)

      INTEGER i, j, k

      DO i = cj, ck
         k = cptr(i + 1) - 1
         DO j = cptr(i), k
            y(ridx(j)) = y(ridx(j)) - x(i - cj + 1) * a(j)
         ENDDO
      ENDDO
      RETURN
   END SUBROUTINE log_ymspax

   SUBROUTINE log_spatx (a, ridx, cptr, n, p, nnz, x, y, cj, ck)
      IMPLICIT NONE
      INTEGER n, p, nnz, cj, ck
      DOUBLE PRECISION, INTENT(in) :: a(nnz)
      INTEGER, INTENT(in) :: ridx(nnz)
      INTEGER, INTENT(in) :: cptr(p+1)
      DOUBLE PRECISION, INTENT(in) :: x(n)
      DOUBLE PRECISION, INTENT(inout) :: y(ck - cj + 1)

      INTEGER i, j, k
      y = 0.0D0

      DO i = cj, ck
         k = i - cj + 1
         DO j = cptr(i), (cptr(i + 1) - 1)
            y(k) = y(k) + x(ridx(j)) * a(j)
         ENDDO
      ENDDO
      RETURN
   END SUBROUTINE log_spatx

   SUBROUTINE log_softthresh(vec, thresh, n)
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
   END SUBROUTINE log_softthresh

END MODULE spmatmul

!---------------------------------------------

MODULE sgl_subfuns

   USE spmatmul
   IMPLICIT NONE

   CONTAINS

   SUBROUTINE log_strong_rule (is_in_E_set, ga, pf, tlam, alsparse)
      IMPLICIT NONE
      INTEGER :: g, k
      INTEGER, DIMENSION (:), INTENT(inout) :: is_in_E_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: ga
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: pf
      DOUBLE PRECISION, INTENT(in) :: tlam, alsparse
      DOUBLE PRECISION :: z
      k = SIZE(is_in_E_set)
      z = tlam * (1 - alsparse)
      
      DO g = 1, k
         IF (is_in_E_set(g) == 1) CYCLE
         IF (ga(g) > pf(g) * z) is_in_E_set(g) = 1
      ENDDO
      RETURN
   END SUBROUTINE log_strong_rule

   SUBROUTINE log_kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama, ga)
      IMPLICIT NONE
      INTEGER :: g, startix, endix
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) :: bs(bn)
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
         IF (is_in_E_set(g) == 1) CYCLE
         startix = ix(g)
         endix = iy(g)
         ALLOCATE(s(bs(g)))
         s = vl(startix:endix)
         CALL log_softthresh(s, lama, bs(g))
         snorm = SQRT(DOT_PRODUCT(s,s))
         ga(g) = snorm
         IF(ga(g) > pf(g) * lam1ma) THEN
            is_in_E_set(g) = 1
            violation = 1
         ENDIF
         DEALLOCATE(s)
      ENDDO
      RETURN
   END SUBROUTINE log_kkt_check


   SUBROUTINE log_update_step(bsg, startix, endix, b, lama, t_for_sg, pfg, lam1ma, x,y, &
         isDifZero, nobs, r, gamg, maxDif,nvars, lb, ub)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: bsg, nobs, nvars
      INTEGER, INTENT(in) :: startix, endix
      DOUBLE PRECISION :: gamg
      DOUBLE PRECISION, INTENT(inout) :: maxDif
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb, s, dd
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: b, r
      DOUBLE PRECISION :: snorm, tea
      DOUBLE PRECISION, INTENT(in) :: lama, t_for_sg, pfg, lam1ma, lb, ub
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: x(nobs,nvars)
      DOUBLE PRECISION :: y(nobs)
      INTEGER, INTENT(inout) :: isDifZero
      INTEGER :: k

      ALLOCATE(s(bsg))
      ALLOCATE(oldb(bsg))
      isDifZero = 0
      oldb = b(startix:endix)
      s = MATMUL(y/(1.0D0+exp(r)), x(:, startix:endix))/nobs
      s = s*t_for_sg + b(startix:endix)
      CALL log_softthresh(s, lama*t_for_sg, bsg)
      snorm = SQRT(DOT_PRODUCT(s,s))
      tea = snorm - t_for_sg * lam1ma * pfg
      IF (tea > 0.0D0) THEN
         b(startix:endix) = s * tea / snorm
         DO k = startix, endix
            b(k) = MIN(MAX(lb, b(k)), ub)
         ENDDO
      ELSE
         b(startix:endix) = 0.0D0
      ENDIF
      ALLOCATE(dd(bsg))
      dd = b(startix:endix) - oldb
      IF (ANY(dd .ne. 0.0D0)) THEN
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd,dd))
         r = r + MATMUL(x(:,startix:endix), dd)
         isDifZero = 1
      ENDIF
      DEALLOCATE(s, oldb, dd)
      RETURN
   END SUBROUTINE log_update_step


   SUBROUTINE log_strong_kkt_check(is_in_E_set,violation,bn,ix,iy,pf,lam1ma,bs,&
         lama,ga,is_in_S_set,x,y,r,nobs,nvars,vl)
      IMPLICIT NONE
      INTEGER, INTENT(in)::nobs
      INTEGER, INTENT(in)::nvars
      DOUBLE PRECISION,INTENT(in):: x(nobs, nvars)
      DOUBLE PRECISION,INTENT(in) :: y(nobs)
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
         IF (is_in_S_set(g) == 1) THEN
            startix = ix(g)
            endix = iy(g)
            ALLOCATE(s(bs(g)))
            s = MATMUL(y/(1.0D0+exp(r)), x(:, startix:endix))/nobs
            vl(startix:endix) = s
            CALL log_softthresh(s, lama, bs(g))
            snorm = SQRT(dot_PRODUCT(s,s))
            ga(g) = snorm
            DEALLOCATE(s)
            IF (is_in_E_set(g) == 1) CYCLE
            IF (ga(g) > pf(g) * lam1ma) THEN
               is_in_E_set(g) = 1
               violation = 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
   END SUBROUTINE log_strong_kkt_check

   SUBROUTINE log_sp_update_step(bsg, startix, endix, b, lama, t_for_sg, pfg, lam1ma, x,&
         xidx, xcptr, nnz, isDifZero, nobs, r, gamg, maxDif, nvars, lb, ub)

      IMPLICIT NONE
      INTEGER, INTENT(in) :: bsg, nobs, nvars, nnz
      INTEGER, INTENT(in) :: startix, endix
      DOUBLE PRECISION :: gamg
      DOUBLE PRECISION, INTENT(inout) :: maxDif
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb, s, dd
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: r
      DOUBLE PRECISION :: snorm, tea
      DOUBLE PRECISION, INTENT(in) :: lama, t_for_sg, pfg, lam1ma, lb, ub
      DOUBLE PRECISION, INTENT(in) :: x(nnz)
      INTEGER, INTENT(in) :: xidx(nnz)
      INTEGER, INTENT(in) :: xcptr(nvars + 1)
      INTEGER, INTENT(inout) :: isDifZero
      INTEGER :: k

      ALLOCATE(s(bsg))
      ALLOCATE(oldb(bsg))
      isDifZero = 0
      s = 0.0D0
      oldb = b(startix:endix)
      ! print *, oldb

      CALL log_spatx(x, xidx, xcptr, nobs, nvars, nnz, r, s, startix, endix)
      s = s * t_for_sg / nobs + b(startix:endix)
      CALL log_softthresh(s, lama * t_for_sg, bsg)
      snorm = SQRT(DOT_PRODUCT(s,s))
      tea = snorm - t_for_sg * lam1ma * pfg
      IF (tea > 0.0D0) THEN
         b(startix:endix) = s * tea / snorm
         DO k = startix, endix
            b(k) = MIN(MAX(b(k), lb), ub)
         ENDDO
      ELSE
         b(startix:endix) = 0.0D0
      ENDIF
      ALLOCATE(dd(bsg))
      dd = b(startix:endix) - oldb
      IF(ANY(dd .ne. 0.0D0)) THEN
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd,dd))
         CALL log_ymspax(x, xidx, xcptr, nobs, nvars, nnz, dd, r, startix, endix, bsg)
         isDifZero = 1
      ENDIF
      DEALLOCATE(s, oldb, dd)
      RETURN
   END SUBROUTINE log_sp_update_step


   SUBROUTINE log_sp_strong_kkt_check(is_in_E_set,violation,bn,ix,iy,pf,lam1ma,bs,&
         lama,ga,is_in_S_set,x,xidx,xcptr,nnz,r,nobs,nvars,vl)

      IMPLICIT NONE
      INTEGER, INTENT(in) :: nobs, nvars, nnz
      DOUBLE PRECISION, INTENT(in) :: x(nnz)
      INTEGER, INTENT(in) :: xidx(nnz)
      INTEGER, INTENT(in) :: xcptr(nvars + 1)
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
            s = 0.0D0
            CALL log_spatx(x, xidx, xcptr, nobs, nvars, nnz, r, s, startix, endix)
            vl(startix:endix) = s / nobs
            CALL log_softthresh(s, lama, bs(g))
            snorm = SQRT(dot_PRODUCT(s,s))
            ! print *, "kkt snorm = ", snorm
            ga(g) = snorm
            DEALLOCATE(s)
            IF(is_in_E_set(g) == 1) CYCLE
            IF(ga(g) > pf(g) * lam1ma) THEN
               is_in_E_set(g) = 1
               violation = 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
   END SUBROUTINE log_sp_strong_kkt_check

END MODULE sgl_subfuns


   
! --------------------------------------------------
SUBROUTINE log_sparse_four (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
   eps,maxit,nalam,beta,activeGroup,nbeta,alam,npass,jerr,alsparse,lb,ub)

   USE sgl_subfuns
   IMPLICIT NONE
   ! - - - arg types - - -
   DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
   INTEGER, PARAMETER :: mnlam = 6
   INTEGER:: isDifZero
   INTEGER:: mnl
   INTEGER:: bn
   INTEGER:: bs(bn)
   INTEGER:: ix(bn)
   INTEGER:: iy(bn)
   INTEGER:: nobs, nvars, dfmax, pmax, nlam, nalam, npass, jerr, maxit
   INTEGER:: activeGroup(pmax)
   INTEGER:: nbeta(nlam)
   DOUBLE PRECISION :: flmin, eps, alsparse, max_gam, maxDif, al, alf, snorm
   DOUBLE PRECISION :: x(nobs,nvars)
   DOUBLE PRECISION :: y(nobs)
   DOUBLE PRECISION :: pf(bn)
   DOUBLE PRECISION :: ulam(nlam)
   DOUBLE PRECISION :: gam(bn)
   DOUBLE PRECISION :: lb(bn), ub(bn)
   DOUBLE PRECISION :: beta(nvars,nlam)
   DOUBLE PRECISION :: alam(nlam)

   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s !need for log_sparse_four
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
   ALLOCATE(b(1:nvars))
   ALLOCATE(oldbeta(1:nvars))
   ALLOCATE(r(1:nobs))
   ALLOCATE(activeGroupIndex(1:bn))
   !    ALLOCATE(al_sparse)
   ! - - - checking pf - ! pf is the relative penalties for each group
   IF(MAXVAL(pf) <= 0.0D0) THEN
      jerr = 10000
      RETURN
   ENDIF
   pf = MAX(0.0D0, pf)
   ! - - - some initial setup - - -
   is_in_E_set = 0
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
   IF (flmin < 1.0D0) THEN ! THIS is the default...
      flmin = MAX(mfl, flmin) ! just sets a threshold above zero
      alf = flmin ** (1.0D0 / (nlam - 1.0D0))
   ENDIF
   ! PRINT *, alf
   vl = MATMUL(y/(1.0D0+exp(r)), x)/nobs
   al0 = 0.0D0
   DO g = 1, bn ! For each group...
      ALLOCATE(u(bs(g)))
      u = vl(ix(g):iy(g))
      ga(g) = SQRT(DOT_PRODUCT(u, u))
      DEALLOCATE(u)
   ENDDO
   DO vl_iter = 1, nvars
      al0 = MAX(al0, ABS(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
   ENDDO
   ! PRINT *, alsparse
   al = al0 ! this value ensures all betas are 0 
   l = 0
   tlam = 0.0D0
   DO WHILE (l < nlam) !! This is the start of the loop over all lambda values...
      ! print *, "l = ", l
      ! print *, "al = ", al
      al0 = al ! store old al value on subsequent loops, first set to al
      IF (flmin >= 1.0D0) THEN ! user supplied lambda value, break out of everything
         l = l+1
         al = ulam(l)
         ! print *, "This is at the flmin step of the while loop"
      ELSE
         IF (l > 1) THEN ! have some active groups
            al = al * alf
            tlam = MAX((2.0 * al - al0), 0.0) ! Here is the strong rule...
            l = l+1
            ! print *, "This is the l>1 step of while loop"
         ELSE IF (l == 0) THEN
            al = al * 0.99
            tlam = al
            ! Trying to find an active group
         ENDIF
      ENDIF
      lama = al * alsparse
      lam1ma = al * (1 - alsparse)
      ! This is the start of the algorithm, for a given lambda...
      CALL log_strong_rule (is_in_S_set, ga, pf, tlam, alsparse) !uses s_set instead of e_set...
      ! --------- outer loop ---------------------------- !
      DO
         ! print *, is_in_E_set
         IF (ni > 0) THEN
            ! print *, "ni > 0"
            DO j = 1, ni
               g = activeGroup(j)
               oldbeta(ix(g):iy(g)) = b(ix(g):iy(g))
            ENDDO
         ENDIF
         ! --inner loop-------------------------------------
         DO
            ! print *, "This is where we enter the inner loop"
            npass = npass + 1
            maxDif = 0.0D0
            isDifZero = 0 !Boolean to check if b-oldb nonzero. Unnec, in fn.
            DO g = 1, bn
               IF (is_in_E_set(g) == 0) CYCLE
               startix = ix(g)
               endix = iy(g)
               CALL log_update_step(bs(g), startix, endix, b, lama, t_for_s(g), pf(g), lam1ma, x,y,& 
                  isDifZero, nobs, r, gam(g), maxDif, nvars, lb(g), ub(g))
               IF (activeGroupIndex(g) == 0 .AND. isDifZero == 1) THEN
                  ni = ni + 1
                  IF (ni > pmax) EXIT
                  activeGroupIndex(g) = ni
                  activeGroup(ni) = g
               ENDIF
            ENDDO
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
         CALL log_strong_kkt_check(is_in_E_set, violation, bn, ix, iy, pf, lam1ma,&
               bs, lama, ga, is_in_S_set, x, y, r, nobs, nvars, vl) ! Step 3
         IF (violation == 1) CYCLE
         ! Need to compute vl/ga for the ones that aren't already updated, before log_kkt_check
         DO g = 1, bn
            IF (is_in_S_set(g) == 0) THEN
               startix = ix(g)
               endix = iy(g)
               ALLOCATE(s(bs(g)))
               s = MATMUL(y/(1.0D0+exp(r)), x(:, startix:endix))/nobs
               vl(startix:endix) = s
               CALL log_softthresh(s, lama, bs(g))
               snorm = SQRT(DOT_PRODUCT(s,s))
               ga(g) = snorm
               DEALLOCATE(s)
            ENDIF
         ENDDO
         CALL log_kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama, ga) ! Step 4
         IF (violation == 1) CYCLE
         EXIT
      ENDDO ! Ends outer loop
      !---------- final update variable and save results------------
      IF (l == 0) THEN
         IF (MAXVAL(is_in_E_set) == 0) THEN
            CYCLE ! don't save anything, we're still decrementing lambda
         ELSE
            l = 2
            alam(1) = al / MAX(alf, .99) ! store previous, larger value
         ENDIF
      ENDIF
      ! PRINT *, "Here is where the final update starts"
      IF(ni > pmax) THEN
         jerr = -10000 - l
         EXIT
      ENDIF
      IF (ni > 0) THEN
         DO j = 1, ni
            g = activeGroup(j)
            beta(ix(g):iy(g),l) = b(ix(g):iy(g))
         ENDDO
      ENDIF
      nbeta(l) = ni
      alam(l) = al
      nalam = l
      IF (l < mnl) CYCLE
      me = 0
      DO j = 1, ni
         g = activeGroup(j)
         IF (ANY(beta(ix(g):iy(g),l) .ne. 0.0D0)) me=me+1
      ENDDO
      IF (me > dfmax) EXIT
   ENDDO ! end lambda loop
   ! print *, is_in_E_set
   DEALLOCATE(b, oldbeta, r, activeGroupIndex)
   RETURN
END SUBROUTINE log_sparse_four


   ! --------------------------------------------------
SUBROUTINE log_spmat_four (bn,bs,ix,iy,gam,nobs,nvars,x,xidx,xcptr,nnz,y,pf,&
   dfmax,pmax,nlam,flmin,ulam,eps,maxit,intr,nalam,b0,beta,&
   activeGroup,nbeta,alam,npass,jerr,alsparse,lb,ub)
   ! --------------------------------------------------
   USE sgl_subfuns
   USE spmatmul
   IMPLICIT NONE
   ! - - - arg types - - -
   DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
   INTEGER, PARAMETER :: mnlam = 6
   INTEGER :: isDifZero, mnl, bn, nobs, nvars, nnz, dfmax, pmax, nlam, nalam
   INTEGER :: npass, jerr, maxit, intr
   INTEGER :: bs(bn)
   INTEGER :: ix(bn)
   INTEGER :: iy(bn)
   INTEGER :: activeGroup(pmax)
   INTEGER :: nbeta(nlam)
   DOUBLE PRECISION :: flmin, eps, max_gam, d, maxDif, al, alf, alsparse, snorm
   DOUBLE PRECISION, INTENT(in) :: x(nnz)
   INTEGER, INTENT(in) :: xidx(nnz)
   INTEGER, INTENT(in) :: xcptr(nvars+1)
   DOUBLE PRECISION, INTENT(in) :: y(nobs)
   DOUBLE PRECISION :: pf(bn)
   DOUBLE PRECISION :: ulam(nlam)
   DOUBLE PRECISION :: gam(bn)
   DOUBLE PRECISION, INTENT(in) :: lb(bn), ub(bn)
   DOUBLE PRECISION :: b0(nlam)
   DOUBLE PRECISION :: beta(nvars,nlam)
   DOUBLE PRECISION :: alam(nlam)
   ! - - - local declarations - - -
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s !need for log_sparse_four
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r ! Residual
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u 
   INTEGER, DIMENSION (:), ALLOCATABLE :: activeGroupIndex
   INTEGER :: g, j, l, ni, me, startix, endix, vl_iter
   ! - - - Aaron's declarations
   DOUBLE PRECISION :: t_for_s(bn) ! this is for now just 1/gamma
   ! - - - begin local declarations - - -
   DOUBLE PRECISION :: tlam, lama, lam1ma, al0
   INTEGER :: violation
   INTEGER :: is_in_E_set(bn)
   INTEGER :: is_in_S_set(bn) !this is for 4-step alg
   DOUBLE PRECISION :: ga(bn) ! What is this for??
   DOUBLE PRECISION :: vl(nvars) ! What is this for?
   ! - - - allocate variables - - -
   ALLOCATE(b(0:nvars))
   ALLOCATE(oldbeta(0:nvars))
   ALLOCATE(r(1:nobs))
   ALLOCATE(activeGroupIndex(1:bn))
   !    ALLOCATE(al_sparse)
   ! - - - checking pf - ! pf is the relative penalties for each group
   IF (MAXVAL(pf) <= 0.0D0) THEN
      jerr = 10000
      RETURN
   ENDIF
   pf = MAX(0.0D0, pf)
   ! - - - some initial setup - - -
   is_in_E_set = 0
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
   t_for_s = 1/gam 
   ! --------- lambda loop ----------------------------
   IF (flmin < 1.0D0) THEN ! THIS is the default...
      flmin = MAX(mfl, flmin) ! just sets a threshold above zero
      alf = flmin**(1.0D0 / (nlam - 1.0D0))
   ENDIF
   vl = 0.0D0
   CALL log_spatx(x, xidx, xcptr, nobs, nvars, nnz, r, vl, 1, nvars)
   vl = vl / nobs
   al0 = 0.0D0
   DO g = 1, bn ! For each group...
      ALLOCATE(u(bs(g)))
      u = vl(ix(g):iy(g))
      ga(g) = SQRT(DOT_PRODUCT(u,u))
      DEALLOCATE(u)
   ENDDO
   DO vl_iter = 1, nvars
      al0 = MAX(al0, ABS(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
   ENDDO
   ! PRINT *, alsparse
   al = al0 !  this value ensures all betas are 0 
   l = 0
   tlam = 0.0D0
   DO WHILE (l < nlam) !! This is the start of the loop over all lambda values...
      al0 = al ! store old al value on subsequent loops, first set to al
      IF (flmin >= 1.0D0) THEN ! user supplied lambda value, break out of everything
         l = l + 1
         al = ulam(l)
         ! print *, "This is at the flmin step of the while loop"
      ELSE
         IF (l > 1) THEN ! have some active groups
            ! print *, "l = ", l
            al = al * alf
            tlam = MAX((2.0 * al - al0), 0.0) ! Here is the strong rule...
            l = l + 1
            ! print *, "This is the l>1 step of while loop"
         ELSE IF (l == 0) THEN
            al = al * .99
            tlam = al
            ! Trying to find an active group
         ENDIF
      ENDIF
      lama = al * alsparse
      lam1ma = al * (1-alsparse)
      ! This is the start of the algorithm, for a given lambda...
      CALL log_strong_rule (is_in_S_set, ga, pf, tlam, alsparse) !uses s_set instead of e_set...
      ! --------- outer loop ---------------------------- !
      DO
         oldbeta(0) = b(0)
         IF (ni > 0) THEN
            DO j = 1, ni
               g = activeGroup(j)
               oldbeta(ix(g):iy(g)) = b(ix(g):iy(g))
            ENDDO
         ENDIF
         ! --inner loop-------------------------------------
         DO
            npass = npass + 1
            maxDif = 0.0D0
            isDifZero = 0 !Boolean to check if b-oldb nonzero. Unnec, in fn.
            DO g = 1, bn
               IF (is_in_E_set(g) == 0) CYCLE
               startix = ix(g)
               endix = iy(g)
               CALL log_sp_update_step(bs(g), startix, endix, b, lama, t_for_s(g),&
                           pf(g), lam1ma, x, xidx, xcptr, nnz, isDifZero, nobs,&
                           r, gam(g), maxDif, nvars, lb(g), ub(g))
               IF (activeGroupIndex(g) == 0 .AND. isDifZero == 1) THEN
                  ni = ni+1
                  IF (ni > pmax) EXIT
                  activeGroupIndex(g) = ni
                  activeGroup(ni) = g
               ENDIF
            ENDDO 
            IF(intr .ne. 0) THEN
               d = sum(r) / nobs
               IF(d .ne. 0.0D0) THEN
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
         violation = 0
         IF (ANY((max_gam * (b - oldbeta) / (1 + ABS(b)))**2 >= eps)) violation = 1 !has beta moved globally
         IF (violation == 1) CYCLE
         CALL log_sp_strong_kkt_check(is_in_E_set, violation, bn, ix, iy, pf,&
                  lam1ma, bs, lama, ga, is_in_S_set, x, xidx, xcptr, nnz,&
                  r,nobs,nvars, vl)
         IF (violation == 1) CYCLE
         ! Need to compute vl/ga for the ones that aren't already updated, before log_kkt_check
         DO g = 1, bn
            IF (is_in_S_set(g) == 0) THEN
               startix = ix(g)
               endix = iy(g)
               ALLOCATE(s(bs(g)))
               s = 0.0D0
               CALL log_spatx(x, xidx, xcptr, nobs, nvars, nnz, r, s, startix, endix)
               vl(startix:endix) = s / nobs
               CALL log_softthresh(s, lama, bs(g))
               snorm = SQRT(DOT_PRODUCT(s,s))
               ga(g) = snorm
               DEALLOCATE(s)
            ENDIF
         ENDDO
         CALL log_kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, lam1ma, bs, lama, ga) ! Step 4
         IF (violation == 1) CYCLE
         EXIT
      ENDDO ! Ends outer loop
      !---------- final update variable and save results------------
      IF (l == 0) THEN
         IF (MAXVAL(is_in_E_set) == 0) THEN
            CYCLE ! don't save anything, we're still decrementing lambda
         ELSE
            l=2
            alam(1) = al / MAX(alf, .99) ! store previous, larger value
         ENDIF
      ENDIF
      ! PRINT *, "Here is where the final update starts"
      IF (ni > pmax) THEN
         jerr = -10000 - l
         EXIT
      ENDIF
      IF (ni > 0) THEN
         DO j = 1, ni
            g = activeGroup(j)
            beta(ix(g):iy(g),l) = b(ix(g):iy(g))
         ENDDO
      ENDIF
      nbeta(l) = ni
      b0(l) = b(0)
      alam(l) = al
      nalam = l
      IF (l < mnl) CYCLE
      me = 0
      DO j = 1, ni
         g = activeGroup(j)
         IF (ANY(beta(ix(g):iy(g),l) .ne. 0.0D0)) me = me + 1
      ENDDO
      IF (me > dfmax) EXIT
   ENDDO ! end lambda loop
   DEALLOCATE(b,oldbeta,r,activeGroupIndex)
   RETURN
END SUBROUTINE log_spmat_four
