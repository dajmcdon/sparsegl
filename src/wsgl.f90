SUBROUTINE wsgl (bn,bs,ix,iy,gam,nobs,nvars,x,r,pf,pfl1,pmax,&
   ulam,eps,maxit,intr,b0,beta,activeGroup,activeGroupIndex,ni,&
   npass,jerr,alsparse,lb,ub,sset,eset,b0old,betaold,al0,findlambda,l,me)

   USE sgl_subfuns
   IMPLICIT NONE
   ! - - - arg types - - -
   DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
   INTEGER:: isDifZero
   INTEGER:: intr
   INTEGER:: bn, findlambda
   INTEGER:: bs(bn)
   INTEGER:: ix(bn)
   INTEGER:: iy(bn)
   INTEGER:: nobs, nvars, pmax, npass, jerr, maxit
   INTEGER:: activeGroup(pmax)
   DOUBLE PRECISION :: eps, alsparse, max_gam, d, maxDif, al, snorm
   DOUBLE PRECISION, INTENT(in) :: x(nobs,nvars)
   ! DOUBLE PRECISION, INTENT(in) :: y(nobs)
   DOUBLE PRECISION, INTENT(in) :: pf(bn)
   DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
   DOUBLE PRECISION :: ulam, al0
   DOUBLE PRECISION :: gam(bn)
   DOUBLE PRECISION, INTENT(in) :: lb(bn), ub(bn)
   DOUBLE PRECISION :: b0, b0old
   DOUBLE PRECISION, INTENT(out) :: beta(nvars) ! our result
   DOUBLE PRECISION, INTENT(in) :: betaold(nvars) ! warm start
   DOUBLE PRECISION, INTENT(inout) :: r(nobs)

   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta ! local, not warm start
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u
   INTEGER :: activeGroupIndex(bn)
   INTEGER:: g, j, l, ni, me, startix, endix, i
   DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma
   DOUBLE PRECISION:: rpr(nvars + 1)

   ! - - - begin local declarations - - -
   DOUBLE PRECISION:: tlam, lama, lam1ma
   INTEGER:: violation
   INTEGER :: sset(bn)
   INTEGER :: eset(bn) ! this is for 4-step alg
   DOUBLE PRECISION:: ga(bn)
   DOUBLE PRECISION:: vl(nvars)
   ! - - - allocate variables - - -
   ALLOCATE(b(0:nvars))
   ALLOCATE(oldbeta(0:nvars))
   ! - - - checking pf - ! pf is the relative penalties for each group
   IF (MAXVAL(pf) <= 0.0D0) THEN
      jerr = 10000
      RETURN
   ENDIF

   ! - - - some initial setup - - -
   DO j = 1, nvars
      b(j) = betaold(j)
   ENDDO
   b(0) = b0old
   oldbeta = b
   npass = 0
   max_gam = MAXVAL(gam)
   t_for_s = 1 / gam
   i = 1

   vl = MATMUL(r, x) / nobs
   DO g = 1, bn ! For each group...
      ALLOCATE(u(bs(g)))
      u = vl(ix(g):iy(g))
      ga(g) = SQRT(DOT_PRODUCT(u, u))
      DEALLOCATE(u)
   ENDDO
   CALL rchkusr()

   ! PRINT *, alsparse
   al = ulam
   tlam = MAX((2.0 * al - al0), 0.0D0)
   lama = al * alsparse
   lam1ma = al * (1 - alsparse)
   ! This is the start of the algorithm, for a given lambda...
   CALL strong_rule(sset, ga, pf, tlam, alsparse) !uses s_set instead of e_set...
   ! --------- outer loop ---------------------------- !
   DO
      ! PRINT *, al
      CALL rchkusr()
      oldbeta(0) = b(0)
      IF (ni > 0) THEN
         DO j = 1, ni
            g = activeGroup(j)
            oldbeta(ix(g):iy(g)) = b(ix(g):iy(g))
         ENDDO
      ENDIF
      ! PRINT *, i
      ! PRINT *, b(0)
      ! PRINT *, b(1)
      ! PRINT *, oldbeta(0)
      ! PRINT *, oldbeta(1)

      ! --inner loop-------------------------------------
      DO
         CALL rchkusr()
         ! print *, "This is where we enter the inner loop"
         npass = npass + 1
         maxDif = 0.0D0
         isDifZero = 0 !Boolean to check if b-oldb nonzero. Unnec, in fn.
         DO g = 1, bn
            IF (eset(g) == 0) CYCLE
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
      ! IF (i > 10) RETURN
      violation = 0
      rpr = (max_gam * (b - oldbeta) / (1 + ABS(b)))**2
      ! print*, rpr
      IF (ANY(rpr >= eps)) violation = 1 !has beta moved globally
      IF (violation == 1) CYCLE
      CALL strong_kkt_check(eset, violation, bn, ix, iy, pf, pfl1, lam1ma,&
            bs, lama, ga, sset, x, r, nobs, nvars, vl) ! Step 3
      IF (violation == 1) CYCLE
      ! Need to compute vl/ga for the ones that aren't already updated, before kkt_check
      DO g = 1, bn
         CALL rchkusr()
         IF (sset(g) == 0) THEN
            startix = ix(g)
            endix = iy(g)
            ALLOCATE(s(bs(g)))
            s = MATMUL(r, x(:,startix:endix)) / nobs
            vl(startix:endix) = s
            CALL softthresh(s, lama * pfl1(startix:endix), bs(g))
            snorm = SQRT(DOT_PRODUCT(s,s))
            ga(g) = snorm
            DEALLOCATE(s)
         ENDIF
      ENDDO
      CALL kkt_check(eset, violation, bn, ix, iy, vl, pf, pfl1,lam1ma, bs, lama, ga, nvars)
      IF (violation == 1) CYCLE
      IF (findlambda == 1 .AND. ni == 0) THEN
         al0 = al
         al = al * 0.99
         tlam = MAX((2.0 * al - al0), 0.0D0)
         lama = al * alsparse
         lam1ma = al * (1 - alsparse)
         CYCLE
      ENDIF
      EXIT
   ENDDO ! Ends outer loop
   findlambda = 0
   ulam = al
   !---------- final update variable and save results------------
   CALL rchkusr()
   ! PRINT *, "Here is where the final update starts"
   IF (ni > pmax) THEN
      jerr = -10000 - l
   ENDIF
   IF (ni > 0) THEN
      DO j = 1, ni
         g = activeGroup(j)
         beta(ix(g):iy(g)) = b(ix(g):iy(g))
      ENDDO
   ENDIF
   b0 = b(0)
   me = 0
   DO j = 1, ni
      g = activeGroup(j)
      IF (ANY(beta(ix(g):iy(g)) .ne. 0.0D0)) me = me+1
   ENDDO
   ! IF (me > dfmax) EXIT
   DEALLOCATE(b, oldbeta)
   RETURN
END SUBROUTINE wsgl


SUBROUTINE spmat_wsgl (bn,bs,ix,iy,gam,nobs,nvars,x,xidx,xcptr,nnz,r,pf,pfl1,pmax,&
   ulam,eps,maxit,intr,b0,beta,activeGroup,activeGroupIndex,ni,&
   npass,jerr,alsparse,lb,ub,sset,eset,b0old,betaold,al0,findlambda,l,me)

   USE sgl_subfuns
   USE spmatmul
   IMPLICIT NONE
   ! - - - arg types - - -
   DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
   INTEGER:: isDifZero
   INTEGER:: intr
   INTEGER:: bn, findlambda
   INTEGER:: bs(bn)
   INTEGER:: ix(bn)
   INTEGER:: iy(bn)
   INTEGER:: nobs, nvars, pmax, npass, jerr, maxit, nnz
   INTEGER:: activeGroup(pmax)
   DOUBLE PRECISION :: eps, alsparse, max_gam, d, maxDif, al, alf, snorm
   DOUBLE PRECISION, INTENT(in) :: x(nnz)
   INTEGER, INTENT(in) :: xidx(nnz)
   INTEGER, INTENT(in) :: xcptr(nvars+1)
   DOUBLE PRECISION, INTENT(in) :: pf(bn)
   DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
   DOUBLE PRECISION :: ulam, al0
   DOUBLE PRECISION :: gam(bn)
   DOUBLE PRECISION, INTENT(in) :: lb(bn), ub(bn)
   DOUBLE PRECISION :: b0, b0old
   DOUBLE PRECISION, INTENT(out) :: beta(nvars) ! our result
   DOUBLE PRECISION, INTENT(in) :: betaold(nvars) ! warm start
   DOUBLE PRECISION, INTENT(inout) :: r(nobs)

   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta ! local, not warm start
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u
   INTEGER :: activeGroupIndex(bn)
   INTEGER:: g, j, l, ni, me, startix, endix, i
   DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma
   DOUBLE PRECISION:: rpr(nvars + 1)

   ! - - - begin local declarations - - -
   DOUBLE PRECISION:: tlam, lama, lam1ma
   INTEGER:: violation
   INTEGER :: sset(bn)
   INTEGER :: eset(bn) ! this is for 4-step alg
   DOUBLE PRECISION:: ga(bn)
   DOUBLE PRECISION:: vl(nvars)
   ! - - - allocate variables - - -
   ALLOCATE(b(0:nvars))
   ALLOCATE(oldbeta(0:nvars))
   ! - - - checking pf - ! pf is the relative penalties for each group
   IF (MAXVAL(pf) <= 0.0D0) THEN
      jerr = 10000
      RETURN
   ENDIF

   ! - - - some initial setup - - -
   DO j = 1, nvars
      b(j) = betaold(j)
   ENDDO
   b(0) = b0old
   oldbeta = b
   npass = 0
   alf = 0.0D0
   max_gam = MAXVAL(gam)
   t_for_s = 1 / gam
   i = 1

   vl = 0.0D0
   CALL spatx(x, xidx, xcptr, nobs, nvars, nnz, r, vl, 1, nvars)
   vl = vl / nobs
   !vl = MATMUL(r, x) / nobs
   DO g = 1, bn ! For each group...
      ALLOCATE(u(bs(g)))
      u = vl(ix(g):iy(g))
      ga(g) = SQRT(DOT_PRODUCT(u, u))
      DEALLOCATE(u)
   ENDDO
   CALL rchkusr()

   ! PRINT *, alsparse
   al = ulam
   tlam = MAX((2.0 * al - al0), 0.0D0)
   lama = al * alsparse
   lam1ma = al * (1 - alsparse)
   ! This is the start of the algorithm, for a given lambda...
   CALL strong_rule(sset, ga, pf, tlam, alsparse) !uses s_set instead of e_set...
   ! --------- outer loop ---------------------------- !
   DO
      ! PRINT *, al
      CALL rchkusr()
      oldbeta(0) = b(0)
      IF (ni > 0) THEN
         DO j = 1, ni
            g = activeGroup(j)
            oldbeta(ix(g):iy(g)) = b(ix(g):iy(g))
         ENDDO
      ENDIF

      ! --inner loop-------------------------------------
      DO
         CALL rchkusr()
         ! print *, "This is where we enter the inner loop"
         npass = npass + 1
         maxDif = 0.0D0
         isDifZero = 0 !Boolean to check if b-oldb nonzero. Unnec, in fn.
         DO g = 1, bn
            IF (eset(g) == 0) CYCLE
            startix = ix(g)
            endix = iy(g)
            CALL sp_update_step(bs(g), startix, endix, b, lama, t_for_s(g),&
                   pf(g), pfl1(startix:endix), lam1ma, x, xidx, xcptr, nnz, isDifZero, nobs,&
                   r, gam(g), maxDif, nvars, lb(g), ub(g))
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
      ! IF (i > 10) RETURN
      violation = 0
      rpr = (max_gam * (b - oldbeta) / (1 + ABS(b)))**2
      ! print*, rpr
      IF (ANY(rpr >= eps)) violation = 1 !has beta moved globally
      IF (violation == 1) CYCLE
      ! change this line
      CALL sp_strong_kkt_check(eset, violation, bn, ix, iy, pf, pfl1,&
            lam1ma, bs, lama, ga, sset, x, xidx, xcptr, nnz,&
            r, nobs, nvars, vl)
      IF (violation == 1) CYCLE
      ! Need to compute vl/ga for the ones that aren't already updated, before kkt_check
      DO g = 1, bn
         CALL rchkusr()
         IF (sset(g) == 0) THEN
            startix = ix(g)
            endix = iy(g)
            ALLOCATE(s(bs(g)))
            ! change this line
            s = 0.0D0
            CALL spatx(x, xidx, xcptr, nobs, nvars, nnz, r, s, startix, endix)
            vl(startix:endix) = s / nobs
            CALL softthresh(s, lama * pfl1(startix:endix), bs(g))
            snorm = SQRT(DOT_PRODUCT(s,s))
            ga(g) = snorm
            DEALLOCATE(s)
         ENDIF
      ENDDO
      CALL kkt_check(eset, violation, bn, ix, iy, vl, pf, pfl1,lam1ma, bs, lama, ga, nvars)
      IF (violation == 1) CYCLE
      IF (findlambda == 1 .AND. ni == 0) THEN
         al0 = al
         al = al * 0.99
         tlam = MAX((2.0 * al - al0), 0.0D0)
         lama = al * alsparse
         lam1ma = al * (1 - alsparse)
         CYCLE
      ENDIF
      EXIT
   ENDDO ! Ends outer loop
   findlambda = 0
   ulam = al
   !---------- final update variable and save results------------
   CALL rchkusr()
   ! PRINT *, "Here is where the final update starts"
   IF (ni > pmax) THEN
      jerr = -10000 - l
   ENDIF
   IF (ni > 0) THEN
      DO j = 1, ni
         g = activeGroup(j)
         beta(ix(g):iy(g)) = b(ix(g):iy(g))
      ENDDO
   ENDIF
   b0 = b(0)
   me = 0
   DO j = 1, ni
      g = activeGroup(j)
      IF (ANY(beta(ix(g):iy(g)) .ne. 0.0D0)) me = me+1
   ENDDO
   ! IF (me > dfmax) EXIT
   DEALLOCATE(b, oldbeta)
   RETURN
END SUBROUTINE spmat_wsgl
