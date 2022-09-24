MODULE irlsgl_subfuns

   USE spmatmul
   IMPLICIT NONE

   CONTAINS
   ! assumes x is already weighted
   
   
   SUBROUTINE irls_update_step(bsg, startix, endix, b, lama, pfg, pfl1, lam1ma, x,&
      isDifZero, nobs, r, gamg, maxDif, nvars, lb, ub, eps)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: bsg, nobs, nvars, startix, endix
      DOUBLE PRECISION :: gamg, maxDif
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: g, oldb, dd, del
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: thll, thl, bl, bll 
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: r
      DOUBLE PRECISION :: snorm, tea
      DOUBLE PRECISION, INTENT(in) :: lama, pfg, lam1ma, lb, ub
      DOUBLE PRECISION, INTENT(in) :: pfl1(bsg)
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: x(nobs,nvars)
      INTEGER :: k, ll, lmax, isDifZero
      DOUBLE PRECISION :: tt, ell_old, ell_new, sc, eps
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: rr

      ALLOCATE(rr(nobs))
      ALLOCATE(g(bsg))
      ALLOCATE(oldb(bsg))
      ALLOCATE(dd(bsg))
      ALLOCATE(del(bsg))
      ALLOCATE(thll(bsg))
      ALLOCATE(thl(bsg))
      ALLOCATE(bl(bsg))
      ALLOCATE(bll(bsg))

      tt = 1.0D0
      lmax = 100
      isDifZero = 0
      oldb = b(startix:endix) ! passed in, for checking convergence
      thll = oldb 
      bll = oldb
      
      DO ll = 1, lmax
         bl = bll
         thl = thll
         dd = oldb - bl
         rr = r - MATMUL(x(:, startix:endix), dd)
         ell_old = DOT_PRODUCT(rr, rr)
         g = - MATMUL(rr, x(:, startix:endix)) / nobs
        
         DO ! backtracking loop
            bll = bl - g*tt
            CALL softthresh(bll, lama*tt*pfl1, bsg)
            snorm = SQRT(DOT_PRODUCT(bll, bll))
            tea = snorm - tt * lam1ma * pfg
            IF (tea > 0.0D0) THEN
               bll = g * tea / snorm
               DO k = startix, endix
                  bll(k) = MIN(MAX(lb, bll(k)), ub)
               ENDDO
            ELSE
               bll = 0.0D0
            ENDIF
            del = bll - bl
            dd = oldb - del
            rr = r - MATMUL(x(:, startix:endix), dd)
            ell_new = DOT_PRODUCT(rr, rr)
            IF (ell_new > ell_old + DOT_PRODUCT(g, del) + 0.5*DOT_PRODUCT(del,del) / tt) THEN
               tt = 0.8 * tt
            ELSE
               thll = bll
               EXIT
            ENDIF
         ENDDO
         ! Nesterov step
         sc = ll / (ll + 3)
         bll = thl + sc * (thll - thl)
         ! Convergence Check
         IF (ABS(ell_new - ell_old) < eps) EXIT
         IF (ALL(((gamg*(bll - bl)/(1+ABS(bll)))**2) < eps)) EXIT
      ENDDO ! finish group update
      
      dd = oldb - bll
      IF (ANY(dd .ne. 0.0D0)) THEN
         b(startix:endix) = bll
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd, dd))
         r = r - MATMUL(x(:,startix:endix), dd)
         isDifZero = 1
      ENDIF
      DEALLOCATE(rr, g, oldb, dd, del, thl, thll, bl, bll)
      RETURN
   END SUBROUTINE irls_update_step

   SUBROUTINE sp_irls_update_step(bsg, startix, endix, b, lama, pfg, pfl1, lam1ma, x,&
      xidx, xcptr, nnz, isDifZero, nobs, r, gamg, maxDif, nvars, lb, ub, eps)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: bsg, nobs, nvars, startix, endix, nnz
      DOUBLE PRECISION :: gamg, maxDif
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: g, oldb, dd, del
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: thll, thl, bl, bll 
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: r
      DOUBLE PRECISION :: snorm, tea
      DOUBLE PRECISION, INTENT(in) :: lama, pfg, lam1ma, lb, ub
      DOUBLE PRECISION, INTENT(in) :: pfl1(bsg)
      DOUBLE PRECISION, INTENT(in) :: x(nnz)
      INTEGER, INTENT(in) :: xidx(nnz)
      INTEGER, INTENT(in) :: xcptr(nvars + 1)
      INTEGER :: k, ll, lmax, isDifZero
      DOUBLE PRECISION :: tt, ell_old, ell_new, sc, eps
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: rr

      ALLOCATE(rr(nobs))
      ALLOCATE(g(bsg))
      ALLOCATE(oldb(bsg))
      ALLOCATE(dd(bsg))
      ALLOCATE(del(bsg))
      ALLOCATE(thll(bsg))
      ALLOCATE(thl(bsg))
      ALLOCATE(bl(bsg))
      ALLOCATE(bll(bsg))

      tt = 1.0D0
      lmax = 100
      isDifZero = 0
      oldb = b(startix:endix) ! passed in, for checking convergence
      thll = oldb 
      bll = oldb
      
      DO ll = 1, lmax
         bl = bll
         thl = thll
         dd = oldb - bl
         rr = r
         CALL ymspax(x, xidx, xcptr, nobs, nvars, nnz, dd, rr, startix, endix, bsg)
         ell_old = DOT_PRODUCT(rr, rr)
         CALL spatx(x, xidx, xcptr, nobs, nvars, nnz, rr, g, startix, endix)
         g = - g / nobs
        
         DO ! backtracking loop
            bll = bl - g*tt
            CALL softthresh(bll, lama*tt*pfl1, bsg)
            snorm = SQRT(DOT_PRODUCT(bll, bll))
            tea = snorm - tt * lam1ma * pfg
            IF (tea > 0.0D0) THEN
               bll = g * tea / snorm
               DO k = startix, endix
                  bll(k) = MIN(MAX(lb, bll(k)), ub)
               ENDDO
            ELSE
               bll = 0.0D0
            ENDIF
            del = bll - bl
            dd = oldb - del
            rr = r
            CALL ymspax(x, xidx, xcptr, nobs, nvars, nnz, dd, rr, startix, endix, bsg)
            ell_new = DOT_PRODUCT(rr, rr)
            IF (ell_new > ell_old + DOT_PRODUCT(g, del) + 0.5*DOT_PRODUCT(del,del) / tt) THEN
               tt = 0.8 * tt
            ELSE
               thll = bll
               EXIT
            ENDIF
         ENDDO
         ! Nesterov step
         sc = ll / (ll + 3)
         bll = thl + sc * (thll - thl)
         ! Convergence Check
         IF (ABS(ell_new - ell_old) < eps) EXIT
         IF (ALL(((gamg*(bll - bl)/(1+ABS(bll)))**2) < eps)) EXIT
      ENDDO ! finish group update
      
      dd = oldb - bll
      IF (ANY(dd .ne. 0.0D0)) THEN
         b(startix:endix) = bll
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd, dd))
         CALL ymspax(x, xidx, xcptr, nobs, nvars, nnz, dd, r, startix, endix, bsg)
         isDifZero = 1
      ENDIF
      DEALLOCATE(rr, g, oldb, dd, del, thl, thll, bl, bll)
      RETURN
   END SUBROUTINE sp_irls_update_step

END MODULE irlsgl_subfuns


