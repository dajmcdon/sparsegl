! Script for testing functions. Comment out all but the function you want to test


!-!----------------------------------------------------
!-subroutine strong_rule (jxx, ga, pf, tlam, alsparse)
!-        !-------------------------------------------
!-        implicit none
!-        integer :: g
!-        integer, dimension (:), intent(inout) :: jxx
!-        double precision, dimension (:), intent(in) :: ga
!-        double precision, dimension (:), intent(in) :: pf
!-        double precision, intent(in) :: tlam, alsparse
!-        !------------------------------------------
!-        do g = 1, size(jxx) 
!-                if(jxx(g) == 1) cycle
!-                if(ga(g) > pf(g)*tlam*(1-alsparse)) jxx(g) = 1
!-        enddo
!-end subroutine strong_rule
!-
!-
!-program test_strong_rule
!-        implicit none
!-        integer :: jxx(4) = (/1,0,0,1/) 
!-        double precision :: ga(4) = (/1,2,3,4/)
!-        double precision :: pf(4) = (/1,1,1,1/)
!-        double precision :: tlam = 2.1, alsparse = 0.5
!-        !-----------------------------
!-        interface
!-                subroutine strong_rule(jxx, ga, pf, tlam, alsparse)
!-                        integer, dimension (:), intent(inout) :: jxx
!-                        double precision, dimension (:), intent(in) :: ga
!-                        double precision, dimension (:), intent(in) :: pf
!-                        double precision, intent(in) :: tlam, alsparse
!-                end subroutine strong_rule
!-        end interface
!-        !------------------------------------------
!-        print *, jxx
!-        call strong_rule(jxx, ga, pf, tlam, alsparse)
!-        print *, jxx
!-end program


!-!--------------------------------------------
!-SUBROUTINE softthresh(vec, thresh)
!-        !------------------------------------
!-  INTEGER :: n, it
!-  double precision :: sg
!-  DOUBLE PRECISION, dimension (:), intent(inout) :: vec
!-  double precision, intent(in) :: thresh
!-  n = size(vec)
!-  DO it=1,n
!-    sg = vec(it)
!-    vec(it) = sign(max(abs(sg) - thresh, 0.0D0), sg)
!-  ENDDO
!-END SUBROUTINE softthresh
!-
!-
!-
!-
!-program test_softthresh
!-        implicit none
!-        integer :: n=10, i
!-        double precision, dimension(:), allocatable :: vec
!-        double precision :: thresh = 2.0
!-        !------------------------------
!-        interface 
!-                subroutine softthresh(vec, thresh)
!-                        double precision, dimension (:), intent(inout) :: vec
!-                        double precision, intent(in) :: thresh
!-                end subroutine softthresh
!-        end interface 
!-        !-----------------------------
!-        allocate(vec(1:n))
!-        do i = 1, n
!-                vec(i) = i
!-        end do
!-        print *, vec
!-        call softthresh(vec, thresh)
!-        print *, vec
!-end program test_softthresh
 

!-!----------------------------------------------
!-subroutine kkt_check(bn, soft_g, jx, startix, endix, jxx, ix, iy, bs, vl, snorm, lama, lam1ma, ga, pf, s)
!-        !--------------------------------------
!-        implicit none
!-        integer :: g, bn, soft_g, jx, startix, endix
!-        integer :: jxx(bn)
!-        integer :: ix(bn)
!-        integer :: iy(bn)
!-        integer :: bs(bn)
!-        double precision :: vl(nvars)
!-        double precision :: snorm, lama, lam1ma
!-        double precision, dimension (:), allocatable :: s
!-        double precision :: ga(bn)
!-        double precision :: pf(bn)
!-        !---------------------------------------
!-        do g = 1, bn
!-                if (jxx(g) == 1) cycle
!-                startix = ix(g)
!-                endix = iy(g)
!-                allocate (s(bs(g)))
!-                s = vl(startix:endix)
!-                do soft_g = 1,bs(g)
!-                        s(soft_g) = sign(max(abs(s(soft_g))-lama, 0.0D0), s(soft_g))
!-                enddo
!-                snorm = sqrt(dot_product(s,s))
!-                ga(g) = snorm
!-                deallocate(s)
!-                if(ga(g) > pf(g)*lam1ma) then
!-                        jxx(g) = 1
!-                        jx = 1
!-                endif
!-        enddo
!-end subroutine
!-


!---------------------------------------
subroutine update_step(bsg, startix, endix, dd, b, lama, t_for_sg, pfg, lam1ma, vl)
        implicit none
        integer, intent(in) :: bsg
        integer, intent(in) :: startix, endix
        double precision, dimension (:), allocatable :: oldb, s
        double precision, dimension (:), intent(inout) :: dd
        double precision, dimension (:), intent(inout) :: b
        double precision :: snorm, tea 
        double precision, intent(in) :: lama, t_for_sg, pfg, lam1ma
        double precision, dimension (:), intent(in) :: vl
        !------------------------------
        allocate(s(bsg))
        allocate(oldb(bsg))
        oldb = b(startix:endix)
        s = vl(startix:endix)
        s = s*t_for_sg + b(startix:endix)
        softthresh(s, lama*t_for_sg)
        snorm = sqrt(dot_product(s,s))
        tea = snorm - t_for_sg*lam1ma*pfg
        if(tea > 0.0D0) then
                b(startix:endix) = s*tea/snorm
        else
                b(startix:endix) = 0.0D0
        endif
        dd = b(startix:endix) - oldb
        deallocate(s, oldb)
end subroutine update_step
