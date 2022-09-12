!
!   My_Function.f90
!   Growth
!
!   Created by Rong Fan on 10/10/21.
!   Copyright 2020 FanRong. All rights reserved.
!
module My_Function
    implicit none

contains

! ======================================================================
! Function approximation
! ======================================================================

subroutine spline(N,xa,ya,yp1,ypn,y2a)
    implicit none
integer, intent(in) :: N                        ! Length of data
real(8), intent(in) :: xa(N), ya(N)             ! Data point y = f(x)
real(8), intent(in) :: yp1, ypn                 ! Derivative at point x1 and xn
real(8), intent(out) :: y2a(N)                  ! Second derivatives
integer, parameter :: Nmax = 500                ! Largest anticipated value of n
integer :: i, k                                 ! Iterator
real(8) :: p, qn, sig, un, u(Nmax)

    y2a(1) = -0.5d0
    u(1) = (3d0/(xa(2)-xa(1)))*((ya(2)-ya(1))/(xa(2)-xa(1))-yp1)
    do i = 2,N-1                                ! Decomposition loop
        sig = (xa(i)-xa(i-1))/(xa(i+1)-xa(i-1))
        p = sig*y2a(i-1)+2d0
        y2a(i) = (sig-1d0)/p
        u(i) = (6d0*((ya(i+1)-ya(i))/(xa(i+1)-xa(i))-(ya(i)-ya(i-1))    &
              /(xa(i)-xa(i-1)))/(xa(i+1)-xa(i-1))-sig*u(i-1))/p
    enddo

    qn = 0.5d0
    un = (3d0/(xa(N)-xa(N-1)))*(ypn-(ya(N)-ya(N-1))/(xa(N)-xa(N-1)))
    y2a(N) = (un-qn*u(N-1))/(qn*y2a(N-1)+1d0)
    do k = N-1,1,-1                             ! Backsubstitution loop
        y2a(k) = y2a(k)*y2a(k+1)+u(k)
    enddo
end subroutine spline


subroutine splint(N,xa,ya,y2a,x,y)
integer, intent(in) :: N
real(8), intent(in) :: xa(N),ya(N),y2a(N)           ! Data point yi = f(xi); y2a from spline
real(8), intent(in) :: x                            ! Point to evaluate
real(8), intent(out) :: y                           ! Interpolated value
integer :: khi, klo
real(8) :: a, b, hx, hy

    klo = 1; khi = N
    call Interpolate(N,xa,x,klo,khi)
    hx = xa(khi)-xa(klo)
    hy = ya(khi)-ya(klo)
    if (hx==0d0 .or. hy==0d0) then
        y = ya(khi)
        return
    endif
    a = (xa(khi)-x)/hx                              ! Cubic spline polynomial is now evaluated.
    b = (x-xa(klo))/hx
    y = a*ya(klo)+b*ya(khi)+((a**3d0-a)*y2a(klo)+(b**3d0-b)*y2a(khi))*(hx**2d0)/6d0
end subroutine splint


subroutine splint_derivative(N,xa,ya,y2a,x,dy)
integer, intent(in) :: N
real(8), intent(in) :: xa(N),ya(N),y2a(N)           ! Data point yi = f(xi); y2a from spline
real(8), intent(in) :: x                            ! Point to evaluate
real(8), intent(out) :: dy                          ! Derivative
integer :: khi, klo
real(8) :: a, b, hx, hy

    klo = 1; khi = N
    call Interpolate(N,xa,x,klo,khi)
    hx = xa(khi)-xa(klo)
    hy = ya(khi)-ya(klo)
    if (hx == 0d0 .or. hy == 0d0) then
        dy = 0d0
        return
    endif
    a = (xa(khi)-x)/hx                              ! Cubic spline polynomial is now evaluated.
    b = (x-xa(klo))/hx
    dy = hy/hx-(3d0*a**2d0-1d0)/6d0*hx*y2a(klo)+(3d0*b**2d0-1d0)/6d0*hx*y2a(khi)
end subroutine splint_derivative


subroutine splie2(M,N,x1a,x2a,ya,y2a)
    implicit none
integer, intent(in) :: M, N                         ! Length of data
real(8), intent(in) :: x1a(M), x2a(N), ya(M,N)      ! Data point y = f(x1,x2)
real(8), intent(out) :: y2a(M,N)                    ! Second derivatives at (x1,x2)
integer :: j                                        ! Iterator
    do j = 1,M
        call spline(N,x2a,ya(j,:),0d0,0d0,y2a(j,:))
    enddo
end subroutine splie2


subroutine splin2(M,N,x1a,x2a,ya,y2a,x1,x2,y)
    implicit none
integer, intent(in) :: M, N                                 ! Length of data
real(8), intent(in) :: x1a(M),x2a(N),y2a(M,N),ya(M,N)       ! Data point y = f(x1,x2); y2a from splie2
real(8), intent(in) :: x1, x2                               ! Point to Interpolate
real(8), intent(out) :: y                                   ! Interpolated value
integer, parameter :: NN = 100                              ! Largest anticipated value of M and N
integer :: j
real(8) :: y2tmp(NN), yytmp(NN)
    do j = 1,M
        call splint(N,x2a,ya(j,:),y2a(j,:),x2,yytmp(j))
    enddo
    call spline(M,x1a,yytmp,0d0,0d0,y2tmp)
    call splint(M,x1a,yytmp,y2tmp,x1,y)
end subroutine splin2


subroutine Interpolate(N,xvec,x,p0,p1)
    implicit none
! Function grid
integer, intent(in) :: N
real(8), dimension(N), intent(in) :: xvec
! Point to Interpolate
real(8), intent(in) :: x
! Interpolate interval
integer, intent(inout) :: p0, p1
! Iterator
integer :: iter, p
    if (x .ge. xvec(p1)) then
        p0 = p1; return
    elseif (x .le. xvec(p0)) then
        p1 = p0; return
    endif
    
    do iter = 1,100
        if (p1-p0<2) exit
        p = (p0+p1)/2
        if (x==xvec(p)) then
            p0 = p; p1 = p
            return
        elseif (xvec(p)>x) then
            p1 = p
        else
            p0 = p
        endif
    enddo
end subroutine Interpolate



real(8) function Brent_Root(func,x1,x2,TOL,flag) result(x0)
    implicit none
! Input of the function
real(8), intent(in) :: x1, x2, TOL          ! Root search interval
integer, intent(out) :: flag
! Others
real(8) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
! Iterator
integer, parameter :: iter_max = 100
integer :: iter
! Tolerance
real(8), parameter :: eps = 3d-8

interface
    real(8) function func(x)
        implicit none
    real(8), intent(in) :: x
    end function func
end interface

    ! Initialize values
    flag = 0
    a = x1; b = x2
    fa = func(a); fb = func(b)
    ! Check if the root is bracketed
    if (bracket(fa,fb) == .false.) return
    c = b; fc = fb
    do iter = 1,iter_max
        if (bracket(fb,fc) == .false.) then
            c = a; fc = fa; d = b-a; e = d
        endif
        ! b must be the best guess
        if (abs(fc)<abs(fb)) then
            a = b; b = c; c = a
            fa = fb; fb = fc; fc = fa
        endif

        ! Check convergence
        tol1 = 2d0*eps*abs(b)+0.5d0*tol
        xm = 0.5d0*(c-b)
        if (abs(xm) .le. tol1 .or. fb==0d0) then
            x0 = b; flag = 1
            return
        endif

        if (abs(e) .ge. tol1 .and. abs(fa)>abs(fb)) then
            s = fb/fa
            if (a==c) then                          ! Secant
                p = 2d0*xm*s; q = 1d0-s
            else                                    ! Inverse quadradic
                q = fa/fc; r = fb/fc
                p = s*(2d0*xm*q*(q-r)-(b-a)*(r-1d0))
                q = (q-1d0)*(r-1d0)*(s-1d0)
            endif

            if (p>0d0) q = -q
            p = abs(p)
            if (2d0*p<min(3d0*xm*q-abs(tol1*q),abs(e*q))) then
                e = d; d = p/q                      ! Secant/Inverse quadradic
            else
                d = xm; e = d                       ! Bisection
            endif
        else
                d = xm; e = d                       ! Bisection
        endif
        a = b; fa = fb
        if (abs(d)>tol1) then
            b = b+d
        else
            b = b+sign(tol1,xm)
        endif
        fb = func(b)
    enddo
    x0 = b

contains

logical function bracket(a,b)
implicit none
real(8), intent(in) :: a, b
    bracket = .TRUE.
    if ((a>0 .and. b>0) .or. (a<0 .and. b<0)) bracket = .false.
end function bracket

end function Brent_Root



subroutine Brent_Min(func,ax,bx,cx,tol,xmin,flag)
    implicit none
! Input and output
real(8), intent(in) :: ax,bx,cx,tol
real(8), intent(out) :: xmin
integer, intent(out) :: flag
! Parameters
integer, parameter :: ITMAX=100
real(8), parameter :: CGOLD=0.3819660d0,ZEPS=1d-3*epsilon(ax)
! Others
integer :: iter
real(8) :: a,b,d,e,etemp,fx,fu,fv,fw,p,q,r,tol1,tol2,u,v,w,x,xm

interface
    real(8) function func(x)
        implicit none
    real(8), intent(in) :: x
    end function func
end interface

    flag = 0
    a = min(ax,cx)
    b = max(ax,cx)
    v = bx
    w = v
    x = v
    e = 0d0
    fx = func(x)        ! here uses the function
    fv = fx
    fw = fx
    do iter = 1,ITMAX
        xm = 0.5d0*(a+b)
        tol1 = tol*abs(x)+ZEPS
        tol2 = 2d0*tol1
        if (abs(x-xm) .le. (tol2-0.5d0*(b-a))) then
            xmin = x
            flag = 1
            RETURN
        endif
        if (abs(e)>tol1) then
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q-(x-w)*r
            q = 2d0*(q-r)
            if (q > 0.0d0) p=-p
            q = abs(q)
            etemp = e
            e = d
            if (abs(p)  .ge.  abs(0.5d0*q*etemp) .or. &
                p  .le.  q*(a-x) .or. p  .ge.  q*(b-x)) then
                e = merge(a-x,b-x, x  .ge.  xm )
                d = CGOLD*e
            else
                d = p/q
                u = x+d
                if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
            endif
        else
            e = merge(a-x,b-x, x  .ge.  xm )
            d = CGOLD*e
        endif
        u = merge(x+d,x+sign(tol1,d),abs(d)  .ge.  tol1)
        fu = func(u)                        ! here uses the function
        if (fu .le. fx) then
            if (u  .ge.  x) then
                a = x
            else
                b = x
            endif
            call shft(v,w,x,u)
            call shft(fv,fw,fx,fu)
        else
            if (u<x) then
                a = u
            else
                b = u
            endif
            if (fu .le. fw .or. w==x) then
                v = w
                fv = fw
                w = u
                fw = fu
            else if (fu .le. fv .or. v==x .or. v==w) then
                v = u
                fv = fu
            endif
        endif
    enddo

contains

!shift the variables
subroutine shft(a,b,c,d)
    implicit none
real(8), intent(out) :: a
real(8), intent(inout) :: b,c
real(8), intent(in) :: d
    a = b
    b = c
    c = d
end subroutine shft

end subroutine Brent_Min


subroutine Amoeba(p,ndim,func,TOL)
    implicit none
real(8), intent(in) :: TOL                                    ! Tolerance
integer, intent(in) :: ndim                                   ! Dimension of x
real(8), intent(inout) :: p(ndim+1,ndim)                      ! Multidimensional minimization
integer, parameter :: nmax = 20, itmax = 5000
real(8), parameter :: alpha = 1d0, beta = 0.5d0, gamma = 2d0
real(8) :: y(ndim+1), pr(nmax), prr(nmax), pbar(nmax), ypr, yprr
integer :: i, j, ihi, inhi, ilo, mpts, iter
real(8) :: ftol, rtol

interface
    real(8) function func(x)
        implicit none
    real(8), dimension(:), intent(in) :: x
    end function func
end interface

    ftol = TOL
    do iter = 1,ndim+1
        y(iter) = func(p(iter,:))
    enddo
    mpts = ndim+1

    do iter = 0,itmax
        ilo = 1
        if (y(1)>y(2)) then
            ihi = 1; inhi = 2
        else
            ihi = 2; inhi = 1
        endif

        do i=1,mpts
            if (y(i)<y(ilo)) ilo = i
            if (y(i)>y(ihi)) then
                inhi = ihi; ihi = i
            elseif (y(i)>y(inhi)) then
                if (i .ne. ihi) inhi = i
            endif
        enddo

        rtol = 2d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
        if (rtol<ftol) return

        do j = 1,ndim
            pbar(j)=0.
        enddo

        do i = 1,mpts
            if (i .ne. ihi) then
                do j = 1,ndim
                    pbar(j) = pbar(j)+p(i,j)
                enddo
            endif
        enddo

        do j = 1,ndim
            pbar(j) = pbar(j)/ndim
            pr(j) = (1d0+alpha)*pbar(j)-alpha*p(ihi,j)
        enddo

        ypr = func(pr)

        if (ypr .le. y(ilo)) then
            do j = 1,ndim
                prr(j) = gamma*pr(j)+(1.-gamma)*pbar(j)
            enddo
            yprr = func(prr)
            if (yprr<y(ilo)) then
                do j = 1,ndim
                    p(ihi,j)=prr(j)
                enddo
                y(ihi) = yprr
            else
                do j = 1,ndim
                    p(ihi,j) = pr(j)
                enddo
                y(ihi) = ypr
            endif
        else if (ypr .ge. y(inhi)) then
            if (ypr<y(ihi)) then
                do j = 1,ndim
                    p(ihi,j) = pr(j)
                enddo
                y(ihi) = ypr
            endif
            do j = 1,ndim
                prr(j) = beta*p(ihi,j)+(1d0-beta)*pbar(j)
            enddo
            yprr = func(prr)
            if (yprr<y(ihi)) then
                do j = 1,ndim
                    p(ihi,j) = prr(j)
                enddo
                y(ihi) = yprr
            else
                do i = 1,mpts
                    if (i .ne. ilo) then
                        do j = 1,ndim
                            pr(j) = 0.5*(p(i,j)+p(ilo,j))
                            p(i,j) = pr(j)
                        enddo
                        y(i) = func(pr)
                    endif
                enddo
            endif
        else
            do j = 1,ndim
                p(ihi,j) = pr(j)
            enddo
            y(ihi) = ypr
        endif
    enddo

end subroutine Amoeba


subroutine mnewt(n,func,x,ntrial,tolx,tolf)
    implicit none
integer, intent(in) :: n                    ! Dimension
integer, intent(in) :: ntrial               ! Number of trial
real(8), intent(in) :: tolf, tolx           ! Tolerance
real(8), intent(inout) :: x(n)              ! Update x
integer, parameter :: NP = 15               ! Maximum dimension
integer :: i, k, indx(NP)
real(8) :: d, errf, errx, fvec(NP), fjac(NP,NP), p(NP)

interface
    subroutine func(n,NP,x,fvec,fjac)
        implicit none
    integer, intent(in) :: n, NP
    real(8), intent(in) :: x(n)
    real(8), intent(out) :: fvec(NP), fjac(NP,NP)
    end subroutine func
end interface

    do k = 1,ntrial
        call func(n,NP,x,fvec,fjac)
        errf = 0d0
        do i = 1,n
            errf = errf+abs(fvec(i))
        enddo
        if (errf .le. tolf) exit
        do i = 1,n
            p(i) = -fvec(i)
        enddo
        call ludcmp(fjac,n,NP,indx,d)
        call lubksb(fjac,n,NP,indx,p)
        errx = 0d0
        do i = 1,n
            errx = errx+abs(p(i))
            x(i) = x(i)+p(i)
        enddo
        if (errx .le. tolx) exit
    enddo
    
end subroutine mnewt


subroutine ludcmp(a,n,NP,indx,d)    ! LU decomposition
    implicit none
real(8), intent(inout) :: a(NP,NP)
integer, intent(in) :: n, NP
integer, intent(out) :: indx(n)
real(8), intent(out) :: d
integer, parameter :: NMAX = 500
real(8), parameter :: TINY = 1d-20
integer :: i, imax, j, k
real(8) :: aamax, dum, sum, vv(NMAX)

    d = 1d0
    do i = 1,n
        aamax = 0d0
        do j = 1,n
            if (abs(a(i,j))>aamax) aamax=abs(a(i,j))
        enddo
        if (aamax==0d0) pause 'singular matrix in ludcmp'
        vv(i) = 1d0/aamax
    enddo
    do j = 1,n
        do i = 1,j-1
            sum = a(i,j)
            do k = 1,i-1
                sum = sum-a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
        enddo
        aamax = 0d0
        do i = j,n
            sum = a(i,j)
            do k = 1,j-1
                sum = sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum .ge. aamax) then
                imax = i
                aamax = dum
            endif
        enddo
        if (j .ne. imax) then
            do k = 1,n
                dum = a(imax,k)
                a(imax,k) = a(j,k)
                a(j,k) = dum
            enddo
            d = -d
            vv(imax) = vv(j)
        endif
        indx(j) = imax
        if (a(j,j)==0d0) a(j,j) = TINY
        if (j .ne. n) then
            dum = 1d0/a(j,j)
            do i = j+1,n
                a(i,j) = a(i,j)*dum
            enddo
        endif
    enddo
    
end subroutine ludcmp
    
    
subroutine lubksb(a,n,NP,indx,b)
! Solve AX = B
    implicit none
real(8), intent(in) :: a(NP,NP)
integer, intent(in) :: n, NP
integer, intent(in) :: indx(n)
real(8), intent(inout) :: b(n)
integer :: i, ii, j, ll
real(8) :: sum

    ii = 0
    do i = 1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if (ii .ne. 0) then
            do j = ii,i-1
                sum = sum-a(i,j)*b(j)
            enddo
        elseif (sum .ne. 0.) then
            ii = i
        endif
        b(i) = sum
    enddo
    do i = n,1,-1
        sum = b(i)
        do j = i+1,n
            sum = sum-a(i,j)*b(j)
        enddo
        b(i) = sum/a(i,i)
    enddo
    
end subroutine lubksb

end module My_Function
