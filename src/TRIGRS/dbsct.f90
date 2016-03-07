subroutine dbsct(mmax, eps, a, b, rtb, dx, x0, imax)
! Refine a root of  f(x) = sin(a*x) + b*x*cos(a*x)  by interval halving.
!
! On entry, rtb is the end of a bracketing interval at which f <= 0 and
! dx is the signed width of that interval (pointing towards the other end).
! Each step halves dx and evaluates f at the trial point x0 = rtb + dx;
! when f(x0) <= 0 the trial point becomes the new rtb.  The search stops
! as soon as |dx| <= eps or f(x0) is exactly zero.
!
! On exit, imax is the number of halving steps taken, or 0 if the
! tolerance was not reached within mmax steps.  rtb, dx and x0 hold the
! state of the search, so a later call can continue where this one ended.
  implicit none
  integer, intent(in) :: mmax
  double precision, intent(in) :: eps, a, b
  double precision, intent(inout) :: rtb, dx, x0
  integer, intent(out) :: imax
  integer :: step
  double precision :: f

  imax = 0
  do step = 1, mmax
    dx = 0.5d0 * dx
    x0 = rtb + dx
    f = sin(a * x0) + b * x0 * cos(a * x0)
    if (f <= 0.d0) rtb = x0
    if (f == 0.d0 .or. abs(dx) <= eps) then
      imax = step
      return
    end if
  end do
end subroutine dbsct
