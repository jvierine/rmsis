! Separate C ABI adapter. NRL's source files are compiled without edits.
module rmsis21_bridge
  use iso_c_binding
  use msis_constants, only: rp, dmissing
  use msis_init, only: msisinit, masswgt
  use msis_calc, only: msiscalc
  implicit none
contains
  subroutine rmsis21_init(path, length, history) bind(C)
    character(c_char), intent(in) :: path(*)
    integer(c_int), value :: length, history
    character(len=128) :: filename
    real(4) :: switches(25)
    real(rp) :: temperature, densities(10)
    integer :: i
    filename = ''
    do i=1,length
      filename(i:i) = path(i)
    end do
    switches = 1.0
    if (history == 1) switches(9) = -1.0
    call msisinit(parmpath='', parmfile=trim(filename), switch_legacy=switches)
    ! Reset the original library's saved evaluation cache after reinitialization,
    ! following the pymsis adapter convention. No NRL source is modified.
    densities = 0.0_rp
    call msiscalc(0.0_rp,0.0_rp,-999.0_rp,-999.0_rp,-1.0_rp,0.0_rp,0.0_rp, &
                 [1.0_rp,1.0_rp,1.0_rp,1.0_rp,1.0_rp,1.0_rp,1.0_rp],temperature,densities)
  end subroutine
  subroutine rmsis21_calc(input, ap, output) bind(C)
    real(c_double), intent(in) :: input(7), ap(7)
    real(c_double), intent(out) :: output(12)
    real(rp) :: tn,tex,dn(10)
    integer :: i
    ! Original mass-density sum includes dn(1) with zero weight. Avoid any
    ! uninitialized NaN in the caller buffer without changing NRL source.
    dn = 0.0_rp
    call msiscalc(real(input(1),rp),real(input(2),rp),real(input(3),rp),real(input(4),rp), &
                  real(input(5),rp),real(input(6),rp),real(input(7),rp),real(ap,rp),tn,dn,tex)
    ! Compute the same weighted mass sum from species only. The original sum
    ! reads its as-yet-unassigned dn(1); 0*NaN can contaminate it. This is an
    ! adapter-level output safeguard, not an edit to the NRL library.
    dn(1) = dot_product(dn(2:),masswgt(2:))
    do i=1,10
      output(i)=real(dn(i),c_double)
      if (dn(i)==dmissing) output(i)=-1.0_c_double
    end do
    output(11)=real(tn,c_double)
    output(12)=real(tex,c_double)
  end subroutine
end module
