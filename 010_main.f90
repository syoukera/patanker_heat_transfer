module variables
    implicit none

    integer, parameter :: max_itr = 10
    integer, parameter :: ni = 12
    integer, parameter :: nim1 = 12-1
    integer, parameter :: ns = 53

    real(16) :: X(ni), Dxep(ni), Dxpw(ni), Sew(ni), Xu(ni)
    real(16) :: U(ni), P(ni), PP(ni), T(ni), W(ni, ns)

    ! Velocity
    integer, parameter :: nswpu = 2
    real(16), parameter :: urfu = 0.4
    real(16) :: resoru, Dxepu(ni), Dxpwu(ni), Sewu(ni)

    ! Pressure
    integer, parameter :: nswpp = 10
    integer, parameter :: ipref = 2
    integer, parameter :: jpref = 2
    real(16), parameter :: urfp = 0.4
    real(16) :: resorm, Du(ni)

    ! Temperature
    integer, parameter :: nswpt = 4
    real(16), parameter :: urft = 0.4
    real(16) :: resort

    ! TDMA coefficients
    real(16) :: Ap(ni), Ae(ni), Aw(ni)
    real(16) :: Su(ni), Sp(ni)

    real(16) :: Viscos(ni), Densit(ni), Sph(ni), Tcn(ni)

    ! ! Temporal initialization
    ! Viscos = 1.0e-3
    ! Densit = 998.2
    ! sph    = 4.1816d03
    ! tcn    = 0.594

    real(16), parameter :: beta   = 0.207d-3
    real(16), parameter :: gravit = 9.8
    
end module variables

subroutine initialize_variables()
    use variables
    implicit none

    U = 0.0
    P = 0.0
    PP = 0.0
    T = 300.0
    T(ni) =600.0
    W = 0.0

end subroutine initialize_variables
    

program main
    use variables
    implicit none

    integer :: i

    ! integer :: max_itr = 10

    ! initialize
    call grid
    call initialize_variables()

    ! Loop of conversion
    do i = 1, max_itr
        call calct
    end do

    print *, 'Solution diverges or needs more iterations'

end program main

subroutine calct()
    use variables
    implicit none

    integer :: i

    do i = 2, nim1
        Ae(i) = 0.0
        Aw(i) = 0.0
        Sp(i) = 0.0
        Su(i) = 0.0
    end do
    
    ! boundary conditions

    ! solve TDMA
    call TDMA(2, T)

end subroutine calct

subroutine grid()
    implicit none

    print *, 'bbb'

end subroutine grid


subroutine tdma(istart, PHI)
    use variables, only: Ap, Ae, Aw, Su, Sp, ni, nim1
    implicit none

    integer, intent(in) :: istart
    real(16), intent(inout) :: PHI(ni)
    integer :: istm1, i, ii
    real(16) :: A(ni), B(ni), C(ni), D(ni), term

    istm1 = istart - 1
    C(istm1) = PHI(istm1)
    ! commence S-N traverse
    do i = istart, nim1
        ! assemble TDMA coefficents
        C(i) = AE(i)*PHI(i+1) + AW(i)*PHI(i-1) + SU(i)
        D(i) = AP(i)
        ! calculate coefficients of recurrence formula
        term = 1./(D(i) - B(i)*A(i-1))
        A(i) = A(i)*term
        C(i) = (C(i) + B(i)*C(i-1))*term
    end do
    ! obtain new PHI's
    do ii = istart, nim1
        i = ni + istm1 - ii
        PHI(i) = A(i)* PHI(i+1) + C(i)
    end do

end subroutine tdma