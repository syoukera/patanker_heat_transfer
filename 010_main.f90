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

    real(16) :: Viscos(ni), Densit(ni), Sph(ni), Tcn(ni)

    ! ! Temporal initialization
    ! Viscos = 1.0e-3
    ! Densit = 998.2
    ! sph    = 4.1816d03
    ! tcn    = 0.594

    real(16), parameter :: beta   = 0.207d-3
    real(16), parameter :: gravit = 9.8
    

end module variables

program main
    use variables
    implicit none

    integer :: i

    ! integer :: max_itr = 10

    ! initialize
    call grid
    call init

    ! Loop of conversion
    do i = 1, max_itr
        call calct
    end do

    print *, 'Solution diverges or needs more iterations'

end program main

subroutine calct()
    implicit none

    print *, 'aaa'

end subroutine calct

subroutine grid()
    implicit none

    print *, 'bbb'

end subroutine grid

subroutine init()
    implicit none

    print *, 'ccc'

end subroutine init
