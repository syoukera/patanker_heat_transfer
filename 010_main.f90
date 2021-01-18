module variables
    implicit none

    integer, parameter :: nx = 12
    real(16), parameter :: xmax = 1.0
    real(16), parameter :: dx = xmax/(nx-2)
    real(16) :: X(nx)
    real(16) :: T(nx)

end module variables


program main
    use variables
    implicit none

    call grid()
    call initialize()
    call solvet()

    print *, 'aaa'
    print *, xmax
    print *, X

end program main

subroutine grid()
    use variables
    implicit none

    integer i

    do i = 1, nx
        X(i) = (i-1.5)*dx
    end do

end subroutine

subroutine initialize()
    implicit none

end subroutine

subroutine solvet()
    implicit none

    call calc_tdma()

end subroutine

subroutine calc_tdma()
    implicit none

end subroutine
