module variables
    implicit none

    integer, parameter :: max_itr = 10
    integer, parameter :: ni = 12
    integer, parameter :: nim1 = ni-1
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
    
    real(16), parameter :: xmax = 1.0
    real(16), parameter :: ymax = 1.0

    real(16) mf_chem(ns)
    !
    data mf_chem /0.00E+00,0.00E+00,0.00E+00,2.20E-01,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00, &
                    0.00E+00,0.00E+00,0.00E+00,0.00E+00,5.51E-02,0.00E+00,0.00E+00,0.00E+00,0.00E+00, &
                    0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00, &
                    0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00, &
                    0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00, &
                    0.00E+00,0.00E+00,7.25E-01,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00/
    
end module variables

subroutine initialize_variables()
    use variables
    implicit none

    integer :: i, j

    U = 0.0
    P = 0.0
    PP = 0.0
    T = 300.0
    T(ni) =600.0
    W = 0.0

    do i = 1, ni
        do j = 1, ns
            W(i, j) = mf_chem(j)
        end do
    end do
    

end subroutine initialize_variables
    

program main
    use variables
    use chemkin_params, only: initialize_chemkin_workarray
    implicit none

    integer :: i

    ! integer :: max_itr = 10

    ! initialize
    call initialize_chemkin_workarray()
    call grid
    call initialize_variables()

    ! Loop of conversion
    do i = 1, max_itr
        call calct
    end do

    print *, T

    print *, 'Solution diverges or needs more iterations'

end program main

subroutine calct()
    use variables
    use chemkin_params, only: get_tranport_data

    implicit none

    integer :: i
    ! output transport data
    ! mixture diffusion coefficient [CM**2/S]
    real(16) :: D_mix(ns) 
    ! mixture thermal conductivity [ERG/CM*K*S]
    real(16) :: Lambda_mix
    !  mean specific heat at constant pressure [ergs/(gm*K)]
    real(16) :: c_p

    do i = 2, nim1
        call get_tranport_data(T(i), P(i), W(:, i), D_mix, Lambda_mix, c_p)

        Ae(i) = Lambda_mix/Dxep(i)
        Aw(i) = Lambda_mix/Dxpw(i)
        Ap(i) = Ae(i) + Aw(i)
        Sp(i) = 0.0
        Su(i) = 0.0
    end do
    
    ! boundary conditions

    ! solve TDMA
    call TDMA(2, T)

end subroutine calct

subroutine grid()
    use variables
    implicit none

    integer :: i, nig
    real(16) :: dx

    nig = ni - 2    ! = 11
    dx = xmax/float(nig) ! 1.0/10 = 0.1

    ! X = [-0.05, 0.05, 0.15, ... 1.05]
    X(1) = -0.5*dx
    X(2) = -X(1)
    do i = 3, ni
        X(i) = X(i-1) + dx
    end do

end subroutine grid


subroutine init()
    use variables
    implicit none
    
    integer :: i
   
    ! Assign dx value
    DXPW(1) = 0.0
    DXEP(ni) = 0.0
    do i = 1, nim1
        DXEP(i) = X(i+1) - X(i)
        DXPW(i+1) = DXEP(i)
    end do

    ! Averaged value of dx
    SEW(1) = 0.0
    SEW(ni) = 0.0
    do i = 2, nim1
        SEW(i) = 0.5*(DXEP(i) + DXPW(i))
    end do

    ! Posion of midpoint
    XU(1) = 0.0
    do i = 2, ni
        XU(i) = 0.5*(X(i) + X(i-1))
    end do
    ! print '(12(f6.2))', XU

    DXPWU(1) = 0.0
    DXPWU(2) = 0.0
    DXEPU(1) = 0.0
    DXEPU(ni) = 0.0
    do i = 1, nim1
        DXEPU(i) = XU(i+1) - XU(i)
        DXPWU(i+1) = DXEPU(i)
    end do
    ! print '(12(f6.2))', DXEPU
    ! print '(12(f6.2))', DXPWU

    SEWU(1) = 0.0
    do i = 2, ni
        SEWU(i) = (X(i) - X(i-1))
    end do
    ! print '(12(f6.2))', SEWU

    do i = 1, ni
        U(i) = 0.0
        P(i) = 0.0
        PP(i) = 0.0
        T(i) = 300
        SU(i) = 0.0
        SP(i) = 0.0
    end do

    do i = 1, ni
        DU(i) = 0.0
    end do   

end subroutine init


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