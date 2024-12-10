module constants
    implicit none
    
    real*8 :: k_bond(6,6), r_bond(6,6), k_angle(6,6,6)
    real*8 :: rho_angle(6,6,6), eps(6), sigma(6)
    real*8 :: eps_m(6,6), sigma_m(6,6)
    real*8 :: A_vdw(6,6), B_vdw(6,6)
    real*8 :: A_torsion, n_torsion, degress_to_rad, pi

contains

    subroutine initialize_parameters()
    integer :: i, j
    
    k_bond = 0.d0; r_bond = 0.d0; k_angle = 0.d0
    rho_angle = 0.d0; eps = 0.d0; sigma = 0.d0
    eps_m = 0.d0; sigma_m = 0.d0; A_vdw = 0.d0
    B_vdw = 0.d0

    pi = dacos(-1.d0)
    degress_to_rad = 2.d0*pi/360.d0
    ! The rows and columns of the matrices are going to be the atomic
    ! numbers of the atoms, that way we can implement the constants
    ! more easily in a matrix
    k_bond(6, 6) = 300.0d0
    k_bond(1, 6) = 350.d0
    k_bond(6, 1) = k_bond(1, 6)

    r_bond(6, 6) = 1.53d0
    r_bond(1, 6) = 1.11d0
    r_bond(6, 1) = r_bond(1, 6)
    
    ! For angles we need an extra dimensions as it depends on 3 atoms
    k_angle(1, 6, 1) = 35.0d0
    k_angle(1, 6, 6) = 35.0d0
    k_angle(6, 6, 1) = k_angle(1, 6, 6)
    k_angle(6, 6, 6) = 60.0d0

    rho_angle(1, 6, 1) = 109.5d0
    rho_angle(1, 6, 6) = 109.5d0
    rho_angle(6, 6, 1) = rho_angle(1, 6, 6)
    rho_angle(6, 6, 6) = 109.5d0

    A_torsion = 0.3d0
    n_torsion = 3.0d0
    
    eps(1)=0.03d0
    eps(6)=0.07d0
    sigma(1)=1.20d0
    sigma(6)=1.75d0
    
    ! Calculation of the Lennard-Jones potential
    do i=1, size(eps)
        do j=1, size(eps)
            A_vdw(i,j) = 4.d0 * (sqrt(eps(i)*eps(j))) * (2.d0*sqrt(sigma(i)*sigma(j)))**12.d0
            B_vdw(i,j) = 4.d0 * (sqrt(eps(i)*eps(j))) * (2.d0*sqrt(sigma(i)*sigma(j)))**6.d0
        end do
    end do
    end subroutine initialize_parameters

end module constants