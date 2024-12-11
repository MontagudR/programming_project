module gradient
    use constants
    use functions
    use energy
    implicit none
    
contains
    
subroutine cart_bond_gradient(bond_list, coordinates, atomic_numbers, cart_bond_grad)
    integer, intent(in) :: bond_list(:,:), atomic_numbers(:)
    real*8, intent(in) :: coordinates(:)
    real*8, intent(out) :: cart_bond_grad(size(atomic_numbers)*3)
    integer :: i, atom1, atom2, at_num1, at_num2, ind1, ind2
    real*8 :: r(3), grad(3)
    ! This function calculates the stretching term of the cartesian gradient for all bonds.
    ! It returns an n_x vector with the cartesian gradient 

    cart_bond_grad = 0.d0
    do i=1, size(bond_list, 1)
        ! Store the atomic indexes and atomic numbers
        atom1 = bond_list(i, 1); atom2 = bond_list(i, 2)
        at_num1 = atomic_numbers(atom1); at_num2 = atomic_numbers(atom2)
        ! Because we have a vector of the form x1, y1, z1, x2 ..., get the starting index
        ! for each atom
        ind1 = atom1*3 - 2; ind2 = atom2*3 - 2
        ! Calculate the distance between atoms
        r(:) = coordinates(ind1:ind1+2) - coordinates(ind2:ind2+2)
        ! Calculate the bond gradient and store it in the vector
        cart_bond_grad(ind1:ind1+2) = cart_bond_grad(ind1:ind1+2) + get_bond_grad(at_num1, at_num2, r)
        cart_bond_grad(ind2:ind2+2) = cart_bond_grad(ind2:ind2+2) - get_bond_grad(at_num1, at_num2, r)
    end do

end subroutine cart_bond_gradient

function get_bond_grad(atom1, atom2, r)
    integer, intent(in) :: atom1, atom2
    real*8, intent(in) :: r(:)
    real*8, dimension(3) :: get_bond_grad
    ! This function contains the bond gradient formula

    get_bond_grad = 2*k_bond(atom1, atom2)*(norm2(r) - r_bond(atom1, atom2))*(r(:)/norm2(r))

end function get_bond_grad

subroutine calculate_cart_angle_gradient(angle_list, coordinates, atomic_numbers, cart_angle_grad)
    integer, intent(in) :: angle_list(:,:), atomic_numbers(:)
    real*8, intent(in) :: coordinates(:)
    real*8, intent(out) :: cart_angle_grad(size(atomic_numbers)*3)
    integer :: i, atom1, atom2, atom3, at_num1, at_num2, at_num3, ind1, ind2, ind3
    real*8 :: r_ba(3), r_bc(3), angle, p(3), g_a(3), g_c(3), g_b(3), tmp
    ! This function calculates the bending term of the cartesian gradient for all angles.
    ! It returns an n_x vector with the cartesian gradient 

    cart_angle_grad = 0.d0
    do i=1, size(angle_list, 1)
        ! Store the atomic indexes and atomic numbers
        atom1 = angle_list(i, 1); atom2 = angle_list(i, 2); atom3 = angle_list(i, 3)
        at_num1 = atomic_numbers(atom1); at_num2 = atomic_numbers(atom2); at_num3 = atomic_numbers(atom3)
        ind1 = atom1*3 - 2; ind2 = atom2*3 - 2; ind3 = atom3*3 - 2
        r_ba(:) = coordinates(ind1:ind1+2) - coordinates(ind2:ind2+2)
        r_bc(:) =  coordinates(ind3:ind3+2) - coordinates(ind2:ind2+2)
        angle = get_angle(r_ba, r_bc)
        p(:) = cross_prod(r_ba, r_bc)
        ! Calculate one part of the bending gradient
        tmp = 2*k_angle(at_num1, at_num2, at_num3)*(angle - degress_to_rad*rho_angle(at_num1, at_num2, at_num3))
        ! Store the gradient in the vector
        cart_angle_grad(ind1:ind1+2) = cart_angle_grad(ind1:ind1+2) + tmp*grad_angle_a(r_ba, p)
        cart_angle_grad(ind2:ind2+2) = cart_angle_grad(ind2:ind2+2) + tmp*grad_angle_b(r_ba, r_bc, p)
        cart_angle_grad(ind3:ind3+2) = cart_angle_grad(ind3:ind3+2) + tmp*grad_angle_c(r_bc, p)
    end do

end subroutine calculate_cart_angle_gradient

function grad_angle_a(r_ba, p)
    real*8, intent(in) :: r_ba(:), p(:)
    real*8 :: grad_angle_a(3)
    ! This function contains the bending gradient term for the atom a

    grad_angle_a = (cross_prod(r_ba, p))/((norm2(r_ba)**2)*norm2(p))

end function grad_angle_a

function grad_angle_b(r_ba, r_bc, p)
    real*8, intent(in) :: r_ba(:), p(:), r_bc(:)
    real*8 :: grad_angle_b(3)
    ! This function contains the bending gradient term for the atom b

    grad_angle_b = (cross_prod(-r_ba, p))/((norm2(r_ba)**2)*norm2(p)) + &
                    (cross_prod(r_bc, p))/((norm2(r_bc)**2)*norm2(p))

end function grad_angle_b

function grad_angle_c(r_bc, p)
    real*8, intent(in) :: r_bc(:), p(:)
    real*8 :: grad_angle_c(3)
    ! This function contains the bending gradient term for the atom c

    grad_angle_c = (cross_prod(-r_bc, p))/((norm2(r_bc)**2)*norm2(p))

end function grad_angle_c

subroutine calculate_cart_torsion_gradient(dihedral_list, coordinates, atomic_numbers, cart_dihedral_grad)
    integer, intent(in) :: dihedral_list(:,:), atomic_numbers(:)
    real*8, intent(in) :: coordinates(:)
    real*8, intent(out) :: cart_dihedral_grad(size(atomic_numbers)*3)
    integer :: i, atom1, atom2, atom3, atom4, ind1, ind2, ind3, ind4
    real*8 :: r_ab(3), r_bc(3), r_cd(3), dihedral, v(3), t(3), u(3), r_ac(3), r_bd(3), tmp1(3), tmp2(3), grad
    ! This function calculates the torsion term of the cartesian gradient for all dihedrals.
    ! It returns an n_x vector with the cartesian gradient 

    cart_dihedral_grad = 0.d0
    do i=1, size(dihedral_list, 1)
        ! Store the atomic indexes
        atom1 = dihedral_list(i, 1); atom2 = dihedral_list(i, 2)
        atom3 = dihedral_list(i, 3); atom4 = dihedral_list(i, 4)
        ind1 = atom1*3 - 2; ind2 = atom2*3 - 2; ind3 = atom3*3 - 2; ind4 = atom4*3 - 2
        r_ab(:) = coordinates(ind2:ind2+2) - coordinates(ind1:ind1+2)
        r_ac(:) = coordinates(ind3:ind3+2) - coordinates(ind1:ind1+2)
        r_bc(:) = coordinates(ind3:ind3+2) - coordinates(ind2:ind2+2)
        r_bd(:) = coordinates(ind4:ind4+2) - coordinates(ind2:ind2+2)
        r_cd(:) = coordinates(ind4:ind4+2) - coordinates(ind3:ind3+2)
        ! Calculate intermediate variables
        t = cross_prod(r_ab, r_bc); u = cross_prod(r_bc, r_cd); v = cross_prod(t, u)
        dihedral = get_dihedral(r_ab, r_bc, r_cd)
        ! Calculate one of the terms of the gradient
        grad = -n_torsion*A_torsion*dsin(dihedral*n_torsion)
        ! Store and calculate the torsional gradient for each atom of the dihedral
        cart_dihedral_grad(ind2:ind2+2) = cart_dihedral_grad(ind2:ind2+2) + grad_torsion_b(t, r_bc, r_ac, r_cd, u)*grad
        cart_dihedral_grad(ind3:ind3+2) = cart_dihedral_grad(ind3:ind3+2) + grad_torsion_c(t, r_bc, r_ab, r_bd, u)*grad
        cart_dihedral_grad(ind1:ind1+2) = cart_dihedral_grad(ind1:ind1+2) + grad_torsion_a(t, r_bc)*grad
        cart_dihedral_grad(ind4:ind4+2) = cart_dihedral_grad(ind4:ind4+2) + grad_torsion_d(u, r_bc)*grad
    end do
    
end subroutine calculate_cart_torsion_gradient

function grad_torsion_a(t, r_bc)
    real*8, intent(in) :: t(:), r_bc(:)
    real*8 :: grad_torsion_a(3)
    ! This function contains the torsion gradient term for the atom a

    grad_torsion_a = cross_prod( (cross_prod(t,r_bc)/((norm2(t)**2)*norm2(r_bc))), r_bc )

end function grad_torsion_a

function grad_torsion_b(t, r_bc, r_ac, r_cd, u)
    real*8, intent(in) :: t(:), r_bc(:), r_ac(:), r_cd(:), u(:)
    real*8 :: grad_torsion_b(3)
    ! This function contains the torsion gradient term for the atom b

    grad_torsion_b = cross_prod( r_ac, (cross_prod(t,r_bc)/((norm2(t)**2)*norm2(r_bc))) ) + &
                    cross_prod( (cross_prod(-u,r_bc)/((norm2(u)**2)*norm2(r_bc))), r_cd )

end function grad_torsion_b

function grad_torsion_c(t, r_bc, r_ab, r_bd, u)
    real*8, intent(in) :: t(:), r_bc(:), r_ab(:), r_bd(:), u(:)
    real*8 :: grad_torsion_c(3)
    ! This function contains the torsion gradient term for the atom c

    grad_torsion_c = cross_prod( (cross_prod(t,r_bc)/((norm2(t)**2)*norm2(r_bc))), r_ab ) + &
                    cross_prod( r_bd, (cross_prod(-u,r_bc)/((norm2(u)**2)*norm2(r_bc))) )

end function grad_torsion_c

function grad_torsion_d(u, r_bc)
    real*8, intent(in) :: u(:), r_bc(:)
    real*8 :: grad_torsion_d(3)
    ! This function contains the torsion gradient term for the atom d

    grad_torsion_d = cross_prod( (cross_prod(-u,r_bc)/((norm2(u)**2)*norm2(r_bc))), r_bc )

end function grad_torsion_d

subroutine calculate_cart_vdw_gradient(vdw_list, coordinates, atomic_numbers, cart_vdw_grad)
    integer, intent(in) :: vdw_list(:,:), atomic_numbers(:)
    real*8, intent(in) :: coordinates(:)
    real*8, intent(out) :: cart_vdw_grad(size(atomic_numbers)*3)
    integer :: i, atom1, atom2, at_num1, at_num2, ind1, ind2
    real*8 :: r_ij(3)
    ! This function calculates the VdW term of the cartesian gradient for the VdW pairs.
    ! It returns an n_x vector with the cartesian gradient

    cart_vdw_grad = 0.d0
    do i=1, size(vdw_list, 1)
        ! Store the atomic indexes and atomic numbers
        atom1 = vdw_list(i, 1); atom2 = vdw_list(i, 2)
        at_num1 = atomic_numbers(atom1); at_num2 = atomic_numbers(atom2)
        ind1 = atom1*3 - 2; ind2 = atom2*3 - 2
        r_ij(:) = coordinates(ind1:ind1+2) - coordinates(ind2:ind2+2)
        ! Store and calculate the VdW gradient
        cart_vdw_grad(ind1:ind1+2) = cart_vdw_grad(ind1:ind1+2) + grad_vdw(at_num1, at_num2, r_ij)
        cart_vdw_grad(ind2:ind2+2) = cart_vdw_grad(ind2:ind2+2) - grad_vdw(at_num1, at_num2, r_ij)
    end do

end subroutine calculate_cart_vdw_gradient

function grad_vdw(atom1, atom2, r)
    integer, intent(in) :: atom1, atom2
    real*8, intent(in) :: r(:)
    real*8 :: grad_vdw(3)
    ! This function calculates the VdW gradient for two atoms

    grad_vdw = r*( ((-12.0d0*A_vdw(atom1, atom2))/(norm2(r)**14.d0)) + ((6.0d0*B_vdw(atom1, atom2))/(norm2(r)**8.d0)) )

end function grad_vdw

subroutine calculate_cart_gradient(bond_list, angle_list, dihedral_list, vdw_list, coordinates, atomic_numbers, total_cart_grad)
    integer, intent(in) :: bond_list(:,:), angle_list(:,:), dihedral_list(:,:), vdw_list(:,:), atomic_numbers(:)
    real*8, intent(in) :: coordinates(:)
    real*8, intent(out) :: total_cart_grad(size(atomic_numbers)*3)
    real*8, dimension(size(atomic_numbers)*3) :: cart_strech_grad, cart_bending_grad, cart_torsion_grad, cart_vdw_grad
    ! This function calculate the total cartesian gradient summing all the contributions

    ! Call all the individual terms that compose the gradient
    call cart_bond_gradient(bond_list, coordinates, atomic_numbers, cart_strech_grad)
    call calculate_cart_angle_gradient(angle_list, coordinates, atomic_numbers, cart_bending_grad)
    call calculate_cart_torsion_gradient(dihedral_list, coordinates, atomic_numbers, cart_torsion_grad)
    call calculate_cart_vdw_gradient(vdw_list, coordinates, atomic_numbers, cart_vdw_grad)

    ! Sum the individual contributions
    total_cart_grad = cart_strech_grad + cart_bending_grad + cart_torsion_grad + cart_vdw_grad

end subroutine calculate_cart_gradient

subroutine cartesian_optimization(n_x, bond_list, angle_list, dihedral_list, vdw_list, coordinates, atomic_numbers, final_coords)
    integer, intent(in) :: bond_list(:,:), angle_list(:,:), dihedral_list(:,:), vdw_list(:,:), atomic_numbers(:), n_x
    real*8, intent(in) :: coordinates(:)
    real*8, intent(out) :: final_coords(n_x) 
    real*8 :: inv_hes(n_x, n_x), gradient(n_x)
    real*8 :: tmp_grad(n_x), p(n_x), tmp_p(n_x)
    real*8 :: total_energy=0., new_coords(n_x), alpha, new_total_energy, s(n_x)
    real*8 :: new_gradient(n_x), y(n_x), v(n_x), tmp, wolfe_rule
    real*8 :: new_inv_hes(n_x, n_x), grms, coords(n_x)
    integer :: i, j, k, num_atoms, count
    ! This function performs the optimization using the cartesian gradient

    num_atoms = size(atomic_numbers)
    coords = coordinates

    ! Generate the initial guess for the inverse hessian
    inv_hes = 0.d0
    do i=1, size(atomic_numbers)*3
        inv_hes(i, i) = 1.d0/300.d0
    end do

    ! Calculate the initial gradient
    call calculate_cart_gradient(bond_list, angle_list, dihedral_list, vdw_list, coords, atomic_numbers, gradient)

    print *, "Starting iterations"
    count = 0
    grms = sqrt(dot_product(gradient, gradient)/dble(n_x))

    ! Start the iterations. The convergence criteria is 1e-3 for the gradient RMS
    do while (grms > 0.001d0) !0.001d0

        p = matmul(inv_hes, -gradient)

        ! Start the linear search along the p direction
        call calculate_total_energy(bond_list, angle_list, dihedral_list, vdw_list, coords, atomic_numbers, total_energy)
        alpha = 0.80d0; tmp = 0.80d0
        new_total_energy = 1.d0 + total_energy
        wolfe_rule = (total_energy + 0.1d0*alpha*dot_product(p, gradient))

        ! Iterate until we satisfy the first wolfe rule
        do while (new_total_energy >= wolfe_rule)

            alpha = tmp
            new_coords = coords + alpha*p
            call calculate_total_energy(bond_list, angle_list, dihedral_list, vdw_list, new_coords, atomic_numbers, new_total_energy)
            wolfe_rule = (total_energy + 0.1d0*alpha*dot_product(p, gradient))
            print *, "Alpha", alpha, "New energy", new_total_energy

            tmp = alpha * 0.80d0

        end do
        ! Finised the linear search

        ! Update s
        s(:) = alpha*p 

        ! Calculate the new gradient using the new coordinates obtained in the 
        ! linear search
        call calculate_cart_gradient(bond_list, angle_list, dihedral_list, vdw_list, new_coords, atomic_numbers, new_gradient)

        ! Calculate y and v vectors
        y(:) = new_gradient - gradient
        v(:) = matmul(inv_hes, y)

        ! Update the hessian
        new_inv_hes = get_new_hessian(n_x, inv_hes, s, y, v)

        ! Use the new gradient to update the GRMS
        grms = sqrt(dot_product(new_gradient, new_gradient)/dble(n_x))
        print *, "Gradient RMS: ", grms
        
        ! Assign the new variables for the new iteration
        inv_hes = new_inv_hes; gradient = new_gradient; coords = new_coords
        count = count + 1

    end do

    print *, "Optimized geometry, number of iterations: ", count
    print *, "Minimized energy", new_total_energy

    final_coords = coords

    do i=1, size(coordinates)/3
        print "(3F12.6)", final_coords(3*i-2:3*i)
    end do

    print *, "---------------------------------------------------------------------------------"
    print *, "---------------------------------------------------------------------------------"

end subroutine cartesian_optimization

subroutine build_B_matrix(n_x, n_q, b_l, a_l, d_l, vdw_l, coords, at_n, B)
    integer, intent(in) :: b_l(:, :), a_l(: ,:), d_l(:, :), vdw_l(:, :), at_n(:), n_x, n_q
    real*8, intent(in) :: coords(:)
    integer :: i, j, count, at_1, at_2, at_3, at_4, ind1, ind2, ind3, ind4
    real*8 :: r(3), p(3), r_ba(3), r_bc(3), r_ab(3), r_cd(3), r_bd(3), r_ac(3), t(3), u(3)
    real*8, intent(out) :: B(n_q, n_x)
    ! This function builds the Wilson B matrix. This matrix has n_q * n_x dimensions
    ! In this function, the rows are each one of the internal coordinates and the 
    ! columns are the x1, y1, z1, x2 .... cartesian terms. First we loop over the 
    ! bonds, then the angles and finally the dihedrals

    count = 1

    B = 0.d0
    ! First, we include the gradient of the stretching terms
    do i=1, size(b_l, 1)
        at_1 = b_l(i, 1); at_2 = b_l(i, 2)
        ind1 = at_1*3 - 2; ind2 = at_2*3 - 2
        r = 0.
        r = coords(ind1:ind1+2) - coords(ind2:ind2+2)
        ! Calculate and store the stretching gradient in the b matrix
        B(count, ind1:ind1+2) =  r(:) / norm2(r)
        B(count, ind2:ind2+2) = -1*r(:) / norm2(r)
        count = count + 1
    end do

    ! Now, we include the gradient of the bending terms
    do i=1, size(a_l, 1)
        at_1 = a_l(i, 1); at_2 = a_l(i, 2); at_3 = a_l(i, 3);
        ind1 = at_1*3 - 2; ind2 = at_2*3 - 2; ind3 = at_3*3 - 2
        r_ba = coords(ind1:ind1+2) - coords(ind2:ind2+2); r_bc = coords(ind3:ind3+2) - coords(ind2:ind2+2)
        ! Calculate intermediate variable
        p = cross_prod(r_ba, r_bc)
        ! Calculate the strectching gradient for each atom in the bending and store it
        ! in the Wilson B matrix
        B(count, ind1:ind1+2) = grad_angle_a(r_ba, p)
        B(count, ind2:ind2+2) = grad_angle_b(r_ba, r_bc, p)
        B(count, ind3:ind3+2) = grad_angle_c(r_bc, p)
        count = count + 1
    end do

    ! Then, calculate the gradient of the torsin terms
    do i=1, size(d_l, 1)
        at_1 = d_l(i, 1); at_2 = d_l(i, 2); at_3 = d_l(i, 3); at_4 = d_l(i, 4)
        ind1 = at_1*3 - 2; ind2 = at_2*3 - 2; ind3 = at_3*3 - 2; ind4 = at_4*3 - 2
        r_ab = coords(ind2:ind2+2) - coords(ind1:ind1+2); r_bc = coords(ind3:ind3+2) - coords(ind2:ind2+2)
        r_bd = coords(ind4:ind4+2) - coords(ind2:ind2+2); r_ac = coords(ind3:ind3+2) - coords(ind1:ind1+2)
        r_cd = coords(ind4:ind4+2) - coords(ind3:ind3+2)
        ! Calculate intermediate variables 
        t = cross_prod(r_ab, r_bc); u = cross_prod(r_bc, r_cd)
        ! Calculate and store the gradient for each atom of the dihedral
        ! and store it in the Wilson B matrix
        B(count, ind1:ind1+2) = grad_torsion_a(t, r_bc)
        B(count, ind2:ind2+2) = grad_torsion_b(t, r_bc, r_ac, r_cd, u)
        B(count, ind3:ind3+2) = grad_torsion_c(t, r_bc, r_ab, r_bd, u)
        B(count, ind4:ind4+2) = grad_torsion_d(u, r_bc)
        count = count + 1
    end do

end subroutine build_B_matrix

subroutine calculate_G_inv(B, n_q, inv_G)
    real*8, intent(in) :: B(:,:)
    integer, intent(in) :: n_q
    real*8 , intent(out) :: inv_G(n_q, n_q)
    real*8 :: G(n_q,n_q), w(n_q), inv_w_m(n_q, n_q)
    real*8, allocatable :: work(:)
    integer :: i, info, lwork
    ! This function calculates the inverse of the G matrix. First, it calculates
    ! the G matrix, then we diagonalize it using dsyev, and we obtain the 
    ! eigenvalues and eigenvectors. With that we can calculate the 
    ! inverse of G

    ! Calculate the matrix G
    G = matmul(B, transpose(B))

    ! We diagonalize the matrix. The first call to dsyev is just to 
    ! get the optimal size for the work array
    allocate(work(1))
    lwork = -1
    call dsyev('V', 'U', n_q, G, n_q, w, work, lwork, info)   
    lwork = int(work(1))
    deallocate(work)

    ! Assignate the optimal size for the work array
    allocate(work(lwork))
    ! This diagonalizes the G matrix.
    call dsyev('V', 'U', n_q, G, n_q, w, work, lwork, info)
    ! On exit, w is the vector of the eigenvalues and G are the 
    ! eigenvectors

    ! Now we create a matrix with the eigenvalues vectors on the diagonal
    inv_w_m = 0.d0
    do i=1, size(w)
        ! Filter the small eigenvalues, so we do not have big values
        ! when doing the inverse
        if (w(i) < 1e-5) then
            inv_w_m(i, i) = 0.d0
        else
            inv_w_m(i, i) = 1.d0 / w(i)
        end if
    end do

    ! With the eigenvectors and the matrix of eigenvalues, we calculate
    ! the inverse of matrix G
    inv_G = matmul(matmul(G, inv_w_m), transpose(G))

end subroutine calculate_G_inv

function initial_inv_hess_internal(b_l, a_l, d_l, n_q) result(inv_hess)
    integer, intent(in) :: b_l(:, :), a_l(:, :), d_l(:, :), n_q
    real*8 :: inv_hess(n_q, n_q)
    integer :: i, count
    ! This function generates the initial inverse hessian for the 
    ! optimization in internal coordinates. We have to generate 
    ! different values for the strectching, bending and torsion
    ! terms

    inv_hess = 0.d0
    count = 1
    ! Values for the stretching terms
    do i=1, size(b_l, 1)
        inv_hess(count, count) = 1.d0 / 600.d0
        count = count + 1
    end do

    ! Values for the bending terms
    do i=1, size(a_l, 1)
        inv_hess(count, count) = 1.d0 / 150.d0
        count = count + 1
    end do

    ! Values for the torsion terms
    do i=1, size(d_l, 1)
        inv_hess(count, count) = 1.d0 / 80.d0
        count = count + 1
    end do

end function initial_inv_hess_internal

subroutine find_opt_cart(b_l, a_l, d_l, vdw_l, at_n, n_q, n_x, q, new_q, s_q, B, inv_G, x)
    real*8, intent(inout) :: x(:), new_q(:)
    real*8, intent(in) :: B(:, :), inv_G(:,:), q(:), s_q(:)
    integer, intent(in) :: n_q, n_x, at_n(:), b_l(:, :), a_l(:, :), d_l(:, :), vdw_l(:, :)
    real*8 :: new_x(n_x), fix_q(n_q), tmp, new_sq(n_q), tmp_sq(n_q), diff
    integer :: count
    ! This function performs the search for the optimal updated cartesian coordinates for 
    ! the optimization in internal coordinates. 

    count = 0
    ! Calculate q_k+1   
    fix_q = s_q + q
    diff = 1
    new_sq = s_q
    
    ! Start the iterations
    do while ( diff > 1.d-5 )

        ! Calculate the nex cartesian coordiantes
        new_x = x + matmul(matmul(transpose(B), inv_G), new_sq)
        ! Calculate the new energy and the new internal coordinates
        call calculate_interanls_and_total_energy(b_l, a_l, d_l, vdw_l, new_x, at_n, n_q, tmp, new_q)
        new_sq = fix_q - new_q
        ! Check that the change in the torsion term is not 
        ! above pi 
        new_sq = scan_torsions(size(b_l, 1) + size(a_l, 1) + 1, new_sq, n_q)

        diff = maxval(abs(new_x - x))
        print *, "Maximum change: ", diff
        x = new_x
        count = count + 1

    end do

    print *, " "
    print *, "Cartersian fiited after: ", count
    print *, " "

end subroutine find_opt_cart

subroutine calculate_internal_gradient(n_x, n_q, b_l, a_l, d_l, vdw_l, i_coords, at_n)
    integer, intent(in) :: b_l(:, :), a_l(: ,:), d_l(:, :), vdw_l(:, :), at_n(:), n_q, n_x
    real*8, intent(in) :: i_coords(:)
    real*8 :: B(n_q, n_x), inv_G(n_q, n_q), inv_hess(n_q, n_q), p_q(n_q), s_q(n_q)
    real*8  :: g_q(n_q), g_x(n_x), tmp_gx(n_x/3, 3), q(n_q), new_q(n_q), new_g_x(n_x), new_g_q(n_q)
    real*8  :: new_B(n_q, n_x), new_inv_G(n_q, n_q), y_q(n_q), v_q(n_q), new_inv_hess(n_q, n_q)
    real*8 :: lamda_max, l_q, tot_e, v_coords(size(i_coords)), coords(n_x)
    real*8 :: new_tot_e, grms
    integer :: i, j, num_at, count
    ! This function performs the optimization in internal coordinates.

    coords = i_coords
    lamda_max = 0.02d0
    count = 0

    ! Construct the B matrix
    call build_B_matrix(n_x, n_q, b_l, a_l, d_l, vdw_l, coords, at_n, B)
    ! Calculate the number of atoms
    num_at = size(at_n)
    
    call calculate_G_inv(B, n_q, inv_G)

    call calculate_cart_gradient(b_l, a_l, d_l, vdw_l, coords, at_n, g_x)
    g_q = matmul(matmul(inv_G, B), g_x)
    grms = sqrt(dot_product(g_x, g_x)/n_x)

    ! Calculate the initial hessian
    inv_hess = initial_inv_hess_internal(b_l, a_l, d_l, n_q)
    call calculate_interanls_and_total_energy(b_l, a_l, d_l, vdw_l, coords, at_n, n_q, tot_e, q)

    ! Start the iterations. The convergence criteria is that the gradient RMS
    ! is below 1e-3.
    do while( grms > 1d-3 ) ! 1e-3

        ! Calculate p
        p_q = -1 * matmul(inv_hess, g_q)
        ! Check that the step is not too big
        l_q = sqrt( (sum(p_q**2)) / (n_q))
        
        ! Rescale the step in case is too big
        if (lamda_max < l_q) then
            p_q = p_q * (lamda_max / l_q)
        end if

        new_q = q + p_q
        s_q = new_q - q
        s_q = scan_torsions(size(b_l, 1) + size(a_l, 1) + 1, s_q, n_q)

        call find_opt_cart(b_l, a_l, d_l, vdw_l, at_n, n_q, n_x, q, new_q, s_q, B, inv_G, coords)
        ! CAREFUL, now v_coords and s_q are the updated because of the search for the optimal
        ! coordinates
        s_q = new_q - q
        s_q = scan_torsions(size(b_l, 1) + size(a_l, 1) + 1, s_q, n_q)

        ! Get the upgraded internal gradient
        call build_B_matrix(n_x, n_q, b_l, a_l, d_l, vdw_l, coords, at_n, new_B)
        call calculate_G_inv(new_B, n_q, new_inv_G)
        call calculate_total_energy(b_l, a_l, d_l, vdw_l, coords, at_n, new_tot_e) ! WARNING NO NEW_Q
        call calculate_cart_gradient(b_l, a_l, d_l, vdw_l, coords, at_n, new_g_x)
        new_g_q = matmul(matmul(new_inv_G, new_B), new_g_x)
        
        print *, "Old e", tot_e, "New e", new_tot_e

        ! Calculate the updated GRMS
        grms = sqrt( dot_product(new_g_x, new_g_x)/n_x )
        print *, " "
        print *, "grms:", grms

        ! This is for debugging purposes
        call print_vector(q, size(q))
        call print_vector(new_q, size(q))
        call print_vector(s_q, size(q))
        call print_vector(g_q, size(q))
        call print_vector(new_g_q, size(q))

        ! Calculate the y and v vectors
        y_q(:)= new_g_q - g_q
        v_q(:) = matmul(inv_hess, y_q)
        ! Update the hessian
        new_inv_hess = get_new_hessian(n_q, inv_hess, s_q, y_q, v_q)

        call print_matrix(new_inv_hess, n_q, n_q)

        ! Assign the variables for the new iteration
        inv_hess = new_inv_hess; g_q = new_g_q; q = new_q; tot_e = new_tot_e
        inv_G = new_inv_G; B = new_B
        count = count + 1

    end do

    !coords = retransform_matrix(v_coords, n_x)

    print *, " "
    print *, "Final optimized energy", tot_e
    print *, " "
    print *, "Number of cycles: ", count

    print*, " "
    print *, "Optimized coordinates"
    do i=1, num_at
        print "(3F12.6)", coords(3*i-2:3*i)
    end do

end subroutine calculate_internal_gradient

end module gradient