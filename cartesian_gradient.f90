module cartesian_gradient
    use constants
    use parse_coordinates
    implicit none
    
contains
    
    subroutine cart_bond_gradient(bond_list, coordinates, atomic_numbers, cart_bond_grad)
        integer, intent(in) :: bond_list(:,:), atomic_numbers(:)
        real*8, intent(in) :: coordinates(:)
        real*8, intent(out) :: cart_bond_grad(size(atomic_numbers)*3)
        integer :: i, atom1, atom2, at_num1, at_num2, ind1, ind2
        real*8 :: r(3), grad(3)

        cart_bond_grad = 0.d0
        do i=1, size(bond_list, 1)
            atom1 = bond_list(i, 1); atom2 = bond_list(i, 2)
            at_num1 = atomic_numbers(atom1); at_num2 = atomic_numbers(atom2)
            ind1 = atom1*3 - 2; ind2 = atom2*3 - 2
            r(:) = coordinates(ind1:ind1+2) - coordinates(ind2:ind2+2)
            cart_bond_grad(ind1:ind1+2) = cart_bond_grad(ind1:ind1+2) + get_bond_grad(at_num1, at_num2, r)
            cart_bond_grad(ind2:ind2+2) = cart_bond_grad(ind2:ind2+2) - get_bond_grad(at_num1, at_num2, r)
            !print *, "Atom 1: ", atom1, "Atom2: ", atom2
            !print *, "Distance: ", norm2(r)
            !print *, "Grad: ", get_bond_grad(at_num1, at_num2, r)
        end do

    end subroutine cart_bond_gradient

    function get_bond_grad(atom1, atom2, r)
        integer, intent(in) :: atom1, atom2
        real*8, intent(in) :: r(:)
        real*8, dimension(3) :: get_bond_grad

        get_bond_grad = 2*k_bond(atom1, atom2)*(norm2(r) - r_bond(atom1, atom2))*(r(:)/norm2(r))

    end function get_bond_grad

    subroutine calculate_cart_angle_gradient(angle_list, coordinates, atomic_numbers, cart_angle_grad)
        integer, intent(in) :: angle_list(:,:), atomic_numbers(:)
        real*8, intent(in) :: coordinates(:)
        real*8, intent(out) :: cart_angle_grad(size(atomic_numbers)*3)
        integer :: i, atom1, atom2, atom3, at_num1, at_num2, at_num3, ind1, ind2, ind3
        real*8 :: r_ba(3), r_bc(3), angle, p(3), g_a(3), g_c(3), g_b(3), tmp

        cart_angle_grad = 0.d0
        do i=1, size(angle_list, 1)
            atom1 = angle_list(i, 1); atom2 = angle_list(i, 2); atom3 = angle_list(i, 3)
            at_num1 = atomic_numbers(atom1); at_num2 = atomic_numbers(atom2); at_num3 = atomic_numbers(atom3)
            ind1 = atom1*3 - 2; ind2 = atom2*3 - 2; ind3 = atom3*3 - 2
            r_ba(:) = coordinates(ind1:ind1+2) - coordinates(ind2:ind2+2)
            r_bc(:) =  coordinates(ind3:ind3+2) - coordinates(ind2:ind2+2)
            angle = get_angle(r_ba, r_bc)
            p(:) = cross_prod(r_ba, r_bc)
            tmp = 2*k_angle(at_num1, at_num2, at_num3)*(angle - degress_to_rad*rho_angle(at_num1, at_num2, at_num3))
            cart_angle_grad(ind1:ind1+2) = cart_angle_grad(ind1:ind1+2) + tmp*grad_angle_a(r_ba, p)
            cart_angle_grad(ind2:ind2+2) = cart_angle_grad(ind2:ind2+2) + tmp*grad_angle_b(r_ba, r_bc, p)
            cart_angle_grad(ind3:ind3+2) = cart_angle_grad(ind3:ind3+2) + tmp*grad_angle_c(r_bc, p)
        end do

    end subroutine calculate_cart_angle_gradient

    function grad_angle_a(r_ba, p)
        real*8, intent(in) :: r_ba(:), p(:)
        real*8 :: grad_angle_a(3)

        grad_angle_a = (cross_prod(r_ba, p))/((norm2(r_ba)**2)*norm2(p))

    end function grad_angle_a

    function grad_angle_b(r_ba, r_bc, p)
        real*8, intent(in) :: r_ba(:), p(:), r_bc(:)
        real*8 :: grad_angle_b(3)

        grad_angle_b = (cross_prod(-r_ba, p))/((norm2(r_ba)**2)*norm2(p)) + &
                        (cross_prod(r_bc, p))/((norm2(r_bc)**2)*norm2(p))

    end function grad_angle_b

    function grad_angle_c(r_bc, p)
        real*8, intent(in) :: r_bc(:), p(:)
        real*8 :: grad_angle_c(3)

        grad_angle_c = (cross_prod(-r_bc, p))/((norm2(r_bc)**2)*norm2(p))

    end function grad_angle_c

    subroutine calculate_cart_torsion_gradient(dihedral_list, coordinates, atomic_numbers, cart_dihedral_grad)
        integer, intent(in) :: dihedral_list(:,:), atomic_numbers(:)
        real*8, intent(in) :: coordinates(:)
        real*8, intent(out) :: cart_dihedral_grad(size(atomic_numbers)*3)
        integer :: i, atom1, atom2, atom3, atom4, ind1, ind2, ind3, ind4
        real*8 :: r_ab(3), r_bc(3), r_cd(3), dihedral, v(3), t(3), u(3), r_ac(3), r_bd(3), tmp1(3), tmp2(3), grad

        cart_dihedral_grad = 0.d0
        do i=1, size(dihedral_list, 1)
            atom1 = dihedral_list(i, 1); atom2 = dihedral_list(i, 2)
            atom3 = dihedral_list(i, 3); atom4 = dihedral_list(i, 4)
            ind1 = atom1*3 - 2; ind2 = atom2*3 - 2; ind3 = atom3*3 - 2; ind4 = atom4*3 - 2
            r_ab(:) = coordinates(ind2:ind2+2) - coordinates(ind1:ind1+2)
            r_ac(:) = coordinates(ind3:ind3+2) - coordinates(ind1:ind1+2)
            r_bc(:) = coordinates(ind3:ind3+2) - coordinates(ind2:ind2+2)
            r_bd(:) = coordinates(ind4:ind4+2) - coordinates(ind2:ind2+2)
            r_cd(:) = coordinates(ind4:ind4+2) - coordinates(ind3:ind3+2)
            t = cross_prod(r_ab, r_bc); u = cross_prod(r_bc, r_cd); v = cross_prod(t, u)
            dihedral = get_dihedral(r_ab, r_bc, r_cd)
            grad = -n_torsion*A_torsion*dsin(dihedral*n_torsion)
            cart_dihedral_grad(ind2:ind2+2) = cart_dihedral_grad(ind2:ind2+2) + grad_torsion_b(t, r_bc, r_ac, r_cd, u)*grad
            cart_dihedral_grad(ind3:ind3+2) = cart_dihedral_grad(ind3:ind3+2) + grad_torsion_c(t, r_bc, r_ab, r_bd, u)*grad
            cart_dihedral_grad(ind1:ind1+2) = cart_dihedral_grad(ind1:ind1+2) + grad_torsion_a(t, r_bc)*grad
            cart_dihedral_grad(ind4:ind4+2) = cart_dihedral_grad(ind4:ind4+2) + grad_torsion_d(u, r_bc)*grad
        end do
        
    end subroutine calculate_cart_torsion_gradient

    function grad_torsion_a(t, r_bc)
        real*8, intent(in) :: t(:), r_bc(:)
        real*8 :: grad_torsion_a(3)

        grad_torsion_a = cross_prod( (cross_prod(t,r_bc)/((norm2(t)**2)*norm2(r_bc))), r_bc )

    end function grad_torsion_a

    function grad_torsion_b(t, r_bc, r_ac, r_cd, u)
        real*8, intent(in) :: t(:), r_bc(:), r_ac(:), r_cd(:), u(:)
        real*8 :: grad_torsion_b(3)

        grad_torsion_b = cross_prod( r_ac, (cross_prod(t,r_bc)/((norm2(t)**2)*norm2(r_bc))) ) + &
                        cross_prod( (cross_prod(-u,r_bc)/((norm2(u)**2)*norm2(r_bc))), r_cd )

    end function grad_torsion_b

    function grad_torsion_c(t, r_bc, r_ab, r_bd, u)
        real*8, intent(in) :: t(:), r_bc(:), r_ab(:), r_bd(:), u(:)
        real*8 :: grad_torsion_c(3)

        grad_torsion_c = cross_prod( (cross_prod(t,r_bc)/((norm2(t)**2)*norm2(r_bc))), r_ab ) + &
                        cross_prod( r_bd, (cross_prod(-u,r_bc)/((norm2(u)**2)*norm2(r_bc))) )

    end function grad_torsion_c

    function grad_torsion_d(u, r_bc)
        real*8, intent(in) :: u(:), r_bc(:)
        real*8 :: grad_torsion_d(3)

        grad_torsion_d = cross_prod( (cross_prod(-u,r_bc)/((norm2(u)**2)*norm2(r_bc))), r_bc )

    end function grad_torsion_d

    subroutine calculate_cart_vdw_gradient(vdw_list, coordinates, atomic_numbers, cart_vdw_grad)
        integer, intent(in) :: vdw_list(:,:), atomic_numbers(:)
        real*8, intent(in) :: coordinates(:)
        real*8, intent(out) :: cart_vdw_grad(size(atomic_numbers)*3)
        integer :: i, atom1, atom2, at_num1, at_num2, ind1, ind2
        real*8 :: r_ij(3)

        cart_vdw_grad = 0.
        do i=1, size(vdw_list, 1)
            atom1 = vdw_list(i, 1); atom2 = vdw_list(i, 2)
            at_num1 = atomic_numbers(atom1); at_num2 = atomic_numbers(atom2)
            ind1 = atom1*3 - 2; ind2 = atom2*3 - 2
            r_ij(:) = coordinates(ind1:ind1+2) - coordinates(ind2:ind2+2)
            cart_vdw_grad(ind1:ind1+2) = cart_vdw_grad(ind1:ind1+2) + grad_vdw(at_num1, at_num2, r_ij)
            cart_vdw_grad(ind2:ind2+2) = cart_vdw_grad(ind2:ind2+2) - grad_vdw(at_num1, at_num2, r_ij)
        end do

    end subroutine calculate_cart_vdw_gradient

    function grad_vdw(atom1, atom2, r)
        integer, intent(in) :: atom1, atom2
        real*8, intent(in) :: r(:)
        real*8 :: grad_vdw(3)

        grad_vdw = r*( ((-12.0d0*A_vdw(atom1, atom2))/(norm2(r)**14.d0)) + ((6.0d0*B_vdw(atom1, atom2))/(norm2(r)**8.d0)) )

    end function grad_vdw

    subroutine calculate_cart_gradient(bond_list, angle_list, dihedral_list, vdw_list, coordinates, atomic_numbers, total_cart_grad)
        integer, intent(in) :: bond_list(:,:), angle_list(:,:), dihedral_list(:,:), vdw_list(:,:), atomic_numbers(:)
        real*8, intent(in) :: coordinates(:)
        real*8, intent(out) :: total_cart_grad(size(atomic_numbers)*3)
        real*8, dimension(size(atomic_numbers)*3) :: cart_strech_grad, cart_bending_grad, cart_torsion_grad, cart_vdw_grad
        integer :: i

        call cart_bond_gradient(bond_list, coordinates, atomic_numbers, cart_strech_grad)
        call calculate_cart_angle_gradient(angle_list, coordinates, atomic_numbers, cart_bending_grad)
        call calculate_cart_torsion_gradient(dihedral_list, coordinates, atomic_numbers, cart_torsion_grad)
        call calculate_cart_vdw_gradient(vdw_list, coordinates, atomic_numbers, cart_vdw_grad)

        !print *, " "
        !print *, "Indivual cartesian gradient: "
        !print *, "Cartesian bonding: "
        !do i=1, size(coordinates)/3
        !    print "(3F12.6)", cart_strech_grad(3*i-2:3*i)
        !end do
        !print *, "Cartesian VdW: "
        !do i=1, size(coordinates)/3
        !    print "(3F12.6)", cart_vdw_grad(3*i-2:3*i)
        !end do

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

        num_atoms = size(atomic_numbers)
        coords = coordinates

        inv_hes = 0.d0
        do i=1, size(atomic_numbers)*3
            inv_hes(i, i) = 1.d0/300.d0
        end do

        call calculate_cart_gradient(bond_list, angle_list, dihedral_list, vdw_list, coords, atomic_numbers, gradient)
        !tmp_grad(:) = transform_matrix(gradient, n_x)

        print *, "Starting iterations"
        count = 0
        grms = sqrt(dot_product(gradient, gradient)/dble(n_x))

        do while (grms > 0.001d0) !0.001d0

            p = matmul(inv_hes, -gradient)

            call calculate_total_energy(bond_list, angle_list, dihedral_list, vdw_list, coords, atomic_numbers, total_energy)
            alpha = 0.80d0; tmp = 0.80d0
            new_total_energy = 1.d0 + total_energy
            wolfe_rule = (total_energy + 0.1d0*alpha*dot_product(p, gradient))

            do while (new_total_energy >= wolfe_rule)

                alpha = tmp
                new_coords = coords + alpha*p
                call calculate_total_energy(bond_list, angle_list, dihedral_list, vdw_list, new_coords, atomic_numbers, new_total_energy)
                wolfe_rule = (total_energy + 0.1d0*alpha*dot_product(p, gradient))
                print *, "Alpha", alpha, "New energy", new_total_energy

                tmp = alpha * 0.80d0

            end do 

            s(:) = alpha*p 

            call calculate_cart_gradient(bond_list, angle_list, dihedral_list, vdw_list, new_coords, atomic_numbers, new_gradient)

            y(:) = new_gradient - gradient

            v(:) = matmul(inv_hes, y)

            new_inv_hes = get_new_hessian(n_x, inv_hes, s, y, v)

            grms = sqrt(dot_product(new_gradient, new_gradient)/dble(n_x))
            print *, "Gradient RMS: ", grms
            
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

    function outer_product(A, B) result(outer)
        real*8, intent(in) :: A(:), B(:)
        real*8 :: outer(size(A), size(B))
        integer :: i, j

        do i=1, size(A)
            do j=1, size(B)
                outer(i, j) = A(i)*B(j)
            end do
        end do

    end function outer_product

    function get_new_hessian(n, inv_hess, s, y, v) result(new_hess)
        implicit none
        real*8, intent(in) :: inv_hess(:, :), s(:), y(:), v(:)
        integer, intent(in) :: n
        real*8 :: new_hess(n, n)
    
        new_hess = 0.d0
        new_hess = inv_hess + (( outer_product((dot_product(s, y) + dot_product(y, v))*s, s) ) / ( (dot_product(s, y))**2.d0 )) - &
                    (( outer_product(v, s) + outer_product(s, v) ) / ( dot_product(s, y) ))
    
    end function get_new_hessian

end module cartesian_gradient