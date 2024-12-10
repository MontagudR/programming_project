module internal_gradient
    use parse_coordinates
    use constants
    use cartesian_gradient
    implicit none
    
contains

subroutine build_B_matrix(n_x, n_q, b_l, a_l, d_l, vdw_l, coords, at_n, B)
    integer, intent(in) :: b_l(:, :), a_l(: ,:), d_l(:, :), vdw_l(:, :), at_n(:), n_x, n_q
    real*8, intent(in) :: coords(:)
    integer :: i, j, count, at_1, at_2, at_3, at_4, ind1, ind2, ind3, ind4
    real*8 :: r(3), p(3), r_ba(3), r_bc(3), r_ab(3), r_cd(3), r_bd(3), r_ac(3), t(3), u(3)
    real*8, intent(out) :: B(n_q, n_x)

    count = 1

    B = 0.
    do i=1, size(b_l, 1)
        at_1 = b_l(i, 1); at_2 = b_l(i, 2)
        ind1 = at_1*3 - 2; ind2 = at_2*3 - 2
        r = 0.
        r = coords(ind1:ind1+2) - coords(ind2:ind2+2)
        B(count, ind1:ind1+2) =  r(:) / norm2(r)
        B(count, ind2:ind2+2) = -1*r(:) / norm2(r)
        count = count + 1
    end do

    do i=1, size(a_l, 1)
        at_1 = a_l(i, 1); at_2 = a_l(i, 2); at_3 = a_l(i, 3);
        ind1 = at_1*3 - 2; ind2 = at_2*3 - 2; ind3 = at_3*3 - 2
        r_ba = coords(ind1:ind1+2) - coords(ind2:ind2+2); r_bc = coords(ind3:ind3+2) - coords(ind2:ind2+2)
        p = cross_prod(r_ba, r_bc)
        B(count, ind1:ind1+2) = grad_angle_a(r_ba, p)
        B(count, ind2:ind2+2) = grad_angle_b(r_ba, r_bc, p)
        B(count, ind3:ind3+2) = grad_angle_c(r_bc, p)
        count = count + 1
    end do

    do i=1, size(d_l, 1)
        at_1 = d_l(i, 1); at_2 = d_l(i, 2); at_3 = d_l(i, 3); at_4 = d_l(i, 4)
        ind1 = at_1*3 - 2; ind2 = at_2*3 - 2; ind3 = at_3*3 - 2; ind4 = at_4*3 - 2
        r_ab = coords(ind2:ind2+2) - coords(ind1:ind1+2); r_bc = coords(ind3:ind3+2) - coords(ind2:ind2+2)
        r_bd = coords(ind4:ind4+2) - coords(ind2:ind2+2); r_ac = coords(ind3:ind3+2) - coords(ind1:ind1+2)
        r_cd = coords(ind4:ind4+2) - coords(ind3:ind3+2)
        t = cross_prod(r_ab, r_bc); u = cross_prod(r_bc, r_cd)
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

    G = matmul(B, transpose(B))

    allocate(work(1))
    lwork = -1
    call dsyev('V', 'U', n_q, G, n_q, w, work, lwork, info)   
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    
    call dsyev('V', 'U', n_q, G, n_q, w, work, lwork, info)

    inv_w_m = 0.d0
    do i=1, size(w)
        if (w(i) < 1e-5) then
            inv_w_m(i, i) = 0.d0
        else
            inv_w_m(i, i) = 1.d0 / w(i)
        end if
    end do

    inv_G = matmul(matmul(G, inv_w_m), transpose(G))

end subroutine calculate_G_inv

function initial_inv_hess_internal(b_l, a_l, d_l, n_q) result(inv_hess)
    integer, intent(in) :: b_l(:, :), a_l(:, :), d_l(:, :), n_q
    real*8 :: inv_hess(n_q, n_q)
    integer :: i, count

    inv_hess = 0.
    count = 1
    do i=1, size(b_l, 1)
        inv_hess(count, count) = 1.d0 / 600.d0
        count = count + 1
    end do

    do i=1, size(a_l, 1)
        inv_hess(count, count) = 1.d0 / 150.d0
        count = count + 1
    end do

    do i=1, size(d_l, 1)
        inv_hess(count, count) = 1.d0 / 80.d0
        count = count + 1
    end do

end function initial_inv_hess_internal

function scan_torsions(index, s_q, n_q) result(new_sq)
    integer, intent(in) :: index, n_q
    real*8, intent(in) :: s_q(:)
    real*8 :: new_sq(n_q)
    integer :: i

    new_sq = s_q

    do i=index, n_q
        if (abs(s_q(i)) > pi ) then
            if (s_q(i) > 0.) then
                
                new_sq(i) = s_q(i) -2.d0*pi
                
            else
                
                new_sq(i) = s_q(i) + 2.d0*pi
                
            end if
        end if
    end do

end function scan_torsions

subroutine find_opt_cart(b_l, a_l, d_l, vdw_l, at_n, n_q, n_x, q, new_q, s_q, B, inv_G, x)
    real*8, intent(inout) :: x(:), new_q(:)
    real*8, intent(in) :: B(:, :), inv_G(:,:), q(:), s_q(:)
    integer, intent(in) :: n_q, n_x, at_n(:), b_l(:, :), a_l(:, :), d_l(:, :), vdw_l(:, :)
    real*8 :: new_x(n_x), fix_q(n_q), tmp, new_sq(n_q), tmp_sq(n_q), diff
    integer :: count
    count = 0

    fix_q = s_q + q
    diff = 1
    new_sq = s_q
    
    do while ( diff > 1.d-5 )

        new_x = x + matmul(matmul(transpose(B), inv_G), new_sq)
        call calculate_interanls_and_total_energy(b_l, a_l, d_l, vdw_l, new_x, at_n, n_q, tmp, new_q)
        new_sq = fix_q - new_q
        new_sq = scan_torsions(size(b_l, 1) + size(a_l, 1) + 1, new_sq, n_q)

        diff = maxval(abs(new_x - x))
        print *, "Maximum change: ", diff
        x = new_x
        count = count + 1

    end do

    print *, " "
    print *, "Cartersian fiited after: ", count
    print *, " "
    !print *, "New q", new_q(:)
    !print *, " "

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

    coords = i_coords
    lamda_max = 0.02d0
    count = 0

    ! Calculation of gradient in internal coordinates
    call build_B_matrix(n_x, n_q, b_l, a_l, d_l, vdw_l, coords, at_n, B)

    num_at = size(at_n)
    
    call calculate_G_inv(B, n_q, inv_G)

    call calculate_cart_gradient(b_l, a_l, d_l, vdw_l, coords, at_n, g_x)
    g_q = matmul(matmul(inv_G, B), g_x)
    grms = sqrt(dot_product(g_x, g_x)/n_x)

    ! Calculate the initial hessian)

    inv_hess = initial_inv_hess_internal(b_l, a_l, d_l, n_q)
    !call print_matrix(inv_hess, n_q, n_q)
    call calculate_interanls_and_total_energy(b_l, a_l, d_l, vdw_l, coords, at_n, n_q, tot_e, q)
    !v_coords = transform_matrix(coords, n_x)

    ! Start the iterations
    do while( grms > 1d-3 ) ! 1e-3

        p_q = -1 * matmul(inv_hess, g_q)
        l_q = sqrt( (sum(p_q**2)) / (n_q))
        
        if (lamda_max < l_q) then
            p_q = p_q * (lamda_max / l_q)
        end if

        new_q = q + p_q
        s_q = new_q - q
        s_q = scan_torsions(size(b_l, 1) + size(a_l, 1) + 1, s_q, n_q)

        call find_opt_cart(b_l, a_l, d_l, vdw_l, at_n, n_q, n_x, q, new_q, s_q, B, inv_G, coords)

        !print *, " "
        !print *, "Updated coordinates"
        !do i=1, num_at
        !    print "(3F12.6)", coords(3*i-2:3*i)
        !end do

        ! CAREFUL, now v_coords and s_q are the updated
        s_q = new_q - q
        s_q = scan_torsions(size(b_l, 1) + size(a_l, 1) + 1, s_q, n_q)
        !print *, " "
        !print *, "New sq: "
        !print *, s_q

        ! Get the upgraded internal gradient
        call build_B_matrix(n_x, n_q, b_l, a_l, d_l, vdw_l, coords, at_n, new_B)
        call calculate_G_inv(new_B, n_q, new_inv_G)
        call calculate_total_energy(b_l, a_l, d_l, vdw_l, coords, at_n, new_tot_e) ! WARNING NO NEW_Q
        call calculate_cart_gradient(b_l, a_l, d_l, vdw_l, coords, at_n, new_g_x)
        new_g_q = matmul(matmul(new_inv_G, new_B), new_g_x)
        
        print *, "Old e", tot_e, "New e", new_tot_e

        grms = sqrt( dot_product(new_g_x, new_g_x)/n_x )
        print *, " "
        !do i=1, num_at
        !    print "(3F12.6)", new_g_x(3*i-2:3*i)
        !end do
        print *, "grms:", grms

        ! Approximate the new hessian

        call print_vector(q, size(q))
        call print_vector(new_q, size(q))
        call print_vector(s_q, size(q))
        call print_vector(g_q, size(q))
        call print_vector(new_g_q, size(q))

        y_q(:)= new_g_q - g_q
        v_q(:) = matmul(inv_hess, y_q)
        new_inv_hess = get_new_hessian(n_q, inv_hess, s_q, y_q, v_q)
        call print_matrix(new_inv_hess, n_q, n_q)

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

subroutine print_matrix(matrix, r, c) 
    real*8, intent(in) :: matrix(:, :)
    integer, intent(in) :: r, c
    character(len=100) :: filename, g
    character(len=5) numb
    integer :: unit, i, j

    filename = "output.out"
    write(numb, '(I0)') r
    g = "(" // trim(numb) // "F12.6)"
    unit = 12
    open(unit=unit, file=filename, status="unknown", position="append")

    do i=1, r
        write(12, g) matrix(i, :)
    end do

    write(12, *) " "

    close(12)

end subroutine print_matrix

subroutine print_vector(vector, r) 
    real*8, intent(in) :: vector(:)
    integer, intent(in) :: r
    character(len=100) :: filename, g
    character(len=5) numb
    integer :: unit, i, j

    filename = "output.out"
    write(numb, '(I0)') r
    g = "(" // trim(numb) // "F12.6)"
    unit = 12
    open(unit=unit, file=filename, status="unknown", position="append")

    write(12, g) vector(:)
    write(12, *) " "

    close(12)

end subroutine print_vector
    
end module internal_gradient