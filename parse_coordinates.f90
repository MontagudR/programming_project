module parse_coordinates
    use constants
    implicit none
    
contains
    
    subroutine read_mol2(filename, coordinates, bonds, atomic_numbers)
        character(len=100), intent(in) :: filename
        real*8, allocatable, intent(out) :: coordinates(:, :)
        integer, allocatable, intent(out) :: bonds(:, :), atomic_numbers(:)
        character(len=3), allocatable :: atom_symbols(:)
        integer :: atom_num, bond_num, i

        open(10, file=filename, status='old')

        read(10, *) atom_num, bond_num
        allocate(coordinates(atom_num, 3), bonds(bond_num, 2), atom_symbols(atom_num), atomic_numbers(atom_num))

        do i=1, atom_num
            read(10, *) coordinates(i, 1), coordinates(i, 2), coordinates(i, 3), atom_symbols(i)
            atomic_numbers(i) = to_atomic_number(atom_symbols(i))
        end do

        do i=1, bond_num
            read(10, *)  bonds(i, 1), bonds(i, 2)
        end do

        close(10)
    end subroutine read_mol2

    integer function to_atomic_number(atom_symbol)
        implicit none
        character(len=3), intent(in) :: atom_symbol

        select case (atom_symbol)
            case ('C')
                to_atomic_number = 6
            case ('H')
                to_atomic_number = 1
            case default
                print *, "Coulnd't find atomic number. Are the symbols properly written?"
        end select

    end function to_atomic_number

    subroutine calculate_bond_properties(coordinates, bond_list, atomic_numbers, bond_lengths, bond_energies)
        real*8, intent(in) :: coordinates(:)
        integer, intent(in) :: bond_list(:, :), atomic_numbers(:)
        real*8, intent(inout) :: bond_lengths(:), bond_energies(:)
        integer :: i, atom1, atom2, at_n1, at_n2
        
        do i=1, size(bond_list, 1)
            atom1 = bond_list(i, 1); atom2 = bond_list(i, 2)
            at_n1=atomic_numbers(atom1); at_n2=atomic_numbers(atom2)
            bond_lengths(i) = norm2(coordinates(3*atom1-2:3*atom1) - coordinates(3*atom2-2:3*atom2))
            bond_energies(i) = bond_energy(at_n1, at_n2, bond_lengths(i))
            !print *, atom_1, atom_2, bond_lengths(i), bond_energies(i)
        end do

    end subroutine calculate_bond_properties

    real*8 function bond_energy(atom1, atom2, bond_length)
        integer, intent(in) :: atom1, atom2
        real*8, intent(in) :: bond_length

        bond_energy = k_bond(atom1, atom2)*(bond_length - r_bond(atom1, atom2))**2

    end function bond_energy

    function build_conectivity_m(bond_list, atom_num) result(conectivity_m)
        integer, intent(in) :: bond_list(:,:), atom_num
        integer :: conectivity_m(atom_num, atom_num)
        integer:: i, atom_1, atom_2

        conectivity_m = 0
        do i=1, size(bond_list, 1)
            atom_1 = bond_list(i, 1)
            atom_2 = bond_list(i, 2)
            conectivity_m(atom_1, atom_2) = 1
            conectivity_m(atom_2, atom_1) = 1
        end do
    end function build_conectivity_m
    
    subroutine quantity_of_angles(conectivity_m, atom_num, num_angles)
        integer, intent(in) :: conectivity_m(:, :), atom_num
        integer, intent(out) :: num_angles
        integer :: i, j, z

        num_angles=0
        do i=1, atom_num
            do j=1, atom_num
                if (conectivity_m(i,j) == 1) then
                    do z=j+1, atom_num
                        if ( conectivity_m(i, z) == 1 ) then
                            num_angles = num_angles + 1
                        end if
                    end do
                end if
            end do
        end do
    end subroutine quantity_of_angles

    subroutine get_angles_list(conectivity_m, atom_num, angles_list)
        integer, intent(in) :: conectivity_m(:, :), atom_num
        integer, allocatable, intent(out) :: angles_list(:, :)
        integer :: i, j, z, angle_num, count

        call quantity_of_angles(conectivity_m, atom_num, angle_num)
        allocate(angles_list(angle_num, 3))

        count=1
        do i=1, atom_num
            do j=1, atom_num
                if (conectivity_m(i,j) == 1) then
                    do z=j+1, atom_num
                        if ( conectivity_m(i, z) == 1 ) then
                            angles_list(count, :) = [j, i, z]
                            count = count + 1
                        end if
                    end do
                end if
            end do
        end do

    end subroutine get_angles_list

    subroutine get_dihedrals_list(conectivity_m, atom_num, atom_list, dihedral_list)
        integer, intent(in) :: conectivity_m(:, :), atom_num, atom_list(:)
        integer, allocatable, intent(out) :: dihedral_list(:, :)
        integer :: i, j, z, k, tmp(atom_num**2, 4), count, num_dihedral
        logical :: not_found

        num_dihedral=0
        do i=1, atom_num
            do j=1, atom_num
                if (conectivity_m(i,j) == 1) then
                    do z=1, atom_num
                        if ( conectivity_m(i, z) == 1 .and. z /= j) then
                            do k=1, atom_num
                                if ( conectivity_m(k, z) == 1 .and. k /= i .and. atom_list(i) == 6 .and. atom_list(z) == 6) then
                                    !print *, j, i, z, k
                                    num_dihedral = num_dihedral + 1
                                    if (num_dihedral > atom_num**2) then
                                        print *, "Too many dihedrals"
                                    end if
                                    tmp(num_dihedral, :) = [j, i, z, k]
                                end if
                            end do
                        end if
                    end do
                end if
            end do
        end do

        ! We have always duplicities
        allocate(dihedral_list(num_dihedral/2, 4))
        dihedral_list = 0
        not_found = .true.

        count = 0
        outer: do i=1, num_dihedral
            inner: do j=1, num_dihedral/2
                if ( all(tmp(i, :) .eq. dihedral_list(j, 4:1:-1)) ) then 
                    not_found = .false.
                    exit inner
                end if
            end do inner
            if (not_found .eqv. .true.) then
                count = count + 1
                if (count > num_dihedral) then
                    print *, "Too many dihedrals, something is wrong"
                end if
                dihedral_list(count, :) = tmp(i, :)
            end if
            not_found = .true.
        end do outer

    end subroutine get_dihedrals_list

    subroutine get_vdw_lists(bonds_list, angles_list, atom_num, vdw_lists)
        integer, intent(in) :: bonds_list(:, :), angles_list(:,:), atom_num
        integer, allocatable, intent(out) :: vdw_lists(:, :)
        integer :: tmp(atom_num**2, 2), i, j, k, l, vdw_count
        logical :: no_bond, no_angle

        vdw_count = 0
        do i=1, atom_num
            do j=1, atom_num
                if (j > i) then
                    no_bond = .true.
                    no_angle = .true.
                    do k=1, size(bonds_list, 1)
                        if (all(bonds_list(k, :) .eq. [i, j]) .or. all(bonds_list(k, :) .eq. [j, i])) then
                            no_bond = .false.
                            exit
                        end if
                    end do
                    if (no_bond .eqv. .true.) then
                        do l=1, size(angles_list, 1)
                            if (any(angles_list(l, :) == i) .and. any(angles_list(l, :) == j)) then
                                no_angle = .false.
                                exit
                            end if
                        end do
                        if (no_angle .eqv. .true.) then
                            vdw_count = vdw_count + 1
                            tmp(vdw_count, :) = [i, j]
                            if ( vdw_count == atom_num**2 ) then
                                print *, "Warning, more vdw interactions than expected!"
                            end if
                            !print *, tmp(vdw_count, :)
                        end if
                    end if
                end if
            end do
        end do

        allocate(vdw_lists(vdw_count, 2))

        do i=1, vdw_count
            vdw_lists(i, :) = tmp(i, :)
        end do

    end subroutine get_vdw_lists

    subroutine calculate_angle_properties(angles_list, coordinates, atomic_numbers, angles, bending_energy)
        integer, intent(in) :: angles_list(:, :), atomic_numbers(:)
        real*8, intent(in) ::  coordinates(:)
        real*8, intent(out) :: angles(size(angles_list, 1)), bending_energy(size(angles_list, 1))
        integer :: i, atom1, atom2, atom3
        real*8 :: r_ba(3), r_bc(3)

        do i=1, size(angles_list, 1)
            atom1 = angles_list(i, 1); atom2 = angles_list(i, 2); atom3 = angles_list(i, 3)
            r_ba(:) = coordinates(3*atom1-2:3*atom1) - coordinates(3*atom2-2:3*atom2)
            r_bc(:) =  coordinates(3*atom3-2:3*atom3) - coordinates(3*atom2-2:3*atom2)
            angles(i) = get_angle(r_ba, r_bc)
            bending_energy(i) = get_bending_energy(atomic_numbers(atom1), atomic_numbers(atom2), atomic_numbers(atom3), &
                                angles(i))
        end do
    end subroutine calculate_angle_properties

    real*8 function get_angle(r_ba, r_bc)
        real*8, intent(in) :: r_ba(3), r_bc(3)

        get_angle = dacos((dot_product(r_ba, r_bc))/(norm2(r_ba)*norm2(r_bc)))
    end function get_angle

    real*8 function get_bending_energy(atom1, atom2, atom3, angle)
        integer, intent(in) :: atom1, atom2, atom3
        real*8, intent(in) :: angle

        get_bending_energy = k_angle(atom1, atom2, atom3)*(angle - degress_to_rad*rho_angle(atom1, atom2, atom3))**2.d0

    end function get_bending_energy

    function cross_prod(a, b) result(cross)
        real*8, intent(in) :: a(3), b(3)
        real*8 :: cross(3)

        cross(1) = a(2)*b(3) - b(2)*a(3)
        cross(2) = a(3)*b(1) - b(3)*a(1)
        cross(3) = a(1)*b(2) - b(1)*a(2)
    end function cross_prod

    subroutine calculate_dihedral_properties(dihedral_list, coordinates, atomic_numbers, dihedrals, torsion_energy)
        integer, intent(in) :: dihedral_list(:, :), atomic_numbers(:)
        real*8, intent(in) ::  coordinates(:)
        real*8, intent(out) :: dihedrals(size(dihedral_list, 1)), torsion_energy(size(dihedral_list, 1))
        integer :: i, atom1, atom2, atom3, atom4
        real*8 :: r_ab(3), r_bc(3), r_cd(3), t(3), u(3), v(3), tmp1, tmp2

        do i=1, size(dihedral_list, 1)
            atom1 = dihedral_list(i, 1); atom2 = dihedral_list(i, 2); atom3 = dihedral_list(i, 3); atom4 = dihedral_list(i, 4)
            r_ab(:) = coordinates(3*atom2-2:3*atom2) - coordinates(3*atom1-2:3*atom1)
            r_bc(:) =  coordinates(3*atom3-2:3*atom3) - coordinates(3*atom2-2:3*atom2)
            r_cd(:) = coordinates(3*atom4-2:3*atom4) - coordinates(3*atom3-2:3*atom3)
            dihedrals(i) = get_dihedral(r_ab, r_bc, r_cd)
            torsion_energy(i) = get_torsion_energy(dihedrals(i))
            !print *, dihedrals(i), torsion_energy(i)
        end do
    end subroutine calculate_dihedral_properties

    real*8 function get_dihedral(r_ab, r_bc, r_cd)
        real*8, intent(in) :: r_ab(:), r_bc(:), r_cd(:)
        real*8 :: t(3), u(3), v(3)

        t = cross_prod(r_ab, r_bc); u = cross_prod(r_bc, r_cd); v = cross_prod(t,u)
        get_dihedral = datan2( (dot_product(r_bc,v))/(norm2(r_bc)*norm2(t)*norm2(u)), &
                        (dot_product(t,u))/(norm2(t)*norm2(u)) )

    end function get_dihedral

    real*8 function get_torsion_energy(dihedral)
        real*8, intent(in) :: dihedral

        get_torsion_energy = A_torsion*(1 + dcos(n_torsion*dihedral))

    end function get_torsion_energy

    subroutine calculate_vdw_energies(vdw_list, coordinates, atomic_numbers, vdw_energy)
        integer, intent(in) :: vdw_list(:, :), atomic_numbers(:)
        real*8, intent(in) :: coordinates(:)
        real*8, intent(out) :: vdw_energy(size(vdw_list, 1))
        integer :: i, atom1, atom2
        real*8 :: r

        do i=1, size(vdw_list, 1)
            atom1 = vdw_list(i, 1); atom2 = vdw_list(i, 2)
            r = norm2(coordinates(3*atom1-2:3*atom1) - coordinates(3*atom2-2:3*atom2))
            vdw_energy(i) = get_vdw_energy(atomic_numbers(atom1), atomic_numbers(atom2), r)
            !print *, atomic_numbers(atom1), atomic_numbers(atom2), vdw_energy(i) 
        end do

    end subroutine calculate_vdw_energies

    real*8 function get_vdw_energy(atom1, atom2, r)
        integer, intent(in) :: atom1, atom2
        real*8, intent(in) :: r

        get_vdw_energy = (A_vdw(atom1, atom2) / (r**12.d0)) - (B_vdw(atom1, atom2) / (r**6.d0))

    end function get_vdw_energy

    subroutine calculate_total_energy(bond_list, angle_list, torsion_list, vdw_list, coordinates, atomic_numbers, total_energy)
        integer, intent(in) :: bond_list(:, :), angle_list(:, :), torsion_list(:, :), vdw_list(:, :), atomic_numbers(:)
        real*8, intent(in) :: coordinates(:)
        real*8, intent(out) :: total_energy
        real*8 :: bonds_lengths(size(bond_list, 1)), bond_energies(size(bond_list, 1))
        real*8 :: angles(size(angle_list, 1)), bending_energies(size(angle_list, 1))
        real*8 :: dihedrals(size(torsion_list, 1)), torsion_energies(size(torsion_list, 1))
        real*8 :: vdw_energies(size(vdw_list, 1))

        total_energy = 0.d0
        call calculate_bond_properties(coordinates, bond_list, atomic_numbers, bonds_lengths, bond_energies)
        call calculate_angle_properties(angle_list, coordinates, atomic_numbers, angles, bending_energies)
        call calculate_dihedral_properties(torsion_list, coordinates, atomic_numbers, dihedrals, torsion_energies)
        call calculate_vdw_energies(vdw_list, coordinates, atomic_numbers, vdw_energies)

        total_energy = sum(bond_energies) + sum(bending_energies) + sum(torsion_energies) + sum(vdw_energies)

    end subroutine calculate_total_energy

    subroutine calculate_interanls_and_total_energy(bond_list, angle_list, torsion_list, vdw_list, coordinates, &
                atomic_numbers, n_q, total_energy, q)
        integer, intent(in) :: bond_list(:, :), angle_list(:, :), torsion_list(:, :), vdw_list(:, :), atomic_numbers(:), n_q
        real*8, intent(in) :: coordinates(:)
        real*8, intent(out) :: total_energy, q(n_q)
        real*8 :: bonds_lengths(size(bond_list, 1)), bond_energies(size(bond_list, 1))
        real*8 :: angles(size(angle_list, 1)), bending_energies(size(angle_list, 1))
        real*8 :: dihedrals(size(torsion_list, 1)), torsion_energies(size(torsion_list, 1))
        real*8 :: vdw_energies(size(vdw_list, 1))

        call calculate_bond_properties(coordinates, bond_list, atomic_numbers, bonds_lengths, bond_energies)
        call calculate_angle_properties(angle_list, coordinates, atomic_numbers, angles, bending_energies)
        call calculate_dihedral_properties(torsion_list, coordinates, atomic_numbers, dihedrals, torsion_energies)
        call calculate_vdw_energies(vdw_list, coordinates, atomic_numbers, vdw_energies)

        total_energy = sum(bond_energies) + sum(bending_energies) + sum(torsion_energies) + sum(vdw_energies)
        q = [bonds_lengths, angles, dihedrals]

    end subroutine calculate_interanls_and_total_energy

    function transform_matrix(matrix, n_x) result(v_coordinates)
        real*8, intent(in) :: matrix(:, :)
        integer, intent(in) :: n_x
        real*8 :: v_coordinates(n_x)

        v_coordinates = reshape(transpose(matrix), (/ n_x /))

    end function transform_matrix

    function retransform_matrix(v_matrix, n_x) result(matrix)
        real*8, intent(in) :: v_matrix(:)
        integer, intent(in) :: n_x
        real*8 :: matrix(n_x/3, 3)

        matrix = reshape(v_matrix, (/n_x/3, 3/), order=[2,1])

    end function retransform_matrix

end module parse_coordinates