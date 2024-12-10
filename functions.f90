module functions
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

function cross_prod(a, b) result(cross)
    real*8, intent(in) :: a(3), b(3)
    real*8 :: cross(3)

    cross(1) = a(2)*b(3) - b(2)*a(3)
    cross(2) = a(3)*b(1) - b(3)*a(1)
    cross(3) = a(1)*b(2) - b(1)*a(2)
end function cross_prod

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

end module functions
    
