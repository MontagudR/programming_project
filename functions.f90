module functions
    use constants
    implicit none
    
contains
    
subroutine read_mol2(filename, coordinates, bonds, atomic_numbers)
    character(len=100), intent(in) :: filename
    real*8, allocatable, intent(out) :: coordinates(:, :)
    integer, allocatable, intent(out) :: bonds(:, :), atomic_numbers(:)
    character(len=3), allocatable :: atom_symbols(:)
    integer :: atom_num, bond_num, i
    ! This subroutine reads the modified mol2 file in which
    ! the first line has the number of atoms and bonds,
    ! then reads the coordinates and after that
    ! reads the bond information

    open(10, file=filename, status='old')

    read(10, *) atom_num, bond_num
    allocate(coordinates(atom_num, 3), bonds(bond_num, 2), atom_symbols(atom_num), atomic_numbers(atom_num))

    ! read the coordinates
    do i=1, atom_num
        read(10, *) coordinates(i, 1), coordinates(i, 2), coordinates(i, 3), atom_symbols(i)
        ! Convert the atomic symbols to atomic numbers
        atomic_numbers(i) = to_atomic_number(atom_symbols(i))
    end do

    ! read the bond information
    do i=1, bond_num
        read(10, *)  bonds(i, 1), bonds(i, 2)
    end do

    close(10)
end subroutine read_mol2

integer function to_atomic_number(atom_symbol)
    implicit none
    character(len=3), intent(in) :: atom_symbol
    ! This function converts the atomic symbol in the atomic number.
    ! Only implemented the C and H, for the rest, the program stops.

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
    ! This function builds the connectivity matrix (n x n). For this,
    ! uses the bonds information, putting 1 if a bond exists or a 0
    ! if there is no bond in between those atoms

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
    ! This function gets the number of angles using the connectivity matrix.
    ! First, we look for a 1 in a row, that will indicate a bond between
    ! atom1 and atom2, then in the same row, we look for another 1 and if it
    ! exists, add the angle to the count.

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
    ! This function gets the atom indexes from the atoms that forms
    ! an angle. First we look for a 1 in a row, that will indicate a bond between
    ! atom1 and atom2, then in the same row, we look for another 1 and if it
    ! exists, that means that we have an angle ([atom2, atom1, atom3]) and 
    ! we append it to the list. 

    ! Get the number of angles, so we can allocate the exact 
    ! quantity of memory.
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
    ! This function gets the atom indexes for the atoms that forms
    ! a dihedral. First we look for a 1 in the matrix, that indicates
    ! a bond between atom1 and atom2. Then we look for another 1 in
    ! the same row, we have a bond between atom1 and atom3. Then,
    ! we look for a bond in the column of atom3. If it exists, 
    ! we add the dihedral [atom2, atom1, atom3, atom4]. We store them
    ! in a temporal array, and because of the method used, we will
    ! have duplicites, so at the end we filter the temporal list

    num_dihedral=0
    do i=1, atom_num
        do j=1, atom_num
            if (conectivity_m(i,j) == 1) then
                do z=1, atom_num
                    ! make sure that z and i are not equal
                    if ( conectivity_m(i, z) == 1 .and. z /= j) then
                        do k=1, atom_num
                            ! here we also add the restrain that the central atoms
                            ! of the dihedral should be carbon.
                            if ( conectivity_m(k, z) == 1 .and. k /= i .and. atom_list(i) == 6 .and. atom_list(z) == 6) then
                                num_dihedral = num_dihedral + 1
                                ! Stop the program if there are more dihedrals
                                ! than expected for the temporal list.
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

    ! We have always duplicities, create another
    ! empty list with half the capacity of the 
    ! total quantity of dihedrals
    allocate(dihedral_list(num_dihedral/2, 4))
    dihedral_list = 0
    not_found = .true.

    ! Two for loops, one iters over the temporal list and 
    ! the other over the empty list
    count = 0
    outer: do i=1, num_dihedral
        inner: do j=1, num_dihedral/2
            ! Check for the reverse array, the duplicate one.
            if ( all(tmp(i, :) .eq. dihedral_list(j, 4:1:-1)) ) then 
                not_found = .false.
                exit inner
            end if
        end do inner
        ! If it did not found a duplicate, store it in the final
        ! list
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
    ! This function gets a list of the index of atoms that forms a 
    ! VdW interaction. We iterate over all the possible paris. 
    ! First we check if there is no bond, and if there is no bond
    ! then, we check if they form angle. If it does not, we add
    ! the indexes to the list. We add them to a temporal list, 
    ! and when we have all the interactions

    vdw_count = 0
    do i=1, atom_num
        do j=1, atom_num
            if (j > i) then
                no_bond = .true.
                no_angle = .true.
                ! Check for bonds between the pairs
                do k=1, size(bonds_list, 1)
                    if (all(bonds_list(k, :) .eq. [i, j]) .or. all(bonds_list(k, :) .eq. [j, i])) then
                        no_bond = .false.
                        exit
                    end if
                end do
                if (no_bond .eqv. .true.) then
                    ! Check for angle between the pairs
                    do l=1, size(angles_list, 1)
                        if (any(angles_list(l, :) == i) .and. any(angles_list(l, :) == j)) then
                            no_angle = .false.
                            exit
                        end if
                    end do
                    if (no_angle .eqv. .true.) then
                        vdw_count = vdw_count + 1
                        ! Failsafe
                        if ( vdw_count > atom_num**2 ) then
                            print *, "Warning, more vdw interactions than expected!"
                        end if
                        tmp(vdw_count, :) = [i, j] 
                    end if
                end if
            end if
        end do
    end do

    allocate(vdw_lists(vdw_count, 2))

    ! Add the VdW to the final list
    do i=1, vdw_count
        vdw_lists(i, :) = tmp(i, :)
    end do


end subroutine get_vdw_lists

function cross_prod(a, b) result(cross)
    real*8, intent(in) :: a(3), b(3)
    real*8 :: cross(3)
    ! This function calculates the cross product between
    ! two vectors of dimension 3

    cross(1) = a(2)*b(3) - b(2)*a(3)
    cross(2) = a(3)*b(1) - b(3)*a(1)
    cross(3) = a(1)*b(2) - b(1)*a(2)
end function cross_prod

function transform_matrix(matrix, n_x) result(v_coordinates)
    real*8, intent(in) :: matrix(:, :)
    integer, intent(in) :: n_x
    real*8 :: v_coordinates(n_x)
    ! This function transform a matrix in vector form.
    ! Mainly used to transform the coordinates

    v_coordinates = reshape(transpose(matrix), (/ n_x /))

end function transform_matrix

function retransform_matrix(v_matrix, n_x) result(matrix)
    real*8, intent(in) :: v_matrix(:)
    integer, intent(in) :: n_x
    real*8 :: matrix(n_x/3, 3)
    ! This functions transforms a vector in a matrix.
    ! Mainly used for the coordinates

    matrix = reshape(v_matrix, (/n_x/3, 3/), order=[2,1])

end function retransform_matrix

function outer_product(A, B) result(outer)
    real*8, intent(in) :: A(:), B(:)
    real*8 :: outer(size(A), size(B))
    integer :: i, j
    ! This function calculatest the outer product
    ! between two vectors

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
    ! This function calculates the new hessian, using
    ! the old hessian and the s, y and v vectors

    new_hess = 0.d0
    new_hess = inv_hess + (( outer_product((dot_product(s, y) + dot_product(y, v))*s, s) ) / ( (dot_product(s, y))**2.d0 )) - &
                (( outer_product(v, s) + outer_product(s, v) ) / ( dot_product(s, y) ))

end function get_new_hessian

function scan_torsions(index, s_q, n_q) result(new_sq)
    integer, intent(in) :: index, n_q
    real*8, intent(in) :: s_q(:)
    real*8 :: new_sq(n_q)
    integer :: i
    ! This function scan the s_q used in the internal coordinates
    ! method. It scans the torsion change and if its above pi,
    ! changes it by the corresponding normalized value

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

subroutine print_matrix(matrix, r) 
    real*8, intent(in) :: matrix(:, :)
    integer, intent(in) :: r
    character(len=100) :: filename, g
    character(len=5) numb
    integer :: unit, i
    ! Auxiliary function to print the matrix in the
    ! output.out file

    filename = "output.out"
    ! Write the dimension r in the numb variable as a string
    write(numb, '(I0)') r
    ! Get the correct format
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
    integer :: unit
    ! Auxiliary function to print the vector in the
    ! output.out file

    filename = "output.out"
    ! Write the dimension r in the numb variable as a string
    write(numb, '(I0)') r
    ! Get the correct format
    g = "(" // trim(numb) // "F12.6)"
    unit = 12
    open(unit=unit, file=filename, status="unknown", position="append")

    write(12, g) vector(:)
    write(12, *) " "

    close(12)

end subroutine print_vector

end module functions
    
