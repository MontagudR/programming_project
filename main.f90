program opt
    use constants
    use functions
    use energy
    use gradient
    implicit none

    character(len=100) :: filename
    real*8, allocatable :: coordinates(:, :), bonds_lengths(:), bond_energies(:)
    real*8, allocatable :: distance_m(:,:), angles(:), bending_energies(:)
    real*8, allocatable :: dihedrals(:), torsion_energies(:), vdw_energies(:)
    real*8, allocatable :: total_cart_grad(:, :), opt_coords(:, :), v_coords(:)
    integer, allocatable :: bonds(:, :), atomic_numbers(:), angles_list(:, :)
    integer, allocatable :: conectivity_m(:,:), dihedrals_list(:, :), vdw_list(:, :)
    integer :: i, j, atom_num, n_q, n_x
    real*8 :: total_energy=0.
    
    print *, "Welcome to the optimizer, please select a file: "
    read(*, *) filename
    !filename = "ethane.mol2"

    ! Read and parse the coordinates and bonds from the mol2 file
    call read_mol2(filename, coordinates, bonds, atomic_numbers)
    ! Initialize all the constants
    call initialize_parameters()
    ! Get the number of atoms and n_x
    atom_num = size(coordinates, 1)
    n_x = 3*atom_num

    allocate(bonds_lengths(size(bonds, 1)), bond_energies(size(bonds, 1)), conectivity_m(atom_num, atom_num))
    allocate(distance_m(atom_num, atom_num), opt_coords(size(coordinates, 1), 3), v_coords(n_x))

    ! Put the coordinates in vector form ej x1, y1, z1, x2, y2 ...
    v_coords = transform_matrix(coordinates, n_x)
    ! Now we buid a connectivity matrix which is an n x n matrix which contains a 1 if there is a bond
    ! between the atoms or 0 otherwise
    conectivity_m = build_conectivity_m(bonds, atom_num)
    ! Whith the connectivity matrix, we get the angles and dihedrals of the molecule
    call get_angles_list(conectivity_m, atom_num, angles_list)
    call get_dihedrals_list(conectivity_m, atom_num, atomic_numbers, dihedrals_list)
    ! We get the vdw interaction list
    call get_vdw_lists(bonds, angles_list, atom_num, vdw_list)

    ! Get the number of internal coordinates from the bonds, angles and dihedrals list
    n_q = size(bonds, 1) + size(angles_list, 1) + size(dihedrals_list, 1)

    print *, "Bonds", size(bonds, 1)
    print *, "Angles", size(angles_list, 1)
    print *, "Torsions", size(dihedrals_list, 1)
    print *, "Vdw", size(vdw_list, 1)

    allocate(angles(size(angles_list)), bending_energies(size(angles_list)), vdw_energies(size(vdw_list)))
    allocate(dihedrals(size(dihedrals_list)), torsion_energies(size(dihedrals_list)))

    ! Calculate the stretching, bending, torsion and VdW energies
    call calculate_bond_properties(v_coords, bonds, atomic_numbers, bonds_lengths, bond_energies)
    call calculate_angle_properties(angles_list, v_coords, atomic_numbers, angles, bending_energies)
    call calculate_dihedral_properties(dihedrals_list, v_coords, atomic_numbers, dihedrals, torsion_energies)
    call calculate_vdw_energies(vdw_list, v_coords, atomic_numbers, vdw_energies)

    ! Sum everything to get the initial total energy of the molecule
    total_energy = sum(bond_energies) + sum(bending_energies) + sum(torsion_energies) + sum(vdw_energies)

    ! Print the total energy and the individual contributions for the starting geometry
    print *, "Total initial energy (kcal/mol)", total_energy
    print *, "Bond energy: ", sum(bond_energies), " Bending energies: ", sum(bending_energies), " Torsion energies", sum(torsion_energies), &
                " Vdw energies: ", sum(vdw_energies)

    allocate(total_cart_grad(atom_num, 3))

    ! Perform the cartesian optimization
    call cartesian_optimization(n_x, bonds, angles_list, dihedrals_list, vdw_list, v_coords, atomic_numbers, opt_coords)

    print *, "Internal gradient part"

    ! Perform the optimization in internal coordinates
    call calculate_internal_gradient(n_x, n_q, bonds, angles_list, dihedrals_list, vdw_list, v_coords, atomic_numbers)

    deallocate(coordinates, bonds, atomic_numbers, bonds_lengths, bond_energies, angles_list, conectivity_m, vdw_list)
    deallocate(distance_m, angles, bending_energies, dihedrals, torsion_energies)
    deallocate(total_cart_grad, v_coords)
    
end program opt