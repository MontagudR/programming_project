module energy
    use constants
    use functions
    implicit none
    
contains
    
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

end module energy