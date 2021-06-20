from ase import Atoms, Atom
from ase.lattice.surface import fcc110
from ase.constraints import FixAtoms, Hookean

def make_atoms():
    constraints = []
    atoms = fcc111('Pt', (2, 2, 1), vacuum=12.)
    # Build and add the CO adsorbates.
    for site in [1, 3]:
        ads = Atoms([Atom('C', (0., 0., 0.)),
                     Atom('C', (0., 0., 1.))])
        ads.translate(atoms[site].position + (0., 0., 1.))
        atoms.extend(ads)
        constraints.append(Hookean(a1=len(atoms) - 1,
                                   a2=len(atoms) - 2,
                                   k=10., rt=1.58))
        constraints.append(Hookean(a1=len(atoms) - 2,
                                   a2=(0., 0., 1.,
                                       -(atoms[site].z + 2.2)),
                                   k=10.))
    constraints.append(FixAtoms(indices=[atom.index for atom in atoms
                                         if atom.symbol == 'Pt']))
    atoms.set_constraint(constraints)
    atoms.set_pbc(True)
    return atoms

if __name__ == '__main__':
    atoms = make_atoms()
    from ase.visualize import view
    view(atoms)