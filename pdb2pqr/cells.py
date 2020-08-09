"""Cell list to facilitate neighbor searching."""
import logging


_LOGGER = logging.getLogger(__name__)


class Cells:
    """Accelerate the search for nearby atoms.

    A pure all versus all search is O(n^2) - for every atom, every other atom
    must be searched.
    This is rather inefficient, especially for large proteins where cells may
    be tens of angstroms apart.
    The cell class breaks down the xyz protein space into several 3-D cells of
    desired size - then by simply examining atoms that fall into the adjacent
    cells one can quickly find nearby cells.
    """

    def __init__(self, cellsize):
        """Initialize the cells.

        :param cellsize:  the size of each cell (in Angstroms)
        :type cellsize:  int
        """
        self.cellmap = {}
        self.cellsize = cellsize

    def assign_cells(self, protein):
        """Place each atom in a virtual cell for easy neighbor comparison.

        :param protein:  protein with atoms to assign to cells
        :type protein:  Protein
        """
        for atom in protein.atoms:
            atom.cell = None
            self.add_cell(atom)

    def add_cell(self, atom):
        """Add an atom to the cell.

        :param atom:  the atom to add
        :type atom:  Atom
        """
        size = self.cellsize
        x = atom.x
        if x < 0:
            x = (int(x) - 1) // size * size
        else:
            x = int(x) // size * size
        y = atom.y
        if y < 0:
            y = (int(y) - 1) // size * size
        else:
            y = int(y) // size * size
        z = atom.z
        if z < 0:
            z = (int(z) - 1) // size * size
        else:
            z = int(z) // size * size
        key = (x, y, z)
        try:
            self.cellmap[key].append(atom)
        except KeyError:
            self.cellmap[key] = [atom]
        atom.cell = key

    def remove_cell(self, atom):
        """Remove an atom from a cell.

        :param atom:  the atom to remove
        :type atom:  Atom
        """
        oldcell = atom.cell
        if oldcell is None:
            return
        atom.cell = None
        self.cellmap[oldcell].remove(atom)

    def get_near_cells(self, atom):
        """Find all atoms in cells bordering an atom.

        :param atom:  the atom to test
        :type atom:  Atom
        :return:  a list of nearby atoms
        :rtype:  [Atom]
        """
        size = self.cellsize
        closeatoms = []
        cell = atom.cell
        if cell is None:
            return closeatoms
        else:
            x = cell[0]
            y = cell[1]
            z = cell[2]
            for i in range(-1 * size, 2 * size, size):
                for j in range(-1 * size, 2 * size, size):
                    for k in range(-1 * size, 2 * size, size):
                        newkey = (x + i, y + j, z + k)
                        try:
                            newatoms = self.cellmap[newkey]
                            for atom2 in newatoms:
                                if atom == atom2:
                                    continue
                                closeatoms.append(atom2)
                        except KeyError:
                            pass
            return closeatoms
