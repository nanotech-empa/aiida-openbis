import numpy as np
from scipy.stats import mode
import re

from IPython.display import clear_output
import ipywidgets as ipw
#import nglview

from traitlets import Instance, default

from ase import Atoms
from ase.data import covalent_radii, chemical_symbols
from ase.neighborlist import NeighborList
import ase.neighborlist

from sklearn.decomposition import PCA

from aiida_openbis.utils import bisutils


class OpenbisMolWidget(ipw.VBox):
    """Conver SMILES into 3D structure."""

    structure = Instance(Atoms, allow_none=True)

    SPINNER = """<i class="fa fa-spinner fa-pulse" style="color:red;" ></i>"""

    def __init__(self, title=""):
        # pylint: disable=unused-import
        self.title = title

        try:
            from openbabel import pybel  # noqa: F401
            from openbabel import openbabel  # noqa: F401
        except ImportError:
            super().__init__([
                ipw.HTML("The SmilesWidget requires the OpenBabel library, "
                         "but the library was not found.")
            ])
            return
        try:
            a = 1
            #from rdkit import Chem  # noqa: F401
            #from rdkit.Chem import AllChem  # noqa: F401
        except ImportError:
            super().__init__([
                ipw.HTML("The SmilesWidget requires the rdkit library, "
                         "but the library was not found.")
            ])
            return
        bisdata = bisutils.log_in()
        mols = bisutils.get_molecules(session=bisdata)
        print(mols)
        bisutils.log_out(session=bisdata)
        self.smiles = ipw.Dropdown(options=[(mol[0], mol[2]) for mol in mols])
        self.create_structure_btn = ipw.Button(description="Generate molecule",
                                               button_style="info")
        self.create_structure_btn.on_click(self._on_button_pressed)
        self.output = ipw.HTML("")

        super().__init__([self.smiles, self.create_structure_btn, self.output])

    def make_ase(self, species, positions):
        """Create ase Atoms object."""
        # Get the principal axes and realign the molecule along z-axis.
        positions = PCA(n_components=3).fit_transform(positions)
        atoms = Atoms(species, positions=positions, pbc=True)
        atoms.cell = np.ptp(atoms.positions, axis=0) + 10
        atoms.center()

        return atoms

    def _pybel_opt(self, smile, steps):
        """Optimize a molecule using force field and pybel (needed for complex SMILES)."""
        from openbabel import pybel as pb
        from openbabel import openbabel as ob

        obconversion = ob.OBConversion()
        obconversion.SetInFormat("smi")
        obmol = ob.OBMol()
        obconversion.ReadString(obmol, smile)

        pbmol = pb.Molecule(obmol)
        pbmol.make3D(forcefield="uff", steps=50)

        pbmol.localopt(forcefield="gaff", steps=200)
        pbmol.localopt(forcefield="mmff94", steps=100)

        f_f = pb._forcefields["uff"]  # pylint: disable=protected-access
        f_f.Setup(pbmol.OBMol)
        f_f.ConjugateGradients(steps, 1.0e-9)
        f_f.GetCoordinates(pbmol.OBMol)
        species = [chemical_symbols[atm.atomicnum] for atm in pbmol.atoms]
        positions = np.asarray([atm.coords for atm in pbmol.atoms])
        return self.make_ase(species, positions)

    def _rdkit_opt(self, smile, steps):
        """Optimize a molecule using force field and rdkit (needed for complex SMILES)."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        smile = smile.replace("[", "").replace("]", "")
        mol = Chem.MolFromSmiles(smile)
        mol = Chem.AddHs(mol)

        AllChem.EmbedMolecule(mol, maxAttempts=20, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol, maxIters=steps)
        positions = mol.GetConformer().GetPositions()
        natoms = mol.GetNumAtoms()
        species = [mol.GetAtomWithIdx(j).GetSymbol() for j in range(natoms)]
        return self.make_ase(species, positions)

    def mol_from_smiles(self, smile, steps=10000):
        """Convert SMILES to ase structure try rdkit then pybel"""
        try:
            #return self._rdkit_opt(smile, steps)
            return self._pybel_opt(smile, steps)
        except ValueError:
            return self._pybel_opt(smile, steps)

    def _on_button_pressed(self, change):  # pylint: disable=unused-argument
        """Convert SMILES to ase structure when button is pressed."""
        self.output.value = ""

        if not self.smiles.value:
            return

        self.output.value = "Screening possible conformers {}".format(
            self.SPINNER)  # font-size:20em;
        self.structure = self.mol_from_smiles(self.smiles.value)
        self.output.value = ""

    @default("structure")
    def _default_structure(self):
        return None
