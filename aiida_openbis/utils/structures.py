import numpy as np

import ipywidgets as ipw
#import nglview

from traitlets import Instance, default

from ase import Atoms
from ase.data import chemical_symbols

from sklearn.decomposition import PCA
from aiida import orm
from aiida_openbis.utils import bisutils

from openbabel import pybel as pb
from openbabel import openbabel as ob
from rdkit import Chem  # pylint: disable=(import-error)
from rdkit.Chem import AllChem  # pylint: disable=(import-error)


class OpenbisMolWidget(ipw.VBox):
    """Conver SMILES into 3D structure."""

    structure = Instance(orm.Node, allow_none=True)

    SPINNER = """<i class="fa fa-spinner fa-pulse" style="color:red;" ></i>"""

    def __init__(self, title=""):
        # pylint: disable=unused-import
        self.title = title

        try:
            from openbabel import pybel  # noqa: F401
            from openbabel import openbabel  # noqa: F401
        except ImportError:
            super().__init__([
                ipw.HTML(
                    "The SmilesWidget requires the OpenBabel library, "
                    "but the library was not found."
                )
            ])
            return
        try:
            from rdkit import Chem  # noqa: F401
            from rdkit.Chem import AllChem  # noqa: F401
        except ImportError:
            super().__init__([
                ipw.HTML(
                    "The SmilesWidget requires the rdkit library, "
                    "but the library was not found."
                )
            ])
            return
        #_output = ipw.Output()
        #with _output:
        self.session = bisutils.log_in(bisurl = "openbis", bisuser = "admin", bispasswd = "123456789")
        mols = bisutils.get_precursors(session=self.session)
        self.create_structure_btn = ipw.Button(description="Generate molecule", button_style="info")
        self.create_structure_btn.on_click(self._on_button_pressed)
        self.output = ipw.HTML("")
        
        
        # Dropdowns
        def observe_spaces(change):
            """Observe OpenBIS space change."""
            self.project_dropdown.options = [(project.code.capitalize(), project.code) for  project in self.session.get_projects(space=change['new'])]
            
        self.space_dropdown = ipw.Dropdown(
            description = "Space:",
            options = [(space.code.capitalize(), space.code) for space in self.session.get_spaces()],
        )
        self.space_dropdown.observe(observe_spaces, names='value')
        
        def observe_projects(change):
            """Observe OpenBIS project change."""
            self.collection_dropdown.options = [(collection.props["$name"].capitalize(), collection.code) for collection in self.session.get_collections(space=self.space_dropdown.value, project=change['new'])]
            
        self.project_dropdown = ipw.Dropdown(
            description="Project:",
            options=[(project.description or project.code, project.code) for project in self.session.get_projects(space=self.space_dropdown.value)],
        )
        
        self.project_dropdown.observe(observe_projects, names='value')
        
        def observe_collections(change):
            """Observe OpenBIS project change."""
            self.objects_dropdown.options = [(obj.props['iupac-name'], {'permId': obj.permId, 'smiles': obj.props["smiles"]}) for obj in self.session.get_objects(space=self.space_dropdown.value, project=self.project_dropdown.value, experiment = change['new'])]
            
        self.collection_dropdown = ipw.Dropdown(
            description="Collection:",
            options=[(collection.description or collection.code, collection.code) for collection in self.session.get_collections(space=self.space_dropdown.value, project=self.project_dropdown.value)],
        )
        
        self.collection_dropdown.observe(observe_collections, names='value')
        
        self.objects_dropdown = ipw.Dropdown(
            description="Object:"
#             options=[(str(mol[1])+": "+str(mol[2]), {'permId': mol[0], 'smiles':mol[3]}) for mol in mols]
        )
        
        #bisutils.log_out(session=self.session)

        super().__init__([ipw.HBox([self.space_dropdown, self.project_dropdown, self.collection_dropdown, self.objects_dropdown]), self.create_structure_btn, self.output])

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
            return self._rdkit_opt(smile, steps)
        except ValueError:
            return self._pybel_opt(smile, steps)

    def _on_button_pressed(self, change):  # pylint: disable=unused-argument
        """Convert SMILES to ase structure when button is pressed."""
        self.output.value = ""

        if not self.objects_dropdown.value:
            return

        self.output.value = "Screening possible conformers {}".format(
            self.SPINNER
        )  # font-size:20em;

        eln_info = {
            "eln_instance": "https://openbis-empa-lab205.labnotebook.ch/",
            "eln_type": "OpenBIS",
            #"sample_uuid": self.objects_dropdown.value.permId,
            "spectrum_type": "molecule",
            "file_name": self.objects_dropdown.label,
        }

        self.structure = orm.StructureData(ase=self.mol_from_smiles(self.objects_dropdown.value['smiles']))
        self.structure.set_extra("eln", eln_info)
        self.output.value = ""

    @default("structure")
    def _default_structure(self):
        return None
