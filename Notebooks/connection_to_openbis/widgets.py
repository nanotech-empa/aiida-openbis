import json
import os
import ipywidgets as ipw
from IPython.display import display, clear_output, Javascript
import sys
sys.path.append('/home/jovyan/aiida-openbis/Notebooks/importer')
import nanonis_importer
import shutil
import utils
import base64
from datetime import datetime
from aiida import orm

CONFIG = utils.read_json("config.json")
CONFIG_ELN = utils.read_json("eln_config.json")
OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis(CONFIG_ELN["url"], CONFIG_ELN["token"])

class AtomisticModelSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Atomistic Model', 
            disabled=False, 
            layout = ipw.Layout(width='982px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "ATOMISTIC_MODEL"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select atomistic model...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )

class SimulationSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Simulation', 
            disabled=False, 
            layout = ipw.Layout(width='982px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'PK', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        qb = orm.QueryBuilder()
        qb.append(orm.WorkChainNode)
        results = qb.all()
        
        # List of calculations that can be exported
        labels = ['QeAppWorkChain','Cp2kGeoOptWorkChain','Cp2kStmWorkChain']
        # Create the QueryBuilder
        qb = orm.QueryBuilder()
        qb.append(
            orm.WorkChainNode,
            filters={
                'attributes.process_label': {'in': labels},  # Filter by process_label
                'extras': {'!has_key': 'exported'}, # Exclude nodes with the 'exported' key in extras
                'attributes.process_state':{'in':['finished']},
            },
            project=['id','uuid', 'attributes.process_label','attributes.metadata_inputs.metadata.description','attributes.metadata_inputs.metadata.label']  # Project the PK (id) and process_label
        )
        # Execute the query
        results = qb.all()
        items_names_pks = []
        for result in results:
            if result[3]:
                name_pk_string = f"{result[3][:20]} - {result[2]} (PK: {result[0]})"
            else:
                name_pk_string = f"{result[2]} ({result[0]})"
                
            name_pk_tuple = (name_pk_string, result[0])
            items_names_pks.append(name_pk_tuple)
            
        items_names_pks.insert(0, (f'Select simulation...', -1))
        
        self.dropdown.options = items_names_pks
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PK"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PK"],
                    [True, False]
                ), 
                names='value'
            )

class AnalysisSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = 'Analysis', disabled = False, 
                                             layout = ipw.Layout(width = '800px'), 
                                             style = {'description_width': "110px"})
        
        self.children = [self.selector]
    
    def load_selector(self, project_id):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "ANALYSIS"
        )
        self.selector.options = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]

class CodeSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = 'Code', disabled = False, 
                                             layout = ipw.Layout(width = '800px'), 
                                             style = {'description_width': "110px"})
        
        self.children = [self.selector]
    
    def load_selector(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "CODE"
        )
        self.selector.options = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]

class ChemistSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Chemist', 
            disabled=False, 
            layout = ipw.Layout(width='982px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "CHEMIST"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select chemist...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )

class StorageSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Storage', 
            disabled=False, 
            layout = ipw.Layout(width='982px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "STORAGE"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select storage...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )

class SupplierSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Supplier', 
            disabled=False, 
            layout = ipw.Layout(width='982px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "SUPPLIER"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select supplier...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )

class PublicationSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Publication', 
            disabled=False, 
            layout = ipw.Layout(width='982px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "PUBLICATION_CUSTOM"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select publication...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )

class MeasurementSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = 'Measurements', disabled = False, 
                                             layout = ipw.Layout(width = '800px'), 
                                             style = {'description_width': "110px"})
        
        self.children = [self.selector]
    
    def load_selector(self, project_id):
        measurements_list = []
        oneD_measurements_items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "1D_MEASUREMENT",
            project = project_id
        )
        measurements_list.extend([(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in oneD_measurements_items])
        
        twoD_measurements_items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "2D_MEASUREMENT",
            project = project_id
        )
        
        measurements_list.extend([(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in twoD_measurements_items])
        self.selector.options = measurements_list
        
class DraftSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(
            description = 'Drafts', disabled = False, 
            layout = ipw.Layout(width = '800px'), 
            style = {'description_width': "110px"}
        )
        
        self.children = [self.selector]
    
    def load_selector(self, project_id):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "DRAFT",
            project = project_id
        )
        self.selector.options = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        
class ResultsSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = 'Results', disabled = False, 
                                             layout = ipw.Layout(width = '800px'), 
                                             style = {'description_width': "110px"})
        
        self.children = [self.selector]
    
    def load_selector(self, project_id):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "RESULTS",
            project = project_id
        )
        self.selector.options = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        
class GrantSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = 'Grants', disabled = False, 
                                             layout = ipw.Layout(width = '800px'), 
                                             style = {'description_width': "110px"})
        
        self.children = [self.selector]
    
    def load_selector(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "GRANT"
        )
        self.selector.options = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]

class AuthorSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = 'Authors', disabled = False, 
                                             layout = ipw.Layout(width = '800px'), 
                                             style = {'description_width': "110px"})
        
        self.children = [self.selector]
    
    def load_selector(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "AUTHOR"
        )
        self.selector.options = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        
class MoleculeSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Molecule', 
            disabled=False, 
            layout = ipw.Layout(width='315px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        # Details textbox (disabled for display purposes)
        self.details_textbox = ipw.Textarea(
            description="", disabled=True, layout=ipw.Layout(width='415px', height='250px')
        )
        
        # Image box (displaying a default image)
        self.image_box = ipw.Image(
            value=utils.read_file(CONFIG["default_image_filepath"]), format='jpg', width='220px', 
            height='250px', layout=ipw.Layout(border='solid 1px #cccccc')
        )
        
        self.children = [self.dropdown_boxes, self.details_textbox, self.image_box]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "MOLECULE"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select molecule...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )
    
    # Function to handle changes in the substances dropdown
    def load_molecule_metadata(self, change):
        if self.dropdown_boxes.children[0].value == -1:
            self.details_textbox.value = ''
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])
            return
        
        # Get metadata
        property_list = CONFIG["objects"]["Molecule"]["properties"]
        molecule_object = OPENBIS_SESSION.get_object(self.dropdown_boxes.children[0].value)
        
        # Get image
        molecule_dataset = molecule_object.get_datasets(type="ELN_PREVIEW")[0]

        if molecule_dataset:
            molecule_dataset.download(destination="images")
            material_image_filepath = molecule_dataset.file_list[0]
            self.image_box.value = utils.read_file(f"images/{molecule_dataset.permId}/{material_image_filepath}")
            shutil.rmtree(f"images/{molecule_dataset.permId}/")
        else:
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])

        molecule_metadata = molecule_object.props.all()
        molecule_metadata_string = ""
        for prop_key in property_list:
            prop_title = CONFIG["properties"][prop_key]["title"]
            prop_value = molecule_metadata.get(prop_key)
            if CONFIG["properties"][prop_key]["property_type"] == "QUANTITY_VALUE" and prop_value is not None:
                prop_value = json.loads(prop_value)
                molecule_metadata_string += f"{prop_title}: {prop_value['has_value']} {prop_value['has_unit']}\n"
            else:
                molecule_metadata_string += f"{prop_title}: {prop_value}\n"

        self.details_textbox.value = molecule_metadata_string

class SubstanceSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Substance', 
            disabled=False, 
            layout = ipw.Layout(width='315px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        # Details textbox (disabled for display purposes)
        self.details_textbox = ipw.Textarea(
            description="", disabled=True, layout=ipw.Layout(width='415px', height='250px')
        )
        
        # Image box (displaying a default image)
        self.image_box = ipw.Image(
            value=utils.read_file(CONFIG["default_image_filepath"]), format='jpg', width='220px', 
            height='250px', layout=ipw.Layout(border='solid 1px #cccccc')
        )
        
        self.children = [self.dropdown_boxes, self.details_textbox, self.image_box]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "SUBSTANCE"
        )
        items_names_permids = [(f"{item.props['empa_number']}{item.props['batch']} ({item.permId})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select substance...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )
    
    # Function to handle changes in the substances dropdown
    def load_substance_metadata(self, change):
        if self.dropdown_boxes.children[0].value == -1:
            self.details_textbox.value = ''
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])
            return
        
        # Get substance metadata
        property_list = CONFIG["objects"]["Substance"]["properties"]
        molecule_property_list = CONFIG["objects"]["Molecule"]["properties"]
        material_object = OPENBIS_SESSION.get_object(self.dropdown_boxes.children[0].value)
        
        # Get substance image
        molecule_object = OPENBIS_SESSION.get_object(material_object.props.all()['has_molecule'])
        molecule_metadata = molecule_object.props.all()
        molecule_dataset = molecule_object.get_datasets(type="ELN_PREVIEW")[0]

        if molecule_dataset:
            molecule_dataset.download(destination="images")
            material_image_filepath = molecule_dataset.file_list[0]
            self.image_box.value = utils.read_file(f"images/{molecule_dataset.permId}/{material_image_filepath}")
            shutil.rmtree(f"images/{molecule_dataset.permId}/")
        else:
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])

        material_metadata = material_object.props.all()
        material_metadata_string = ""
        for prop_key in property_list:
            prop_title = CONFIG["properties"][prop_key]["title"]
            prop_value = material_metadata.get(prop_key)
            if CONFIG["properties"][prop_key]["property_type"] == "QUANTITY_VALUE" and prop_value is not None:
                prop_value = json.loads(prop_value)
                material_metadata_string += f"{prop_title}: {prop_value['has_value']} {prop_value['has_unit']}\n"
            elif prop_key == "has_molecule":
                material_metadata_string += f"{prop_title}:\n"
                for mol_prop_key in molecule_property_list:
                    mol_prop_title = CONFIG["properties"][mol_prop_key]["title"]
                    mol_prop_value = molecule_metadata.get(mol_prop_key)
                    material_metadata_string += f"- {mol_prop_title}: {mol_prop_value}\n"
            else:
                material_metadata_string += f"{prop_title}: {prop_value}\n"

        self.details_textbox.value = material_metadata_string

class ReactionProductSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Reaction Product', 
            disabled=False, 
            layout = ipw.Layout(width='315px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        # Details textbox (disabled for display purposes)
        self.details_textbox = ipw.Textarea(
            description="", disabled=True, layout=ipw.Layout(width='415px', height='250px')
        )
        
        # Image box (displaying a default image)
        self.image_box = ipw.Image(
            value=utils.read_file(CONFIG["default_image_filepath"]), format='jpg', width='220px', 
            height='250px', layout=ipw.Layout(border='solid 1px #cccccc')
        )
        
        self.children = [self.dropdown_boxes, self.details_textbox, self.image_box]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "REACTION_PRODUCT"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select product...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )
    
    # Function to handle changes in the products dropdown
    def load_product_metadata(self, change):
        if self.dropdown_boxes.children[0].value == -1:
            self.details_textbox.value = ''
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])
            return
        
        # Get product metadata
        property_list = CONFIG["objects"]["Reaction Product"]["properties"]
        
        material_object = OPENBIS_SESSION.get_object(self.dropdown.value)
        material_metadata = material_object.props.all()
        material_metadata_string = ""
        for prop_key in property_list:
            prop_title = CONFIG["properties"][prop_key]["title"]
            prop_value = material_metadata.get(prop_key)
            if CONFIG["properties"][prop_key]["property_type"] == "QUANTITY_VALUE" and prop_value is not None:
                prop_value = json.loads(prop_value)
                material_metadata_string += f"{prop_title}: {prop_value['has_value']} {prop_value['has_unit']}\n"

        self.details_textbox.value = material_metadata_string
            
class MaterialSelectionWidget(ipw.Output):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='', 
            disabled=False, 
            layout = ipw.Layout(width='315px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.HBox([self.dropdown, self.sorting_checkboxes_list])
        
        # Details textbox (disabled for display purposes)
        self.details_textbox = ipw.Textarea(
            description="", disabled=True, layout=ipw.Layout(width='415px', height='250px')
        )
        
        # Image box (displaying a default image)
        self.image_box = ipw.Image(
            value=utils.read_file(CONFIG["default_image_filepath"]), format='jpg', width='220px', 
            height='250px', layout=ipw.Layout(border='solid 1px #cccccc')
        )
    
    def load_dropdown_box(self, object_type, placeholder):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = object_type
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select {placeholder}...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )

class ProjectSelectionWidget(ipw.VBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Project', 
            disabled=False, 
            layout = ipw.Layout(width='982px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_projects(
            OPENBIS_SESSION
        )
        items_names_permids = [(f"{item.description} ({item.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select project...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )

class ExperimentSelectionWidget(ipw.VBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.create_new_experiment_button = utils.Button(
            description = '', disabled = False, button_style = '', 
            tooltip = 'Add', icon = 'plus', 
            layout = ipw.Layout(width = '50px', height = '25px')
        )
        self.dropdown = utils.Dropdown(
            description='Experiment', 
            disabled=False, 
            layout = ipw.Layout(width = '993px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        self.dropdown_details = ipw.HBox([self.dropdown, self.create_new_experiment_button])
        self.add_experiment_output = ipw.Output()
        self.children = [self.dropdown_details, self.sorting_checkboxes_list, self.add_experiment_output]
        
        self.save_new_experiment_button = utils.Button(
            description = '', disabled = False, button_style = '', 
            tooltip = 'Save', icon = 'save', 
            layout = ipw.Layout(width = '50px', height = '35px', margin = '0 0 0 90px')
        )
        
        self.cancel_new_experiment_button = utils.Button(
            description = '', disabled = False, button_style = '', tooltip = 'Cancel', 
            icon = 'times', layout = ipw.Layout(width = '50px', height = '35px', margin = '0 0 0 5px')
        )
        
        self.new_experiment_name_textbox = utils.Text(
            description = "Name", disabled = False, layout = ipw.Layout(width = '400px'), 
            placeholder = f"Write experiment name here...", style = {'description_width': "110px"}
        )
        
        self.projects_dropdown_boxes = ProjectSelectionWidget()
        self.projects_dropdown_boxes.load_dropdown_box()
        
        self.create_new_experiment_button.on_click(self.create_new_experiment_button_on_click)
        self.cancel_new_experiment_button.on_click(self.cancel_new_experiment_button_on_click)
        self.save_new_experiment_button.on_click(self.save_new_experiment_button_on_click)
    
    def load_dropdown_box(self):
        items = utils.get_openbis_collections(
            OPENBIS_SESSION,
            type = "EXPERIMENT"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select experiment...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )
    
    def create_new_experiment_button_on_click(self, b):
        display_list = [
            self.new_experiment_name_textbox, self.projects_dropdown_boxes,
            ipw.HBox([self.save_new_experiment_button, self.cancel_new_experiment_button])
        ]
        
        with self.add_experiment_output:
            display(ipw.VBox(display_list))
    
    def cancel_new_experiment_button_on_click(self, b):
        with self.add_experiment_output:
            clear_output()
    
    def save_new_experiment_button_on_click(self, b):
        utils.create_experiment_in_openbis(OPENBIS_SESSION, self.projects_dropdown_boxes.children[0].children[0].value, self.new_experiment_name_textbox.value)
        self.load_dropdown_box()
        with self.add_experiment_output:
            clear_output()

class SampleSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Sample', 
            disabled=False, 
            layout = ipw.Layout(width='385px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        # Details textbox (disabled for display purposes)
        self.details_textbox = ipw.Textarea(
            description="", disabled=True, layout=ipw.Layout(width='589px', height='250px')
        )
        
        self.children = [self.dropdown_boxes, self.details_textbox]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "SAMPLE"
        )
        items = utils.filter_samples(OPENBIS_SESSION, items)
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select sample...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )

class InstrumentSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown = utils.Dropdown(
            description='Instrument', 
            disabled=False, 
            layout = ipw.Layout(width='982px'), 
            style = {'description_width': "110px"}, 
            options = [-1]
        )
        
        self.sorting_checkboxes_list = ipw.HBox(
            [
                ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                utils.Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                utils.Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
            ]
        )
        
        self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = "INSTRUMENT"
        )
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        items_names_permids.insert(0, (f'Select instrument...', -1))
        self.dropdown.options = items_names_permids
        self.dropdown.value = -1
        
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            ["Name", "PermID"],
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    ["Name", "PermID"],
                    [True, False]
                ), 
                names='value'
            )
            
class SamplePreparationSelectionWidget(ipw.SelectMultiple):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        self.sample_preparation_options = [object_key for object_key, object_info in CONFIG["objects"].items() if object_info["object_type"] == "sample_preparation"]
        self.description = "Processes"
        self.disabled = False
        self.style = {'description_width': "65px"}
        self.options = self.sample_preparation_options
        selector_height = len(self.sample_preparation_options) * 18
        selector_height = f"{selector_height}px"
        self.layout = ipw.Layout(width = '220px', height = selector_height)

class SamplePreparationAccordionWidget(ipw.Accordion):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        self.tasks_properties_widgets = []

class MultipleSelectorWidget(ipw.VBox):
    def __init__(self, selector_type):
        # Initialize the parent VBox
        super().__init__()
        self.selectors = []
        self.selector_type = selector_type
        
    def add_selector(self):
        if self.selector_type == "molecule":
            selector = MoleculeSelectionWidget()
            selector.load_dropdown_box()
            selector.dropdown.observe(selector.load_molecule_metadata, names = "value")
            
        elif self.selector_type == "product":
            selector = ReactionProductSelectionWidget()
            selector.load_dropdown_box()
            selector.dropdown.observe(selector.load_product_metadata, names = "value")
            
        self.selectors.append(selector)
        self.children = self.selectors
        
    def remove_selector(self):
        if self.selectors:
            self.selectors.pop()
            self.children = self.selectors
            
    def reset_selector_list(self):
        self.selectors = []
        self.children = self.selectors
        
class ObjectPropertiesWidgets(ipw.VBox):
    def __init__(self, task):
        super().__init__()
        
        self.properties_widgets = {}
        properties = CONFIG["objects"][task]["properties"]
        for prop_key in properties:
            property = CONFIG["properties"][prop_key]
            if property["property_widget"] == "TEXT":
                prop_widget = utils.Text(
                    description = property["title"], disabled = property["disabled"], 
                    layout = ipw.Layout(width = property["box_layout"]["width"]), 
                    placeholder = property["placeholder"], 
                    style = {'description_width': property["box_layout"]["description_width"]}
                )
                
            elif property["property_widget"] == "TEXTAREA":
                widget_args = {
                    "description": property["title"], "disabled": property["disabled"], 
                    "layout": ipw.Layout(width = property["box_layout"]["width"]), 
                    "placeholder": property["placeholder"], 
                    "style": {'description_width': property["box_layout"]["description_width"]}
                }
                if property["property_type"] == "JSON":
                    widget_args["value"] = property["default_value"]
                    
                prop_widget = utils.Textarea(**widget_args)
                    
            
            elif property["property_widget"] == "CHECKBOX":
                prop_widget = utils.Checkbox(
                    description = property["title"], disabled = property["disabled"], 
                    layout = ipw.Layout(width = property["box_layout"]["width"]),
                    style = {'description_width': property["box_layout"]["description_width"]},
                    indent = True, value = False
                )
            elif property["property_widget"] == "DATE":
                prop_widget = ipw.DatePicker(
                    description = property["title"], disabled = property["disabled"],
                    layout = ipw.Layout(width = property["box_layout"]["width"]),
                    style = {'description_width': property["box_layout"]["description_width"]}
                )
            elif property["property_widget"] == "INTTEXT" or property["property_widget"] == "FLOATTEXT":
                prop_widget = utils.IntText(
                    description = property["title"], disabled = property["disabled"], 
                    layout = ipw.Layout(width = property["box_layout"]["width"]), 
                    placeholder = property["placeholder"], 
                    style = {'description_width': property["box_layout"]["description_width"]}
                )
            
            elif property["property_widget"] == "DROPDOWN":
                prop_widget = utils.Dropdown(
                    description = property["title"], disabled = property["disabled"], 
                    layout = ipw.Layout(width = property["box_layout"]["width"]), 
                    options = property["options"],
                    value = property["value"],
                    style = {'description_width': property["box_layout"]["description_width"]}
                )
                
            elif property["property_widget"] == "FLOAT_TEXT_W_DROPDOWN":
                prop_widget = utils.FloatTextwithDropdownWidget(
                    property["title"], ipw.Layout(width = property["box_layout"]["width"]), 
                    property["default_value"], {'description_width': property["box_layout"]["description_width"]}, 
                    ipw.Layout(width = property["dropdown_layout"]["width"]), property["units"], property["default_unit"]
                )
            
            elif property["property_widget"] == "INTEGER_SLIDER_W_DETAILS":
                prop_widget = utils.IntSliderwithTextWidget(
                    property["default_value"], property["title"], property["values_list"], 
                    ipw.Layout(width = property["slider_layout"]["width"]), 
                    {'description_width': property["slider_layout"]["description_width"]}, 
                    property["placeholder"], ipw.Layout(width = property["box_layout"]["width"])
                )
            
            elif property["property_widget"] == "FLOAT_SLIDER":
                prop_widget = utils.FloatSlider(
                    description = property["title"], disabled = property["disabled"],
                    layout = ipw.Layout(width = property["slider_layout"]["width"]),
                    style = {'description_width': property["slider_layout"]["description_width"]},
                    value = property["default_value"], min = min(property["values_list"]),
                    max = max(property["values_list"]), step = property["slider_step"],
                    readout_format = '.2f'
                )
                
            self.properties_widgets[prop_key] = prop_widget
            
        self.children = list(self.properties_widgets.values())
        