import ipywidgets as ipw
import utils
from IPython.display import display
import pandas as pd
import json
import re
from collections import Counter
import subprocess
import shutil

OPENBIS_SAMPLES_CACHE = {}
ACTIONS_CODES = {
    "ANNEALING": "HEAT",
    "COOLDOWN": "COOL",
    "SPUTTERING": "IONB",
    "DOSING": "GASD",
    "DEPOSITION": "DEPO",
    "MECHANICAL_PRESSING": "MECH",
    "FIELD_EMISSION": "FIEM",
    "ETCHING": "ETCH",
    "LIGHT_IRRADIATION": "LITE",
    "DELAMINATION": "DELA",
    "FISHING": "FISH",
    "COATING": "COAT",
    "RINSE": "RINS"
}

ACTIONS_LABELS = [
    ("Annealing", "ANNEALING"),
    ("Coating", "COATING"),
    ("Cooldown", "COOLDOWN"),
    ("Delamination", "DELAMINATION"),
    ("Deposition", "DEPOSITION"),
    ("Dosing", "DOSING"),
    ("Etching", "ETCHING"),
    ("Field emission", "FIELD_EMISSION"),
    ("Fishing", "FISHING"),
    ("Light irradiation", "LIGHT_IRRADIATION"),
    ("Mechanical Pressing", "MECHANICAL_PRESSING"),
    ("Rinse", "RINSE"),
    ("Sputtering", "SPUTTERING")
]

OBSERVABLES_LABELS = [
    ("Current", "CURRENT_OBSERVABLE"),
    ("Elemental Composition", "ELEMENTAL_COMPOSITION_OBSERVABLE"),
    ("Flux", "FLUX_OBSERVABLE"),
    ("Force", "FORCE_OBSERVABLE"),
    ("Inductance", "INDUCTANCE_OBSERVABLE"),
    ("Observable", "OBSERVABLE"),
    ("pH", "PH_VALUE_OBSERVABLE"),
    ("Pressure", "PRESSURE_OBSERVABLE"),
    ("Resistance", "RESISTANCE_OBSERVABLE"),
    ("Speed", "SPEED_OBSERVABLE"),
    ("Temperature", "TEMPERATURE_OBSERVABLE"),
    ("Voltage", "VOLTAGE_OBSERVABLE")
]

MATERIALS_LABELS = [
    ("Crystal", "CRYSTAL"),
    ("2D layer material", "2D_LAYER_MATERIAL"),
    ("Wafer substrate", "WAFER_SUBSTRATE")
]

MATERIALS_CONCEPTS_LABELS = [
    ("Crystal concept", "CRYSTAL_CONCEPT"),
    ("2D layer material", "2D_LAYER_MATERIAL"),
    ("Wafer substrate", "WAFER_SUBSTRATE")
]

SIMULATION_OBJECT_TYPES = [
    "BAND_STRUCTURE",
    "GEOMETRY_OPTIMISATION",
    "PDOS",
    "MEASUREMENT_SESSION",
    "UNCLASSIFIED_SIMULATION",
    "VIBRATIONAL_SPECTROSCOPY"
]

WORKCHAIN_VIEWERS = {
    'QeAppWorkChain': "/apps/apps/quantum-espresso/qe.ipynb",
    'Cp2kGeoOptWorkChain': "/apps/apps/surfaces/view_geometry_optimization.ipynb",
    'Cp2kStmWorkChain': "/apps/apps/surfaces/view_stm.ipynb",
}

def get_cached_object(ob_session, obj_id):
    if obj_id in OPENBIS_SAMPLES_CACHE:
        sample_object = OPENBIS_SAMPLES_CACHE[obj_id]
    else:
        sample_object = utils.get_openbis_object(ob_session, sample_ident = obj_id)
        # Save object information
        OPENBIS_SAMPLES_CACHE[obj_id] = sample_object
        
    return sample_object

def stringify_quantity_value(json_obj, unit_type):
    if isinstance(json_obj, str):
        json_obj = json.loads(json_obj)
        
    value = json_obj["value"]
    unit = json_obj[unit_type]
    quantity_value_str = ""
    if value:
        quantity_value_str = f"{value} {unit}"
        
    return quantity_value_str

def find_openbis_simulations(ob_session, obj):
    simulation_objects = set()
    children = obj.children
    if children is not None:
        for child in children:
            child_object = get_cached_object(ob_session, child)
            if child_object.type in SIMULATION_OBJECT_TYPES:
                simulation_objects.add(child_object)
            simulation_objects.update(find_openbis_simulations(ob_session, child_object))
    return simulation_objects

# Import/export simulations widgets
class ImportSimulationsWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_molecules_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select molecules</span>"
        )
        
        self.select_reacprod_concepts_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select reaction product concepts</span>"
        )
        
        self.select_slab_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select slab</span>"
        )
        
        self.search_simulations_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Search simulations</span>"
        )
        
        self.molecules_accordion = ipw.Accordion()
        self.add_molecule_button = ipw.Button(
            description = 'Add', disabled = False, 
            button_style = 'success', tooltip = 'Add molecule', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.reacprod_concepts_accordion = ipw.Accordion()
        self.add_reacprod_concept_button = ipw.Button(
            description = 'Add', disabled = False, 
            button_style = 'success', tooltip = 'Add reaction product concept', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.select_material_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select material</span>"
        )
        
        material_type_options = MATERIALS_CONCEPTS_LABELS
        material_type_options.insert(0, ("Select material type...", "-1"))
        
        self.material_type_dropdown = ipw.Dropdown(
            options = material_type_options,
            value = material_type_options[0][1]
        )
        
        self.material_details_vbox = ipw.VBox()
        
        self.search_logical_operator_label = ipw.Label(value = "Search logical operator")
        self.search_logical_operator_dropdown = ipw.Dropdown(
            value = "AND", options = ["AND", "OR"],
            layout = ipw.Layout(width = '150px')
        )
        self.search_logical_operator_hbox = ipw.HBox(
            children = [self.search_logical_operator_label, self.search_logical_operator_dropdown]
        )
        
        self.search_button = ipw.Button(
            disabled = False, icon = 'search',
            tooltip = 'Search simulations in openBIS', 
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        self.found_simulations_label = ipw.Label(value = "Found simulations")
        self.found_simulations_select_multiple = ipw.SelectMultiple(
            description = '', disabled = False, 
            layout = ipw.Layout(width = '500px'), 
            style = {'description_width': "110px"}
        )
        self.found_simulations_hbox = ipw.HBox(
            children = [self.found_simulations_label, self.found_simulations_select_multiple]
        )
        
        self.import_simulations_button = utils.Button(
            tooltip = 'Import simulations', icon = 'download', 
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        self.import_simulations_message_html = ipw.HTML()
        
        # Increase search button icon size
        increase_search_button = ipw.HTML(
            """<style>
            .fa-search {font-size: 2em !important;}
            .fa-download {font-size: 2em !important;}
            </style>
            """
        )
        
        # Add functionality to the widgets
        self.material_type_dropdown.observe(self.load_material_type_widgets, names = "value")
        self.add_molecule_button.on_click(self.add_molecule)
        self.add_reacprod_concept_button.on_click(self.add_reacprod_concept)
        self.search_button.on_click(self.search_simulations)
        
        self.children = [
            self.select_molecules_title,
            self.molecules_accordion,
            self.add_molecule_button,
            self.select_reacprod_concepts_title,
            self.reacprod_concepts_accordion,
            self.add_reacprod_concept_button,
            self.select_slab_title,
            self.material_type_dropdown,
            self.material_details_vbox,
            self.search_simulations_title,
            self.search_logical_operator_hbox,
            increase_search_button,
            self.search_button,
            self.found_simulations_hbox,
            self.import_simulations_button,
            self.import_simulations_message_html
        ]
    
    def search_simulations(self, b):
        parents_permid_list = []
        for molecule_widget in self.molecules_accordion.children:
            molecule_permid = molecule_widget.dropdown.value
            if molecule_permid != "-1":
                parents_permid_list.append(molecule_permid)
        
        if self.material_type_dropdown.value != "-1":
            material_permid = self.material_details_vbox.children[0].children[0].value
            if material_permid != "-1":
                parents_permid_list.append(material_permid)
        
        for reacprod_concept_widget in self.reacprod_concepts_accordion.children:
            reacprod_concept_permid = reacprod_concept_widget.dropdown.value
            if reacprod_concept_permid != "-1":
                parents_permid_list.append(reacprod_concept_permid)
        
        simulation_permid_set = set()
        logical_operator = self.search_logical_operator_dropdown.value
        
        if logical_operator == "OR": # In OR, all the simulations found are added to the list
            for parent in parents_permid_list:
                parent_object = get_cached_object(self.openbis_session, parent)
                simulation_objects_children = find_openbis_simulations(self.openbis_session, parent_object)
                
                for simulation_object in simulation_objects_children:
                    simulation_permid = simulation_object.permId
                    # Measurement Session is used for both simulation and experiments
                    if simulation_object.type == "MEASUREMENT_SESSION":
                        simulation_measurement = False
                        for parent in simulation_object.parents:
                            parent_object = get_cached_object(self.openbis_session, parent)
                            if parent_object.type == "ATOMISTIC_MODEL":
                                simulation_measurement = True
                                break
                            
                        if simulation_measurement:
                            simulation_permid_set.add(simulation_permid)
                    else:
                        simulation_permid_set.add(simulation_permid)
                        
        else:  # In AND, only the simulations that appear in all selected materials are added to the list
            for idx, parent in enumerate(parents_permid_list):
                parent_object = get_cached_object(self.openbis_session, parent)
                simulation_objects_children = find_openbis_simulations(self.openbis_session, parent_object)
                
                parent_simulation_permid_list = []
                for simulation_object in simulation_objects_children:
                    simulation_permid = simulation_object.permId
                    if simulation_object.type == "MEASUREMENT_SESSION": # 2D Measurement is used for both simulation and experiments
                        simulation_measurement = False
                        for parent in simulation_object.parents:
                            parent_object = get_cached_object(self.openbis_session, parent)
                            if parent_object.type == "ATOMISTIC_MODEL":
                                simulation_measurement = True
                                break
                            
                        if simulation_measurement:
                            parent_simulation_permid_list.append(simulation_permid)
                    else:
                        parent_simulation_permid_list.append(simulation_permid)
                
                if idx == 0:   
                    simulation_permid_set = set(parent_simulation_permid_list)
                else:
                    simulation_permid_set.intersection_update(parent_simulation_permid_list)

        simulation_permid_list = list(simulation_permid_set)
        aiida_node_permid_list = []
        for simulation_permid in simulation_permid_list:
            simulation_object = get_cached_object(self.openbis_session, simulation_permid)
            aiida_node_found = False
            for parent in simulation_object.parents:
                parent_object = get_cached_object(self.openbis_session, parent)
                if parent_object.type == "AIIDA_NODE":
                    aiida_node_found = True
                    aiida_node_permid_list.append(parent)
                    break
            
            if aiida_node_found == False:
                aiida_node_permid_list.append(None)
        
        simulation_aiida_node_list = []
        for idx, simulation_permid in enumerate(simulation_permid_list):
            simulation_object = get_cached_object(self.openbis_session, simulation_permid)
            simulation_info = f"{simulation_object.props['name']}"
            aiida_node_permid = aiida_node_permid_list[idx]
            simulation_aiida_node_list.append([simulation_info, aiida_node_permid])
        
        self.found_simulations_select_multiple.options = simulation_aiida_node_list
    
    def import_aiida_nodes(self, b):
        selected_simulations = self.found_simulations_select_multiple.value
        selected_labels = [label for label, value in self.found_simulations_select_multiple.options if value in selected_simulations]
        
        selected_simulations_messages = ""
        for idx, aiida_node_permid in enumerate(selected_simulations):
            selected_label = selected_labels[idx]
            
            if aiida_node_permid:
                aiida_node_object = get_cached_object(self.openbis_session, aiida_node_permid)
                object_datasets = aiida_node_object.get_datasets()
                
                for dataset in object_datasets:
                    dataset_filenames = dataset.file_list
                    is_aiida_file = False
                    if len(dataset_filenames) == 1:
                        for filename in dataset_filenames:
                            if ".aiida" in filename:
                                is_aiida_file = True
                    
                    if is_aiida_file:
                        dataset.download(destination = 'aiida_nodes')
                        aiida_node_filename = dataset.file_list[0]
                        aiida_node_filepath = f"aiida_nodes/{dataset.permId}/{aiida_node_filename}"
                        command = ["verdi", "archive", "import", aiida_node_filepath]
                        
                        # Execute the command
                        result = subprocess.run(command, capture_output=True, text=True)
                        if result.returncode != 0:
                            print(f"An error occurred: {result.stderr}")
                        else:
                            workchain = load_node(aiida_node_object.props["wfms_uuid"])
                            workchain_viewer_link = WORKCHAIN_VIEWERS[workchain.process_label]
                            notebook_link = f"{workchain_viewer_link}?pk={workchain.pk}"
                            
                            simulation_message = f'<a href="{notebook_link}">Workchain {selected_label} successfully imported.</a>\n'
                            selected_simulations_messages += simulation_message
                        
                        shutil.rmtree(f"aiida_nodes/{dataset.permId}/")
            else:
                simulation_message = f"Workchain {selected_label} cannot be imported because it was done manually.\n"
                selected_simulations_messages += simulation_message
    
    def load_material_type_widgets(self, change):
        if self.material_type_dropdown.value == "-1":
            self.material_details_vbox.children = []
            return
        else:
            material_options = [
                ("Select material...", "-1")
            ]
            
            material_dropdown = ipw.Dropdown(
                options = material_options,
                value = material_options[0][1]
            )
            
            sort_material_label = ipw.Label(
                value = "Sort by:", 
                layout=ipw.Layout(margin='0px', width='50px'),
                style = {'description_width': 'initial'}
            )
            
            name_checkbox = ipw.Checkbox(
                indent = False,
                layout=ipw.Layout(margin='2px', width='20px')
            )
            
            name_label = ipw.Label(
                value = "Name", 
                layout=ipw.Layout(margin='0px', width='50px'),
                style = {'description_width': 'initial'}
            )
            
            registration_date_checkbox = ipw.Checkbox(
                indent = False,
                layout=ipw.Layout(margin='2px', width='20px')
            )
            
            registration_date_label = ipw.Label(
                value = "Registration date",
                layout=ipw.Layout(margin='0px', width='110px'),
                style = {'description_width': 'initial'}
            )
        
            select_material_box = ipw.HBox(
                children = [
                    material_dropdown,
                    sort_material_label,
                    name_checkbox,
                    name_label,
                    registration_date_checkbox,
                    registration_date_label
                ]
            )
            
            material_details_html = ipw.HTML()
            
            self.material_details_vbox.children = [
                select_material_box,
                material_details_html
            ]
            
            material_type = self.material_type_dropdown.value
            material_objects = utils.get_openbis_objects(
                self.openbis_session,
                type = material_type
            )
            materials_objects_names_permids = [(obj.props["$name"], obj.permId) for obj in material_objects]
            material_options += materials_objects_names_permids
            material_dropdown.options = material_options
            
            def sort_material_dropdown(change):
                options = material_options[1:]
                
                df = pd.DataFrame(options, columns=["$name", "registration_date"])
                if name_checkbox.value and not registration_date_checkbox.value:
                    df = df.sort_values(by="$name", ascending=True)
                elif not name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by="registration_date", ascending=False)
                elif name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by=["$name", "registration_date"], ascending=[True, False])

                options = list(df.itertuples(index=False, name=None))
                options.insert(0, material_options[0])
                material_dropdown.options = options
            
            def load_material_details(change):
                obj_permid = material_dropdown.value
                if obj_permid == "-1":
                    return
                else:
                    obj = self.openbis_session.get_object(obj_permid)
                    obj_props = obj.props.all()
                    obj_name = obj_props.get("$name", "")
                    obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
                    for key, value in obj_props.items():
                        if value:
                            prop_type = self.openbis_session.get_property_type(key)
                            prop_label = prop_type.label
                            obj_details_string += f"<p><b>{prop_label}:</b> {value}</p>"
                    
                    obj_details_string += "</div>"
                    
                    material_details_html.value = obj_details_string
                    
            name_checkbox.observe(sort_material_dropdown, names = "value")
            registration_date_checkbox.observe(sort_material_dropdown, names = "value")
            material_dropdown.observe(load_material_details, names = "value")
    
    def add_molecule(self, b):
        molecules_accordion_children = list(self.molecules_accordion.children)
        molecule_index = len(molecules_accordion_children)
        molecule_widget = MoleculeWidget(self.openbis_session, self.molecules_accordion, molecule_index)
        molecules_accordion_children.append(molecule_widget)
        self.molecules_accordion.children = molecules_accordion_children
    
    def add_reacprod_concept(self, b):
        reacprod_concepts_accordion_children = list(self.reacprod_concepts_accordion.children)
        reacprod_concept_index = len(reacprod_concepts_accordion_children)
        reacprod_concept_widget = ReacProdConceptWidget(self.openbis_session, self.reacprod_concepts_accordion, reacprod_concept_index)
        reacprod_concepts_accordion_children.append(reacprod_concept_widget)
        self.reacprod_concepts_accordion.children = reacprod_concepts_accordion_children

class ExportSimulationsWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_experiment_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select experiment</span>"
        )
        
        self.simulation_details_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Simulation details</span>"
        )
        
        self.used_aiida_checkbox = ipw.Checkbox(
            value = False,
            description = "Simulation developed using AiiDA",
            indent = False
        )
        
        self.select_molecules_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select molecules</span>"
        )
        
        self.select_reacprod_concepts_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select reaction product concepts</span>"
        )
        
        self.select_slab_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select slab</span>"
        )
        
        self.simulation_details_vbox = ipw.VBox()
        
        self.children = [
            self.select_experiment_title,
            self.simulation_details_title,
            self.used_aiida_checkbox,
            self.select_molecules_title,
            self.select_reacprod_concepts_title,
            self.select_slab_title,
            self.simulation_details_vbox
        ]

class MoleculeWidget(ipw.VBox):
    def __init__(self, openbis_session, parent_accordion, object_index):
        super().__init__()
        self.openbis_session = openbis_session
        self.parent_accordion = parent_accordion
        self.object_index = object_index
        self.title = ""
        
        molecules_objects = utils.get_openbis_objects(
            self.openbis_session, 
            collection = "/MATERIALS/MOLECULES/PRECURSOR_COLLECTION",
            type = "MOLECULE"
        )
        dropdown_list = [(obj.props["$name"], obj.permId) for obj in molecules_objects]
        dropdown_list.insert(0, ("Select a molecule...", "-1"))
        self.dropdown = ipw.Dropdown(value = "-1",options = dropdown_list)
        self.details_vbox = ipw.VBox()
        
        self.remove_molecule_button = ipw.Button(
            description = 'Remove', disabled = False, 
            button_style = 'danger', tooltip = 'Remove molecule', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.dropdown.observe(self.load_details, names = "value")
        self.remove_molecule_button.on_click(self.remove_molecule)
        self.children = [self.dropdown, self.details_vbox, self.remove_molecule_button]
    
    def load_details(self, change):
        obj_permid = self.dropdown.value
        if obj_permid == "-1":
            return
        else:
            obj = get_cached_object(self.openbis_session, obj_permid)
            obj_props = obj.props.all()
            obj_name = obj_props.get("$name", "")
            self.parent_accordion.set_title(self.object_index, obj_name)
            self.title = obj_name
            
            obj_details_html = ipw.HTML()
            obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
            for key, value in obj_props.items():
                if value:
                    prop_type = self.openbis_session.get_property_type(key)
                    prop_label = prop_type.label
                    obj_details_string += f"<p><b>{prop_label}:</b> {value}</p>"
            
            obj_details_string += "</div>"
            obj_details_html.value = obj_details_string
            self.details_vbox.children = [obj_details_html]
    
    def remove_molecule(self, b):
        molecules_accordion_children = list(self.parent_accordion.children)
        num_molecules = len(molecules_accordion_children)
        molecules_accordion_children.pop(self.object_index)
        
        for index, molecule in enumerate(molecules_accordion_children):
            if index >= self.object_index:
                molecule.object_index -= 1
                self.parent_accordion.set_title(molecule.object_index, molecule.title)

        self.parent_accordion.set_title(num_molecules - 1, "")
        self.parent_accordion.children = molecules_accordion_children

class ReacProdConceptWidget(ipw.VBox):
    def __init__(self, openbis_session, parent_accordion, object_index):
        super().__init__()
        self.openbis_session = openbis_session
        self.parent_accordion = parent_accordion
        self.object_index = object_index
        self.title = ""
        
        molecules_objects = utils.get_openbis_objects(
            self.openbis_session, 
            collection = "/MATERIALS/MOLECULES/PRODUCT_COLLECTION",
            type = "REACTION_PRODUCT_CONCEPT"
        )
        dropdown_list = [(obj.props["$name"], obj.permId) for obj in molecules_objects]
        dropdown_list.insert(0, ("Select a reaction product concept...", "-1"))
        self.dropdown = ipw.Dropdown(value = "-1",options = dropdown_list)
        self.details_vbox = ipw.VBox()
    
        self.remove_reacprod_concept_button = ipw.Button(
            description = 'Remove', disabled = False, 
            button_style = 'danger', tooltip = 'Remove reaction product concept', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.dropdown.observe(self.load_details, names = "value")
        self.remove_reacprod_concept_button.on_click(self.remove_reacprod_concept)
        self.children = [self.dropdown, self.details_vbox, self.remove_reacprod_concept_button]
    
    def load_details(self, change):
        obj_permid = self.dropdown.value
        if obj_permid == "-1":
            return
        else:
            obj = get_cached_object(self.openbis_session, obj_permid)
            obj_props = obj.props.all()
            obj_name = obj_props.get("$name", "")
            self.parent_accordion.set_title(self.object_index, obj_name)
            self.title = obj_name
            
            obj_details_html = ipw.HTML()
            obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
            for key, value in obj_props.items():
                if value:
                    prop_type = self.openbis_session.get_property_type(key)
                    prop_label = prop_type.label
                    obj_details_string += f"<p><b>{prop_label}:</b> {value}</p>"
            
            obj_details_string += "</div>"
            obj_details_html.value = obj_details_string
            self.details_vbox.children = [obj_details_html]
    
    def remove_reacprod_concept(self, b):
        reacprod_concepts_accordion_children = list(self.parent_accordion.children)
        num_reacprod_concepts = len(reacprod_concepts_accordion_children)
        reacprod_concepts_accordion_children.pop(self.object_index)
        
        for index, reacprod_concept in enumerate(reacprod_concepts_accordion_children):
            if index >= self.object_index:
                reacprod_concept.object_index -= 1
                self.parent_accordion.set_title(reacprod_concept.object_index, reacprod_concept.title)

        self.parent_accordion.set_title(num_reacprod_concepts - 1, "")
        self.parent_accordion.children = reacprod_concepts_accordion_children

# Sample preparation widgets
class CreateSampleWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_material_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select material</span>"
        )
        
        material_type_options = MATERIALS_LABELS
        material_type_options.insert(0, ("Select material type...", "-1"))
        
        self.material_type_dropdown = ipw.Dropdown(
            options = material_type_options,
            value = material_type_options[0][1]
        )
        
        self.material_details_vbox = ipw.VBox()
        
        self.sample_name_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Sample name</span>"
        )
        
        self.sample_name_textbox = ipw.Text(
            placeholder = "Write sample name..."
        )
        
        self.save_button = ipw.Button(
            description = '', 
            disabled = False, 
            button_style = '', 
            tooltip = 'Save', 
            icon = 'save', 
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        self.children = [
            self.select_material_title,
            self.material_type_dropdown,
            self.material_details_vbox,
            self.sample_name_title,
            self.sample_name_textbox,
            self.save_button
        ]
        
        self.material_type_dropdown.observe(
            self.load_material_type_widgets,
            names = "value"
        )
        
        self.save_button.on_click(self.save_sample)
    
    def load_material_type_widgets(self, change):
        if self.material_type_dropdown.value == "-1":
            self.material_details_vbox.children = []
            return
        else:
            material_options = [
                ("Select material...", "-1")
            ]
            
            material_dropdown = ipw.Dropdown(
                options = material_options,
                value = material_options[0][1]
            )
            
            sort_material_label = ipw.Label(
                value = "Sort by:", 
                layout=ipw.Layout(margin='0px', width='50px'),
                style = {'description_width': 'initial'}
            )
            
            name_checkbox = ipw.Checkbox(
                indent = False,
                layout=ipw.Layout(margin='2px', width='20px')
            )
            
            name_label = ipw.Label(
                value = "Name", 
                layout=ipw.Layout(margin='0px', width='50px'),
                style = {'description_width': 'initial'}
            )
            
            registration_date_checkbox = ipw.Checkbox(
                indent = False,
                layout=ipw.Layout(margin='2px', width='20px')
            )
            
            registration_date_label = ipw.Label(
                value = "Registration date",
                layout=ipw.Layout(margin='0px', width='110px'),
                style = {'description_width': 'initial'}
            )
        
            select_material_box = ipw.HBox(
                children = [
                    material_dropdown,
                    sort_material_label,
                    name_checkbox,
                    name_label,
                    registration_date_checkbox,
                    registration_date_label
                ]
            )
            
            material_details_html = ipw.HTML()
            
            self.material_details_vbox.children = [
                select_material_box,
                material_details_html
            ]
            
            material_type = self.material_type_dropdown.value
            material_objects = utils.get_openbis_objects(self.openbis_session, type = material_type)
            materials_objects_names_permids = [(obj.props["$name"], obj.permId) for obj in material_objects]
            material_options += materials_objects_names_permids
            material_dropdown.options = material_options
            
            def sort_material_dropdown(change):
                options = material_options[1:]
                
                df = pd.DataFrame(options, columns=["$name", "registration_date"])
                if name_checkbox.value and not registration_date_checkbox.value:
                    df = df.sort_values(by="$name", ascending=True)
                elif not name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by="registration_date", ascending=False)
                elif name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by=["$name", "registration_date"], ascending=[True, False])

                options = list(df.itertuples(index=False, name=None))
                options.insert(0, material_options[0])
                material_dropdown.options = options
            
            def load_material_details(change):
                obj_permid = material_dropdown.value
                if obj_permid == "-1":
                    return
                else:
                    obj = get_cached_object(self.openbis_session, obj_permid)
                    obj_props = obj.props.all()
                    obj_name = obj_props.get("$name", "")
                    obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
                    for key, value in obj_props.items():
                        if value:
                            prop_type = self.openbis_session.get_property_type(key)
                            prop_label = prop_type.label
                            obj_details_string += f"<p><b>{prop_label}:</b> {value}</p>"
                    
                    obj_details_string += "</div>"
                    
                    material_details_html.value = obj_details_string
                    
                    current_datetime = utils.get_current_datetime()
                    current_datetime_str = utils.convert_datetime_to_string(current_datetime)
                    self.sample_name_textbox.value = f"{current_datetime_str}_{obj_name}"
                    
            name_checkbox.observe(sort_material_dropdown, names = "value")
            registration_date_checkbox.observe(sort_material_dropdown, names = "value")
            material_dropdown.observe(load_material_details, names = "value")
    
    def save_sample(self, b):
        if self.material_type_dropdown.value == "-1":
            return
        else:
            if self.material_details_vbox.children[0].children[0].value == "-1":
                return
            else:
                material_object = self.material_details_vbox.children[0].children[0].value
                sample_name = self.sample_name_textbox.value
                sample_props = {
                    "$name": sample_name,
                    "exists": True
                }
                utils.create_openbis_object(
                    self.openbis_session,
                    type = "SAMPLE",
                    collection = "/MATERIALS/SAMPLES/SAMPLE_COLLECTION",
                    props = sample_props,
                    parents = [material_object]
                )
                
                display(utils.Javascript(data = "alert('Sample created successfully!')"))
                
                # Clear interface
                self.material_type_dropdown.value = "-1"
                self.sample_name_textbox.value = ""     

class RegisterProcessWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        self.sample_preparation_object = None
        
        self.select_experiment_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select experiment</span>"
        )
        
        self.select_sample_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select sample</span>"
        )
        
        self.sample_history_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Sample history</span>"
        )
        
        self.new_processes_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Register new steps</span>"
        )
        
        self.select_experiment_dropdown = SelectExperimentWidget(self.openbis_session)
        self.select_sample_dropdown = SelectSampleWidget(self.openbis_session)
        self.sample_history_vbox = SampleHistoryWidget(self.openbis_session)
        self.new_processes_accordion = ipw.Accordion()
        
        self.load_process_button = ipw.Button(
            description = 'Load process', 
            disabled = False, 
            button_style = 'success', 
            tooltip = 'Add process step', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.add_process_step_button = ipw.Button(
            description = 'Add process step', 
            disabled = False, 
            button_style = 'success', 
            tooltip = 'Add process step', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.process_buttons_hbox = ipw.HBox(
            children = [
                self.load_process_button,
                self.add_process_step_button
            ]
        )
        
        self.load_processes_vbox = ipw.VBox()
        self.processes_dropdown = ipw.Dropdown()
        
        self.save_button = ipw.Button(
            description = '', 
            disabled = False, 
            button_style = '', 
            tooltip = 'Save', 
            icon = 'save', 
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        self.children = [
            self.select_experiment_title,
            self.select_experiment_dropdown,
            self.select_sample_title,
            self.select_sample_dropdown,
            self.sample_history_title,
            self.sample_history_vbox,
            self.new_processes_title,
            self.load_processes_vbox,
            self.new_processes_accordion,
            self.process_buttons_hbox,
            self.save_button
        ]
        
        self.select_sample_dropdown.sample_dropdown.observe(
            self.load_sample_data,
            names = "value"
        )
        
        self.load_process_button.on_click(self.load_process)
        self.processes_dropdown.observe(self.load_process_settings, names = "value")
        self.add_process_step_button.on_click(self.add_process_step)
        self.save_button.on_click(self.save_processes)
    
    def load_sample_data(self, change):
        if self.select_sample_dropdown.sample_dropdown.value == "-1":
            self.sample_history_vbox.sample_history.children = []
            return
        
        sample_identifier = self.select_sample_dropdown.sample_dropdown.value
        sample_object = get_cached_object(self.openbis_session, sample_identifier)
        
        sample_object_parents = sample_object.parents
        most_recent_parent = None
        
        for parent_id in sample_object_parents:
            parent_object = get_cached_object(self.openbis_session, parent_id)
            
            parent_type = parent_object.type
            if parent_type == "PROCESS_STEP":
                if most_recent_parent:
                    if parent_object.registrationDate > most_recent_parent.registrationDate:
                        most_recent_parent = parent_object
                else:
                    most_recent_parent = parent_object
        
        if most_recent_parent:
            if most_recent_parent.experiment.permId != self.select_experiment_dropdown.experiment_dropdown.value:
                self.select_experiment_dropdown.experiment_dropdown.value = most_recent_parent.experiment.permId
                display(utils.Javascript(data = "alert('Experiment was changed!')"))
            
            for parent in most_recent_parent.parents:
                parent_object = get_cached_object(self.openbis_session, parent)
                
                if parent_object.type == "PREPARATION":
                    self.sample_preparation_object = parent_object
                    break

        # If sample was used in a measurement session, a new preparation should start
        sample_object_children = sample_object.children
        for child_id in sample_object_children:
            child_object = get_cached_object(self.openbis_session, child_id)
            if child_object.type == "MEASUREMENT_SESSION":
                self.sample_preparation_object = None
                break
        
        # Load sample history
        self.sample_history_vbox.load_sample_history(sample_object)
            
    def add_process_step(self, b):
        processes_accordion_children = list(self.new_processes_accordion.children)
        process_step_index = len(processes_accordion_children)
        new_process_step_widget = RegisterProcessStepWidget(self.openbis_session, self.new_processes_accordion, process_step_index)
        processes_accordion_children.append(new_process_step_widget)
        self.new_processes_accordion.children = processes_accordion_children
    
    def load_process(self, b):
        openbis_processes = utils.get_openbis_objects(self.openbis_session, type = "PROCESS")
        processes_options = [(obj.props["$name"], obj.permId) for obj in openbis_processes]
        processes_options.insert(0, ("Select a process...", "-1"))
        self.processes_dropdown.options = processes_options
        self.processes_dropdown.value = "-1"
        self.load_processes_vbox.children = [self.processes_dropdown]
    
    def load_process_settings(self, change):
        process_id = self.processes_dropdown.value
        if process_id == "-1":
            return
        else:
            process_object = get_cached_object(self.openbis_session, process_id)
            
            instrument_id = process_object.props["instrument"]
            process_steps_settings = process_object.props["process_steps_settings"]
            if process_steps_settings:
                process_steps_settings = json.loads(process_steps_settings)
                for step_settings in process_steps_settings:
                    step_settings["instrument"] = instrument_id
                    processes_accordion_children = list(self.new_processes_accordion.children)
                    process_step_index = len(processes_accordion_children)
                    new_process_step_widget = RegisterProcessStepWidget(
                        self.openbis_session, 
                        self.new_processes_accordion, 
                        process_step_index,
                        step_settings
                    )
                    processes_accordion_children.append(new_process_step_widget)
                    self.new_processes_accordion.children = processes_accordion_children

            self.load_processes_vbox.children = []
            

    def save_processes(self, b):
        processes_widgets = self.new_processes_accordion.children
        if processes_widgets:
            openbis_transaction_objects = []
            experiment_object = self.openbis_session.get_experiment(
                self.select_experiment_dropdown.experiment_dropdown.value
            )
            experiment_project_code = experiment_object.project.identifier
            
            current_sample_id = self.select_sample_dropdown.sample_dropdown.value
            current_sample = get_cached_object(self.openbis_session, current_sample_id)
            
            # Create preparation object when it does not exist
            if self.sample_preparation_object is None:
                self.sample_preparation_object = self.openbis_session.new_object(
                    type = "PREPARATION",
                    experiment = experiment_object.identifier,
                    props = {"$name": current_sample.props["$name"]}
                )
                self.sample_preparation_object.save()
            
            sample_preparation_id = self.sample_preparation_object.permId
            for process_widget in processes_widgets:
                # Reload sample preparation object to load children that was added in the cycle (e.g. process steps)
                self.sample_preparation_object = utils.get_openbis_object(self.openbis_session, sample_ident = sample_preparation_id)
                
                process_code = ""
                current_sample.props["exists"] = False
                current_sample_name = current_sample.props["$name"]
                current_sample.save()
                
                new_process_object = self.openbis_session.new_sample(
                    type = "PROCESS_STEP",
                    experiment = experiment_object.identifier
                )
                
                process_properties = {
                    "$name": process_widget.name_textbox.value,
                    "description": process_widget.description_textbox.value,
                    "instrument": process_widget.instrument_dropdown.value,
                    "comments": process_widget.comments_textarea.value,
                }
                
                actions_widgets = process_widget.actions_accordion.children
                observables_widgets = process_widget.observables_accordion.children
                actions = []
                observables = []
                
                actions_codes = []
                if actions_widgets:
                    for action_widget in actions_widgets:
                        action_properties = {}
                        action_type = action_widget.action_type_dropdown.value
                        duration_days = action_widget.duration_days_intbox.value
                        duration_hours = action_widget.duration_hours_intbox.value
                        duration_minutes = action_widget.duration_minutes_intbox.value
                        duration_seconds = action_widget.duration_seconds_intbox.value
                        
                        if action_type == "ANNEALING":
                            target_temperature = {
                                "value": action_widget.target_temperature_value_textbox.value,
                                "temperature_unit": action_widget.target_temperature_unit_dropdown.value,
                            }
                            action_properties["target_temperature"] = json.dumps(target_temperature)
                            
                        elif action_type == "COOLDOWN":
                            target_temperature = {
                                "value": action_widget.target_temperature_value_textbox.value,
                                "temperature_unit": action_widget.target_temperature_unit_dropdown.value,
                            }
                            action_properties["target_temperature"] = json.dumps(target_temperature)
                            action_properties["cryogen"] = action_widget.cryogen_textbox.value
                            
                        elif action_type == "DEPOSITION":
                            substrate_temperature = {
                                "value": action_widget.substrate_temperature_value_textbox.value,
                                "temperature_unit": action_widget.substrate_temperature_unit_dropdown.value,
                            }
                            substance = action_widget.substance_dropdown.value
                            if substance != "-1":
                                action_properties["substance"] = substance
                                
                            action_properties["substrate_temperature"] = json.dumps(substrate_temperature)
                            
                        elif action_type == "DOSING":
                            substrate_temperature = {
                                "value": action_widget.substrate_temperature_value_textbox.value,
                                "temperature_unit": action_widget.substrate_temperature_unit_dropdown.value,
                            }
                            pressure = {
                                "value": action_widget.pressure_value_textbox.value,
                                "pressure_unit": action_widget.pressure_unit_dropdown.value,
                            }
                            gas = action_widget.gas_dropdown.value
                            action_properties["dosing_gas"] = gas
                            action_properties["substrate_temperature"] = json.dumps(substrate_temperature)
                            action_properties["pressure"] = json.dumps(pressure)
                            
                        elif action_type == "SPUTTERING":
                            pressure = {
                                "value": action_widget.pressure_value_textbox.value,
                                "pressure_unit": action_widget.pressure_unit_dropdown.value,
                            }
                            current = {
                                "value": action_widget.current_value_textbox.value,
                                "current_unit": action_widget.current_unit_dropdown.value,
                            }
                            angle = {
                                "value": action_widget.angle_value_textbox.value,
                                "angle_unit": action_widget.angle_unit_dropdown.value,
                            }
                            substrate_temperature = {
                                "value": action_widget.substrate_temperature_value_textbox.value,
                                "temperature_unit": action_widget.substrate_temperature_unit_dropdown.value,
                            }
                            action_properties["pressure"] = json.dumps(pressure)
                            action_properties["current"] = json.dumps(current)
                            action_properties["angle"] = json.dumps(angle)
                            action_properties["substrate_temperature"] = json.dumps(substrate_temperature)
                            action_properties["sputter_ion"] = action_widget.sputter_ion_textbox.value
                        
                        component_permid = action_widget.component_dropdown.value
                        if component_permid != "-1":
                            component_object = get_cached_object(self.openbis_session, component_permid)
                            component_object_settings = component_object.props["actions_settings"]
                            
                            if component_object_settings:
                                component_settings = {}
                                component_object_settings = json.loads(component_object_settings)
                                action_component_settings = []
                                for action_settings in component_object_settings:
                                    if action_settings["action_type"] == action_type:
                                        action_component_settings = action_settings["action_properties"]
                                        break
                                    
                                for setting in action_component_settings:
                                    if setting == "target_temperature":
                                        component_settings["target_temperature"] = {
                                            "value": action_widget.target_temperature_value_comp_textbox.value,
                                            "temperature_unit": action_widget.target_temperature_unit_comp_dropdown.value,
                                        }
                                    elif setting == "bias_voltage":
                                        component_settings["bias_voltage"] = {
                                            "value": action_widget.bias_voltage_value_textbox.value,
                                            "voltage_unit": action_widget.bias_voltage_unit_dropdown.value,
                                        }
                                    elif setting == "discharge_voltage":
                                        component_settings["discharge_voltage"] = {
                                            "value": action_widget.discharge_voltage_value_textbox.value,
                                            "voltage_unit": action_widget.discharge_voltage_unit_dropdown.value,
                                        }
                                    elif setting == "discharge_current":
                                        component_settings["discharge_current"] = {
                                            "value": action_widget.discharge_current_value_textbox.value,
                                            "discharge_current_unit": action_widget.discharge_current_unit_dropdown.value,
                                        }
                                    elif setting == "evaporator_p_value":
                                        component_settings["evaporator_p_value"] = float(action_widget.evaporator_p_value_textbox.value)
                                        
                                    elif setting == "evaporator_i_value":
                                        component_settings["evaporator_i_value"] = float(action_widget.evaporator_i_value_textbox.value)
                                        
                                    elif setting == "ep_percentage":
                                        component_settings["ep_percentage"] = float(action_widget.ep_percentage_textbox.value)
                                
                                # Update current component settings
                                for setting_key, setting_value in component_settings.items():
                                    if isinstance(setting_value, dict):
                                        component_object.props[setting_key] = json.dumps(setting_value)
                                    else:
                                        component_object.props[setting_key] = setting_value
                                component_object.save()
                                
                                action_properties["component_settings"] = json.dumps(component_settings)

                            action_properties["$name"] = action_widget.name_textbox.value
                            action_properties["description"] = action_widget.description_textbox.value
                            action_properties["duration"] = f"{duration_days} days {duration_hours:02}:{duration_minutes:02}:{duration_seconds:02}"
                            action_properties["component"] = component_permid
                            action_properties["comments"] = action_widget.comments_textarea.value
                            
                            action_collection_code = "ACTIONS_COLLECTION"
                            openbis_experiments = self.openbis_session.get_experiments(code = action_collection_code, project = experiment_project_code)
                            
                            if openbis_experiments.df.empty:
                                actions_collection_object = self.openbis_session.new_collection(
                                    type = "COLLECTION",
                                    code = action_collection_code,
                                    project = experiment_project_code,
                                    props = {"$name": "Actions"}
                                )
                                actions_collection_object.save()
                            
                            new_action_object = self.openbis_session.new_sample(
                                type = action_type,
                                experiment = f"{experiment_project_code}/{action_collection_code}",
                                props = action_properties
                            )
                            new_action_object.save()
                            
                            # Append action to list of actions
                            actions.append(new_action_object.permId)
                            
                            # Get action code
                            actions_codes.append(ACTIONS_CODES[action_type])
                
                if observables_widgets:
                    for observable_widget in observables_widgets:
                        observable_properties = {}
                        observable_type = observable_widget.observable_type_dropdown.value
                        
                        component_permid = observable_widget.component_dropdown.value
                        if component_permid != "-1":
                            component_object = get_cached_object(self.openbis_session, component_permid)
                            component_object_settings = component_object.props["observables_settings"]
                            
                            if component_object_settings:
                                component_settings = {}
                                component_object_settings = json.loads(component_object_settings)
                                observable_component_settings = []
                                for observable_settings in component_object_settings:
                                    if observable_settings["observable_type"] == observable_type:
                                        observable_component_settings = observable_settings["observable_properties"]
                                        break
                                    
                                for setting in observable_component_settings:
                                    if setting == "density":
                                        component_settings["density"] = {
                                            "value": observable_widget.density_value_textbox.value,
                                            "density_unit": observable_widget.density_unit_dropdown.value,
                                        }
                                        
                                    elif setting == "filament":
                                        component_settings["filament"] = observable_widget.filament_textbox.value
                                        
                                    elif setting == "filament_current":
                                        component_settings["filament_current"] = {
                                            "value": observable_widget.filament_current_value_textbox.value,
                                            "current_unit": observable_widget.filament_current_unit_dropdown.value,
                                        }
                                
                                # Update current component settings
                                for setting_key, setting_value in component_settings.items():
                                    if isinstance(setting_value, dict):
                                        component_object.props[setting_key] = json.dumps(setting_value)
                                    else:
                                        component_object.props[setting_key] = setting_value
                                component_object.save()
                                
                                observable_properties["component_settings"] = json.dumps(component_settings)
                        
                            observable_properties["$name"] = observable_widget.name_textbox.value
                            observable_properties["description"] = observable_widget.description_textbox.value
                            observable_properties["channel_name"] = observable_widget.ch_name_textbox.value
                            observable_properties["component"] = component_permid
                            observable_properties["comments"] = observable_widget.comments_textarea.value
                            
                            observable_collection_code = "OBSERVABLES_COLLECTION"
                            openbis_experiments = self.openbis_session.get_experiments(code = observable_collection_code, project = experiment_project_code)
                            
                            if openbis_experiments.df.empty:
                                observables_collection_object = self.openbis_session.new_collection(
                                    type = "COLLECTION",
                                    code = observable_collection_code,
                                    project = experiment_project_code,
                                    props = {"$name": "Observables"}
                                )
                                observables_collection_object.save()
                            
                            new_observable_object = self.openbis_session.new_sample(
                                type = observable_type,
                                experiment = f"{experiment_project_code}/{observable_collection_code}",
                                props = observable_properties
                            )
                            new_observable_object.save()
                            
                            openbis_transaction_objects.append(new_observable_object)
                            
                            # Append observable to list of observables
                            observables.append(new_observable_object.permId)
                        
                process_properties["actions"] = actions
                process_properties["observables"] = observables
                
                if actions_codes:
                    # Compute process code based on the selected actions
                    counts = Counter(actions_codes)
                    unique_codes = list(counts.keys())
                    num_repeats = counts[unique_codes[0]] if len(counts) == 1 or all(v == next(iter(counts.values())) for v in counts.values()) else 1

                    if len(actions_codes) == 1:
                        process_code = actions_codes[0]
                    elif num_repeats > 1 and all(v == num_repeats for v in counts.values()):
                        process_code = f"({':'.join(unique_codes)}){num_repeats}"
                    else:
                        process_code = f"[{':'.join(actions_codes)}]"
                
                new_sample_name = f"{current_sample_name}:{process_code}"
                self.sample_preparation_object.props["$name"] = f"PREP_{new_sample_name}"
                self.sample_preparation_object.save()
                
                new_process_object.props = process_properties
                new_process_object.parents = [self.sample_preparation_object, current_sample]
                new_process_object.save()
                
                new_sample = self.openbis_session.new_sample(
                    type = "SAMPLE",
                    experiment = "/MATERIALS/SAMPLES/SAMPLE_COLLECTION",
                    parents = [new_process_object],
                    props = {"$name": new_sample_name, "exists": True}
                )
                new_sample.save()
                
                # After a process step, the current sample is now the new one
                current_sample = new_sample
            
            # Refresh sample dropdown and sample history
            self.select_sample_dropdown.load_samples()
            self.select_sample_dropdown.sample_dropdown.value = new_sample.permId
            
            # Reset new processes accordion
            processes_accordion_children = list(self.new_processes_accordion.children)
            for index, process_step in enumerate(processes_accordion_children):
                self.new_processes_accordion.set_title(index, "")
            
            self.new_processes_accordion.children = []
                
class SelectExperimentWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.create_experiment_button = ipw.Button(
            tooltip = 'Add', 
            icon = 'plus', 
            layout = ipw.Layout(width = '50px', height = '25px')
        )
        
        self.experiment_label = ipw.Label(
            value = "Experiment"
        )
        
        self.experiment_options = self.load_experiments()
        self.experiment_dropdown = ipw.Dropdown(
            options = self.experiment_options,
            layout=ipw.Layout(width='500px'),
            value = "-1"
        )
        
        
        self.sort_experiment_label = ipw.Label(
            value = "Sort by:"
        )
        
        self.sort_name_label = ipw.Label(
            value = "Name", 
            layout=ipw.Layout(margin='2px', width='50px'),
            style = {'description_width': 'initial'}
        )
        
        self.sort_name_checkbox = ipw.Checkbox(
            indent = False,
            layout=ipw.Layout(margin='2px', width='20px')
        )
        
        self.sort_registration_date_label = ipw.Label(
            value = "Registration date", 
            layout=ipw.Layout(margin='2px', width='110px'),
            style = {'description_width': 'initial'}
        )
        
        self.sort_registration_date_checkbox = ipw.Checkbox(
            indent = False,
            layout=ipw.Layout(margin='2px', width='20px')
        )
        
        self.sort_experiment_widgets = ipw.HBox(
            children = [
                self.sort_experiment_label,
                self.sort_name_checkbox,
                self.sort_name_label,
                self.sort_registration_date_checkbox,
                self.sort_registration_date_label
            ]
        )
        
        self.experiment_dropdown_boxes = ipw.HBox(
            children = [
                self.experiment_label,
                self.experiment_dropdown,
                self.create_experiment_button
            ]
        )
        
        self.create_new_experiment_widgets = ipw.VBox()
        
        self.children = [
            self.experiment_dropdown_boxes,
            self.sort_experiment_widgets,
            self.create_new_experiment_widgets
        ]
        
        self.create_experiment_button.on_click(self.create_new_experiment)
        self.sort_name_checkbox.observe(self.sort_experiment_dropdown)
        self.sort_registration_date_checkbox.observe(self.sort_experiment_dropdown)
    
    def load_experiments(self):
        experiments = utils.get_openbis_collections(
            self.openbis_session,
            type = "EXPERIMENT"
        )
        experiment_options = []
        for exp in experiments:
            if "$name" in exp.props.all():
                exp_option = (f"{exp.props['$name']} from Project {exp.project.code} and Space {exp.project.space}", exp.permId)
            else:
                exp_option = (f"{exp.code} from Project {exp.project.code} and Space {exp.project.space}", exp.permId)
            experiment_options.append(exp_option)
        experiment_options.insert(0, ("Select experiment...", "-1"))
        return experiment_options

    def sort_experiment_dropdown(self, change):
        options = self.experiment_options[1:]
        
        df = pd.DataFrame(options, columns=["$name", "registration_date"])
        if self.sort_name_checkbox.value and not self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="$name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="registration_date", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by=["$name", "registration_date"], ascending=[True, False])

        options = list(df.itertuples(index=False, name=None))
        options.insert(0, self.experiment_options[0])
        self.experiment_dropdown.options = options
    
    def create_new_experiment(self, b):
        new_experiment_name_label = ipw.Label(
            value = "Name"
        )
        
        new_experiment_name_textbox = ipw.Text(
            placeholder = "Write experiment name..."
        )
        
        new_experiment_name_hbox = ipw.HBox(
            children = [
                new_experiment_name_label,
                new_experiment_name_textbox
            ]
        )
        
        project_label = ipw.Label(
            value = "Project"
        )
        
        projects = utils.get_openbis_projects(self.openbis_session)
        project_dropdown_options = [(f"{proj.code} from Space {proj.space}", proj.permId) for proj in projects]
        project_dropdown_options.insert(0, ("Select project...", "-1"))
        
        project_dropdown = ipw.Dropdown(
            options = project_dropdown_options,
            value = "-1"
        )
        
        project_hbox = ipw.HBox(
            children = [
                project_label,
                project_dropdown
            ]
        )
        
        sort_project_label = ipw.Label(
            value = "Sort by:"
        )
        
        sort_project_name_label = ipw.Label(
            value = "Name", 
            layout=ipw.Layout(margin='2px', width='50px'),
            style = {'description_width': 'initial'}
        )
        
        sort_project_name_checkbox = ipw.Checkbox(
            indent = False,
            layout=ipw.Layout(margin='2px', width='20px')
        )
        
        sort_project_registration_date_label = ipw.Label(
            value = "Registration date", 
            layout=ipw.Layout(margin='2px', width='110px'),
            style = {'description_width': 'initial'}
        )
        
        sort_project_registration_date_checkbox = ipw.Checkbox(
            indent = False,
            layout=ipw.Layout(margin='2px', width='20px')
        )
        
        sort_project_hbox = ipw.HBox(
            children = [
                sort_project_label,
                sort_project_name_checkbox,
                sort_project_name_label,
                sort_project_registration_date_checkbox,
                sort_project_registration_date_label
            ]
        )
        
        save_button = ipw.Button(
            description = '', 
            disabled = False, 
            button_style = '', 
            tooltip = 'Save', 
            icon = 'save', 
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        cancel_button = ipw.Button(
            description = '', 
            disabled = False, 
            button_style = '', 
            tooltip = 'Cancel', 
            icon = 'times', 
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        buttons_hbox = ipw.HBox(
            children = [
                save_button,
                cancel_button
            ]
        )
        
        def save_new_experiment(b):
            experiment_props = {
                "$name": new_experiment_name_textbox.value,
                "$default_collection_view": "LIST_VIEW"
            }
            
            project_id = project_dropdown.value
            
            if project_id == "-1":
                display(utils.Javascript(data = "alert('Select a project.')"))
                return
            else:
                try:
                    utils.create_experiment_in_openbis(
                        self.openbis_session, 
                        project_dropdown.value, 
                        experiment_props
                    )
                    self.create_new_experiment_widgets.children = []
                    self.experiment_dropdown.options = self.load_experiments()
                    display(utils.Javascript(data = "alert('Experiment successfully created!')"))
                except ValueError as e:
                    display(utils.Javascript(data = "alert('Error! Check if experiment already exists (either in ELN or Trash).')"))
        
        def cancel_new_experiment(b):
            self.create_new_experiment_widgets.children = []
            
        def sort_project_dropdown(change):
            options = project_dropdown_options[1:]
            
            df = pd.DataFrame(options, columns=["$name", "registration_date"])
            if sort_project_name_checkbox.value and not sort_project_registration_date_checkbox.value:
                df = df.sort_values(by="$name", ascending=True)
            elif not sort_project_name_checkbox.value and sort_project_registration_date_checkbox.value:
                df = df.sort_values(by="registration_date", ascending=False)
            elif sort_project_name_checkbox.value and sort_project_registration_date_checkbox.value:
                df = df.sort_values(by=["$name", "registration_date"], ascending=[True, False])

            options = list(df.itertuples(index=False, name=None))
            options.insert(0, project_dropdown_options[0])
            project_dropdown.options = options
        
        save_button.on_click(save_new_experiment)
        cancel_button.on_click(cancel_new_experiment)
        sort_project_name_checkbox.observe(sort_project_dropdown)
        sort_project_registration_date_checkbox.observe(sort_project_dropdown)
        
        self.create_new_experiment_widgets.children = [
            new_experiment_name_hbox,
            project_hbox,
            sort_project_hbox,
            buttons_hbox
        ]
        
class SelectSampleWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.sample_label = ipw.Label(
            value = "Sample"
        )
        
        self.sample_dropdown = ipw.Dropdown()
        self.load_samples()
        
        
        self.sort_sample_label = ipw.Label(
            value = "Sort by:"
        )
        
        self.sort_name_label = ipw.Label(
            value = "Name", 
            layout=ipw.Layout(margin='2px', width='50px'),
            style = {'description_width': 'initial'}
        )
        
        self.sort_name_checkbox = ipw.Checkbox(
            indent = False,
            layout=ipw.Layout(margin='2px', width='20px')
        )
        
        self.sort_registration_date_label = ipw.Label(
            value = "Registration date", 
            layout=ipw.Layout(margin='2px', width='110px'),
            style = {'description_width': 'initial'}
        )
        
        self.sort_registration_date_checkbox = ipw.Checkbox(
            indent = False,
            layout=ipw.Layout(margin='2px', width='20px')
        )
        
        self.sort_sample_hbox = ipw.HBox(
            children = [
                self.sort_sample_label,
                self.sort_name_checkbox,
                self.sort_name_label,
                self.sort_registration_date_checkbox,
                self.sort_registration_date_label
            ]
        )
        
        self.sample_dropdown_hbox = ipw.HBox(
            children = [
                self.sample_label,
                self.sample_dropdown,
            ]
        )
        
        self.children = [
            self.sample_dropdown_hbox,
            self.sort_sample_hbox
        ]
    
    def load_samples(self):
        samples = utils.get_openbis_objects(
            self.openbis_session,
            type = "SAMPLE"
        )
        sample_options = [(f"{obj.props['$name']}", obj.permId) for obj in samples if obj.props["exists"] == "true"]
        sample_options.insert(0, ("Select sample...", "-1"))
        self.sample_dropdown.options =  sample_options
        self.sample_dropdown.value = "-1"
        
        # TODO: Apply functions to the sort buttons

class SampleHistoryWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        self.sample_history_objects = {}
        self.sample_history = ipw.Accordion()
        
        self.children = [
            self.sample_history
        ]
    
    def load_sample_history(self, sample_object):
        process_steps = []
        sample_parents = sample_object.parents

        while sample_parents:
            next_parents = []
            for parent in sample_parents:
                parent_code = parent.split("/")[-1]
                if "PRST" in parent_code or "SAMP" in parent_code:
                    parent_object = get_cached_object(self.openbis_session, parent)
                    next_parents.extend(parent_object.parents)
                    
                    if "PRST" in parent_code:
                        process_steps.append(parent_object)
                        # Save object information
                        OPENBIS_SAMPLES_CACHE[parent] = parent_object
            sample_parents = next_parents
            
        sample_history_children = []
        for i, process_step in enumerate(process_steps):
            process_step_widget = ProcessStepHistoryWidget(self.openbis_session, process_step)
            sample_history_children.append(process_step_widget)
            process_step_title = process_step_widget.name_html.value + " (" + process_step_widget.registration_date + ")"
            self.sample_history.set_title(i, process_step_title)
        
        self.sample_history.children = sample_history_children

class ProcessStepHistoryWidget(ipw.VBox):
    def __init__(self, openbis_session, openbis_object):
        super().__init__()
        self.openbis_session = openbis_session
        self.openbis_object = openbis_object
        self.registration_date = None
        
        self.name_label = ipw.Label(value = "Name:")
        self.name_html = ipw.HTML()
        self.name_hbox = ipw.HBox(children = [self.name_label, self.name_html])
        
        self.description_label = ipw.Label(value = "Description:")
        self.description_html = ipw.HTML()
        self.description_hbox = ipw.HBox(children = [self.description_label, self.description_html])
        
        self.comments_label = ipw.Label(value = "Comments:")
        self.comments_html = ipw.HTML()
        self.comments_hbox = ipw.HBox(children = [self.comments_label, self.comments_html])
        
        self.instrument_label = ipw.Label(value = "Instrument:")
        self.instrument_html = ipw.HTML()
        self.instrument_hbox = ipw.HBox(children = [self.instrument_label, self.instrument_html])
        
        self.actions_label = ipw.Label(value = "Actions:")
        self.actions_accordion = ipw.Accordion()
        self.actions_vbox = ipw.VBox(children = [self.actions_label, self.actions_accordion])
        
        self.observables_label = ipw.Label(value = "Observables:")
        self.observables_accordion = ipw.Accordion()
        self.observables_vbox = ipw.VBox(children = [self.observables_label, self.observables_accordion])
        
        self.load_process_step_data()
        
        self.children = [
            self.name_hbox,
            self.description_hbox,
            self.comments_hbox,
            self.instrument_hbox,
            self.actions_vbox,
            self.observables_vbox
        ]
    
    def load_process_step_data(self):
        openbis_object_props = self.openbis_object.props.all()
        if openbis_object_props["$name"]:
            self.name_html.value = openbis_object_props["$name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
            
        self.registration_date = self.openbis_object.registrationDate
        
        instrument_id = openbis_object_props["instrument"]
        
        if instrument_id:
            instrument_object = get_cached_object(self.openbis_session, instrument_id)
            self.instrument_html.value = instrument_object.props["$name"]
        
        self.load_actions()
        self.load_observables()
    
    def load_actions(self):
        actions_ids = self.openbis_object.props["actions"]
        if actions_ids:
            actions_accordion_children = []
            for i, act_id in enumerate(actions_ids):
                act_object = get_cached_object(self.openbis_session, act_id)
                act_widget = ActionHistoryWidget(self.openbis_session, act_object)
                actions_accordion_children.append(act_widget)
                act_title = act_widget.name_html.value
                self.actions_accordion.set_title(i, act_title)
            
            self.actions_accordion.children = actions_accordion_children

    def load_observables(self):
        observables_ids = self.openbis_object.props["observables"]
        if observables_ids:
            observables_accordion_children = []
            for i, obs_id in enumerate(observables_ids):
                obs_object = get_cached_object(self.openbis_session, obs_id)
                obs_widget = ObservableHistoryWidget(self.openbis_session, obs_object)
                observables_accordion_children.append(obs_widget)
                obs_title = obs_widget.name_html.value
                self.observables_accordion.set_title(i, obs_title)
            
            self.observables_accordion.children = observables_accordion_children

class ActionHistoryWidget(ipw.VBox):
    def __init__(self, openbis_session, openbis_object):
        super().__init__()
        self.openbis_session = openbis_session
        self.openbis_object = openbis_object
        self.object_type = self.openbis_object.type
        
        self.name_label = ipw.Label(value = "Name:")
        self.name_html = ipw.HTML()
        self.name_hbox = ipw.HBox(children = [self.name_label, self.name_html])
        
        self.duration_label = ipw.Label(value = "Duration:")
        self.duration_html = ipw.HTML()
        self.duration_hbox = ipw.HBox(children = [self.duration_label, self.duration_html])
        
        self.description_label = ipw.Label(value = "Description:")
        self.description_html = ipw.HTML()
        self.description_hbox = ipw.HBox(children = [self.description_label, self.description_html])
        
        self.component_label = ipw.Label(value = "Component:")
        self.component_html = ipw.HTML()
        self.component_hbox = ipw.HBox(children = [self.component_label, self.component_html])
        
        self.component_settings_label = ipw.Label(value = "Component settings:")
        self.component_settings_html = ipw.HTML()
        self.component_settings_hbox = ipw.HBox(children = [self.component_settings_label, self.component_settings_html])
        
        self.comments_label = ipw.Label(value = "Comments:")
        self.comments_html = ipw.HTML()
        self.comments_hbox = ipw.HBox(children = [self.comments_label, self.comments_html])
        
        widget_children = [
            self.name_hbox,
            self.description_hbox,
            self.duration_hbox
        ]
        
        if self.object_type == "ANNEALING":
            self.target_temperature_label = ipw.Label(value = "Name:")
            self.target_temperature_html = ipw.HTML()
            self.target_temperature_hbox = ipw.HBox(children = [self.target_temperature_label, self.target_temperature_html])
            widget_children.append(self.target_temperature_hbox)
            
        elif self.object_type == "COOLDOWN":
            self.target_temperature_label = ipw.Label(value = "Target temperature:")
            self.target_temperature_html = ipw.HTML()
            self.target_temperature_hbox = ipw.HBox(children = [self.target_temperature_label, self.target_temperature_html])
            
            self.cryogen_label = ipw.Label(value = "Cryogen:")
            self.cryogen_html = ipw.HTML()
            self.cryogen_hbox = ipw.HBox(children = [self.cryogen_label, self.cryogen_html])
            
            widget_children.append(self.cryogen_hbox)
            widget_children.append(self.target_temperature_hbox)
        
        elif self.object_type == "DEPOSITION":
            self.substrate_temperature_label = ipw.Label(value = "Substrate temperature:")
            self.substrate_temperature_html = ipw.HTML()
            self.substrate_temperature_hbox = ipw.HBox(children = [self.substrate_temperature_label, self.substrate_temperature_html])
            
            self.substance_label = ipw.Label(value = "Substance:")
            self.substance_html = ipw.HTML()
            self.substance_hbox = ipw.HBox(children = [self.substance_label, self.substance_html])
            
            widget_children.append(self.substance_hbox)
            widget_children.append(self.substrate_temperature_hbox)
        
        elif self.object_type == "DOSING":
            self.substrate_temperature_label = ipw.Label(value = "Substrate temperature:")
            self.substrate_temperature_html = ipw.HTML()
            self.substrate_temperature_hbox = ipw.HBox(children = [self.substrate_temperature_label, self.substrate_temperature_html])
            
            self.pressure_label = ipw.Label(value = "Pressure:")
            self.pressure_html = ipw.HTML()
            self.pressure_hbox = ipw.HBox(children = [self.pressure_label, self.pressure_html])
            
            self.dosing_gas_label = ipw.Label(value = "Dosing gas:")
            self.dosing_gas_html = ipw.HTML()
            self.dosing_gas_hbox = ipw.HBox(children = [self.dosing_gas_label, self.dosing_gas_html])
            
            widget_children.append(self.dosing_gas_hbox)
            widget_children.append(self.pressure_hbox)
            widget_children.append(self.substrate_temperature_hbox)

        elif self.object_type == "SPUTTERING":
            self.sputter_ion_label = ipw.Label(value = "Sputter ion:")
            self.sputter_ion_html = ipw.HTML()
            self.sputter_ion_hbox = ipw.HBox(children = [self.sputter_ion_label, self.sputter_ion_html])
            
            self.pressure_label = ipw.Label(value = "Pressure:")
            self.pressure_html = ipw.HTML()
            self.pressure_hbox = ipw.HBox(children = [self.pressure_label, self.pressure_html])
            
            self.current_label = ipw.Label(value = "Current:")
            self.current_html = ipw.HTML()
            self.current_hbox = ipw.HBox(children = [self.current_label, self.current_html])
            
            self.angle_label = ipw.Label(value = "Angle:")
            self.angle_html = ipw.HTML()
            self.angle_hbox = ipw.HBox(children = [self.angle_label, self.angle_html])
            
            self.substrate_temperature_label = ipw.Label(value = "Substrate temperature:")
            self.substrate_temperature_html = ipw.HTML()
            self.substrate_temperature_hbox = ipw.HBox(children = [self.substrate_temperature_label, self.substrate_temperature_html])
            
            widget_children.append(self.dosing_gas_hbox)
            widget_children.append(self.pressure_hbox)
            widget_children.append(self.substrate_temperature_hbox)
            
        widget_children.append(self.comments_hbox)
        widget_children.append(self.component_hbox)
        widget_children.append(self.component_settings_hbox)
        
        self.load_action_data()
        self.children = widget_children
    
    def load_action_data(self):
        openbis_object_props = self.openbis_object.props.all()
        if openbis_object_props["$name"]:
            self.name_html.value = openbis_object_props["$name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["duration"]:
            self.duration_html.value = openbis_object_props["duration"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
        
        if "target_temperature" in openbis_object_props:
            target_temperature = openbis_object_props["target_temperature"]
            if target_temperature:
                self.target_temperature_html.value = stringify_quantity_value(target_temperature, "temperature_unit")
        
        if "cryogen" in openbis_object_props:
            if openbis_object_props["cryogen"]:
                self.cryogen_html.value = openbis_object_props["cryogen"]
        
        if "substrate_temperature" in openbis_object_props:
            substrate_temperature = openbis_object_props["substrate_temperature"]
            if substrate_temperature:
                self.substrate_temperature_html.value = stringify_quantity_value(substrate_temperature, "temperature_unit")
        
        if "dosing_gas" in openbis_object_props:
            if openbis_object_props["dosing_gas"]:
                self.dosing_gas_html.value = openbis_object_props["dosing_gas"]
        
        if "pressure" in openbis_object_props:
            pressure = openbis_object_props["pressure"]
            if pressure:
                self.pressure_html.value = stringify_quantity_value(pressure, "pressure_unit")
        
        if "current" in openbis_object_props:
            current = openbis_object_props["current"]
            if current:
                self.current_html.value = stringify_quantity_value(current, "current_unit")
        
        if "angle" in openbis_object_props:
            angle = openbis_object_props["angle"]
            if angle:
                self.angle_html.value = stringify_quantity_value(angle, "angle_unit")

        if "substance" in openbis_object_props:
            substance_id = openbis_object_props["substance"]
            if substance_id:
                substance_object = get_cached_object(self.openbis_session, substance_id)
                self.substance_html.value = substance_object.props["$name"]
        
        component_id = openbis_object_props["component"]
        if component_id:
            component_object = get_cached_object(self.openbis_session, component_id)
            self.component_html.value = component_object.props["$name"]
            
        if openbis_object_props["component_settings"]:
            component_settings = json.loads(openbis_object_props["component_settings"])
            component_settings_string = ""
            for prop_key, prop_value in component_settings.items():
                prop_type = self.openbis_session.get_property_type(prop_key)
                prop_label = prop_type.label
                
                # Convert quantity values from json objects to strings
                if prop_value and isinstance(prop_value, dict):
                    unit_type = ""
                    for key in prop_value.keys():
                        if "_unit" in key:
                            unit_type = key
                            break
                    prop_value = stringify_quantity_value(prop_value, unit_type)
                    
                component_settings_string += f"<p>&bull; {prop_label}: {prop_value}</p>"
            
            self.component_settings_html.value = component_settings_string
        
class ObservableHistoryWidget(ipw.VBox):
    def __init__(self, openbis_session, openbis_object):
        super().__init__()
        self.openbis_session = openbis_session
        self.openbis_object = openbis_object
        
        self.name_label = ipw.Label(value = "Name:")
        self.name_html = ipw.HTML()
        self.name_hbox = ipw.HBox(children = [self.name_label, self.name_html])
        
        self.description_label = ipw.Label(value = "Description:")
        self.description_html = ipw.HTML()
        self.description_hbox = ipw.HBox(children = [self.description_label, self.description_html])
        
        self.ch_name_label = ipw.Label(value = "Channel name:")
        self.ch_name_html = ipw.HTML()
        self.ch_name_hbox = ipw.HBox(children = [self.ch_name_label, self.ch_name_html])
        
        self.component_label = ipw.Label(value = "Component:")
        self.component_html = ipw.HTML()
        self.component_hbox = ipw.HBox(children = [self.component_label, self.component_html])
        
        self.component_settings_label = ipw.Label(value = "Component settings:")
        self.component_settings_html = ipw.HTML()
        self.component_settings_hbox = ipw.HBox(children = [self.component_settings_label, self.component_settings_html])
        
        self.comments_label = ipw.Label(value = "Comments:")
        self.comments_html = ipw.HTML()
        self.comments_hbox = ipw.HBox(children = [self.comments_label, self.comments_html])
        
        self.load_observable_data()
        
        self.children = [
            self.name_hbox,
            self.description_hbox,
            self.comments_hbox,
            self.ch_name_hbox,
            self.component_hbox,
            self.component_settings_hbox
        ]
    
    def load_observable_data(self):
        openbis_object_props = self.openbis_object.props.all()
        if openbis_object_props["$name"]:
            self.name_html.value = openbis_object_props["$name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["channel_name"]:
            self.ch_name_html.value = openbis_object_props["channel_name"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
            
        component_id = openbis_object_props["component"]
            
        if component_id:
            component_object = get_cached_object(self.openbis_session, component_id)
            self.component_html.value = component_object.props["$name"]
            
        if openbis_object_props["component_settings"]:
            component_settings = json.loads(openbis_object_props["component_settings"])
            component_settings_string = ""
            for prop_key, prop_value in component_settings.items():
                prop_type = self.openbis_session.get_property_type(prop_key)
                prop_label = prop_type.label
                
                # Convert quantity values from json objects to strings
                if prop_value and isinstance(prop_value, dict):
                    unit_type = ""
                    for key in prop_value.keys():
                        if "_unit" in key:
                            unit_type = key
                            break
                    prop_value = stringify_quantity_value(prop_value, unit_type)
                
                component_settings_string += f"<p>&bull; {prop_label}: {prop_value}</p>"
            
            self.component_settings_html.value = component_settings_string

class RegisterProcessStepWidget(ipw.VBox):
    def __init__(self, openbis_session, processes_accordion, process_step_index, step_settings = None):
        super().__init__()
        self.openbis_session = openbis_session
        self.processes_accordion = processes_accordion
        self.process_step_index = process_step_index
        
        self.name_label = ipw.Label(value = "Name")
        self.name_textbox = ipw.Text()
        self.name_hbox = ipw.HBox(children = [self.name_label, self.name_textbox])
        
        self.description_label = ipw.Label(value = "Description")
        self.description_textbox = ipw.Text()
        self.description_hbox = ipw.HBox(children = [self.description_label, self.description_textbox])
        
        self.instrument_label = ipw.Label(value = "Instrument")
        instrument_objects = utils.get_openbis_objects(self.openbis_session, collection = "/EQUIPMENT/ILOG/INSTRUMENT_COLLECTION")
        instrument_options = [(obj.props["$name"], obj.permId) for obj in instrument_objects]
        instrument_options.insert(0, ("Select an instrument...", "-1"))
        self.instrument_dropdown = ipw.Dropdown(options = instrument_options, value = "-1")
        self.instrument_hbox = ipw.HBox(children = [self.instrument_label, self.instrument_dropdown])
        
        self.comments_label = ipw.Label(value = "Comments")
        self.comments_textarea = ipw.Textarea()
        self.comments_hbox = ipw.HBox(children = [self.comments_label, self.comments_textarea])
        
        self.actions_label = ipw.Label(value = "Actions")
        self.actions_accordion = ipw.Accordion()
        self.add_action_button = ipw.Button(
            description = 'Add action', 
            disabled = False, 
            button_style = 'success', 
            tooltip = 'Add action', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        self.actions_vbox = ipw.VBox(children = [self.actions_label, self.actions_accordion, self.add_action_button])
        
        self.observables_label = ipw.Label(value = "Observables")
        self.observables_accordion = ipw.Accordion()
        self.add_observable_button = ipw.Button(
            description = 'Add observable', 
            disabled = False, 
            button_style = 'success', 
            tooltip = 'Add observable', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        self.observables_vbox = ipw.VBox(children = [self.observables_label, self.observables_accordion, self.add_observable_button])
        
        self.remove_process_step_button = ipw.Button(
            description = 'Remove', 
            disabled = False, 
            button_style = 'danger', 
            tooltip = 'Remove process step', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.remove_process_step_button.on_click(self.remove_process_step)
        self.name_textbox.observe(self.change_process_step_title, names = "value")
        self.add_action_button.on_click(self.add_action)
        self.add_observable_button.on_click(self.add_observable)
        
        if step_settings:
            self.load_process_step(step_settings)
        
        self.children = [
            self.name_hbox,
            self.description_hbox,
            self.instrument_hbox,
            self.comments_hbox,
            self.actions_vbox,
            self.observables_vbox,
            self.remove_process_step_button
        ]
    
    def load_process_step(self, settings):
        self.name_textbox.value = settings.get("name", "")
        self.description_textbox.value = settings.get("description", "")
        self.instrument_dropdown.value = settings.get("instrument", "-1")
        self.comments_textarea.value = settings.get("comments", "")
        actions = settings.get("actions", {})
        observables = settings.get("observables", {})
        for action in actions:
            actions_accordion_children = list(self.actions_accordion.children)
            action_index = len(actions_accordion_children)
            actions[action]["type"] = action
            new_action_widget = RegisterActionWidget(
                self.openbis_session, 
                self.actions_accordion, 
                action_index, 
                self.instrument_dropdown.value,
                actions[action]
            )
            actions_accordion_children.append(new_action_widget)
            self.actions_accordion.children = actions_accordion_children
        
        for observable in observables:
            observables_accordion_children = list(self.observables_accordion.children)
            observable_index = len(observables_accordion_children)
            observables[observable]["type"] = observable
            new_observable_widget = RegisterObservableWidget(
                self.openbis_session, 
                self.observables_accordion, 
                observable_index, 
                self.instrument_dropdown.value,
                observables[observable]
            )
            observables_accordion_children.append(new_observable_widget)
            self.observables_accordion.children = observables_accordion_children
    
    def change_process_step_title(self, change):
        title = self.name_textbox.value
        self.processes_accordion.set_title(self.process_step_index, title)
    
    def remove_process_step(self, b):
        processes_accordion_children = list(self.processes_accordion.children)
        num_process_steps = len(processes_accordion_children)
        processes_accordion_children.pop(self.process_step_index)
        
        for index, process_step in enumerate(processes_accordion_children):
            if index >= self.process_step_index:
                process_step.process_step_index -= 1
                self.processes_accordion.set_title(process_step.process_step_index, process_step.name_textbox.value)

        self.processes_accordion.set_title(num_process_steps - 1, "")
        self.processes_accordion.children = processes_accordion_children
    
    def add_action(self, b):
        instrument_permid = self.instrument_dropdown.value
        if instrument_permid != "-1":
            actions_accordion_children = list(self.actions_accordion.children)
            action_index = len(actions_accordion_children)
            new_action_widget = RegisterActionWidget(self.openbis_session, self.actions_accordion, action_index, instrument_permid)
            actions_accordion_children.append(new_action_widget)
            self.actions_accordion.children = actions_accordion_children
    
    def add_observable(self, b):
        instrument_permid = self.instrument_dropdown.value
        if instrument_permid != "-1":
            observables_accordion_children = list(self.observables_accordion.children)
            observable_index = len(observables_accordion_children)
            new_observable_widget = RegisterObservableWidget(self.openbis_session, self.observables_accordion, observable_index, instrument_permid)
            observables_accordion_children.append(new_observable_widget)
            self.observables_accordion.children = observables_accordion_children

class RegisterActionWidget(ipw.VBox):
    def __init__(self, openbis_session, actions_accordion, action_index, instrument_permid, action_settings = None):
        super().__init__()
        self.openbis_session = openbis_session
        self.actions_accordion = actions_accordion
        self.action_index = action_index
        self.instrument_permid = instrument_permid
        
        self.action_type_label = ipw.Label(value = "Action type")
        action_type_options = ACTIONS_LABELS
        action_type_options.insert(0, ("Select an action type...", "-1"))
        
        self.action_type_dropdown = ipw.Dropdown(
            options = action_type_options,
            value = "-1"
        )
        
        self.action_type_hbox = ipw.HBox(children = [self.action_type_label, self.action_type_dropdown])
        self.action_properties_widgets = ipw.VBox()
        
        self.name_label = ipw.Label(value = "Name")
        self.name_textbox = ipw.Text()
        self.name_hbox = ipw.HBox(children = [self.name_label, self.name_textbox])
        
        self.description_label = ipw.Label(value = "Description")
        self.description_textbox = ipw.Text()
        self.description_hbox = ipw.HBox(children = [self.description_label, self.description_textbox])
        
        self.duration_label = ipw.Label(value = "Duration")
        self.duration_days_intbox = ipw.BoundedIntText(value = 0, min = 0, layout = ipw.Layout(width = "40px"))
        self.duration_days_label = ipw.Label("days")
        self.duration_hours_intbox = ipw.BoundedIntText(value = 0, max = 23, layout = ipw.Layout(width = "40px"))
        self.duration_hours_label = ipw.Label(":")
        self.duration_minutes_intbox = ipw.BoundedIntText(value = 0, max = 59, layout = ipw.Layout(width = "40px"))
        self.duration_minutes_label = ipw.Label(":")
        self.duration_seconds_intbox = ipw.BoundedIntText(value = 0, max = 59, layout = ipw.Layout(width = "40px"))
        self.duration_hbox = ipw.HBox(
            children = [
                self.duration_label, self.duration_days_intbox, self.duration_days_label,
                self.duration_hours_intbox, self.duration_hours_label,
                self.duration_minutes_intbox, self.duration_minutes_label,
                self.duration_seconds_intbox
            ]
        )
        
        self.substance_label = ipw.Label(value = "Substance")
        substances_list = utils.get_openbis_objects(
            self.openbis_session, 
            collection = "/MATERIALS/MOLECULES/PRECURSOR_COLLECTION",
            type = "SUBSTANCE"
        )
        substances_names_ids = []
        for obj in substances_list:
            obj_props = obj.props.all()
            name = obj_props["empa_number"] + obj_props["batch"]
            if "vial" in obj_props:
                if obj_props["vial"]:
                    name += obj_props["vial"]
            substances_names_ids.append((name,obj.permId))
                
        substance_options = substances_names_ids
        substance_options.insert(0, ("Select a substance...", "-1"))
        self.substance_dropdown = ipw.Dropdown(
            options = substance_options,
            value = "-1"
        )
        self.substance_hbox = ipw.HBox(children = [self.substance_label, self.substance_dropdown])
        
        self.gas_label = ipw.Label(value = "Dosing gas")
        gas_list = utils.get_openbis_objects(self.openbis_session, collection = "/MATERIALS/RAW_MATERIALS/CHEMICAL_COLLECTION")
        gas_options = [(obj.props["$name"], obj.permId) for obj in gas_list]
        gas_options.insert(0, ("Select a dosing gas...", "-1"))
        self.gas_dropdown = ipw.Dropdown(
            options = gas_options,
            value = "-1"
        )
        self.gas_hbox = ipw.HBox(children = [self.gas_label, self.gas_dropdown])
        
        self.target_temperature_label = ipw.Label("Target temperature")
        self.target_temperature_value_textbox = ipw.Text()
        self.target_temperature_unit_dropdown = ipw.Dropdown(options = ["K", "C"], value = "C")
        self.target_temperature_hbox = ipw.HBox(
            children = [self.target_temperature_label, self.target_temperature_value_textbox, self.target_temperature_unit_dropdown]
        )
        
        self.cryogen_label = ipw.Label("Cryogen")
        self.cryogen_textbox = ipw.Text()
        self.cryogen_hbox = ipw.HBox(children = [self.cryogen_label, self.cryogen_textbox])
        
        self.substrate_temperature_label = ipw.Label("Substrate temperature")
        self.substrate_temperature_value_textbox = ipw.Text()
        self.substrate_temperature_unit_dropdown = ipw.Dropdown(options = ["K", "C"], value = "C")
        self.substrate_temperature_hbox = ipw.HBox(
            children = [self.substrate_temperature_label, self.substrate_temperature_value_textbox, self.substrate_temperature_unit_dropdown]
        )
        
        self.pressure_label = ipw.Label("Pressure")
        self.pressure_value_textbox = ipw.Text()
        self.pressure_unit_dropdown = ipw.Dropdown(options = ["mBar", "Pa"], value = "mBar")
        self.pressure_hbox = ipw.HBox(
            children = [self.pressure_label, self.pressure_value_textbox, self.pressure_unit_dropdown]
        )
        
        self.sputter_ion_label = ipw.Label("Sputter Ion")
        self.sputter_ion_textbox = ipw.Text()
        self.sputter_ion_hbox = ipw.HBox(children = [self.sputter_ion_label, self.sputter_ion_textbox])
        
        self.current_label = ipw.Label("Current")
        self.current_value_textbox = ipw.Text()
        self.current_unit_dropdown = ipw.Dropdown(options = ["A"], value = "A")
        self.current_hbox = ipw.HBox(
            children = [self.current_label, self.current_value_textbox, self.current_unit_dropdown]
        )
        
        self.angle_label = ipw.Label("Angle")
        self.angle_value_textbox = ipw.Text()
        self.angle_unit_dropdown = ipw.Dropdown(options = ["deg"], value = "deg")
        self.angle_hbox = ipw.HBox(
            children = [self.angle_label, self.angle_value_textbox, self.angle_unit_dropdown]
        )
        
        self.instrument_type_components_dictionary = {
            "INSTRUMENT.STM": [
                "pumps", "gauges", "vacuum_chambers", "ports_valves",
                "preparation_tools", "analysers", "mechanical_components",
                "stm_components", "control_data_acquisition", "temperature_environment_control",
                "auxiliary_components", "tips_sensors", "accessories"
            ]
        }
        
        self.component_label = ipw.Label(value = "Component")
        self.component_dropdown = ipw.Dropdown()
        self.component_hbox = ipw.HBox(children = [self.component_label, self.component_dropdown])
        self.component_settings_label = ipw.Label(value = "Component settings:")
        self.component_settings_vbox = ipw.VBox()
        self.component_settings_hbox = ipw.HBox(children = [self.component_settings_label, self.component_settings_vbox])
        self.component_vbox = ipw.VBox(children = [self.component_hbox, self.component_settings_hbox])
        
        # BEGIN - Widgets for component properties
        self.bias_voltage_label = ipw.Label("Bias voltage")
        self.bias_voltage_value_textbox = ipw.Text()
        self.bias_voltage_unit_dropdown = ipw.Dropdown(options = ["V"], value = "V")
        self.bias_voltage_hbox = ipw.HBox(
            children = [self.bias_voltage_label, self.bias_voltage_value_textbox, self.bias_voltage_unit_dropdown]
        )
        
        self.discharge_voltage_label = ipw.Label("Discharge voltage")
        self.discharge_voltage_value_textbox = ipw.Text()
        self.discharge_voltage_unit_dropdown = ipw.Dropdown(options = ["V", "kV"], value = "V")
        self.discharge_voltage_hbox = ipw.HBox(
            children = [self.discharge_voltage_label, self.discharge_voltage_value_textbox, self.discharge_voltage_unit_dropdown]
        )
        
        self.discharge_current_label = ipw.Label("Discharge current")
        self.discharge_current_value_textbox = ipw.Text()
        self.discharge_current_unit_dropdown = ipw.Dropdown(options = ["A"], value = "A")
        self.discharge_current_hbox = ipw.HBox(
            children = [self.discharge_current_label, self.discharge_current_value_textbox, self.discharge_current_unit_dropdown]
        )
        
        # Target temperature set in the component
        self.target_temperature_comp_label = ipw.Label("Target temperature")
        self.target_temperature_value_comp_textbox = ipw.Text()
        self.target_temperature_unit_comp_dropdown = ipw.Dropdown(options = ["K", "C"], value = "C")
        self.target_temperature_comp_hbox = ipw.HBox(
            children = [self.target_temperature_comp_label, 
                        self.target_temperature_value_comp_textbox, 
                        self.target_temperature_unit_comp_dropdown]
        )
        
        self.evaporator_p_value_label = ipw.Label("P-value")
        self.evaporator_p_value_textbox = ipw.Text()
        self.evaporator_p_value_hbox = ipw.HBox(children = [self.evaporator_p_value_label, self.evaporator_p_value_textbox])
        
        self.evaporator_i_value_label = ipw.Label("I-value")
        self.evaporator_i_value_textbox = ipw.Text()
        self.evaporator_i_value_hbox = ipw.HBox(children = [self.evaporator_i_value_label, self.evaporator_i_value_textbox])
        
        self.ep_percentage_label = ipw.Label("EP (%)")
        self.ep_percentage_textbox = ipw.Text()
        self.ep_percentage_hbox = ipw.HBox(children = [self.ep_percentage_label, self.ep_percentage_textbox])
        
        self.components_properties_widgets = {
            "target_temperature": self.target_temperature_comp_hbox,
            "bias_voltage": self.bias_voltage_hbox,
            "discharge_voltage": self.discharge_voltage_hbox,
            "discharge_current": self.discharge_current_hbox,
            "evaporator_p_value": self.evaporator_p_value_hbox,
            "evaporator_i_value": self.evaporator_i_value_hbox,
            "ep_percentage": self.ep_percentage_hbox
        }
        # END - Widgets for component properties
        
        self.comments_label = ipw.Label(value = "Comments")
        self.comments_textarea = ipw.Textarea()
        self.comments_hbox = ipw.HBox(children = [self.comments_label, self.comments_textarea])
        
        self.remove_action_button = ipw.Button(
            description = 'Remove', 
            disabled = False, 
            button_style = 'danger', 
            tooltip = 'Remove action', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.name_textbox.observe(self.change_action_title, names = "value")
        self.action_type_dropdown.observe(self.load_action_properties, names = "value")
        self.component_dropdown.observe(self.load_component_settings_list, names = "value")
        self.remove_action_button.on_click(self.remove_action)

        if action_settings:
            self.load_action(action_settings)
        
        self.children = [
            self.action_type_hbox,
            self.action_properties_widgets,
            self.remove_action_button
        ]
    
    def load_action(self, settings):
        self.action_type_dropdown.value = settings["type"]
        self.name_textbox.value = settings.get("name", "")
        duration_str = settings.get("duration", "")
        
        # Split into days and time
        if "days" in duration_str:
            days_part, time_part = duration_str.split(" days ")
            self.duration_days_intbox.value = int(days_part)
        else:
            self.duration_days_intbox.value = 0
            time_part = duration_str
        
        # Split time into hours, minutes, seconds
        hours, minutes, seconds = map(int, time_part.split(":"))
        self.duration_hours_intbox.value = hours
        self.duration_minutes_intbox.value = minutes
        self.duration_seconds_intbox.value = seconds
        
        self.description_textbox.value = settings.get("description", "")
        action_target_temperature = settings.get("target_temperature", "")
        action_cryogen = settings.get("cryogen", "")
        action_substance = settings.get("substance", "")
        
        if action_cryogen:
            self.cryogen_textbox.value = action_cryogen
        
        if action_target_temperature:
            self.target_temperature_value_textbox.value = str(action_target_temperature["value"])
            self.target_temperature_unit_dropdown.value = action_target_temperature["temperature_unit"]
        
        if action_substance:
            self.substance_dropdown.value = action_substance
        
        try:
            component_settings = settings.get("component", {})
            if component_settings:
                self.component_dropdown.value = component_settings["identifier"]
                if "target_temperature" in component_settings:
                    self.target_temperature_value_comp_textbox.value = str(component_settings["target_temperature"]["value"])
                    self.target_temperature_unit_comp_dropdown.value = component_settings["target_temperature"]["temperature_unit"]
                    
                elif "bias_voltage" in component_settings:
                    self.bias_voltage_value_textbox.value = str(component_settings["bias_voltage"]["value"])
                    self.bias_voltage_unit_dropdown.value = component_settings["bias_voltage"]["voltage_unit"]
                
                elif "discharge_voltage" in component_settings:
                    self.discharge_voltage_value_textbox.value = str(component_settings["discharge_voltage"]["value"])
                    self.discharge_voltage_unit_dropdown.value = component_settings["discharge_voltage"]["voltage_unit"]
                
                elif "discharge_current" in component_settings:
                    self.discharge_current_value_textbox.value = str(component_settings["discharge_current"]["value"])
                    self.discharge_current_unit_dropdown.value = component_settings["discharge_current"]["current_unit"]
                
                self.evaporator_p_value_textbox.value = str(component_settings.get("evaporator_p_value", ""))
                self.evaporator_i_value_textbox.value = str(component_settings.get("evaporator_i_value", ""))
                self.ep_percentage_textbox.value = str(component_settings.get("ep_percentage", ""))
        except:
            pass

        action_substrate_temperature = settings.get("substrate_temperature", "")
        if action_substrate_temperature:
            self.substrate_temperature_value_textbox.value = str(action_substrate_temperature["value"])
            self.substrate_temperature_unit_dropdown.value = action_substrate_temperature["temperature_unit"]
        
        self.comments_textarea.value = settings.get("comments", "")
    
    def load_action_properties(self, change):
        action_type = self.action_type_dropdown.value
        if action_type == "-1":
            self.action_properties_widgets.children = []
        else:
            action_properties = [
                self.name_hbox,
                self.description_hbox,
                self.duration_hbox
            ]
            
            if action_type == "ANNEALING":
                action_properties.append(self.target_temperature_hbox)
                
            elif action_type == "COOLDOWN":
                action_properties.append(self.target_temperature_hbox)
                action_properties.append(self.cryogen_hbox)
                
            elif action_type == "DEPOSITION":
                action_properties.append(self.substance_hbox)
                action_properties.append(self.substrate_temperature_hbox)
                
            elif action_type == "DOSING":
                action_properties.append(self.substrate_temperature_hbox)
                action_properties.append(self.pressure_hbox)
                action_properties.append(self.gas_hbox)
                
            elif action_type == "SPUTTERING":
                action_properties.append(self.pressure_hbox)
                action_properties.append(self.current_hbox)
                action_properties.append(self.angle_hbox)
                action_properties.append(self.substrate_temperature_hbox)
                action_properties.append(self.sputter_ion_hbox)
            
            action_properties.append(self.component_hbox)
            action_properties.append(self.component_settings_hbox)
            action_properties.append(self.comments_hbox)

            component_options = self.get_instrument_components(self.instrument_permid, action_type)
            component_options.insert(0, ("Select a component...", "-1"))
            self.component_dropdown.options = component_options
            self.component_dropdown.value = "-1"
            
            self.action_properties_widgets.children = action_properties
    
    def get_instrument_components(self, instrument_permid, action_type):
        component_list = []
        if instrument_permid != "-1":
            instrument_object = get_cached_object(self.openbis_session, instrument_permid)
            instrument_type = str(instrument_object.type)
            instrument_components_properties = self.instrument_type_components_dictionary[instrument_type]
            all_instrument_components = []
            for prop in instrument_components_properties:
                prop_value = instrument_object.props[prop]
                if prop_value:
                    for component_id in prop_value:
                        if component_id in OPENBIS_SAMPLES_CACHE:
                            component_object = OPENBIS_SAMPLES_CACHE[component_id]
                        else:
                            component_object = utils.get_openbis_object(self.openbis_session, sample_ident = component_id)
                            # Save object information
                            OPENBIS_SAMPLES_CACHE[component_id] = component_object
                        
                        all_instrument_components.append(component_object)
                        
                        # Evaporators contain sub-components (evaporator slots)
                        if component_object.type == "EVAPORATOR":
                            evaporator_slots = component_object.props["evaporator_slots"]
                            for evaporator_slot_id in evaporator_slots:
                                evaporator_slot_object = get_cached_object(self.openbis_session, evaporator_slot_id)
                                all_instrument_components.append(evaporator_slot_object)
                
            for component_object in all_instrument_components:
                component_actions_settings_prop = component_object.props["actions_settings"]
                if component_actions_settings_prop:
                    component_actions_settings = json.loads(component_actions_settings_prop)
                    for component_action_settings in component_actions_settings:
                        component_action_type = component_action_settings["action_type"]
                        if action_type == component_action_type:
                            component_list.append((component_object.props["$name"], component_object.permId))

        return component_list
    
    def load_component_settings_list(self, change):
        action_type = self.action_type_dropdown.value
        if action_type != "-1":
            component_permid = self.component_dropdown.value
            if component_permid != "-1":
                if component_permid in OPENBIS_SAMPLES_CACHE:
                    component_object = OPENBIS_SAMPLES_CACHE[component_permid]
                else:
                    component_object = utils.get_openbis_object(self.openbis_session, sample_ident = component_permid)
                    # Save object information
                    OPENBIS_SAMPLES_CACHE[component_permid] = component_object
                
                component_settings_property = component_object.props["actions_settings"]
                component_settings_widgets = []
                if component_settings_property:
                    component_settings = json.loads(component_settings_property)
                    for settings in component_settings:
                        if action_type == settings["action_type"]:
                            settings_properties = settings.get("action_properties", {})
                            for prop in settings_properties:
                                setting_widget = self.components_properties_widgets.get(prop, None)
                                if setting_widget:
                                    if component_object.props[prop]:
                                        if prop == "target_temperature":
                                            target_temperature_comp = json.loads(component_object.props[prop])
                                            self.target_temperature_value_comp_textbox.value = target_temperature_comp["value"]
                                            self.target_temperature_unit_comp_dropdown.value = target_temperature_comp["temperature_unit"]
                                        
                                        elif prop == "bias_voltage":
                                            bias_voltage_comp = json.loads(component_object.props[prop])
                                            self.bias_voltage_value_comp_textbox.value = bias_voltage_comp["value"]
                                            self.bias_voltage_unit_comp_dropdown.value = bias_voltage_comp["voltage_unit"]
                                            
                                        elif prop == "discharge_voltage":
                                            discharge_voltage_comp = json.loads(component_object.props[prop])
                                            self.discharge_voltage_value_comp_textbox.value = discharge_voltage_comp["value"]
                                            self.discharge_voltage_unit_comp_dropdown.value = discharge_voltage_comp["voltage_unit"]
                                        
                                        elif prop == "discharge_current":
                                            discharge_current_comp = json.loads(component_object.props[prop])
                                            self.discharge_current_value_comp_textbox.value = discharge_current_comp["value"]
                                            self.discharge_current_unit_comp_dropdown.value = discharge_current_comp["current_unit"]
                                        
                                        elif prop == "evaporator_p_value":
                                            self.evaporator_p_value_textbox.value = str(component_object.props[prop])
                                        
                                        elif prop == "evaporator_i_value":
                                            self.evaporator_i_value_textbox.value = str(component_object.props[prop])  
                                        
                                        elif prop == "ep_percentage":
                                            self.ep_percentage_textbox.value = str(component_object.props[prop])  
                                        
                                    component_settings_widgets.append(self.components_properties_widgets[prop])
                            
                self.component_settings_vbox.children = component_settings_widgets

    def change_action_title(self, change):
        title = self.name_textbox.value
        self.actions_accordion.set_title(self.action_index, title)
                    
    def remove_action(self, b):
        actions_accordion_children = list(self.actions_accordion.children)
        num_actions = len(actions_accordion_children)
        actions_accordion_children.pop(self.action_index)
        
        for index, action in enumerate(actions_accordion_children):
            if index >= self.action_index:
                action.action_index -= 1
                self.actions_accordion.set_title(action.action_index, action.name_textbox.value)

        self.actions_accordion.set_title(num_actions - 1, "")
        self.actions_accordion.children = actions_accordion_children

class RegisterObservableWidget(ipw.VBox):
    def __init__(self, openbis_session, observables_accordion, observable_index, instrument_permid, observable_settings = None):
        super().__init__()
        self.openbis_session = openbis_session
        self.observables_accordion = observables_accordion
        self.observable_index = observable_index
        self.instrument_permid = instrument_permid
        
        self.observable_type_label = ipw.Label(value = "Observable type")
        observable_type_options = OBSERVABLES_LABELS
        observable_type_options.insert(0, ("Select an observable type...", "-1"))
        self.observable_type_dropdown = ipw.Dropdown(
            options = observable_type_options,
            value = "-1"
        )
        
        self.observable_type_hbox = ipw.HBox(children = [self.observable_type_label, self.observable_type_dropdown])
        self.observable_properties_widgets = ipw.VBox()
        
        self.name_label = ipw.Label(value = "Name")
        self.name_textbox = ipw.Text()
        self.name_hbox = ipw.HBox(children = [self.name_label, self.name_textbox])
        
        self.description_label = ipw.Label(value = "Description")
        self.description_textbox = ipw.Text()
        self.description_hbox = ipw.HBox(children = [self.description_label, self.description_textbox])
        
        self.ch_name_label = ipw.Label(value = "Channel name")
        self.ch_name_textbox = ipw.Text()
        self.ch_name_hbox = ipw.HBox(children = [self.ch_name_label, self.ch_name_textbox])
        
        self.instrument_type_components_dictionary = {
            "INSTRUMENT.STM": [
                "pumps", "gauges", "vacuum_chambers", "ports_valves",
                "preparation_tools", "analysers", "mechanical_components",
                "stm_components", "control_data_acquisition", "temperature_environment_control",
                "auxiliary_components", "tips_sensors", "accessories"
            ]
        }
        
        self.component_label = ipw.Label(value = "Component")
        self.component_dropdown = ipw.Dropdown()
        self.component_hbox = ipw.HBox(children = [self.component_label, self.component_dropdown])
        self.component_settings_label = ipw.Label(value = "Component settings:")
        self.component_settings_vbox = ipw.VBox()
        self.component_settings_hbox = ipw.HBox(children = [self.component_settings_label, self.component_settings_vbox])
        self.component_vbox = ipw.VBox(children = [self.component_hbox, self.component_settings_hbox])
        
        # BEGIN - Widgets for component properties
        self.filament_label = ipw.Label("Filament")
        self.filament_textbox = ipw.Text()
        self.filament_hbox = ipw.HBox(
            children = [self.filament_label, self.filament_textbox]
        )
        
        self.filament_current_label = ipw.Label("Filament current")
        self.filament_current_value_textbox = ipw.Text()
        self.filament_current_unit_dropdown = ipw.Dropdown(options = ["A", "mA"], value = "A")
        self.filament_current_hbox = ipw.HBox(
            children = [self.filament_current_label, self.filament_current_value_textbox, self.filament_current_unit_dropdown]
        )
        
        self.density_label = ipw.Label("Density")
        self.density_value_textbox = ipw.Text()
        self.density_unit_dropdown = ipw.Dropdown(options = ["g/m3", "g/cm3"], value = "g/m3")
        self.density_hbox = ipw.HBox(
            children = [self.density_label, self.density_value_textbox, self.density_unit_dropdown]
        )
        # END
        
        self.components_properties_widgets = {
            "filament": self.filament_hbox,
            "filament_current": self.filament_current_hbox,
            "density": self.density_hbox
        }
        
        self.comments_label = ipw.Label(value = "Comments")
        self.comments_textarea = ipw.Textarea()
        self.comments_hbox = ipw.HBox(children = [self.comments_label, self.comments_textarea])
        
        self.upload_readings_label = ipw.Label(value = "Upload readings")
        self.upload_readings_widget = ipw.FileUpload(multiple = True)
        self.upload_readings_hbox = ipw.HBox(children = [self.upload_readings_label, self.upload_readings_widget])
        
        self.remove_observable_button = ipw.Button(
            description = 'Remove', 
            disabled = False, 
            button_style = 'danger', 
            tooltip = 'Remove observable', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.name_textbox.observe(self.change_observable_title, names = "value")
        self.observable_type_dropdown.observe(self.load_observable_properties, names = "value")
        self.component_dropdown.observe(self.load_component_settings_list, names = "value")
        self.remove_observable_button.on_click(self.remove_observable)
        
        if observable_settings:
            self.load_observable(observable_settings)
        
        self.children = [
            self.observable_type_hbox,
            self.observable_properties_widgets,
            self.upload_readings_hbox,
            self.remove_observable_button
        ]
    
    def load_observable(self, settings):
        self.observable_type_dropdown.value = settings["type"]
        self.name_textbox.value = settings.get("name", "")
        self.description_textbox.value = settings.get("description", "")
        self.ch_name_textbox.value = settings.get("channel_name", "")
        self.comments_textarea.value = settings.get("comments", "")
        
        try:
            component_settings = settings.get("component", {})
            if component_settings:
                self.component_dropdown.value = component_settings["identifier"]
                if "density" in component_settings:
                    self.density_value_textbox.value = str(component_settings["density"]["value"])
                    self.density_unit_dropdown.value = component_settings["density"]["density_unit"]
                
                elif "filament_current" in component_settings:
                    self.filament_current_value_textbox.value = str(component_settings["filament_current"]["value"])
                    self.filament_current_unit_dropdown.value = component_settings["filament_current"]["current_unit"]
                
                self.filament_textbox.value = component_settings.get("filament", "")
        except:
            pass
        
    def load_observable_properties(self, change):
        observable_type = self.observable_type_dropdown.value
        if observable_type == "-1":
            self.observable_properties_widgets.children = []
        else:
            observable_properties = [
                self.name_hbox,
                self.description_hbox,
                self.ch_name_hbox,
                self.component_hbox,
                self.component_settings_hbox,
                self.comments_hbox
            ]

            component_options = self.get_instrument_components(self.instrument_permid, observable_type)
            component_options.insert(0, ("Select a component...", "-1"))
            self.component_dropdown.options = component_options
            self.component_dropdown.value = "-1"
            
            self.observable_properties_widgets.children = observable_properties
    
    def get_instrument_components(self, instrument_permid, observable_type):
        component_list = []
        if instrument_permid != "-1":
            instrument_object = get_cached_object(self.openbis_session, instrument_permid)
            instrument_type = str(instrument_object.type)
            instrument_components_properties = self.instrument_type_components_dictionary[instrument_type]
            all_instrument_components = []
            for prop in instrument_components_properties:
                prop_value = instrument_object.props[prop]
                if prop_value:
                    for component_id in prop_value:
                        component_object = get_cached_object(self.openbis_session, component_id)
                        all_instrument_components.append(component_object)
                        
                        # Evaporators contain sub-components (evaporator slots)
                        if component_object.type == "EVAPORATOR":
                            evaporator_slots = component_object.props["evaporator_slots"]
                            for evaporator_slot_id in evaporator_slots:
                                evaporator_slot_object = get_cached_object(self.openbis_session, evaporator_slot_id)
                                all_instrument_components.append(evaporator_slot_object)
                
            for component_object in all_instrument_components:
                component_observables_settings_prop = component_object.props["observables_settings"]
                if component_observables_settings_prop:
                    component_observables_settings = json.loads(component_observables_settings_prop)
                    for component_observable_settings in component_observables_settings:
                        component_observable_type = component_observable_settings["observable_type"]
                        if observable_type == component_observable_type:
                            component_list.append((component_object.props["$name"], component_object.permId))
                            
            for prop in instrument_components_properties:
                prop_value = instrument_object.props[prop]
                if prop_value:
                    for component_id in prop_value:
                        component_object = get_cached_object(self.openbis_session, component_id)
                        component_observables_settings_prop = component_object.props["observables_settings"]
                        if component_observables_settings_prop:
                            component_observables_settings = json.loads(component_observables_settings_prop)
                            if observable_type in component_observables_settings:
                                component_list.append((component_object.props["$name"], component_id))
        return component_list
    
    def load_component_settings_list(self, change):
        observable_type = self.observable_type_dropdown.value
        if observable_type != "-1":
            component_permid = self.component_dropdown.value
            if component_permid != "-1":
                component_object = get_cached_object(self.openbis_session, component_permid)
                component_settings_property = component_object.props["observables_settings"]
                component_settings_widgets = []
                if component_settings_property:
                    component_settings = json.loads(component_settings_property)
                    for settings in component_settings:
                        if observable_type == settings["observable_type"]:
                            settings_properties = settings.get("observable_properties", {})
                            for prop in settings_properties:
                                setting_widget = self.components_properties_widgets.get(prop, None)
                                if setting_widget:
                                    if component_object.props[prop]:
                                        if prop == "density":
                                            density_comp = json.loads(component_object.props[prop])
                                            self.density_value_textbox.value = density_comp["value"]
                                            self.density_unit_dropdown.value = density_comp["temperature_unit"]
                                        
                                        elif prop == "filament_current":
                                            filament_current_comp = json.loads(component_object.props[prop])
                                            self.filament_current_value_textbox.value = filament_current_comp["value"]
                                            self.filament_current_unit_dropdown.value = filament_current_comp["voltage_unit"]
                                        
                                        elif prop == "filament":
                                            self.filament_textbox.value = component_object.props[prop]
                                        
                                    component_settings_widgets.append(self.components_properties_widgets[prop])
                            
                self.component_settings_vbox.children = component_settings_widgets
                
    def change_observable_title(self, change):
        title = self.name_textbox.value
        self.observables_accordion.set_title(self.observable_index, title)
    
    def remove_observable(self, b):
        observables_accordion_children = list(self.observables_accordion.children)
        num_observables = len(observables_accordion_children)
        observables_accordion_children.pop(self.observable_index)
        
        for index, observable in enumerate(observables_accordion_children):
            if index >= self.observable_index:
                observable.observable_index -= 1
                self.observables_accordion.set_title(observable.observable_index, observable.name_textbox.value)

        self.observables_accordion.set_title(num_observables - 1, "")
        self.observables_accordion.children = observables_accordion_children