import ipywidgets as ipw
import src.utils as utils
from IPython.display import display, Javascript
import pandas as pd
import json
import re
import os
from collections import Counter
import subprocess
import shutil
from aiida import orm
import atexit
import ipyfilechooser
import src.aiida_utils as aiida_utils
import signal

OPENBIS_SAMPLES_CACHE = {}

ACTIONS_CODES = {
    "ANNEALING": "HEAT",
    "COOLDOWN": "COOL",
    "SPUTTERING": "IONB",
    "DOSING": "GASD",
    "DEPOSITION": "DEPO",
    "MECHANICAL_PRESSING": "MEPR",
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

SIMULATION_TYPES = [
    ("Geometry optimisation", "GEOMETRY_OPTIMISATION"),
    ("Band structure", "BAND_STRUCTURE"),
    ("PDOS", "PDOS"),
    ("Vibrational spectroscopy", "VIBRATIONAL_SPECTROSCOPY"),
    ("Unclassified simulation", "UNCLASSIFIED_SIMULATION")
]

INSTRUMENTS_TYPES = [
    "INSTRUMENT.STM"
]

OPENBIS_OBJECT_TYPES = {
    "Process": "PROCESS",
    "Process Step": "PROCESS_STEP",
    "Measurement Session": "MEASUREMENT_SESSION",
    "Atomistic Model": "ATOMISTIC_MODEL",
    "Aiida Node": "AIIDA_NODE",
    "Band Structure": "BAND_STRUCTURE",
    "Geometry Optimisation": "GEOMETRY_OPTIMISATION",
    "PDOS": "PDOS",
    "Vibrational Spectroscopy": "VIBRATIONAL_SPECTROSCOPY",
    "Unclassified Simulation": "UNCLASSIFIED_SIMULATION",
    "Code": "CODE",
    "Molecule": "MOLECULE",
    "Reaction Product Concept": "REACTION_PRODUCT_CONCEPT",
    "Reaction Product": "REACTION_PRODUCT",
    "Analysis": "ANALYSIS",
    "Sample": "SAMPLE",
    "Preparation": "PREPARATION",
    "Annealing": "ANNEALING",
    "Cooldown": "COOLDOWN",
    "Deposition": "DEPOSITION",
    "Dosing": "DOSING",
    "Sputtering": "SPUTTERING",
    "Substance": "SUBSTANCE",
    "Instrument STM": "INSTRUMENT.STM",
    "Evaporator": "EVAPORATOR",
    "Evaporator Slot": "EVAPORATOR_SLOT",
    "PBN Stage": "PBN_STAGE",
    "Sputter Gun": "SPUTTER_GUN",
    "Ion Gauge": "ION_GAUGE",
    "Analyser": "ANALYSER",
    "Result": "RESULT"
    
}

OPENBIS_OBJECT_CODES = {
    "Process Step": "PRST",
    "Sample": "SAMP",
}

OPENBIS_COLLECTIONS_PATHS = {
    "Atomistic Model": "/MATERIALS/ATOMISTIC_MODELS/ATOMISTIC_MODEL_COLLECTION",
    "Precursor Molecule": "/MATERIALS/MOLECULES/PRECURSOR_COLLECTION",
    "Precursor Substance": "/MATERIALS/MOLECULES/PRECURSOR_COLLECTION",
    "Reaction Product": "/MATERIALS/MOLECULES/PRODUCT_COLLECTION",
    "Reaction Product Concept": "/MATERIALS/MOLECULES/PRODUCT_COLLECTION",
    "Sample": "/MATERIALS/SAMPLES/SAMPLE_COLLECTION",
    "Instrument": "/EQUIPMENT/ILOG/INSTRUMENT_COLLECTION",
    "Chemical": "/MATERIALS/RAW_MATERIALS/CHEMICAL_COLLECTION"
}

def is_quantity_value(d):
    d_keys = set(d.keys())
    if d_keys == {"value", "unit"}:
        return True
    else:
        return False

def get_cached_object(ob_session, obj_id):
    if obj_id in OPENBIS_SAMPLES_CACHE:
        sample_object = OPENBIS_SAMPLES_CACHE[obj_id]
    else:
        sample_object = utils.get_openbis_object(ob_session, sample_ident = obj_id)
        # Save object information
        OPENBIS_SAMPLES_CACHE[obj_id] = sample_object
        
    return sample_object

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

class RunningMeasurementWatchdogsWidget(ipw.VBox):
    def __init__(self, openbis_session, session_data):
        super().__init__()
        self.openbis_session = openbis_session
        self.session_data = session_data
        
        self.running_watchdogs_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Running watchdogs</span>"
        )
        
        self.running_watchdogs_widget = ipw.SelectMultiple(
            layout = ipw.Layout(width = "500px")
        )
        
        self.stop_watchdog_button = ipw.Button(
            icon = "ban",
            tooltip = "Stop watchdog",
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        self.stop_watchdog_button.on_click(self.stop_watchdog)
        
        self.children = [
            self.running_watchdogs_title,
            self.running_watchdogs_widget,
            self.stop_watchdog_button
        ]
    
    def stop_watchdog(self, b):
        selected_pids = self.running_watchdogs_widget.value

        if selected_pids:
            for pid in selected_pids:
                os.kill(pid, signal.SIGTERM)
                print("Watchdog stopped.")

            self.running_watchdogs_widget.options = [
                (directory, pid) for directory, pid in self.running_watchdogs_widget.options
                if pid not in selected_pids
            ]
            
        else:
            display(Javascript(data="alert('Select at least one directory.')"))

class GenerateMeasurementsWatchdogWidget(ipw.VBox):
    def __init__(self, openbis_session, session_data, running_watchdogs_widget):
        super().__init__()
        self.openbis_session = openbis_session
        self.session_data = session_data
        self.running_watchdogs_widget = running_watchdogs_widget
        
        self.select_experiment_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select experiment</span>"
        )
        
        self.select_experiment_widget = SelectExperimentWidget(self.openbis_session)
        
        self.select_sample_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select sample</span>"
        )
        
        self.select_sample_widget = SelectSampleWidget(self.openbis_session)
        
        self.select_instrument_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select instrument</span>"
        )
        
        self.select_instrument_widget = SelectInstrumentWidget(self.openbis_session)
        
        self.select_measurements_folder_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select measurements directory</span>"
        )
        
        self.select_measurements_folder_widget = ipyfilechooser.FileChooser(
            path = '.', select_default=True, use_dir_icons=True, show_only_dirs = True)
        
        self.generate_watchdog_button = ipw.Button(
            description = '', disabled = False, button_style = '', tooltip = 'Save', 
            icon = 'save', layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        self.select_sample_widget.sample_dropdown.observe(self.load_sample_data, names = "value")
        self.generate_watchdog_button.on_click(self.generate_watchdog)
        
        self.watchdog_processes = []
        
        # Ensure process is killed on notebook shutdown / kernel restart
        atexit.register(self.cleanup_watchdog)
        
        self.children = [
            self.select_experiment_title,
            self.select_experiment_widget,
            self.select_sample_title,
            self.select_sample_widget,
            self.select_instrument_title,
            self.select_instrument_widget,
            self.select_measurements_folder_title,
            self.select_measurements_folder_widget,
            self.generate_watchdog_button
        ]
    
    def load_sample_data(self, change):
        sample_id = self.select_sample_widget.sample_dropdown.value
        if sample_id == "-1":
            return

        sample_object = get_cached_object(self.openbis_session, sample_id)
        
        sample_object_parents = sample_object.parents
        most_recent_parent = None
        
        for parent_id in sample_object_parents:
            parent_object = get_cached_object(self.openbis_session, parent_id)
            
            parent_type = parent_object.type
            if parent_type == OPENBIS_OBJECT_TYPES["Process Step"]:
                if most_recent_parent:
                    if parent_object.registrationDate > most_recent_parent.registrationDate:
                        most_recent_parent = parent_object
                else:
                    most_recent_parent = parent_object
        
        if most_recent_parent:
            experiment_id = self.select_experiment_widget.experiment_dropdown.value
            if most_recent_parent.experiment.permId != experiment_id:
                self.select_experiment_widget.experiment_dropdown.value = most_recent_parent.experiment.permId
                display(Javascript(data = "alert('Experiment was changed!')"))
    
    def generate_watchdog(self, b):
        experiment_id = self.select_experiment_widget.experiment_dropdown.value
        if experiment_id == "-1":
            return
        
        sample_id = self.select_sample_widget.sample_dropdown.value
        if sample_id == "-1":
            return
        
        sample_object = utils.get_openbis_object(self.openbis_session, sample_ident = sample_id)
        sample_name = sample_object.props["name"]

        instrument_id = self.select_instrument_widget.instrument_dropdown.value
        if instrument_id == "-1":
            return

        measurement_session_object = utils.create_openbis_object(
            self.openbis_session,
            type = OPENBIS_OBJECT_TYPES["Measurement Session"],
            collection = experiment_id,
            parents = [sample_id, instrument_id],
            props = {"name": f"Measurement Session on Sample {sample_name}", "default_object_view": "IMAGING_GALLERY_VIEW"}
        )
        measurement_session_id = measurement_session_object.permId
        measurements_directory = self.select_measurements_folder_widget.selected_path
        watchdog_file = f"/home/jovyan/aiida-openbis/Notebooks/openbis_app/src/measurements_uploader.py"
        shutil.copy(watchdog_file, measurements_directory)
        
        watchdog_process = subprocess.Popen(
            [
                "python", 
                f"{measurements_directory}/measurements_uploader.py", 
                "--openbis_url", self.session_data["url"],
                "--openbis_token", self.session_data["token"],
                "--measurement_session_id", measurement_session_id,
                "--data_folder", measurements_directory
            ]
        )
        print("Watchdog process started with PID:", watchdog_process.pid)
        
        self.watchdog_processes.append(watchdog_process)
        
        running_watchdogs = self.running_watchdogs_widget.running_watchdogs_widget.options
        running_watchdogs = list(running_watchdogs)
        running_watchdogs.append((measurements_directory, watchdog_process.pid))
        self.running_watchdogs_widget.running_watchdogs_widget.options = running_watchdogs
    
    def cleanup_watchdog(self):
        if self.watchdog_processes:
            for process in self.watchdog_processes:
                print(f"Terminating watchdog process with PID: {process.pid}")
                process.terminate()
                self.watchdog_processes = []

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
        
        material_type_options = MATERIALS_CONCEPTS_LABELS.copy()
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
        self.search_button = ipw.Button(
            disabled = False, icon = 'search',
            tooltip = 'Search simulations in openBIS', 
            layout = ipw.Layout(width = '50px', height = '25px')
        )
        self.search_logical_operator_hbox = ipw.HBox(
            children = [self.search_logical_operator_label, self.search_logical_operator_dropdown, self.search_button]
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
        
        self.import_simulations_button = ipw.Button(
            tooltip = 'Import simulations', icon = 'download', 
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        self.import_simulations_message_html = ipw.HTML()
        
        # Increase search button icon size
        increase_search_button = ipw.HTML(
            """<style>
            .fa-search {font-size: 1.5em !important;}
            .fa-download {font-size: 2em !important;}
            </style>
            """
        )
        
        # Add functionality to the widgets
        self.material_type_dropdown.observe(self.load_material_type_widgets, names = "value")
        self.add_molecule_button.on_click(self.add_molecule)
        self.add_reacprod_concept_button.on_click(self.add_reacprod_concept)
        self.search_button.on_click(self.search_simulations)
        self.import_simulations_button.on_click(self.import_aiida_nodes)
        
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
                    if simulation_object.type == OPENBIS_OBJECT_TYPES["Measurement Session"]:
                        simulation_measurement = False
                        for parent in simulation_object.parents:
                            parent_object = get_cached_object(self.openbis_session, parent)
                            if parent_object.type == OPENBIS_OBJECT_TYPES["Atomistic Model"]:
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
                    if simulation_object.type == OPENBIS_OBJECT_TYPES["Measurement Session"]: # 2D Measurement is used for both simulation and experiments
                        simulation_measurement = False
                        for parent in simulation_object.parents:
                            parent_object = get_cached_object(self.openbis_session, parent)
                            if parent_object.type == OPENBIS_OBJECT_TYPES["Atomistic Model"]:
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
                if parent_object.type == OPENBIS_OBJECT_TYPES["Aiida Node"]:
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
                aiida_node_object_type = aiida_node_object.type.code
                aiida_node_object_type_lower = aiida_node_object_type.lower()
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
                            workchain = orm.load_node(aiida_node_object.props["wfms_uuid"])
                            workchain_viewer_link = WORKCHAIN_VIEWERS[workchain.process_label]
                            notebook_link = f"{workchain_viewer_link}?pk={workchain.pk}"
                            
                            simulation_message = f'<a href="{notebook_link}">Workchain {selected_label} successfully imported.</a>\n'
                            selected_simulations_messages += simulation_message
                        
                        shutil.rmtree(f"aiida_nodes/{dataset.permId}/")
            else:
                simulation_message = f"Workchain {selected_label} cannot be imported because it was done manually.\n"
                selected_simulations_messages += simulation_message
    
        self.import_simulations_message_html.value = selected_simulations_messages
    
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
            materials_objects_names_permids = [(obj.props["name"], obj.permId) for obj in material_objects]
            material_options += materials_objects_names_permids
            material_dropdown.options = material_options
            
            def sort_material_dropdown(change):
                options = material_options[1:]
                
                df = pd.DataFrame(options, columns=["name", "registration_date"])
                if name_checkbox.value and not registration_date_checkbox.value:
                    df = df.sort_values(by="name", ascending=True)
                elif not name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by="registration_date", ascending=False)
                elif name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

                options = list(df.itertuples(index=False, name=None))
                options.insert(0, material_options[0])
                material_dropdown.options = options
            
            def load_material_details(change):
                obj_permid = material_dropdown.value
                if obj_permid == "-1":
                    return
                else:
                    obj = utils.get_openbis_object(obj_permid)
                    obj_props = obj.props.all()
                    obj_name = obj_props.get("name", "")
                    obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
                    for key, value in obj_props.items():
                        if value:
                            prop_type = utils.get_property_type(self.openbis_session, code = key)
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
        
        self.select_experiment_widget = SelectExperimentWidget(self.openbis_session)
        
        self.used_aiida_checkbox = ipw.Checkbox(
            value = False,
            description = "Simulation developed using AiiDA",
            indent = False
        )
        
        self.simulation_details_vbox = SimulationDetailsWidget(self.openbis_session, True)
        
        self.save_simulations_button = ipw.Button(
            tooltip = 'Import simulations', icon = 'save', 
            layout = ipw.Layout(width = '100px', height = '50px')
        )
        
        # Increase search button icon size
        increase_search_button = ipw.HTML(
            """<style>
            .fa-search {font-size: 1.5em !important;}
            .fa-download {font-size: 2em !important;}
            .fa-save {font-size: 2em !important;}
            </style>
            """
        )
        
        # Add functionality to the widgets
        self.used_aiida_checkbox.observe(self.load_simulations_details_widgets, names = "value")
        self.save_simulations_button.on_click(self.export_simulation_to_openbis)
        self.used_aiida_checkbox.value = True
        
        self.children = [
            self.select_experiment_title,
            self.select_experiment_widget,
            self.simulation_details_title,
            self.used_aiida_checkbox,
            self.simulation_details_vbox,
            increase_search_button,
            self.save_simulations_button
        ]
    
    def load_simulations_details_widgets(self, change):
        used_aiida = self.used_aiida_checkbox.value
        self.simulation_details_vbox.load_widgets(used_aiida)
        
    def export_simulation_to_openbis(self, b):
        selected_experiment_id = self.select_experiment_widget.experiment_dropdown.value
        if selected_experiment_id == "-1":
            display(Javascript(data = "alert('Select an experiment.')"))
        else:
            selected_molecules_widgets = self.simulation_details_vbox.molecules_accordion.children
            selected_reac_prods_widgets = self.simulation_details_vbox.reacprod_concepts_accordion.children
            
            # Get material
            if self.simulation_details_vbox.material_type_dropdown.value == "-1":
                selected_slab = []
            else:
                selected_material_id = self.simulation_details_vbox.material_details_vbox.children[0].children[0].value
                if selected_material_id == "-1":
                    selected_slab = []
                else:
                    selected_slab = [selected_material_id]
            
            # Get molecules
            selected_molecules_ids = []
            for mol_widget in selected_molecules_widgets:
                mol_id = mol_widget.dropdown.value
                if mol_id == "-1":
                    continue
                else:
                    selected_molecules_ids.append(mol_id)
            
            # Get reaction products concepts
            selected_reac_prods_ids = []
            for reac_prod_widget in selected_reac_prods_widgets:
                reac_prod_id = reac_prod_widget.dropdown.value
                if reac_prod_id == "-1":
                    continue
                else:
                    selected_reac_prods_ids.append(reac_prod_id)
                    
            if self.used_aiida_checkbox.value:
                selected_simulation_id = self.simulation_details_vbox.simulations_dropdown.value
                if selected_simulation_id == "-1":
                    display(Javascript(data = "alert('Select a simulation.')"))
                else:
                    atom_model_parents = selected_slab + selected_molecules_ids + selected_reac_prods_ids
                    last_export = aiida_utils.export_workchain(
                        self.openbis_session,
                        selected_experiment_id,
                        selected_simulation_id
                    )
                    
                    if last_export:
                        first_atom_model = utils.find_first_atomistic_model(
                            self.openbis_session,
                            last_export,
                            OPENBIS_OBJECT_TYPES["Atomistic Model"]
                        )
                        
                        if len(first_atom_model.parents) == 0:
                            first_atom_model.parents = atom_model_parents
                            first_atom_model.save()
                        display(Javascript(data = "alert('Upload successful!')"))
                    else:
                        display(Javascript(data = "alert('Simulation is already in openBIS!')"))
            
            else:
                simulation_type = self.simulation_details_vbox.simulation_type_dropdown.value
                if simulation_type == "-1":
                    return
                else:
                    # Get atomistic models
                    atom_model_widget = self.simulation_details_vbox.atom_model_widget
                    selected_atom_model_id = atom_model_widget.atom_model_dropdown.value
                    if selected_atom_model_id == "-1":
                        selected_atom_model_id = []
                    else:
                        selected_atom_model_id = [selected_atom_model_id]
                    
                    selected_codes_ids = self.simulation_details_vbox.codes_multi_selector.value
                    selected_codes_ids = list(selected_codes_ids)
                    
                    simulation_props_widget = self.simulation_details_vbox.simulation_properties_widget
                    
                    level_theory_params = simulation_props_widget.level_theory_parameters_textbox.value
                    if utils.is_valid_json(level_theory_params) == False:
                        level_theory_params = ""
                        
                    input_parameters = simulation_props_widget.method_input_parameters_textbox.value
                    if utils.is_valid_json(input_parameters) == False:
                        input_parameters = ""
                        
                    output_parameters = simulation_props_widget.method_output_parameters_textbox.value
                    if utils.is_valid_json(output_parameters) == False:
                        output_parameters = ""
                    
                    simulation_props = {
                        "name": simulation_props_widget.name_textbox.value,
                        "wfms_uuid": simulation_props_widget.wfms_uuid_textbox.value,
                        "input_parameters": input_parameters,
                        "output_parameters": output_parameters,
                        "comments": simulation_props_widget.comments_textbox.value
                    }
                    
                    simulation_types_with_level_theory = [
                        OPENBIS_OBJECT_TYPES["Band Structure"],
                        OPENBIS_OBJECT_TYPES["Geometry Optimisation"],
                        OPENBIS_OBJECT_TYPES["PDOS"],
                        OPENBIS_OBJECT_TYPES["Vibrational Spectroscopy"]
                    ]
                    if simulation_type in simulation_types_with_level_theory:
                        simulation_props["level_theory_method"] = simulation_props_widget.level_theory_method_textbox.value
                        simulation_props["level_theory_parameters"] = level_theory_params
                        
                        if simulation_type == OPENBIS_OBJECT_TYPES["Band Structure"]:
                            band_gap = {
                                "value": simulation_props_widget.band_gap_value_textbox.value,
                                "unit": simulation_props_widget.band_gap_unit_textbox.value,
                            }
                            simulation_props["band_gap"] = band_gap
                    
                        elif simulation_type == OPENBIS_OBJECT_TYPES["Geometry Optimisation"]:
                            simulation_props["cell_opt_constraints"] = simulation_props_widget.cell_opt_constraints_textbox.value
                            simulation_props["cell_optimised"] = simulation_props_widget.cell_optimised_checkbox.value
                            simulation_props["driver_code"] = simulation_props_widget.driver_code_textbox.value
                            simulation_props["constrained"] = simulation_props_widget.constrained_checkbox.value
                            
                            force_convergence_threshold = {
                                "value": simulation_props_widget.force_convergence_threshold_value_textbox.value,
                                "unit": simulation_props_widget.force_convergence_threshold_unit_textbox.value,
                            }
                            simulation_props["force_convergence_threshold"] = force_convergence_threshold
                            
                    elif simulation_type == OPENBIS_OBJECT_TYPES["Unclassified Simulation"]:
                        simulation_props["description"] = simulation_props_widget.description_textbox.value
                    
                    simulation_props["codes"] = selected_codes_ids
                    
                    simulation_parents = selected_atom_model_id
                    
                    simulation_obj = utils.create_openbis_object(
                        self.openbis_session,
                        type = simulation_type,
                        collection = selected_experiment_id,
                        parents = simulation_parents,
                        props = simulation_props
                    )
                    
                    # Simulation preview
                    for filename in self.simulation_details_vbox.upload_image_preview_uploader.value:
                        file_info = self.simulation_details_vbox.upload_image_preview_uploader.value[filename]
                        utils.write_file(file_info['content'], filename)
                        utils.create_openbis_dataset(
                            self.openbis_session,
                            type = "ELN_PREVIEW", 
                            sample = simulation_obj, 
                            files = [filename]
                        )
                        os.remove(filename)
                    
                    # Simulation datasets
                    for filename in self.simulation_details_vbox.upload_datasets_uploader.value:
                        file_info = self.simulation_details_vbox.upload_datasets_uploader.value[filename]
                        utils.write_file(file_info['content'], filename)
                        utils.create_openbis_dataset(
                            self.openbis_session,
                            type = "ATTACHMENT", 
                            sample = simulation_obj, 
                            files = [filename]
                        )
                        os.remove(filename)
                
class SimulationDetailsWidget(ipw.VBox):
    def __init__(self, openbis_session, used_aiida):
        super().__init__()
        self.openbis_session = openbis_session
        self.used_aiida = used_aiida
        
        self.select_molecules_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select molecules</span>"
        )
        
        self.select_reacprod_concepts_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select reaction product concepts</span>"
        )
        
        self.select_slab_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select slab</span>"
        )
        
        self.select_simulation_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select simulation</span>"
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
        
        material_type_options = MATERIALS_CONCEPTS_LABELS.copy()
        material_type_options.insert(0, ("Select material type...", "-1"))
        
        self.material_type_dropdown = ipw.Dropdown(
            options = material_type_options,
            value = material_type_options[0][1]
        )
        
        self.material_details_vbox = ipw.VBox()
        
        self.simulations_label = ipw.Label(value = "Simulation")
        self.simulations_dropdown = ipw.Dropdown()
        self.load_aiida_simulations()
        self.sort_simulations_label = ipw.Label(value = "Sort by:")
        
        self.sort_name_label = ipw.Label(
            value = "Name", 
            layout=ipw.Layout(margin='2px', width='50px'),
            style = {'description_width': 'initial'}
        )
        
        self.sort_name_checkbox = ipw.Checkbox(
            indent = False,
            layout=ipw.Layout(margin='2px', width='20px')
        )
        
        self.sort_pk_label = ipw.Label(
            value = "PK", 
            layout=ipw.Layout(margin='2px', width='110px'),
            style = {'description_width': 'initial'}
        )
        
        self.sort_pk_checkbox = ipw.Checkbox(
            indent = False,
            layout=ipw.Layout(margin='2px', width='20px')
        )
        
        self.sort_simulations_hbox = ipw.HBox(
            children = [
                self.sort_simulations_label,
                self.sort_name_checkbox,
                self.sort_name_label,
                self.sort_pk_checkbox,
                self.sort_pk_label
            ]
        )
        
        self.simulations_dropdown_hbox = ipw.HBox(
            children = [
                self.simulations_label,
                self.simulations_dropdown,
            ]
        )
        
        self.select_simulation_type_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select simulation type</span>"
        )
        
        self.simulation_type_label = ipw.Label(value = "Simulation type")
        simulation_types = SIMULATION_TYPES.copy()
        simulation_types.insert(0, ("Select simulation type...", "-1"))
        
        self.simulation_type_dropdown = ipw.Dropdown(
            options = simulation_types,
            value = "-1"
        )
        
        self.simulation_type_hbox = ipw.HBox(
            children = [
                self.simulation_type_label,
                self.simulation_type_dropdown,
            ]
        )
        
        self.simulation_properties_widget = SimulationPropertiesWidget(self.openbis_session)
        
        self.select_atom_model_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select atomistic model</span>"
        )
        self.atom_model_widget = AtomModelWidget(self.openbis_session)
        
        self.select_codes_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select codes</span>"
        )
        
        self.codes_label = ipw.Label(value = "Codes")
        codes_objects = utils.get_openbis_objects(self.openbis_session, type = OPENBIS_OBJECT_TYPES["Code"])
        codes_options = [(obj.props["name"], obj.permId) for obj in codes_objects]
        self.codes_multi_selector = ipw.SelectMultiple(options = codes_options)
        self.codes_hbox = ipw.HBox(children = [self.codes_label, self.codes_multi_selector])
        
        self.upload_image_preview_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Upload image preview</span>"
        )
        
        self.upload_image_preview_uploader = ipw.FileUpload(
            multiple = False,
            accept = '.jpg, .jpeg, .png,'
        )
        
        self.upload_datasets_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Upload datasets</span>"
        )
        
        self.upload_datasets_uploader = ipw.FileUpload(multiple = True)
        
        self.simulation_type_dropdown.observe(self.load_simulation_type_properties, names = "value")
        self.material_type_dropdown.observe(self.load_material_type_widgets, names = "value")
        self.add_molecule_button.on_click(self.add_molecule)
        self.add_reacprod_concept_button.on_click(self.add_reacprod_concept)
        
    def load_simulation_type_properties(self, change):
        simulation_type = self.simulation_type_dropdown.value
        self.simulation_properties_widget.load_widgets(simulation_type)
    
    def load_aiida_simulations(self):
        qb = orm.QueryBuilder()
        qb.append(orm.WorkChainNode)
        results = qb.all()
        
        # List of calculations that can be exported
        labels = list(WORKCHAIN_VIEWERS.keys())
        
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
        options = []
        for result in results:
            if result[3]:
                name_pk_string = f"{result[3][:20]} - {result[2]} (PK: {result[0]})"
            else:
                name_pk_string = f"{result[2]} (PK: {result[0]})"
                
            name_pk_tuple = (name_pk_string, result[0])
            options.append(name_pk_tuple)
        
        options.insert(0, (f'Select a simulation...', "-1"))
        self.simulations_dropdown.options = options
        self.simulations_dropdown.value = "-1"

    def sort_simulations_dropdown(self, change):
        options = self.simulations_dropdown.options[1:]
        
        df = pd.DataFrame(options, columns=["name", "PK"])
        if self.sort_name_checkbox.value and not self.sort_pk_checkbox.value:
            df = df.sort_values(by="name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_pk_checkbox.value:
            df = df.sort_values(by="PK", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_pk_checkbox.value:
            df = df.sort_values(by=["name", "PK"], ascending=[True, False])

        options = list(df.itertuples(index=False, name=None))
        options.insert(0, self.simulations_dropdown.options[0])
        self.simulations_dropdown.options = options

    def load_widgets(self, used_aiida):
        self.used_aiida = used_aiida
        if self.used_aiida:
            self.children = [
                self.select_molecules_title,
                self.molecules_accordion,
                self.add_molecule_button,
                self.select_reacprod_concepts_title,
                self.reacprod_concepts_accordion,
                self.add_reacprod_concept_button,
                self.select_material_title,
                self.material_type_dropdown,
                self.material_details_vbox,
                self.select_simulation_title,
                self.simulations_dropdown_hbox,
                self.sort_simulations_hbox
            ]
            
        else:
            self.children = [
                self.select_simulation_type_title,
                self.simulation_type_hbox,
                self.simulation_properties_widget,
                self.select_atom_model_title,
                self.atom_model_widget,
                self.select_codes_title,
                self.codes_hbox,
                self.upload_image_preview_title,
                self.upload_image_preview_uploader,
                self.upload_datasets_title,
                self.upload_datasets_uploader
            ]

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
            materials_objects_names_permids = [(obj.props["name"], obj.permId) for obj in material_objects]
            material_options += materials_objects_names_permids
            material_dropdown.options = material_options
            
            def sort_material_dropdown(change):
                options = material_options[1:]
                
                df = pd.DataFrame(options, columns=["name", "registration_date"])
                if name_checkbox.value and not registration_date_checkbox.value:
                    df = df.sort_values(by="name", ascending=True)
                elif not name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by="registration_date", ascending=False)
                elif name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

                options = list(df.itertuples(index=False, name=None))
                options.insert(0, material_options[0])
                material_dropdown.options = options
            
            def load_material_details(change):
                obj_permid = material_dropdown.value
                if obj_permid == "-1":
                    return
                else:
                    obj = utils.get_openbis_object(self.openbis_session, sample_ident = obj_permid)
                    obj_props = obj.props.all()
                    obj_name = obj_props.get("name", "")
                    obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
                    for key, value in obj_props.items():
                        if value:
                            prop_type = utils.get_openbis_property_type(self.openbis_session, code = key)
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
   
class SimulationPropertiesWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        self.simulation_type = ""
        
        self.simulation_properties_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Simulation properties</span>"
        )
        
        self.name_label = ipw.Label(value = "Name")
        self.name_textbox = ipw.Text()
        self.name_hbox = ipw.HBox(children = [self.name_label, self.name_textbox])
        
        self.description_label = ipw.Label(value = "Description")
        self.description_textbox = ipw.Textarea()
        self.description_hbox = ipw.HBox(children = [self.description_label, self.description_textbox])
        
        self.wfms_uuid_label = ipw.Label(value = "WFMS UUID")
        self.wfms_uuid_textbox = ipw.Text()
        self.wfms_uuid_hbox = ipw.HBox(children = [self.wfms_uuid_label, self.wfms_uuid_textbox])
        
        self.band_gap_label = ipw.Label(value = "Band gap:")
        self.band_gap_value_textbox = ipw.Text()
        self.band_gap_unit_textbox = ipw.Dropdown(
            options = ["eV", "J", "Ha"],
            value = "Ha"
        )
        self.band_gap_hbox = ipw.HBox(
            children = [
                self.band_gap_label,
                self.band_gap_value_textbox, 
                self.band_gap_unit_textbox
            ]
        )
        
        self.cell_opt_constraints_label = ipw.Label(value = "Cell opt constraints")
        self.cell_opt_constraints_textbox = ipw.Text()
        self.cell_opt_constraints_hbox = ipw.HBox(children = [self.cell_opt_constraints_label, self.cell_opt_constraints_textbox])
        
        self.cell_optimised_label = ipw.Label(value = "Cell optimised")
        self.cell_optimised_checkbox = ipw.Checkbox(indent = False)
        self.cell_optimised_hbox = ipw.HBox(children = [self.cell_optimised_label, self.cell_optimised_checkbox])
        
        self.driver_code_label = ipw.Label(value = "Driver code")
        self.driver_code_textbox = ipw.Text()
        self.driver_code_hbox = ipw.HBox(children = [self.driver_code_label, self.driver_code_textbox])
        
        self.constrained_label = ipw.Label(value = "Constrained")
        self.constrained_checkbox = ipw.Checkbox(indent = False)
        self.constrained_hbox = ipw.HBox(children = [self.constrained_label, self.constrained_checkbox])
        
        self.force_convergence_threshold_label = ipw.Label(value = "Force convergence threshold: ")
        self.force_convergence_threshold_value_label = ipw.Label(value = "Value")
        self.force_convergence_threshold_value_textbox = ipw.Text(layout = ipw.Layout(width = "100px"))
        self.force_convergence_threshold_unit_label = ipw.Label(value = "Unit")
        self.force_convergence_threshold_unit_textbox = ipw.Text(layout = ipw.Layout(width = "100px"))
        self.force_convergence_threshold_hbox = ipw.HBox(
            children = [
                self.force_convergence_threshold_label,
                self.force_convergence_threshold_value_label,
                self.force_convergence_threshold_value_textbox,
                self.force_convergence_threshold_unit_label,
                self.force_convergence_threshold_unit_textbox
            ]
        )
        
        self.level_theory_method_label = ipw.Label(value = "Level of theory (method)")
        self.level_theory_method_textbox = ipw.Text()
        self.level_theory_method_hbox = ipw.HBox(children = [self.level_theory_method_label, self.level_theory_method_textbox])
        
        self.level_theory_parameters_label = ipw.Label(value = "Level of theory (parameters)")
        self.level_theory_parameters_textbox = ipw.Textarea(
            placeholder = '{"uks": false, "charge": 0.0, "plus_u": false, "xc_functional": "PBESOL"}'
        )
        self.level_theory_parameters_hbox = ipw.HBox(
            children = [
                self.level_theory_parameters_label, 
                self.level_theory_parameters_textbox
            ]
        )
        
        self.method_input_parameters_label = ipw.Label(value = "Method input parameters")
        self.method_input_parameters_textbox = ipw.Textarea(
            placeholder = '{"degauss": 0.0, "volume": 0.0}'
        )
        self.method_input_parameters_hbox = ipw.HBox(
            children = [
                self.method_input_parameters_label, 
                self.method_input_parameters_textbox
            ]
        )
        
        self.method_output_parameters_label = ipw.Label(value = "Method output parameters")
        self.method_output_parameters_textbox = ipw.Textarea(
            placeholder = '{"energy": 0.0, "has_electric_field": false}'
        )
        self.method_output_parameters_hbox = ipw.HBox(
            children = [
                self.method_output_parameters_label, 
                self.method_output_parameters_textbox
            ]
        )
        
        self.comments_label = ipw.Label(value = "Comments")
        self.comments_textbox = ipw.Textarea()
        self.comments_hbox = ipw.HBox(children = [self.comments_label, self.comments_textbox])
    
    def load_widgets(self, simulation_type):
        self.simulation_type = simulation_type
        
        if simulation_type == OPENBIS_OBJECT_TYPES["Band Structure"]:
            self.children = [
                self.simulation_properties_title,
                self.name_hbox,
                self.wfms_uuid_hbox,
                self.band_gap_hbox,
                self.level_theory_method_hbox,
                self.level_theory_parameters_hbox,
                self.method_input_parameters_hbox,
                self.method_output_parameters_hbox,
                self.comments_hbox
            ]

        elif simulation_type == OPENBIS_OBJECT_TYPES["Geometry Optimisation"]:
            self.children = [
                self.simulation_properties_title,
                self.name_hbox,
                self.wfms_uuid_hbox,
                self.cell_opt_constraints_hbox,
                self.cell_optimised_hbox,
                self.driver_code_hbox,
                self.constrained_hbox,
                self.force_convergence_threshold_hbox,
                self.level_theory_method_hbox,
                self.level_theory_parameters_hbox,
                self.method_input_parameters_hbox,
                self.method_output_parameters_hbox,
                self.comments_hbox
            ]

        elif simulation_type == OPENBIS_OBJECT_TYPES["PDOS"]:
            self.children = [
                self.simulation_properties_title,
                self.name_hbox,
                self.wfms_uuid_hbox,
                self.level_theory_method_hbox,
                self.level_theory_parameters_hbox,
                self.method_input_parameters_hbox,
                self.method_output_parameters_hbox,
                self.comments_hbox
            ]

        elif simulation_type == OPENBIS_OBJECT_TYPES["Vibrational Spectroscopy"]:
            self.children = [
                self.simulation_properties_title,
                self.name_hbox,
                self.wfms_uuid_hbox,
                self.level_theory_method_hbox,
                self.level_theory_parameters_hbox,
                self.method_input_parameters_hbox,
                self.method_output_parameters_hbox,
                self.comments_hbox
            ]

        elif simulation_type == OPENBIS_OBJECT_TYPES["Unclassified Simulation"]:
            self.children = [
                self.simulation_properties_title,
                self.name_hbox,
                self.description_hbox,
                self.method_input_parameters_hbox,
                self.method_output_parameters_hbox,
                self.comments_hbox
            ]

class SelectInstrumentWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.instrument_label = ipw.Label(
            value = "Instrument"
        )

        self.instrument_dropdown = ipw.Dropdown()
        self.load_instruments()


        self.sort_instrument_label = ipw.Label(
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
        
        self.sort_instrument_hbox = ipw.HBox(
            children = [
                self.sort_instrument_label,
                self.sort_name_checkbox,
                self.sort_name_label,
                self.sort_registration_date_checkbox,
                self.sort_registration_date_label
            ]
        )
        
        self.instrument_dropdown_hbox = ipw.HBox(
            children = [
                self.instrument_label,
                self.instrument_dropdown,
            ]
        )
        
        self.sort_name_checkbox.observe(self.sort_instrument_dropdown, names = "value")
        self.sort_registration_date_checkbox.observe(self.sort_instrument_dropdown, names = "value")
        
        self.children = [
            self.instrument_dropdown_hbox,
            self.sort_instrument_hbox
        ]

    def load_instruments(self):
        instruments = []
        for instrument_type in INSTRUMENTS_TYPES:
            instruments_objects = utils.get_openbis_objects(
                self.openbis_session,
                type = instrument_type
            )
            instruments.extend(instruments_objects)
            
        instrument_options = [(f"{obj.props['name']}", obj.permId) for obj in instruments]
        instrument_options.insert(0, ("Select instrument...", "-1"))
        self.instrument_dropdown.options = instrument_options
        self.instrument_dropdown.value = "-1"

    def sort_instrument_dropdown(self, change):
        options = self.instrument_dropdown.options[1:]
        
        df = pd.DataFrame(options, columns=["name", "registration_date"])
        if self.sort_name_checkbox.value and not self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="registration_date", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

        options = list(df.itertuples(index=False, name=None))
        options.insert(0, self.instrument_dropdown.options[0])
        self.instrument_dropdown.options = options

class AtomModelWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.atom_model_label = ipw.Label(
            value = "Atomistic model"
        )
        
        self.atom_model_dropdown = ipw.Dropdown()
        self.load_atom_models()
        
        
        self.sort_atom_model_label = ipw.Label(
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
        
        self.sort_atom_model_hbox = ipw.HBox(
            children = [
                self.sort_atom_model_label,
                self.sort_name_checkbox,
                self.sort_name_label,
                self.sort_registration_date_checkbox,
                self.sort_registration_date_label
            ]
        )
        
        self.create_atom_model_button = ipw.Button(
            tooltip = 'Add', 
            icon = 'plus', 
            layout = ipw.Layout(width = '50px', height = '25px')
        )
        
        self.atom_model_hbox = ipw.HBox(
            children = [
                self.atom_model_label,
                self.atom_model_dropdown,
                self.create_atom_model_button
            ]
        )
        
        self.create_new_atom_model_widgets = ipw.VBox()
        
        self.children = [
            self.atom_model_hbox,
            self.sort_atom_model_hbox,
            self.create_new_atom_model_widgets
        ]
        
        self.sort_name_checkbox.observe(self.sort_atom_model_dropdown, names = "value")
        self.sort_registration_date_checkbox.observe(self.sort_atom_model_dropdown, names = "value")
        self.create_atom_model_button.on_click(self.create_atom_model)
    
    def load_atom_models(self):
        atom_models = utils.get_openbis_objects(
            self.openbis_session,
            type = OPENBIS_OBJECT_TYPES["Atomistic Model"]
        )
        atom_model_options = [(f"{obj.props['name']}", obj.permId) for obj in atom_models]
        atom_model_options.insert(0, ("Select atomistic model...", "-1"))
        self.atom_model_dropdown.options = atom_model_options
        self.atom_model_dropdown.value = "-1"
    
    def sort_atom_model_dropdown(self, change):
        options = self.atom_model_dropdown.options[1:]
        
        df = pd.DataFrame(options, columns=["name", "registration_date"])
        if self.sort_name_checkbox.value and not self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="registration_date", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

        options = list(df.itertuples(index=False, name=None))
        options.insert(0, self.atom_model_dropdown.options[0])
        self.atom_model_dropdown.options = options
    
    def create_atom_model(self, b):
        select_molecules_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select molecules</span>"
        )
        
        select_reacprod_concepts_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select reaction product concepts</span>"
        )
        
        select_slab_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Select slab</span>"
        )
        
        molecules_accordion = ipw.Accordion()
        add_molecule_button = ipw.Button(
            description = 'Add', disabled = False, 
            button_style = 'success', tooltip = 'Add molecule', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        reacprod_concepts_accordion = ipw.Accordion()
        add_reacprod_concept_button = ipw.Button(
            description = 'Add', disabled = False, 
            button_style = 'success', tooltip = 'Add reaction product concept', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        material_type_options = MATERIALS_CONCEPTS_LABELS.copy()
        material_type_options.insert(0, ("Select material type...", "-1"))
        
        material_type_dropdown = ipw.Dropdown(
            options = material_type_options,
            value = material_type_options[0][1]
        )
        
        material_details_vbox = ipw.VBox()
        
        atom_model_props_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Atomistic model properties</span>"
        )
        
        name_label = ipw.Label(value = "Name")
        name_textbox = ipw.Text()
        name_hbox = ipw.HBox([name_label, name_textbox])
        
        wfms_uuid_label = ipw.Label(value = "WFMS UUID")
        wfms_uuid_textbox = ipw.Text()
        wfms_uuid_hbox = ipw.HBox([wfms_uuid_label, wfms_uuid_textbox])
        
        cell_label = ipw.Label(value = "Cell")
        cell_textbox = ipw.Textarea(placeholder = "{[[1,1,1], [2,2,2], [3,3,3]]}")
        cell_hbox = ipw.HBox([cell_label, cell_textbox])
        
        dimensionality_label = ipw.Label(value = "Dimensionality")
        dimensionality_intbox = ipw.IntText()
        dimensionality_hbox = ipw.HBox([dimensionality_label, dimensionality_intbox])
        
        pbc_label = ipw.Label(value = "PBC")
        pbc_x_checkbox = ipw.Checkbox(indent = False, layout = ipw.Layout(width = "20px"))
        pbc_y_checkbox = ipw.Checkbox(indent = False, layout = ipw.Layout(width = "20px"))
        pbc_z_checkbox = ipw.Checkbox(indent = False, layout = ipw.Layout(width = "20px"))
        pbc_hbox = ipw.HBox([pbc_label, pbc_x_checkbox, pbc_y_checkbox, pbc_z_checkbox])
        
        volume_label = ipw.Label(value = "Volume")
        volume_floatbox = ipw.FloatText()
        volume_hbox = ipw.HBox([volume_label, volume_floatbox])
        
        comments_label = ipw.Label(value = "Comments")
        comments_textbox = ipw.Textarea()
        comments_hbox = ipw.HBox([comments_label, comments_textbox])
        
        atom_model_preview_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Upload image preview</span>"
        )
        
        atom_model_preview_uploader = ipw.FileUpload(multiple = False)
        
        atom_model_datasets_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 18px;'>Upload datasets</span>"
        )
        
        atom_model_datasets_uploader = ipw.FileUpload()
        
        atom_model_properties_widgets = ipw.VBox(
            children = [
                atom_model_props_title,
                name_hbox,
                wfms_uuid_hbox,
                cell_hbox,
                dimensionality_hbox,
                pbc_hbox,
                volume_hbox,
                comments_hbox
            ]
        )
        
        atom_model_datasets_widgets = ipw.VBox(
            children = [
                atom_model_preview_title,
                atom_model_preview_uploader,
                atom_model_datasets_title,
                atom_model_datasets_uploader
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
        
        def load_material_type_widgets(change):
            if material_type_dropdown.value == "-1":
                material_details_vbox.children = []
                return
            else:
                material_options = [("Select material...", "-1")]
                material_dropdown = ipw.Dropdown(options = material_options, value = material_options[0][1])
                
                sort_material_label = ipw.Label(
                    value = "Sort by:", 
                    layout=ipw.Layout(margin='0px', width='50px'),
                    style = {'description_width': 'initial'}
                )
                
                name_checkbox = ipw.Checkbox(indent = False, layout=ipw.Layout(margin='2px', width='20px'))
                
                name_label = ipw.Label(
                    value = "Name", 
                    layout=ipw.Layout(margin='0px', width='50px'),
                    style = {'description_width': 'initial'}
                )
                
                registration_date_checkbox = ipw.Checkbox(indent = False,layout=ipw.Layout(margin='2px', width='20px'))
                
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
                
                material_details_vbox.children = [
                    select_material_box,
                    material_details_html
                ]
                
                material_type = material_type_dropdown.value
                material_objects = utils.get_openbis_objects(
                    self.openbis_session,
                    type = material_type
                )
                materials_objects_names_permids = [(obj.props["name"], obj.permId) for obj in material_objects]
                material_options += materials_objects_names_permids
                material_dropdown.options = material_options
                
                def sort_material_dropdown(change):
                    options = material_options[1:]
                    
                    df = pd.DataFrame(options, columns=["name", "registration_date"])
                    if name_checkbox.value and not registration_date_checkbox.value:
                        df = df.sort_values(by="name", ascending=True)
                    elif not name_checkbox.value and registration_date_checkbox.value:
                        df = df.sort_values(by="registration_date", ascending=False)
                    elif name_checkbox.value and registration_date_checkbox.value:
                        df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

                    options = list(df.itertuples(index=False, name=None))
                    options.insert(0, material_options[0])
                    material_dropdown.options = options
                
                def load_material_details(change):
                    obj_permid = material_dropdown.value
                    if obj_permid == "-1":
                        return
                    else:
                        obj = utils.get_openbis_object(self.openbis_session, sample_ident = obj_permid)
                        obj_props = obj.props.all()
                        obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
                        for key, value in obj_props.items():
                            if value:
                                prop_type = utils.get_openbis_property_type(self.openbis_session, code = key)
                                prop_label = prop_type.label
                                obj_details_string += f"<p><b>{prop_label}:</b> {value}</p>"
                        
                        obj_details_string += "</div>"
                        material_details_html.value = obj_details_string
                        
                name_checkbox.observe(sort_material_dropdown, names = "value")
                registration_date_checkbox.observe(sort_material_dropdown, names = "value")
                material_dropdown.observe(load_material_details, names = "value")
        
        def add_molecule(b):
            molecules_accordion_children = list(molecules_accordion.children)
            molecule_index = len(molecules_accordion_children)
            molecule_widget = MoleculeWidget(self.openbis_session, molecules_accordion, molecule_index)
            molecules_accordion_children.append(molecule_widget)
            molecules_accordion.children = molecules_accordion_children
        
        def add_reacprod_concept(b):
            reacprod_concepts_accordion_children = list(reacprod_concepts_accordion.children)
            reacprod_concept_index = len(reacprod_concepts_accordion_children)
            reacprod_concept_widget = ReacProdConceptWidget(self.openbis_session, reacprod_concepts_accordion, reacprod_concept_index)
            reacprod_concepts_accordion_children.append(reacprod_concept_widget)
            reacprod_concepts_accordion.children = reacprod_concepts_accordion_children
        
        def save_new_atom_model(b):
            atom_model_type = OPENBIS_OBJECT_TYPES["Atomistic Model"]
            atom_model_type_lower = atom_model_type.lower()
            atom_models_objs = utils.get_openbis_objects(self.openbis_session, type = atom_model_type)
            wfms_uuid = wfms_uuid_textbox.value
            atom_models_uuids = [obj.props["wfms_uuid"] for obj in atom_models_objs]
            if wfms_uuid in atom_models_uuids:
                display(Javascript(data = "alert('Atomistic model already in openBIS!')"))
            else:
                pbc = [pbc_x_checkbox.value, pbc_y_checkbox.value, pbc_z_checkbox.value]
                cell_json = cell_textbox.value
                if utils.is_valid_json(cell_json) == False:
                    cell_json = ""
                    
                atom_model_props = {
                    "name": name_textbox.value,
                    "wfms_uuid": wfms_uuid,
                    "cell": cell_json,
                    "dimensionality": dimensionality_intbox.value,
                    "periodic_boundary_conditions": pbc,
                    "volume": volume_floatbox.value,
                    "comments": comments_textbox.value
                }
                
                selected_molecules_widgets = molecules_accordion.children
                selected_reac_prods_widgets = reacprod_concepts_accordion.children
                
                # Get material
                if material_type_dropdown.value == "-1":
                    selected_slab = []
                else:
                    selected_material_id = material_details_vbox.children[0].children[0].value
                    if selected_material_id == "-1":
                        selected_slab = []
                    else:
                        selected_slab = [selected_material_id]
                
                # Get molecules
                selected_molecules_ids = []
                for mol_widget in selected_molecules_widgets:
                    mol_id = mol_widget.dropdown.value
                    if mol_id == "-1":
                        continue
                    else:
                        selected_molecules_ids.append(mol_id)
                
                # Get reaction products concepts
                selected_reac_prods_ids = []
                for reac_prod_widget in selected_reac_prods_widgets:
                    reac_prod_id = reac_prod_widget.dropdown.value
                    if reac_prod_id == "-1":
                        continue
                    else:
                        selected_reac_prods_ids.append(reac_prod_id)
                
                atom_model_parents = selected_slab + selected_molecules_ids + selected_reac_prods_ids
                
                atom_model_obj = utils.create_openbis_object(
                    self.openbis_session, 
                    type = atom_model_type,
                    collection = OPENBIS_COLLECTIONS_PATHS["Atomistic Model"],
                    props = atom_model_props,
                    parents = atom_model_parents
                )
                
                # Atomistic model preview
                for filename in atom_model_preview_uploader.value:
                    file_info = atom_model_preview_uploader.value[filename]
                    utils.write_file(file_info['content'], filename)
                    utils.create_openbis_dataset(
                        self.openbis_session,
                        type = "ELN_PREVIEW", 
                        sample = atom_model_obj, 
                        files = [filename]
                    )
                    os.remove(filename)
                
                # Atomistic model datasets
                for filename in atom_model_datasets_uploader.value:
                    file_info = atom_model_datasets_uploader.value[filename]
                    utils.write_file(file_info['content'], filename)
                    utils.create_openbis_dataset(
                        self.openbis_session,
                        type = "ATTACHMENT", 
                        sample = atom_model_obj, 
                        files = [filename]
                    )
                    os.remove(filename)
                
                self.create_new_atom_model_widgets.children = []
                self.atom_model_dropdown.options = self.load_atom_models()
                display(Javascript(data = "alert('Atomistic model successfully created!')"))
        
        def cancel_new_atom_model(b):
            self.create_new_atom_model_widgets.children = []
        
        save_button.on_click(save_new_atom_model)
        cancel_button.on_click(cancel_new_atom_model)
        material_type_dropdown.observe(load_material_type_widgets, names = "value")
        add_molecule_button.on_click(add_molecule)
        add_reacprod_concept_button.on_click(add_reacprod_concept)

        self.create_new_atom_model_widgets.children = [
            select_molecules_title,
            molecules_accordion,
            add_molecule_button,
            select_reacprod_concepts_title,
            reacprod_concepts_accordion,
            add_reacprod_concept_button,
            select_slab_title,
            material_type_dropdown,
            material_details_vbox,
            atom_model_properties_widgets,
            atom_model_datasets_widgets,
            buttons_hbox
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
            collection = OPENBIS_COLLECTIONS_PATHS["Precursor Molecule"],
            type = OPENBIS_OBJECT_TYPES["Molecule"]
        )
        dropdown_list = [(obj.props["name"], obj.permId) for obj in molecules_objects]
        dropdown_list.insert(0, ("Select a molecule...", "-1"))
        self.dropdown = ipw.Dropdown(value = "-1",options = dropdown_list)
        self.details_vbox = ipw.VBox()
        
        self.remove_molecule_button = ipw.Button(
            description = 'Remove', disabled = False, 
            button_style = 'danger', tooltip = 'Remove molecule', 
            layout = ipw.Layout(width = '150px', height = '25px')
        )
        
        self.molecule_sketch = ipw.Image(layout = ipw.Layout(width = '300px', height = '300px'))
        
        self.dropdown.observe(self.load_details, names = "value")
        self.remove_molecule_button.on_click(self.remove_molecule)
        self.children = [self.dropdown, self.details_vbox, self.molecule_sketch, self.remove_molecule_button]
    
    def load_details(self, change):
        obj_permid = self.dropdown.value
        if obj_permid == "-1":
            return
        else:
            obj = get_cached_object(self.openbis_session, obj_permid)
            obj_datasets = obj.get_datasets(type = "ELN_PREVIEW")
            obj_props = obj.props.all()
            obj_name = obj_props.get("name", "")
            self.parent_accordion.set_title(self.object_index, obj_name)
            self.title = obj_name
            
            obj_details_html = ipw.HTML()
            obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
            for key, value in obj_props.items():
                if value:
                    prop_type = utils.get_openbis_property_type(self.openbis_session, code = key)
                    prop_label = prop_type.label
                    obj_details_string += f"<p><b>{prop_label}:</b> {value}</p>"
            
            obj_details_string += "</div>"
            obj_details_html.value = obj_details_string
            
            if obj_datasets:
                object_dataset = obj_datasets[0]
                object_dataset.download(destination="images")
                object_image_filepath = object_dataset.file_list[0]
                self.molecule_sketch.value = utils.read_file(f"images/{object_dataset.permId}/{object_image_filepath}")
                
                # Erase file after downloading it
                shutil.rmtree(f"images/{object_dataset.permId}")
            
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
            collection = OPENBIS_COLLECTIONS_PATHS["Reaction Product"],
            type = OPENBIS_OBJECT_TYPES["Reaction Product Concept"]
        )
        dropdown_list = [(obj.props["name"], obj.permId) for obj in molecules_objects]
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
            obj_name = obj_props.get("name", "")
            self.parent_accordion.set_title(self.object_index, obj_name)
            self.title = obj_name
            
            obj_details_html = ipw.HTML()
            obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
            for key, value in obj_props.items():
                if value:
                    prop_type = utils.get_openbis_property_type(self.openbis_session, code = key)
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

# Create analysis/results/drafts widgets
class CreateDraftsWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_project_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select project</span>"
        )
        
        self.select_project_widget = SelectProjectWidget(self.openbis_session)
        
        self.select_results_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select results</span>"
        )
        
        self.select_results_label = ipw.Label(value = "Results")
        self.select_results_selector = ipw.SelectMultiple()
        self.select_results_hbox = ipw.HBox([self.select_results_label, self.select_results_selector])
        
        self.drafts_properties_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Properties</span>"
        )
        
        self.name_label = ipw.Label(value = "Name")
        self.name_textbox = ipw.Text()
        self.name_hbox = ipw.HBox([self.name_label, self.name_textbox])
        
        self.draft_type_label = ipw.Label(value = "Draft type")
        self.draft_type_dropdown = ipw.Dropdown(value = "Preprint", options = ["Preprint", "Postprint"])
        self.draft_type_hbox = ipw.HBox([self.draft_type_label, self.draft_type_dropdown])
        
        self.description_label = ipw.Label(value = "Description")
        self.description_textbox = ipw.Textarea()
        self.description_hbox = ipw.HBox([self.description_label, self.description_textbox])
        
        self.comments_label = ipw.Label(value = "Comments")
        self.comments_textbox = ipw.Textarea()
        self.comments_hbox = ipw.HBox([self.comments_label, self.comments_textbox])
        
        self.support_files_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Support files</span>"
        )
        self.support_files_uploader = ipw.FileUpload()
        
        self.select_project_widget.project_dropdown.observe(self.load_results)
        
        self.children = [
            self.select_project_title,
            self.select_project_widget,
            self.select_results_title,
            self.select_results_hbox,
            self.drafts_properties_title,
            self.name_hbox,
            self.draft_type_hbox,
            self.description_hbox,
            self.comments_hbox,
            self.support_files_title,
            self.support_files_uploader
        ]
    
    def load_results(self, change):
        project_id = self.select_project_widget.project_dropdown.value
        results_objects = []
        
        if project_id != "-1":
            results = utils.get_openbis_objects(
                self.openbis_session,
                type = OPENBIS_OBJECT_TYPES["Result"],
                project = project_id
            )
            
            results_objects = [(obj.props["name"], obj.permId) for obj in results]
            
        self.select_results_selector.options = results_objects

class CreateAnalysisWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_project_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select project</span>"
        )
        
        self.select_project_widget = SelectProjectWidget(self.openbis_session)
        
        self.select_simulations_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select simulations</span>"
        )
        
        self.select_simulations_label = ipw.Label(value = "Simulations")
        self.select_simulations_selector = ipw.SelectMultiple()
        self.select_simulations_hbox = ipw.HBox([self.select_simulations_label, self.select_simulations_selector])
        
        self.select_measurements_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select measurements</span>"
        )
        
        self.select_measurements_label = ipw.Label(value = "Measurements")
        self.select_measurements_selector = ipw.SelectMultiple()
        self.select_measurements_hbox = ipw.HBox([self.select_measurements_label, self.select_measurements_selector])
        
        self.select_software_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select software</span>"
        )
        
        self.select_software_label = ipw.Label(value = "Software")
        openbis_software = utils.get_openbis_objects(self.openbis_session, type = "SOFTWARE")
        software_options = [(obj.props["name"], obj.permId) for obj in openbis_software]
        self.select_software_selector = ipw.SelectMultiple(options = software_options)
        self.select_software_hbox = ipw.HBox([self.select_software_label, self.select_software_selector])
        
        self.select_code_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select code</span>"
        )
        
        self.select_code_label = ipw.Label(value = "Code")
        openbis_code = utils.get_openbis_objects(self.openbis_session, type = "CODE")
        code_options = [(obj.props["name"], obj.permId) for obj in openbis_code]
        self.select_code_selector = ipw.SelectMultiple(options = code_options)
        self.select_code_hbox = ipw.HBox([self.select_code_label, self.select_code_selector])
        
        self.analysis_properties_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Properties</span>"
        )
        
        self.name_label = ipw.Label(value = "Name")
        self.name_textbox = ipw.Text()
        self.name_hbox = ipw.HBox([self.name_label, self.name_textbox])
        
        self.description_label = ipw.Label(value = "Description")
        self.description_textbox = ipw.Textarea()
        self.description_hbox = ipw.HBox([self.description_label, self.description_textbox])
        
        self.comments_label = ipw.Label(value = "Comments")
        self.comments_textbox = ipw.Textarea()
        self.comments_hbox = ipw.HBox([self.comments_label, self.comments_textbox])
        
        self.support_files_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Support files</span>"
        )
        self.support_files_uploader = ipw.FileUpload()
        
        self.select_project_widget.project_dropdown.observe(self.load_measurements_and_simulations)
        
        self.children = [
            self.select_project_title,
            self.select_project_widget,
            self.select_simulations_title,
            self.select_simulations_hbox,
            self.select_measurements_title,
            self.select_measurements_hbox,
            self.select_software_title,
            self.select_software_hbox,
            self.select_code_title,
            self.select_code_hbox,
            self.analysis_properties_title,
            self.name_hbox,
            self.description_hbox,
            self.comments_hbox,
            self.support_files_title,
            self.support_files_uploader
        ]
    
    def load_measurements_and_simulations(self, change):
        project_id = self.select_project_widget.project_dropdown.value
        simulations_objects = []
        measurements_objects = []
        
        if project_id != "-1":
            for _, simulation_type in SIMULATION_TYPES:
                openbis_objs = utils.get_openbis_objects(
                    self.openbis_session,
                    type = simulation_type,
                    project = project_id
                )
                objs = [(obj.props["name"], obj.permId) for obj in openbis_objs]
                simulations_objects += objs
            
            measurements = utils.get_openbis_objects(
                self.openbis_session,
                type = OPENBIS_OBJECT_TYPES["Measurement Session"],
                project = project_id
            )
            
            for measurement in measurements:
                measurement_details = (measurement.props["name"], measurement.permId)
                if measurement.props["wfms_uuid"]:
                    simulations_objects.append(measurement_details)
                else:
                    measurements_objects.append(measurement_details)
            
        self.select_simulations_selector.options = simulations_objects
        self.select_measurements_selector.options = measurements_objects     

class CreateResultsWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_project_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select project</span>"
        )
        
        self.select_project_widget = SelectProjectWidget(self.openbis_session)
        
        self.select_simulations_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select simulations</span>"
        )
        
        self.select_simulations_label = ipw.Label(value = "Simulations")
        self.select_simulations_selector = ipw.SelectMultiple()
        self.select_simulations_hbox = ipw.HBox([self.select_simulations_label, self.select_simulations_selector])
        
        self.select_measurements_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select measurements</span>"
        )
        
        self.select_measurements_label = ipw.Label(value = "Measurements")
        self.select_measurements_selector = ipw.SelectMultiple()
        self.select_measurements_hbox = ipw.HBox([self.select_measurements_label, self.select_measurements_selector])
        
        self.select_analysis_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select analysis</span>"
        )
        self.select_analysis_label = ipw.Label(value = "Analysis")
        self.select_analysis_selector = ipw.SelectMultiple()
        self.select_analysis_hbox = ipw.HBox([self.select_analysis_label, self.select_analysis_selector])
        
        self.results_properties_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Properties</span>"
        )
        
        self.name_label = ipw.Label(value = "Name")
        self.name_textbox = ipw.Text()
        self.name_hbox = ipw.HBox([self.name_label, self.name_textbox])
        
        self.description_label = ipw.Label(value = "Description")
        self.description_textbox = ipw.Textarea()
        self.description_hbox = ipw.HBox([self.description_label, self.description_textbox])
        
        self.comments_label = ipw.Label(value = "Comments")
        self.comments_textbox = ipw.Textarea()
        self.comments_hbox = ipw.HBox([self.comments_label, self.comments_textbox])
        
        self.support_files_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Support files</span>"
        )
        self.support_files_uploader = ipw.FileUpload()
        
        self.select_project_widget.project_dropdown.observe(self.load_analysis_measurements_and_simulations)
        
        self.children = [
            self.select_project_title,
            self.select_project_widget,
            self.select_simulations_title,
            self.select_simulations_hbox,
            self.select_measurements_title,
            self.select_measurements_hbox,
            self.select_analysis_title,
            self.select_analysis_hbox,
            self.results_properties_title,
            self.name_hbox,
            self.description_hbox,
            self.comments_hbox,
            self.support_files_title,
            self.support_files_uploader
        ]
    
    def load_analysis_measurements_and_simulations(self, change):
        project_id = self.select_project_widget.project_dropdown.value
        analysis_objects = []
        simulations_objects = []
        measurements_objects = []
        
        if project_id != "-1":
            for _, simulation_type in SIMULATION_TYPES:
                openbis_objs = utils.get_openbis_objects(
                    self.openbis_session,
                    type = simulation_type,
                    project = project_id
                )
                objs = [(obj.props["name"], obj.permId) for obj in openbis_objs]
                simulations_objects += objs
            
            measurements = utils.get_openbis_objects(
                self.openbis_session,
                type = OPENBIS_OBJECT_TYPES["Measurement Session"],
                project = project_id
            )
            
            for measurement in measurements:
                measurement_details = (measurement.props["name"], measurement.permId)
                if measurement.props["wfms_uuid"]:
                    simulations_objects.append(measurement_details)
                else:
                    measurements_objects.append(measurement_details)
            
            analysis = utils.get_openbis_objects(
                self.openbis_session,
                type = OPENBIS_OBJECT_TYPES["Analysis"],
                project = project_id
            )
            
            analysis_objects = [(obj.props["name"], obj.permId) for obj in analysis]
            
        self.select_simulations_selector.options = simulations_objects
        self.select_measurements_selector.options = measurements_objects
        self.select_analysis_selector.options = analysis_objects  

class SelectProjectWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.project_label = ipw.Label(
            value = "Project"
        )
        
        self.project_dropdown = ipw.Dropdown(
            layout=ipw.Layout(width='500px')
        )
        self.load_projects()
        
        self.sort_project_label = ipw.Label(
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
        
        self.sort_project_widgets = ipw.HBox(
            children = [
                self.sort_project_label,
                self.sort_name_checkbox,
                self.sort_name_label,
                self.sort_registration_date_checkbox,
                self.sort_registration_date_label
            ]
        )
        
        self.project_dropdown_boxes = ipw.HBox(
            children = [
                self.project_label,
                self.project_dropdown
            ]
        )
        
        self.children = [
            self.project_dropdown_boxes,
            self.sort_project_widgets
        ]
        
        self.sort_name_checkbox.observe(self.sort_project_dropdown, names = "value")
        self.sort_registration_date_checkbox.observe(self.sort_project_dropdown, names = "value")
    
    def load_projects(self):
        projects = utils.get_openbis_projects(self.openbis_session)
        project_options = []
        for prj in projects:
            prj_option = (f"{prj.code} from Space {prj.space.code}", prj.permId)
            project_options.append(prj_option)
            
        project_options.insert(0, ("Select project...", "-1"))
        self.project_dropdown.options = project_options
        self.project_dropdown.value = "-1"

    def sort_project_dropdown(self, change):
        options = self.project_dropdown.options[1:]
        
        df = pd.DataFrame(options, columns=["name", "registration_date"])
        if self.sort_name_checkbox.value and not self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="registration_date", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

        options = list(df.itertuples(index=False, name=None))
        options.insert(0, self.project_dropdown.options[0])
        self.project_dropdown.options = options

# Sample preparation widgets
class CreateSampleWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_material_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select material</span>"
        )
        
        material_type_options = MATERIALS_LABELS.copy()
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
            materials_objects_names_permids = [(obj.props["name"], obj.permId) for obj in material_objects]
            material_options += materials_objects_names_permids
            material_dropdown.options = material_options
            
            def sort_material_dropdown(change):
                options = material_options[1:]
                
                df = pd.DataFrame(options, columns=["name", "registration_date"])
                if name_checkbox.value and not registration_date_checkbox.value:
                    df = df.sort_values(by="name", ascending=True)
                elif not name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by="registration_date", ascending=False)
                elif name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

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
                    obj_name = obj_props.get("name", "")
                    obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
                    for key, value in obj_props.items():
                        if value:
                            prop_type = utils.get_openbis_property_type(self.openbis_session, code = key)
                            prop_label = prop_type.label
                            prop_datatype = prop_type.dataType
                            if prop_datatype == OPENBIS_OBJECT_TYPES["Sample"]:
                                if isinstance(value, list):
                                    prop_obj_names = []
                                    for id in value:
                                        prop_obj = get_cached_object(self.openbis_session, id)
                                        prop_obj_name = prop_obj.props["name"]
                                        prop_obj_names.append(prop_obj_name)
                                    value = ", ".join(prop_obj_names)
                                else:
                                    obj = get_cached_object(self.openbis_session, value)
                                    value = obj.props["name"]
                                    
                            elif prop_datatype == "JSON":
                                json_content = json.loads(value)
                                if is_quantity_value(json_content):
                                    value = f"<p>{json_content['value']} {json_content['unit']}</p>"
                                else:
                                    value = "<ul>"
                                    for k, v in json_content.items():
                                        if isinstance(v, dict):
                                            if is_quantity_value(v):
                                                value += f"<li><b>{k}:</b> {v['value']} {v['unit']}</li>"
                                            else:
                                                value += f"<li><b>{k}:</b> {v}</li>"
                                        else:
                                            value += f"<li><b>{k}:</b> {v}</li>"
                                    
                                    value += "</ul>"
                            
                            elif prop_datatype == "XML" and prop_type.metaData["custom_widget"] == "Spreadsheet":
                                table_headers = value.headers
                                table_data = value.data
                                
                                # Build table header
                                table_html = "<table style='width:100%; border-collapse:collapse;'>"
                                table_html += "<thead><tr>"
                                for h in table_headers:
                                    table_html += f"<th style='padding:0; text-align:left; font-weight:bold;'>{h}</th>"
                                table_html += "</tr></thead>"

                                # Build table body
                                table_html += "<tbody>"
                                for row in table_data:
                                    table_html += "<tr>"
                                    for cell in row:
                                        table_html += f"<td style='padding:0;'>{cell}</td>"
                                    table_html += "</tr>"
                                table_html += "</tbody></table>"
                                value = table_html
                                
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
                sample_type = OPENBIS_OBJECT_TYPES["Sample"]
                sample_props = {
                    "name": sample_name,
                    "exists": True
                }
                utils.create_openbis_object(
                    self.openbis_session,
                    type = sample_type,
                    collection = OPENBIS_COLLECTIONS_PATHS["Sample"],
                    props = sample_props,
                    parents = [material_object]
                )
                
                display(Javascript(data = "alert('Sample created successfully!')"))
                
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
            if parent_type == OPENBIS_OBJECT_TYPES["Process Step"]:
                if most_recent_parent:
                    if parent_object.registrationDate > most_recent_parent.registrationDate:
                        most_recent_parent = parent_object
                else:
                    most_recent_parent = parent_object
        
        if most_recent_parent:
            if most_recent_parent.experiment.permId != self.select_experiment_dropdown.experiment_dropdown.value:
                self.select_experiment_dropdown.experiment_dropdown.value = most_recent_parent.experiment.permId
                display(Javascript(data = "alert('Experiment was changed!')"))
            
            for parent in most_recent_parent.parents:
                parent_object = get_cached_object(self.openbis_session, parent)

                if parent_object.type == OPENBIS_OBJECT_TYPES["Preparation"]:
                    self.sample_preparation_object = parent_object
                    break

        # If sample was used in a measurement session, a new preparation should start
        sample_object_children = sample_object.children
        for child_id in sample_object_children:
            child_object = get_cached_object(self.openbis_session, child_id)
            if child_object.type == OPENBIS_OBJECT_TYPES["Measurement Session"]:
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
        processes_options = [(obj.props["name"], obj.permId) for obj in openbis_processes]
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
            process_type = OPENBIS_OBJECT_TYPES["Process"]
            process_type_lower = process_type.lower()
            process_step_type = OPENBIS_OBJECT_TYPES["Process Step"]
            process_step_type_lower = process_step_type.lower()
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
        experiment_id = self.select_experiment_dropdown.experiment_dropdown.value
        if experiment_id == "-1":
            display(Javascript(data = "alert('Select an experiment.')"))
            return
        
        processes_widgets = self.new_processes_accordion.children
        if processes_widgets:
            openbis_transaction_objects = []
            experiment_object = utils.get_openbis_collection(
                self.openbis_session,
                code = experiment_id
            )
            experiment_project_code = experiment_object.project.identifier
            
            current_sample_id = self.select_sample_dropdown.sample_dropdown.value
            current_sample = get_cached_object(self.openbis_session, current_sample_id)
            
            # Create preparation object when it does not exist
            if self.sample_preparation_object is None:
                self.sample_preparation_object = utils.create_openbis_object(
                    self.openbis_session,
                    type = OPENBIS_OBJECT_TYPES["Preparation"],
                    experiment = experiment_object.identifier,
                    props = {"name": current_sample.props["name"]}
                )
            
            sample_preparation_id = self.sample_preparation_object.permId
            for process_widget in processes_widgets:
                # Reload sample preparation object to load children that was added in the cycle (e.g. process steps)
                self.sample_preparation_object = utils.get_openbis_object(self.openbis_session, sample_ident = sample_preparation_id)
                sample_type = OPENBIS_OBJECT_TYPES["Sample"]
                sample_type_lower = sample_type.lower()
                process_code = ""
                current_sample.props["exists"] = False
                current_sample_name = current_sample.props["name"]
                current_sample.save()
                
                process_step_type = OPENBIS_OBJECT_TYPES["Process Step"]
                process_step_type_lower = process_step_type.lower()
                new_process_object = utils.create_openbis_object(
                    self.openbis_session,
                    type = process_step_type,
                    experiment = experiment_object.identifier
                )
                
                process_properties = {
                    "name": process_widget.name_textbox.value,
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
                        action_type_lower = action_type.lower()
                        duration_days = action_widget.duration_days_intbox.value
                        duration_hours = action_widget.duration_hours_intbox.value
                        duration_minutes = action_widget.duration_minutes_intbox.value
                        duration_seconds = action_widget.duration_seconds_intbox.value

                        if action_type == OPENBIS_OBJECT_TYPES["Annealing"]:
                            target_temperature = {
                                "value": action_widget.target_temperature_value_textbox.value,
                                "unit": action_widget.target_temperature_unit_dropdown.value,
                            }
                            action_properties["target_temperature"] = json.dumps(target_temperature)

                        elif action_type == OPENBIS_OBJECT_TYPES["Cooldown"]:
                            target_temperature = {
                                "value": action_widget.target_temperature_value_textbox.value,
                                "unit": action_widget.target_temperature_unit_dropdown.value,
                            }
                            action_properties["target_temperature"] = json.dumps(target_temperature)
                            action_properties["cryogen"] = action_widget.cryogen_textbox.value

                        elif action_type == OPENBIS_OBJECT_TYPES["Deposition"]:
                            substrate_temperature = {
                                "value": action_widget.substrate_temperature_value_textbox.value,
                                "unit": action_widget.substrate_temperature_unit_dropdown.value,
                            }
                            substance = action_widget.substance_dropdown.value
                            if substance != "-1":
                                action_properties["substance"] = substance
                                
                            action_properties["substrate_temperature"] = json.dumps(substrate_temperature)

                        elif action_type == OPENBIS_OBJECT_TYPES["Dosing"]:
                            substrate_temperature = {
                                "value": action_widget.substrate_temperature_value_textbox.value,
                                "unit": action_widget.substrate_temperature_unit_dropdown.value,
                            }
                            pressure = {
                                "value": action_widget.pressure_value_textbox.value,
                                "unit": action_widget.pressure_unit_dropdown.value,
                            }
                            gas = action_widget.gas_dropdown.value
                            action_properties["gas"] = gas
                            action_properties["substrate_temperature"] = json.dumps(substrate_temperature)
                            action_properties["pressure"] = json.dumps(pressure)

                        elif action_type == OPENBIS_OBJECT_TYPES["Sputtering"]:
                            pressure = {
                                "value": action_widget.pressure_value_textbox.value,
                                "unit": action_widget.pressure_unit_dropdown.value,
                            }
                            current = {
                                "value": action_widget.current_value_textbox.value,
                                "unit": action_widget.current_unit_dropdown.value,
                            }
                            angle = {
                                "value": action_widget.angle_value_textbox.value,
                                "unit": action_widget.angle_unit_dropdown.value,
                            }
                            substrate_temperature = {
                                "value": action_widget.substrate_temperature_value_textbox.value,
                                "unit": action_widget.substrate_temperature_unit_dropdown.value,
                            }
                            action_properties["pressure"] = json.dumps(pressure)
                            action_properties["current"] = json.dumps(current)
                            action_properties["angle"] = json.dumps(angle)
                            action_properties["substrate_temperature"] = json.dumps(substrate_temperature)
                            action_properties["ion"] = action_widget.sputter_ion_textbox.value
                        
                        component_permid = action_widget.component_dropdown.value
                        if component_permid != "-1":
                            component_object = get_cached_object(self.openbis_session, component_permid)
                            component_type = component_object.type.code
                            component_type_lower = component_type.lower()
                            component_object_settings = component_object.props["actions_settings"]
                            
                            if component_object_settings:
                                component_settings = {}
                                component_object_settings = json.loads(component_object_settings)
                                action_component_settings = []
                                for action_settings in component_object_settings:
                                    if action_settings["action_type"] == action_type:
                                        action_component_settings = action_settings["component_properties"]
                                        break
                                    
                                for setting in action_component_settings:
                                    if setting == "target_temperature":
                                        component_settings["target_temperature"] = {
                                            "value": action_widget.target_temperature_value_comp_textbox.value,
                                            "unit": action_widget.target_temperature_unit_comp_dropdown.value,
                                        }
                                    elif setting == "bias_voltage":
                                        component_settings["bias_voltage"] = {
                                            "value": action_widget.bias_voltage_value_textbox.value,
                                            "unit": action_widget.bias_voltage_unit_dropdown.value,
                                        }
                                    elif setting == "discharge_voltage":
                                        component_settings["discharge_voltage"] = {
                                            "value": action_widget.discharge_voltage_value_textbox.value,
                                            "unit": action_widget.discharge_voltage_unit_dropdown.value,
                                        }
                                    elif setting == "discharge_current":
                                        component_settings["discharge_current"] = {
                                            "value": action_widget.discharge_current_value_textbox.value,
                                            "unit": action_widget.discharge_current_unit_dropdown.value,
                                        }
                                    elif setting == "p_value":
                                        component_settings["p_value"] = float(action_widget.evaporator_p_value_textbox.value)
                                        
                                    elif setting == "i_value":
                                        component_settings["i_value"] = float(action_widget.evaporator_i_value_textbox.value)
                                        
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

                            action_properties["name"] = action_widget.name_textbox.value
                            action_properties["description"] = action_widget.description_textbox.value
                            action_properties["duration"] = f"{duration_days} days {duration_hours:02}:{duration_minutes:02}:{duration_seconds:02}"
                            action_properties["component"] = component_permid
                            action_properties["comments"] = action_widget.comments_textarea.value
                            
                            action_collection_code = "ACTIONS_COLLECTION"
                            openbis_experiments = utils.get_openbis_collections(
                                self.openbis_session,
                                code = action_collection_code, 
                                project = experiment_project_code
                            )
                            
                            if openbis_experiments.df.empty:
                                actions_collection_object = utils.create_openbis_collection(
                                    self.openbis_session,
                                    type = "COLLECTION",
                                    code = action_collection_code,
                                    project = experiment_project_code,
                                    props = {"name": "Actions"}
                                )
                            
                            new_action_object = utils.create_openbis_object(
                                self.openbis_session,
                                type = action_type,
                                experiment = f"{experiment_project_code}/{action_collection_code}",
                                props = action_properties
                            )
                            
                            # Append action to list of actions
                            actions.append(new_action_object.permId)
                            
                            # Get action code
                            actions_codes.append(ACTIONS_CODES[action_type])
                
                if observables_widgets:
                    for observable_widget in observables_widgets:
                        observable_properties = {}
                        observable_type = observable_widget.observable_type_dropdown.value
                        observable_type_lower = observable_type.lower()
                        
                        component_permid = observable_widget.component_dropdown.value
                        if component_permid != "-1":
                            component_object = get_cached_object(self.openbis_session, component_permid)
                            component_type = component_object.type.code
                            component_type_lower = component_type.lower()
                            component_object_settings = component_object.props["observables_settings"]
                            
                            if component_object_settings:
                                component_settings = {}
                                component_object_settings = json.loads(component_object_settings)
                                observable_component_settings = []
                                for observable_settings in component_object_settings:
                                    if observable_settings["observable_type"] == observable_type:
                                        observable_component_settings = observable_settings["component_properties"]
                                        break
                                    
                                for setting in observable_component_settings:
                                    if setting == "density":
                                        component_settings["density"] = {
                                            "value": observable_widget.density_value_textbox.value,
                                            "unit": observable_widget.density_unit_dropdown.value,
                                        }
                                        
                                    elif setting == "filament":
                                        component_settings["filament"] = observable_widget.filament_textbox.value
                                        
                                    elif setting == "filament_current":
                                        component_settings["filament_current"] = {
                                            "value": observable_widget.filament_current_value_textbox.value,
                                            "unit": observable_widget.filament_current_unit_dropdown.value,
                                        }
                                
                                # Update current component settings
                                for setting_key, setting_value in component_settings.items():
                                    if isinstance(setting_value, dict):
                                        component_object.props[setting_key] = json.dumps(setting_value)
                                    else:
                                        component_object.props[setting_key] = setting_value
                                component_object.save()
                                
                                observable_properties["component_settings"] = json.dumps(component_settings)
                        
                            observable_properties["name"] = observable_widget.name_textbox.value
                            observable_properties["description"] = observable_widget.description_textbox.value
                            observable_properties["channel_name"] = observable_widget.ch_name_textbox.value
                            observable_properties["component"] = component_permid
                            observable_properties["comments"] = observable_widget.comments_textarea.value
                            
                            observable_collection_code = "OBSERVABLES_COLLECTION"
                            openbis_experiments = utils.get_openbis_collections(
                                self.openbis_session,
                                code = observable_collection_code, 
                                project = experiment_project_code
                            )
                            
                            if openbis_experiments.df.empty:
                                observables_collection_object = utils.create_openbis_collection(
                                    self.openbis_session,
                                    type = "COLLECTION",
                                    code = observable_collection_code,
                                    project = experiment_project_code,
                                    props = {"name": "Observables"}
                                )
                            
                            new_observable_object = utils.create_openbis_object(
                                self.openbis_session,
                                type = observable_type,
                                experiment = f"{experiment_project_code}/{observable_collection_code}",
                                props = observable_properties
                            )
                            
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
                self.sample_preparation_object.props["name"] = f"PREP_{new_sample_name}"
                self.sample_preparation_object.save()
                
                new_process_object.props = process_properties
                new_process_object.parents = [self.sample_preparation_object, current_sample]
                new_process_object.save()
                
                new_sample = utils.create_openbis_object(
                    self.openbis_session,
                    type = sample_type,
                    experiment = OPENBIS_COLLECTIONS_PATHS["Sample"],
                    parents = [new_process_object],
                    props = {"name": new_sample_name, "exists": True}
                )
                
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
        
        self.experiment_dropdown = ipw.Dropdown(
            layout=ipw.Layout(width='500px')
        )
        self.load_experiments()
        
        
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
        self.sort_name_checkbox.observe(self.sort_experiment_dropdown, names = "value")
        self.sort_registration_date_checkbox.observe(self.sort_experiment_dropdown, names = "value")
    
    def load_experiments(self):
        experiments = utils.get_openbis_collections(
            self.openbis_session,
            type = "EXPERIMENT"
        )
        experiment_options = []
        for exp in experiments:
            if "name" in exp.props.all():
                exp_option = (f"{exp.props['name']} from Project {exp.project.code} and Space {exp.project.space}", exp.permId)
            else:
                exp_option = (f"{exp.code} from Project {exp.project.code} and Space {exp.project.space}", exp.permId)
            experiment_options.append(exp_option)
        experiment_options.insert(0, ("Select experiment...", "-1"))
        self.experiment_dropdown.options = experiment_options
        self.experiment_dropdown.value = "-1"

    def sort_experiment_dropdown(self, change):
        options = self.experiment_dropdown.options[1:]
        
        df = pd.DataFrame(options, columns=["name", "registration_date"])
        if self.sort_name_checkbox.value and not self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="registration_date", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

        options = list(df.itertuples(index=False, name=None))
        options.insert(0, self.experiment_dropdown.options[0])
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
                "name": new_experiment_name_textbox.value,
                "default_collection_view": "IMAGING_GALLERY_VIEW"
            }
            
            project_id = project_dropdown.value
            
            if project_id == "-1":
                display(Javascript(data = "alert('Select a project.')"))
                return
            else:
                try:
                    utils.create_openbis_collection(
                        self.openbis_session,
                        type = "EXPERIMENT",
                        project = project_dropdown.value, 
                        props = experiment_props
                    )
                    self.create_new_experiment_widgets.children = []
                    self.experiment_dropdown.options = self.load_experiments()
                    display(Javascript(data = "alert('Experiment successfully created!')"))
                except ValueError as e:
                    display(Javascript(data = "alert('Error! Check if experiment already exists (either in ELN or Trash).')"))
        
        def cancel_new_experiment(b):
            self.create_new_experiment_widgets.children = []
            
        def sort_project_dropdown(change):
            options = project_dropdown_options[1:]
            
            df = pd.DataFrame(options, columns=["name", "registration_date"])
            if sort_project_name_checkbox.value and not sort_project_registration_date_checkbox.value:
                df = df.sort_values(by="name", ascending=True)
            elif not sort_project_name_checkbox.value and sort_project_registration_date_checkbox.value:
                df = df.sort_values(by="registration_date", ascending=False)
            elif sort_project_name_checkbox.value and sort_project_registration_date_checkbox.value:
                df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

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
        
        self.sort_name_checkbox.observe(self.sort_sample_dropdown, names = "value")
        self.sort_registration_date_checkbox.observe(self.sort_sample_dropdown, names = "value")
        
        self.children = [
            self.sample_dropdown_hbox,
            self.sort_sample_hbox
        ]
    
    def load_samples(self):
        sample_type = OPENBIS_OBJECT_TYPES["Sample"]
        sample_type_lower = sample_type.lower()
        samples = utils.get_openbis_objects(
            self.openbis_session,
            type = sample_type
        )
        sample_options = [(f"{obj.props['name']}", obj.permId) for obj in samples if obj.props["exists"] == "true"]
        sample_options.insert(0, ("Select sample...", "-1"))
        self.sample_dropdown.options = sample_options
        self.sample_dropdown.value = "-1"
    
    def sort_sample_dropdown(self, change):
        options = self.sample_dropdown.options[1:]
        
        df = pd.DataFrame(options, columns=["name", "registration_date"])
        if self.sort_name_checkbox.value and not self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="registration_date", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

        options = list(df.itertuples(index=False, name=None))
        options.insert(0, self.sample_dropdown.options[0])
        self.sample_dropdown.options = options

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
                if OPENBIS_OBJECT_CODES["Process Step"] in parent_code or OPENBIS_OBJECT_CODES["Sample"] in parent_code:
                    parent_object = get_cached_object(self.openbis_session, parent)
                    next_parents.extend(parent_object.parents)

                    if OPENBIS_OBJECT_CODES["Process Step"] in parent_code:
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
        self.process_step_type = OPENBIS_OBJECT_TYPES["Process Step"]
        self.process_step_type_lower = self.process_step_type.lower()
        
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
        if openbis_object_props["name"]:
            self.name_html.value = openbis_object_props["name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
            
        self.registration_date = self.openbis_object.registrationDate
        
        instrument_id = openbis_object_props["instrument"]
        
        if instrument_id:
            instrument_object = get_cached_object(self.openbis_session, instrument_id)
            self.instrument_html.value = instrument_object.props["name"]
        
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
        self.object_type = self.openbis_object.type.code
        self.object_type_lower = self.object_type.lower()
        
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
        
        if self.object_type == OPENBIS_OBJECT_TYPES["Annealing"]:
            self.target_temperature_label = ipw.Label(value = "Name:")
            self.target_temperature_html = ipw.HTML()
            self.target_temperature_hbox = ipw.HBox(children = [self.target_temperature_label, self.target_temperature_html])
            widget_children.append(self.target_temperature_hbox)

        elif self.object_type == OPENBIS_OBJECT_TYPES["Cooldown"]:
            self.target_temperature_label = ipw.Label(value = "Target temperature:")
            self.target_temperature_html = ipw.HTML()
            self.target_temperature_hbox = ipw.HBox(children = [self.target_temperature_label, self.target_temperature_html])
            
            self.cryogen_label = ipw.Label(value = "Cryogen:")
            self.cryogen_html = ipw.HTML()
            self.cryogen_hbox = ipw.HBox(children = [self.cryogen_label, self.cryogen_html])
            
            widget_children.append(self.cryogen_hbox)
            widget_children.append(self.target_temperature_hbox)

        elif self.object_type == OPENBIS_OBJECT_TYPES["Deposition"]:
            self.substrate_temperature_label = ipw.Label(value = "Substrate temperature:")
            self.substrate_temperature_html = ipw.HTML()
            self.substrate_temperature_hbox = ipw.HBox(children = [self.substrate_temperature_label, self.substrate_temperature_html])
            
            self.substance_label = ipw.Label(value = "Substance:")
            self.substance_html = ipw.HTML()
            self.substance_hbox = ipw.HBox(children = [self.substance_label, self.substance_html])
            
            widget_children.append(self.substance_hbox)
            widget_children.append(self.substrate_temperature_hbox)

        elif self.object_type == OPENBIS_OBJECT_TYPES["Dosing"]:
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

        elif self.object_type == OPENBIS_OBJECT_TYPES["Sputtering"]:
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
        if openbis_object_props["name"]:
            self.name_html.value = openbis_object_props["name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["duration"]:
            self.duration_html.value = openbis_object_props["duration"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
        
        if "target_temperature" in openbis_object_props:
            target_temperature = openbis_object_props["target_temperature"]
            if target_temperature:
                self.target_temperature_html.value = utils.stringify_quantity_value(target_temperature, "unit")
        
        if "cryogen" in openbis_object_props:
            if openbis_object_props["cryogen"]:
                self.cryogen_html.value = openbis_object_props["cryogen"]
        
        if "substrate_temperature" in openbis_object_props:
            substrate_temperature = openbis_object_props["substrate_temperature"]
            if substrate_temperature:
                self.substrate_temperature_html.value = utils.stringify_quantity_value(substrate_temperature, "unit")
        
        if "dosing_gas" in openbis_object_props:
            if openbis_object_props["gas"]:
                self.dosing_gas_html.value = openbis_object_props["gas"]
        
        if "pressure" in openbis_object_props:
            pressure = openbis_object_props["pressure"]
            if pressure:
                self.pressure_html.value = utils.stringify_quantity_value(pressure, "unit")
        
        if "current" in openbis_object_props:
            current = openbis_object_props["current"]
            if current:
                self.current_html.value = utils.stringify_quantity_value(current, "unit")
        
        if "angle" in openbis_object_props:
            angle = openbis_object_props["angle"]
            if angle:
                self.angle_html.value = utils.stringify_quantity_value(angle, "unit")

        if "substance" in openbis_object_props:
            substance_id = openbis_object_props["substance"]
            if substance_id:
                substance_object = get_cached_object(self.openbis_session, substance_id)
                self.substance_html.value = substance_object.props["name"]
        
        component_id = openbis_object_props["component"]
        if component_id:
            component_object = get_cached_object(self.openbis_session, component_id)
            self.component_html.value = component_object.props["name"]
            
        if openbis_object_props["component_settings"]:
            component_settings = json.loads(openbis_object_props["component_settings"])
            component_settings_string = ""
            for prop_key, prop_value in component_settings.items():
                prop_type = utils.get_openbis_property_type(self.openbis_session, code = prop_key)
                prop_label = prop_type.label
                
                # Convert quantity values from json objects to strings
                if prop_value and isinstance(prop_value, dict):
                    unit_type = ""
                    for key in prop_value.keys():
                        if key == "unit":
                            unit_type = key
                            break
                    prop_value = utils.stringify_quantity_value(prop_value, unit_type)
                    
                component_settings_string += f"<p>&bull; {prop_label}: {prop_value}</p>"
            
            self.component_settings_html.value = component_settings_string
        
class ObservableHistoryWidget(ipw.VBox):
    def __init__(self, openbis_session, openbis_object):
        super().__init__()
        self.openbis_session = openbis_session
        self.openbis_object = openbis_object
        self.observable_type = self.openbis_object.type.code
        self.observable_type_lower = self.observable_type.lower()
        
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
        if openbis_object_props["name"]:
            self.name_html.value = openbis_object_props["name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["channel_name"]:
            self.ch_name_html.value = openbis_object_props["channel_name"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
            
        component_id = openbis_object_props["component"]
            
        if component_id:
            component_object = get_cached_object(self.openbis_session, component_id)
            self.component_html.value = component_object.props["name"]
            
        if openbis_object_props["component_settings"]:
            component_settings = json.loads(openbis_object_props["component_settings"])
            component_settings_string = ""
            for prop_key, prop_value in component_settings.items():
                prop_type = utils.get_openbis_property_type(self.openbis_session, code = prop_key)
                prop_label = prop_type.label
                
                # Convert quantity values from json objects to strings
                if prop_value and isinstance(prop_value, dict):
                    unit_type = ""
                    for key in prop_value.keys():
                        if key == "unit":
                            unit_type = key
                            break
                    prop_value = utils.stringify_quantity_value(prop_value, unit_type)
                
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
        instrument_objects = utils.get_openbis_objects(self.openbis_session, collection = OPENBIS_COLLECTIONS_PATHS["Instrument"])
        instrument_options = [(obj.props["name"], obj.permId) for obj in instrument_objects]
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
        process_step_code = OPENBIS_OBJECT_TYPES["Process Step"]
        process_step_code_lower = process_step_code.lower()
        self.name_textbox.value = settings.get("name", "")
        self.description_textbox.value = settings.get("description", "")
        self.instrument_dropdown.value = settings.get("instrument", "-1")
        self.comments_textarea.value = settings.get("comments", "")
        actions_settings = settings.get("actions_settings", {})
        observables_settings = settings.get("observables_settings", {})
        for action_settings in actions_settings:
            actions_accordion_children = list(self.actions_accordion.children)
            action_index = len(actions_accordion_children)
            new_action_widget = RegisterActionWidget(
                self.openbis_session, 
                self.actions_accordion, 
                action_index, 
                self.instrument_dropdown.value,
                action_settings
            )
            actions_accordion_children.append(new_action_widget)
            self.actions_accordion.children = actions_accordion_children
        
        for observable_settings in observables_settings:
            observables_accordion_children = list(self.observables_accordion.children)
            observable_index = len(observables_accordion_children)
            new_observable_widget = RegisterObservableWidget(
                self.openbis_session, 
                self.observables_accordion, 
                observable_index, 
                self.instrument_dropdown.value,
                observable_settings
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
        action_type_options = ACTIONS_LABELS.copy()
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
            collection = OPENBIS_COLLECTIONS_PATHS["Precursor Substance"],
            type = OPENBIS_OBJECT_TYPES["Substance"]
        )
        
        substance_code = OPENBIS_OBJECT_TYPES["Substance"]
        substance_code_lower = substance_code.lower()
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
        gas_list = utils.get_openbis_objects(self.openbis_session, collection = OPENBIS_COLLECTIONS_PATHS["Chemical"])
        gas_options = [(obj.props["name"], obj.permId) for obj in gas_list]
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
        
        instrument_stm_type = OPENBIS_OBJECT_TYPES["Instrument STM"]
        instrument_stm_type_lower = instrument_stm_type.lower()
        self.instrument_type_components_dictionary = {
            OPENBIS_OBJECT_TYPES["Instrument STM"]: [
                "pumps", "gauges", 
                "vacuum_chambers", "ports_valves",
                "preparation_tools", "analysers", 
                "mechanical_components", "stm_components", 
                "control_data_acquisition", 
                "temperature_environment_control",
                "auxiliary_components", "tips_sensors", 
                "accessories"
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
        
        evaporator_slot_type_lower = OPENBIS_OBJECT_TYPES["Evaporator Slot"].lower()
        pbn_stage_type_lower = OPENBIS_OBJECT_TYPES["PBN Stage"].lower()
        sputter_gun_type_lower = OPENBIS_OBJECT_TYPES["Sputter Gun"].lower()
        self.components_properties_widgets = {
            "target_temperature": self.target_temperature_comp_hbox,
            "bias_voltage": self.bias_voltage_hbox,
            "discharge_voltage": self.discharge_voltage_hbox,
            "discharge_current": self.discharge_current_hbox,
            "p_value": self.evaporator_p_value_hbox,
            "i_value": self.evaporator_i_value_hbox,
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
        action_permid = settings["action"]
        action_object = get_cached_object(self.openbis_session, action_permid)
        action_props = action_object.props.all()
        action_type = action_object.type.code
        action_type_lower = action_type.lower()
        self.action_type_dropdown.value = action_type
        self.name_textbox.value = action_props["name"] or ""
        duration_str = action_props["duration"]
        
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
        
        self.description_textbox.value = action_props["description"] or ""
        action_target_temperature = action_props.get("target_temperature", "")
        action_cryogen = action_props.get("cryogen", "")
        action_substance = action_props.get("substance", "")
        
        if action_cryogen:
            self.cryogen_textbox.value = action_cryogen
        
        if action_target_temperature:
            action_target_temperature = json.loads(action_target_temperature)
            self.target_temperature_value_textbox.value = str(action_target_temperature["value"])
            self.target_temperature_unit_dropdown.value = action_target_temperature["unit"]
        
        if action_substance:
            self.substance_dropdown.value = action_substance
        
        try:
            component_settings = action_props.get("component_settings", {})
            component_settings = json.loads(component_settings)
            component_permid = action_props.get("component", {})
            if component_settings:
                component_object = get_cached_object(self.openbis_session, component_permid)
                component_type = component_object.type.code
                component_type_lower = component_type.lower()
                self.component_dropdown.value = component_permid
                if "target_temperature" in component_settings:
                    self.target_temperature_value_comp_textbox.value = str(component_settings["target_temperature"]["value"])
                    self.target_temperature_unit_comp_dropdown.value = component_settings["target_temperature"]["unit"]
                    
                if "bias_voltage" in component_settings:
                    self.bias_voltage_value_textbox.value = str(component_settings["bias_voltage"]["value"])
                    self.bias_voltage_unit_dropdown.value = component_settings["bias_voltage"]["unit"]
                
                if "discharge_voltage" in component_settings:
                    self.discharge_voltage_value_textbox.value = str(component_settings["discharge_voltage"]["value"])
                    self.discharge_voltage_unit_dropdown.value = component_settings["discharge_voltage"]["unit"]
                
                if "discharge_current" in component_settings:
                    self.discharge_current_value_textbox.value = str(component_settings["discharge_current"]["value"])
                    self.discharge_current_unit_dropdown.value = component_settings["discharge_current"]["unit"]
                
                self.evaporator_p_value_textbox.value = str(component_settings.get("p_value", ""))
                self.evaporator_i_value_textbox.value = str(component_settings.get("i_value", ""))
                self.ep_percentage_textbox.value = str(component_settings.get("ep_percentage", ""))
        except:
            pass
        
        action_substrate_temperature = action_props.get("substrate_temperature", "")
        if action_substrate_temperature:
            action_substrate_temperature = json.loads(action_substrate_temperature)
            self.substrate_temperature_value_textbox.value = str(action_substrate_temperature["value"])
            self.substrate_temperature_unit_dropdown.value = action_substrate_temperature["unit"]
        
        self.comments_textarea.value = action_props["comments"] or ""
    
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

            if action_type == OPENBIS_OBJECT_TYPES["Annealing"]:
                action_properties.append(self.target_temperature_hbox)

            elif action_type == OPENBIS_OBJECT_TYPES["Cooldown"]:
                action_properties.append(self.target_temperature_hbox)
                action_properties.append(self.cryogen_hbox)

            elif action_type == OPENBIS_OBJECT_TYPES["Deposition"]:
                action_properties.append(self.substance_hbox)
                action_properties.append(self.substrate_temperature_hbox)

            elif action_type == OPENBIS_OBJECT_TYPES["Dosing"]:
                action_properties.append(self.substrate_temperature_hbox)
                action_properties.append(self.pressure_hbox)
                action_properties.append(self.gas_hbox)

            elif action_type == OPENBIS_OBJECT_TYPES["Sputtering"]:
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
                        component_object = get_cached_object(self.openbis_session, component_id)
                        all_instrument_components.append(component_object)
                        
                        # Evaporators contain sub-components (evaporator slots)
                        evaporator_type = OPENBIS_OBJECT_TYPES["Evaporator"]
                        evaporator_type_lower = evaporator_type.lower()
                        if component_object.type == evaporator_type:
                            evaporator_slots = component_object.props["evaporator_slots"]
                            for evaporator_slot_id in evaporator_slots:
                                evaporator_slot_object = get_cached_object(self.openbis_session, evaporator_slot_id)
                                all_instrument_components.append(evaporator_slot_object)
            for component_object in all_instrument_components:
                component_type = component_object.type.code
                component_type_lower = component_type.lower()
                component_actions_settings_prop = component_object.props["actions_settings"]
                if component_actions_settings_prop:
                    component_actions_settings = json.loads(component_actions_settings_prop)
                    for component_action_settings in component_actions_settings:
                        component_action_type = component_action_settings["action_type"]
                        if action_type == component_action_type:
                            component_list.append((component_object.props["name"], component_object.permId))
        
        return component_list
    
    def load_component_settings_list(self, change):
        action_type = self.action_type_dropdown.value
        if action_type != "-1":
            component_permid = self.component_dropdown.value
            if component_permid != "-1":
                component_object = get_cached_object(self.openbis_session, component_permid)
                component_type = component_object.type.code
                component_type_lower = component_type.lower()
                component_settings_property = component_object.props["actions_settings"]
                component_settings_widgets = []
                if component_settings_property:
                    component_settings = json.loads(component_settings_property)
                    for settings in component_settings:
                        if action_type == settings["action_type"]:
                            settings_properties = settings.get("component_properties", {})
                            for prop in settings_properties:
                                setting_widget = self.components_properties_widgets.get(prop, None)
                                if setting_widget:
                                    if component_object.props[prop]:
                                        if prop == "target_temperature":
                                            target_temperature_comp = json.loads(component_object.props[prop])
                                            self.target_temperature_value_comp_textbox.value = target_temperature_comp["value"]
                                            self.target_temperature_unit_comp_dropdown.value = target_temperature_comp["unit"]
                                        
                                        elif prop == "bias_voltage":
                                            bias_voltage_comp = json.loads(component_object.props[prop])
                                            self.bias_voltage_value_comp_textbox.value = bias_voltage_comp["value"]
                                            self.bias_voltage_unit_comp_dropdown.value = bias_voltage_comp["unit"]
                                            
                                        elif prop == "discharge_voltage":
                                            discharge_voltage_comp = json.loads(component_object.props[prop])
                                            self.discharge_voltage_value_comp_textbox.value = discharge_voltage_comp["value"]
                                            self.discharge_voltage_unit_comp_dropdown.value = discharge_voltage_comp["unit"]
                                        
                                        elif prop == "discharge_current":
                                            discharge_current_comp = json.loads(component_object.props[prop])
                                            self.discharge_current_value_comp_textbox.value = discharge_current_comp["value"]
                                            self.discharge_current_unit_comp_dropdown.value = discharge_current_comp["unit"]
                                        
                                        elif prop == "p_value":
                                            self.evaporator_p_value_textbox.value = str(component_object.props[prop])
                                        
                                        elif prop == "i_value":
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
        observable_type_options = OBSERVABLES_LABELS.copy()
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
        
        instrument_stm_type = OPENBIS_OBJECT_TYPES["Instrument STM"]
        instrument_stm_type_lower = instrument_stm_type.lower()
        self.instrument_type_components_dictionary = {
            OPENBIS_OBJECT_TYPES["Instrument STM"]: [
                "pumps", "gauges", 
                "vacuum_chambers", "ports_valves",
                "preparation_tools", "analysers", 
                "mechanical_components", "stm_components", 
                "control_data_acquisition", 
                "temperature_environment_control",
                "auxiliary_components", "tips_sensors", 
                "accessories"
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
        
        ion_gauge_type_lower = OPENBIS_OBJECT_TYPES["Ion Gauge"].lower()
        analyser_type_lower = OPENBIS_OBJECT_TYPES["Analyser"].lower()
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
        observable_permid = settings["observable"]
        observable_object = get_cached_object(self.openbis_session, observable_permid)
        observable_props = observable_object.props.all()
        observable_type = observable_object.type.code
        observable_type_lower = observable_type.lower()
        self.observable_type_dropdown.value = observable_type
        self.name_textbox.value = observable_props["name"] or ""
        self.description_textbox.value = observable_props["description"] or ""
        self.ch_name_textbox.value = observable_props["channel_name"] or ""
        self.comments_textarea.value = observable_props["comments"] or ""
        
        try:
            component_settings = observable_props.get("component_settings", {})
            component_settings = json.loads(component_settings)
            component_permid = observable_props.get("component", {})
            if component_settings:
                component_object = get_cached_object(self.openbis_session, component_permid)
                component_type = component_object.type.code
                component_type_lower = component_type.lower()
                self.component_dropdown.value = component_permid
                if "density" in component_settings:
                    self.density_value_textbox.value = str(component_settings["density"]["value"])
                    self.density_unit_dropdown.value = component_settings["density"]["unit"]
                if "filament_current" in component_settings:
                    self.filament_current_value_textbox.value = str(component_settings["filament_current"]["value"])
                    self.filament_current_unit_dropdown.value = component_settings["filament_current"]["unit"]
                
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
                        evaporator_type = OPENBIS_OBJECT_TYPES["Evaporator"]
                        evaporator_type_lower = evaporator_type.lower()
                        if component_object.type == evaporator_type:
                            evaporator_slots = component_object.props["evaporator_slots"]
                            for evaporator_slot_id in evaporator_slots:
                                evaporator_slot_object = get_cached_object(self.openbis_session, evaporator_slot_id)
                                all_instrument_components.append(evaporator_slot_object)
                
            for component_object in all_instrument_components:
                component_type = component_object.type.code
                component_type_lower = component_type.lower()
                component_observables_settings_prop = component_object.props["observables_settings"]
                if component_observables_settings_prop:
                    component_observables_settings = json.loads(component_observables_settings_prop)
                    for component_observable_settings in component_observables_settings:
                        component_observable_type = component_observable_settings["observable_type"]
                        if observable_type == component_observable_type:
                            component_list.append((component_object.props["name"], component_object.permId))

        return component_list
    
    def load_component_settings_list(self, change):
        observable_type = self.observable_type_dropdown.value
        if observable_type != "-1":
            component_permid = self.component_dropdown.value
            if component_permid != "-1":
                component_object = get_cached_object(self.openbis_session, component_permid)
                component_type = component_object.type.code
                component_type_lower = component_type.lower()
                component_settings_property = component_object.props["observables_settings"]
                component_settings_widgets = []
                if component_settings_property:
                    component_settings = json.loads(component_settings_property)
                    for settings in component_settings:
                        if observable_type == settings["observable_type"]:
                            settings_properties = settings.get("component_properties", {})
                            for prop in settings_properties:
                                setting_widget = self.components_properties_widgets.get(prop, None)
                                if setting_widget:
                                    if component_object.props[prop]:
                                        if prop == "density":
                                            density_comp = json.loads(component_object.props[prop])
                                            self.density_value_textbox.value = density_comp["value"]
                                            self.density_unit_dropdown.value = density_comp["unit"]
                                        
                                        elif prop == "filament_current":
                                            filament_current_comp = json.loads(component_object.props[prop])
                                            self.filament_current_value_textbox.value = filament_current_comp["value"]
                                            self.filament_current_unit_dropdown.value = filament_current_comp["unit"]
                                        
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