import json
from pathlib import Path
import os
import ipywidgets as widgets
from IPython.display import display, Javascript, HTML, clear_output
from pybis import Openbis
import ipyfilechooser
import sys
sys.path.append('/home/jovyan/aiida-openbis/Notebooks/importer')
import nanonis_importer
import pandas as pd
import numpy as np
import yaml
import shutil

def read_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)

def read_yaml(filepath: str):
    return yaml.safe_load(open(filepath, 'r'))

def create_json(json_content, filename):
    with open(filename, 'w') as file:
        json.dump(json_content, file, indent=4)

def close_notebook(b):
    display(Javascript('window.location.replace("home.ipynb")'))

class AppWidgets():
    def __init__(self, method_type, config_filename):
        self.method_type = method_type
        self.openbis_session, self.session_data = self.connect_openbis()
        self.config = read_json(config_filename)
        
        if self.method_type == "home":
            # Home page configuration
            self.openbis_connection_status_htmlbox = self.HTML(value = '')
            if self.openbis_session:
                self.openbis_connection_status_htmlbox.value = "<span style='color: green; font-weight: bold;'>Connected to openBIS!</span>"
                self.open_notebooks_html_enable_code = ''.join(self.config["home_page"]["enable_links"])
                self.open_notebooks_htmlbox = self.HTML(value = self.open_notebooks_html_enable_code)
            else:
                self.open_notebooks_html_disable_code = ''.join(self.config["home_page"]["disable_links"])
                self.openbis_connection_status_htmlbox.value = "<span style='color: red; font-weight: bold;'>Please connect to openBIS using Configure ELN available in start page!</span>"
                self.open_notebooks_htmlbox = self.HTML(value = self.open_notebooks_html_disable_code)
        else:
            self.raw_materials_types = self.config["raw_materials_types"]
            self.process_sample_types = self.config["process_sample_types"]
            self.samples_collection_openbis_path = self.config["samples_collection_openbis_path"]
            self.measurement_file_extensions = self.config["measurement_file_extensions"]
            
            self.openbis_schema = read_yaml("/home/jovyan/aiida-openbis/Notebooks/Metadata_Schemas_LinkML/materialMLinfo.yaml")
            
            self.support_files_uploader = widgets.FileUpload(multiple = True)
            
            self.materials_dropdown = self.Dropdown(description='', disabled=False, layout = widgets.Layout(width = '350px'))
            self.material_details_textbox = self.Textarea(description = "", disabled = True, layout = widgets.Layout(width = '425px', height = '200px'))
            self.material_image_box = self.Image(value = open("images/white_screen.jpg", "rb").read(), format = 'jpg', width = '200px', height = '300px', layout=widgets.Layout(border='solid 1px #cccccc'))
            self.material_metadata_boxes = widgets.HBox([self.material_details_textbox, self.material_image_box])
            self.material_selection_radio_button = self.Radiobuttons(description = '', options = self.config["materials_options"], disabled = False, layout = widgets.Layout(width = '300px'), style = {'description_width': "100px"})
            material_sorting_checkboxes_list = self.SortingCheckboxes("50px", "60px", "200px")
            self.material_sorting_checkboxes = widgets.HBox([e for e in material_sorting_checkboxes_list])
            
            self.samples_dropdown = self.Dropdown(description='Sample', disabled=False, layout = widgets.Layout(width = '400px'), style = {'description_width': "110px"})
            self.sample_details_textbox = self.Textarea(description = "", disabled = True, layout = widgets.Layout(width = '589px', height = '300px'))
            sample_sorting_checkboxes_list = self.SortingCheckboxes("130px", "60px", "200px")
            self.sample_sorting_checkboxes = widgets.HBox([e for e in sample_sorting_checkboxes_list])
            self.sample_metadata_boxes = widgets.HBox([widgets.VBox([self.samples_dropdown, self.sample_sorting_checkboxes]), self.sample_details_textbox])

            self.instruments_dropdown = self.Dropdown(description='Instrument', disabled=False, layout = widgets.Layout(width = '993px'), style = {'description_width': "110px"})
            instrument_sorting_checkboxes_list = self.SortingCheckboxes("130px", "60px", "200px")
            self.instrument_sorting_checkboxes = widgets.HBox([e for e in instrument_sorting_checkboxes_list])
            self.instruments_dropdown_boxes = widgets.VBox([self.instruments_dropdown, self.instrument_sorting_checkboxes])
            
            self.create_new_experiment_button = self.Button(description = '', disabled = False, button_style = '', tooltip = 'Add', icon = 'plus', layout = widgets.Layout(width = '50px', height = '25px'))
            self.save_new_experiment_button = self.Button(description = '', disabled = False, button_style = '', tooltip = 'Save', icon = 'save', layout = widgets.Layout(width = '50px', height = '35px', margin = '0 0 0 90px'))
            self.cancel_new_experiment_button = self.Button(description = '', disabled = False, button_style = '', tooltip = 'Cancel', icon = 'times', layout = widgets.Layout(width = '50px', height = '35px', margin = '0 0 0 5px'))
            self.new_experiment_name_textbox = self.Text(description = "Name", disabled = False, layout = widgets.Layout(width = '400px'), placeholder = f"Write experiment name here...", style = {'description_width': "110px"})
            
            self.projects_dropdown = self.Dropdown(description='Project', disabled=False, layout = widgets.Layout(width = '993px'), style = {'description_width': "110px"})
            project_sorting_checkboxes_list = self.SortingCheckboxes("130px", "60px", "200px")
            self.project_sorting_checkboxes = widgets.HBox([e for e in project_sorting_checkboxes_list])
            self.projects_dropdown_boxes = widgets.VBox([self.projects_dropdown, self.project_sorting_checkboxes])
            
            self.experiments_dropdown = self.Dropdown(description='Experiment', disabled=False, layout = widgets.Layout(width = '993px'), style = {'description_width': "110px"})
            experiment_sorting_checkboxes_list = self.SortingCheckboxes("130px", "60px", "200px")
            self.experiment_sorting_checkboxes = widgets.HBox([e for e in experiment_sorting_checkboxes_list])
            self.experiments_dropdown_details = widgets.HBox([self.experiments_dropdown, self.create_new_experiment_button])
            self.experiments_dropdown_boxes = widgets.VBox([self.experiments_dropdown_details, self.experiment_sorting_checkboxes])
            
            self.molecules_dropdown = self.Dropdown(description='Substance', disabled=False, layout = widgets.Layout(width = '350px'), style = {'description_width': "110px"})
            self.molecule_details_textbox = self.Textarea(description = "", disabled = True, layout = widgets.Layout(width = '415px', height = '250px'))
            self.molecule_image_box = self.Image(value = open("images/white_screen.jpg", "rb").read(), format = 'jpg', width = '220px', height = '250px', layout=widgets.Layout(border='solid 1px #cccccc'))
            molecule_sorting_checkboxes_list = self.SortingCheckboxes("170px", "60px", "200px")
            self.molecule_sorting_checkboxes = widgets.HBox([e for e in molecule_sorting_checkboxes_list])
            self.molecule_metadata_boxes = widgets.HBox([widgets.VBox([self.molecules_dropdown, self.molecule_sorting_checkboxes]), self.molecule_details_textbox, self.molecule_image_box])

            self.method_name_textbox = self.Text(description = "Name", disabled = False, layout = widgets.Layout(width = '400px'), placeholder = f"Write {method_type} task name here...", style = {'description_width': "150px"})
            self.comments_textbox = self.Textarea(description = "Comments", disabled = False, layout = widgets.Layout(width = '993px', height = '200px'), placeholder = "Write comments here...", style = {'description_width': "150px"})
            
            # Quantity properties widgets
            self.properties = ["duration", "pressure", "discharge_voltage", "voltage", "temperature",
                               "current", "angle", "stabilisation_time", "deposition_time", "substrate_temperature",
                               "molecule_temperature"]
            self.property_widgets = {}
            
            for prop in self.properties:
                if self.config["properties"][prop]["property_type"] == "quantity_value":
                    self.property_widgets[prop] = self.FloatTextwithDropdownWidget(
                        self.config["properties"][prop]["description"], 
                        widgets.Layout(width = self.config["properties"][prop]["box_layout"]["width"]), 
                        0, 
                        {'description_width': self.config["properties"][prop]["box_layout"]["description_width"]}, 
                        widgets.Layout(width = self.config["properties"][prop]["dropdown_layout"]["width"]), 
                        self.config["properties"][prop]["units"], 
                        self.config["properties"][prop]["units"][0]
                    )
            
            self.property_widgets["sum_formula"] = self.Text(
                description = "Substance sum formula", 
                disabled = False, 
                layout = widgets.Layout(
                    width = self.config["properties"]["sum_formula"]["box_layout"]["width"]), 
                    placeholder = f"Write substance formula here...",
                    style = {'description_width': self.config["properties"]["sum_formula"]["box_layout"]["description_width"]}
                )
            
            self.property_widgets["evaporator_slot"] = self.IntSliderwithTextWidget(
                1, 
                self.config["properties"]["evaporator_slot"]["description"], 
                [1,6], 
                widgets.Layout(width = self.config["properties"]["evaporator_slot"]["slider_layout"]["width"]), 
                {'description_width': self.config["properties"]["evaporator_slot"]["slider_layout"]["description_width"]}, 
                "Write evaporator slot details...", 
                widgets.Layout(width = self.config["properties"]["evaporator_slot"]["box_layout"]["width"])
            )
            
            if self.method_type in self.config["objects"].keys():
                items = [self.property_widgets[prop] for prop in self.config["objects"][self.method_type]["properties"]]
                self.method_properties = widgets.VBox(items)

            self.sample_out_name_textbox = self.Text(description = "Name", disabled = False, layout = widgets.Layout(width = '400px'), placeholder = f"Write sample name here...", style = {'description_width': "110px"})
            self.create_button = self.Button(description = '', disabled = False, button_style = '', tooltip = 'Save', icon = 'save', layout = widgets.Layout(width = '100px', height = '50px'))
            self.quit_button = self.Button(description = '', disabled = False, button_style = '', tooltip = 'Main menu', icon = 'home', layout = widgets.Layout(width = '100px', height = '50px'))
            self.bottom_buttons_hbox = widgets.HBox([self.create_button, self.quit_button])
            self.folder_selector = self.FileChooser(path = '.', select_default=True, use_dir_icons=True, show_only_dirs = True)
            
            # Assign functions to widgets
            self.create_new_experiment_button.on_click(self.create_new_experiment_button_on_click)
            self.cancel_new_experiment_button.on_click(self.cancel_new_experiment_button_on_click)
            self.save_new_experiment_button.on_click(self.save_new_experiment_button_on_click)
            
            # Increase buttons icons' size
            self.increase_buttons_size = HTML("""
            <style>
                .fa-save {
                    font-size: 2em !important; /* Increase icon size */
                }
                .fa-home {
                    font-size: 2em !important; /* Increase icon size */
                }
                .fa-times {
                    font-size: 2em !important; /* Increase icon size */
                }
            </style>
            """)
            display(self.increase_buttons_size)
    
    # @staticmethod
    # def DropdownwithSortingCheckboxesWidget(*args):
    #     dropdown = AppWidgets.Dropdown(description='Instrument', disabled=False, layout = widgets.Layout(width = '993px'), style = {'description_width': "110px"})
    #     sorting_checkboxes_list = AppWidgets.SortingCheckboxes("130px", "60px", "200px")
    #     widgets = [dropdown]
    #     widgets.extend(sorting_checkboxes_list)
    #     return widgets.VBox(widgets)
    
    @staticmethod
    def IntSliderwithTextWidget(*args):
        intslider = AppWidgets.IntSlider(value = args[0], description = args[1], min = args[2][0], max = args[2][1], disabled = False, layout = args[3], style = args[4])
        textbox = AppWidgets.Text(value = '', description = '', placeholder= args[5], disabled = False, layout = args[6])
        return widgets.HBox([intslider, textbox])
    
    @staticmethod
    def SortingCheckboxes(*args):
        return [
            widgets.Label(value = "Sort by:", layout = widgets.Layout(width = args[0], justify_content='flex-end')),
            AppWidgets.Checkbox(description = 'Name', value = False, disabled = False, layout = widgets.Layout(width = args[1]), indent = False),
            AppWidgets.Checkbox(description = 'Registration date', value = False, disabled = False, layout = widgets.Layout(width = args[2]), indent = False)
        ]
    
    @staticmethod
    def FloatTextwithDropdownWidget(*args):
        floatbox = AppWidgets.FloatText(description=args[0], disabled=False, layout = args[1], value = args[2], style = args[3])
        dropdown = AppWidgets.Dropdown(description='', disabled=False, layout = args[4], options = args[5], value = args[6])
        return widgets.HBox([floatbox, dropdown])
    
    @staticmethod
    def IntSlider(**kwargs):
        return widgets.IntSlider(**kwargs)
    
    @staticmethod
    def HTML(**kwargs):
        return widgets.HTML(**kwargs)
    
    @staticmethod
    def FileChooser(**kwargs):
        return ipyfilechooser.FileChooser(**kwargs)
    
    @staticmethod
    def Image(**kwargs):
        return widgets.Image(**kwargs)
    
    @staticmethod
    def Button(**kwargs):
        return widgets.Button(**kwargs)
    
    @staticmethod
    def Checkbox(**kwargs):
        return widgets.Checkbox(**kwargs)
    
    @staticmethod
    def Textarea(**kwargs):
        textarea_box = widgets.Textarea(**kwargs)
        return textarea_box

    @staticmethod
    def Text(**kwargs):
        text_box = widgets.Text(**kwargs)
        return text_box
        
    @staticmethod
    def Radiobuttons(**kwargs):
        radio_buttons = widgets.RadioButtons(**kwargs)
        return radio_buttons

    @staticmethod
    def Dropdown(**kwargs):
        dropdown_box = widgets.Dropdown(**kwargs)
        return dropdown_box

    @staticmethod
    def FloatText(**kwargs):
        floattext_box = widgets.FloatText(**kwargs)
        return floattext_box

    @staticmethod
    def read_file(filename):
        file = open(filename, "rb")
        return file.read()
    
    @staticmethod
    def sort_dataframe(df, columns, ascending_columns):
        df = df.sort_values(by = columns, ascending = ascending_columns)
        return AppWidgets.dataframe_to_list_of_tuples(df)

    @staticmethod
    def dataframe_to_list_of_tuples(df):
        return list(df.itertuples(index = False, name = None))

    @staticmethod
    def connect_openbis():
        try:
            eln_config = Path.home() / ".aiidalab" / "aiidalab-eln-config.json"
            eln_config.parent.mkdir(parents=True, exist_ok=True)  # making sure that the folder exists.
            config = read_json(eln_config)
            eln_url = "https://local.openbis.ch"
            eln_token = config[eln_url]["token"]
            session_data = {"url": eln_url, "token": eln_token}
            openbis_session = Openbis(eln_url, verify_certificates = False)
            openbis_session.set_token(eln_token)
        except ValueError:
            openbis_session = None
            session_data = {}
        
        return openbis_session, session_data
    
    def create_new_experiment_button_on_click(self, b):
        clear_output()
        self.load_list("PROJECT", self.projects_dropdown, self.project_sorting_checkboxes, "project")
        display(self.experiments_dropdown_boxes, self.new_experiment_name_textbox, self.projects_dropdown_boxes,
                widgets.HBox([self.save_new_experiment_button, self.cancel_new_experiment_button]), 
                self.sample_metadata_boxes, self.instruments_dropdown_boxes)
    
    def cancel_new_experiment_button_on_click(self, b):
        clear_output()
        display(self.experiments_dropdown_boxes, self.sample_metadata_boxes, self.instruments_dropdown_boxes)
    
    def save_new_experiment_button_on_click(self, b):
        self.create_experiment_in_openbis(self.projects_dropdown.value, self.new_experiment_name_textbox.value)
        clear_output()
        self.load_list("EXPERIMENT", self.experiments_dropdown, self.experiment_sorting_checkboxes, "experiment")
        display(self.experiments_dropdown_boxes, self.sample_metadata_boxes, self.instruments_dropdown_boxes)
        
    def sort_dropdown(self, sorting_checkboxes, dropdown_box):
        dropdown_list = list(dropdown_box.options[1:]) # Default -1 message should not be sorted.
        df = pd.DataFrame(dropdown_list, columns = ['Name', 'PermID'])
        
        # Determine sort columns and order based on checkboxes
        columns, ascending = [], []
        if sorting_checkboxes.children[1].value:
            columns.append('Name')
            ascending.append(True)
            
        if sorting_checkboxes.children[2].value:
            columns.append('PermID')
            ascending.append(False)
        
        if columns:
            dropdown_list = self.sort_dataframe(df, columns, ascending)

        dropdown_list.insert(0, dropdown_box.options[0])
        dropdown_box.options = dropdown_list
    
    def get_next_experiment_code(self):
        experiments = self.openbis_session.get_experiments(type = "EXPERIMENT")
        experiment_number = max(int(exp.code.rsplit('_')[-1]) for exp in experiments) + 1
        return f"{experiments[0].code.rsplit('_')[0]}_{experiment_number}"
    
    def create_experiment_in_openbis(self, project_id, experiment_name):
        experiment_code = self.get_next_experiment_code()
        experiment = self.openbis_session.new_experiment(code = experiment_code, type = "EXPERIMENT", project = project_id, props = {"$name": experiment_name})
        experiment.save()
    
    def create_openbis_object(self, **kwargs):
        openbis_object = self.openbis_session.new_object(**kwargs)
        openbis_object.save()
        return openbis_object

    def create_openbis_collection(self, **kwargs):
        measurements_collection = self.openbis_session.new_collection(**kwargs)
        measurements_collection.save()
        return measurements_collection
    
    # Function to create sample object inside openBIS using information selected in the app
    def create_sample_action(self, b):
        samples_names = [sample.props["$name"] for sample in self.openbis_session.get_objects(type = "SAMPLE")]
        if self.sample_out_name_textbox.value in samples_names:
            display(Javascript(f"alert('{'Sample name already exists!'}')"))
            return
        else:
            sample_parents = [] if self.materials_dropdown.value == -1 or self.material_selection_radio_button.value == "No material" else [self.materials_dropdown.value]
            sample_props = {"$name": self.sample_out_name_textbox.value, "exists": True}
            self.create_openbis_object(type="SAMPLE", collection=self.samples_collection_openbis_path, props=sample_props, parents=sample_parents)
            print("Upload successful!")
    
    # Function to handle changes in the materials dropdown
    def load_material_metadata(self, change):
        if self.materials_dropdown.value == -1:
            self.material_details_textbox.value = ''
            self.material_image_box.value = self.read_file("images/white_screen.jpg")
            return
        
        property_lists = self.config["property_lists"]["materials"]
        property_list = property_lists.get(self.material_selection_radio_button.value, [])
        
        material_object = self.openbis_session.get_object(self.materials_dropdown.value)
        material_dataset = material_object.get_datasets()[0]
        
        if material_dataset:
            material_dataset.download(destination="images")
            self.material_image_box.value = self.read_file(f"images/{material_dataset.permId}/{material_dataset.file_list[0]}")
            shutil.rmtree(f"images/{material_dataset.permId}")
        else:
            self.material_image_box.value = self.read_file("images/white_screen.jpg")
        
        material_metadata = material_object.props.all()
        material_metadata_string = ""
        for prop_name, prop_key in property_list:
            if prop_key in self.config["property_lists"]["qunit_properties"]:
                value = material_metadata.get(prop_key)
                if value:
                    prop_dict = json.loads(value)
                    material_metadata_string += f"{prop_name}: {prop_dict['value']} {prop_dict['unit']}\n"
                else:
                    material_metadata_string += f"{prop_name}: {value}\n"
            else:
                material_metadata_string += f"{prop_name}: {material_metadata.get(prop_key)}\n"

        self.material_details_textbox.value = material_metadata_string
    
    def select_material_radio_change(self, change):
        self.material_details_textbox.value = ''
        clear_output()
        display(self.material_selection_radio_button)

        material_types = {
            "Crystal": ("CRYSTAL", "Select crystal..."),
            "Wafer substrate": ("WAFER_SUBSTRATE", "Select wafer substrate..."),
            "2D-layer material": ("2D_LAYERED_MATERIAL", "Select 2D-layer material...")
        }

        material_type = self.material_selection_radio_button.value
        if material_type == "No material":
            return

        material_class, placeholder = material_types.get(material_type, (None, None))
        if material_class:
            materials = self.openbis_session.get_objects(type = material_class)
            materials_names_permids = [(f"{mat.props['$name']} ({mat.permId})", mat.permId) for mat in materials]
            self.materials_dropdown.options = [(placeholder, -1)] + materials_names_permids
            self.materials_dropdown.value = -1

            self.sort_dropdown(self.material_sorting_checkboxes, self.materials_dropdown)

            self.materials_dropdown.observe(self.load_material_metadata, names='value')
            for checkbox in self.material_sorting_checkboxes.children[1:3]:
                checkbox.observe(lambda change: self.sort_dropdown(self.material_sorting_checkboxes, self.materials_dropdown), names='value')

            display(widgets.HBox([self.materials_dropdown, self.material_sorting_checkboxes]))
            display(self.material_metadata_boxes)
    
    def upload_datasets(self, method_object):
        for file_info in self.support_files_uploader.value:
            filename = file_info['name']
            with open(filename, 'wb') as f:  # Save file content
                f.write(file_info['content'])
            self.openbis_session.new_dataset(type = 'RAW_DATA', sample = method_object, files = [filename]).save()
            os.remove(filename)
    
    def create_process_action(self, b):
        samples_names = [sample.props["$name"] for sample in self.openbis_session.get_objects(type="SAMPLE")]
        
        if self.sample_out_name_textbox.value in samples_names:
            display(Javascript("alert('Sample name already exists!')"))
            return
        
        if self.experiments_dropdown.value == -1:
            print("Select an experiment.")
            return
        
        if self.samples_dropdown.value == -1:
            print("Select a sample.")
            return
        
        if self.instruments_dropdown.value == -1:
            print("Select an instrument.")
            return
        
        if self.molecules_dropdown.value == -1 and self.method_type == "deposition":
            print("Select a substance.")
            return

        # Prepare sample parents based on method type
        sample_parents = [self.samples_dropdown.value, self.instruments_dropdown.value]
        if self.method_type == "deposition":
            sample_parents.append(self.molecules_dropdown.value)

        object_properties = {
            "$name": self.method_name_textbox.value,
            "comments": self.comments_textbox.value,
            }
        
        for prop in self.config["objects"][self.method_type]["properties"]:
            if prop == "evaporator_slot":
                object_properties[prop] = {"evaporator_number": self.property_widgets[prop].children[0].value, "details": self.property_widgets[prop].children[1].value}
            elif prop == "sum_formula":
                object_properties[prop] = self.property_widgets[prop]
            else:
                object_properties[prop] = {"has_value": self.property_widgets[prop].children[0].value, "has_unit": self.property_widgets[prop].children[1].value}

        method_object = self.create_openbis_object(type = self.method_type.upper(), collection = self.experiments_dropdown.value, props = object_properties, parents = sample_parents)
        self.upload_datasets(method_object)

        # Turn off sample visibility
        parent_sample = self.openbis_session.get_object(self.samples_dropdown.value)
        parent_sample.props["exists"] = False
        parent_sample.save()

        sample_props = {"$name": self.sample_out_name_textbox.value, "exists": True}
        self.create_openbis_object(type = "SAMPLE", collection = self.samples_collection_openbis_path, props = sample_props, parents = [method_object])
        print("Upload successful!")
    
    # Function to handle changes in the substances dropdown
    def load_molecule_metadata(self, change):
        if self.molecules_dropdown.value == -1:
            self.molecule_details_textbox.value = ''
            self.molecule_image_box.value = self.read_file("images/white_screen.jpg")
            return
        
        # Get substance metadata
        property_list = self.config["property_lists"]["materials"]["Substance"]
        molecule_property_list = self.config["property_lists"]["materials"]["Molecule"]
        material_object = self.openbis_session.get_object(self.molecules_dropdown.value)
        
        # Get molecule image
        molecule_object = self.openbis_session.get_object(material_object.props.all()['has_molecule'])
        molecule_metadata = molecule_object.props.all()
        molecule_dataset = molecule_object.get_datasets(type="ELN_PREVIEW")[0]

        if molecule_dataset:
            molecule_dataset.download(destination="images")
            material_image_filepath = molecule_dataset.file_list[0]
            self.molecule_image_box.value = self.read_file(f"images/{molecule_dataset.permId}/{material_image_filepath}")
            shutil.rmtree(f"images/{molecule_dataset.permId}/")
        else:
            self.molecule_image_box.value = self.read_file("images/white_screen.jpg")

        material_metadata = material_object.props.all()
        material_metadata_string = ""
        for name, key in property_list:
            value = material_metadata.get(key)
            if key in self.config["property_lists"]["qunit_properties"] and value is not None:
                value = json.loads(value)
                material_metadata_string += f"{name}: {value['value']} {value['unit']}\n"
            elif key == "has_molecule":
                material_metadata_string += f"Molecule:\n"
                for prop_name, key in molecule_property_list:
                    prop_value = molecule_metadata.get(key)
                    material_metadata_string += f"- {prop_name}: {prop_value}\n"
            else:
                material_metadata_string += f"{name}: {value}\n"

        self.molecule_details_textbox.value = material_metadata_string

    def load_sample_metadata(self, change):
        if self.samples_dropdown.value == -1:
            self.experiments_dropdown.value = -1
            self.instruments_dropdown.value = -1
            self.sample_details_textbox.value = ''
            
            if self.method_type.upper() in self.process_sample_types:
                self.sample_out_name_textbox.value = ''
            return
        
        sample_object = self.openbis_session.get_object(self.samples_dropdown.value)
        sample_parents_metadata = self.get_parents_recursive(sample_object, [])
        
        last_sample_process = None
        sample_strings = {"processes": [], "materials": []}
        number_parents = len(sample_parents_metadata)
        parent_idx = 0
        while parent_idx < number_parents:
            parent_metadata = sample_parents_metadata[parent_idx]
            if parent_metadata[0] == "DEPOSITION":
                deposition_object = self.openbis_session.get_object(parent_metadata[1])
                
                # Get substance used in the specific deposition
                substance_metadata = []
                for parent_id in deposition_object.parents:
                    deposition_parent_object = self.openbis_session.get_object(parent_id)
                    if deposition_parent_object.type == "SUBSTANCE":
                        deposition_parent_object_metadata = deposition_parent_object.props.all()
                        substance_metadata = [deposition_parent_object_metadata["$name"], deposition_parent_object.permId]
                        
                if substance_metadata: # If deposition does not contain any substance, the app must not try to display it
                    sample_metadata_string = f"> {parent_metadata[0]} ({parent_metadata[3]}, {parent_metadata[1]}, {parent_metadata[2]}) [{substance_metadata[0]} ({substance_metadata[1]})]"
                else:
                    sample_metadata_string = f"> {parent_metadata[0]} ({parent_metadata[3]}, {parent_metadata[1]}, {parent_metadata[2]})"
                parent_idx += 1
            else:
                sample_metadata_string = f"> {parent_metadata[0]} ({parent_metadata[3]}, {parent_metadata[1]}, {parent_metadata[2]})"
            
            if parent_metadata[0] in self.process_sample_types:
                sample_strings["processes"].append(sample_metadata_string)
                # Get the last sample preparation method performed on the sample in order to search the correct experiment where the sample is being used.
                if last_sample_process is None:
                    last_sample_process = parent_metadata[1]
            elif parent_metadata[0] in self.raw_materials_types:
                sample_strings["materials"].append(sample_metadata_string)
                
            parent_idx += 1
        
        sample_metadata_string = (f"PermId: {sample_object.attrs.permId}\nMaterial:\n" +
                              "\n".join(sample_strings["materials"]) + 
                              "\nProcesses:\n" + 
                              "\n".join(sample_strings["processes"]))

        if last_sample_process:
            last_sample_process_object = self.openbis_session.get_object(last_sample_process)
            
            if self.method_type.upper() in self.process_sample_types:
                # Automatically select the experiment where the last sample process task was saved
                last_sample_process_experiment = self.openbis_session.get_experiment(last_sample_process_object.attrs.experiment)
                if self.experiments_dropdown.value != -1:
                    display(Javascript(f"alert('{'Experiment was changed!'}')"))
                self.experiments_dropdown.value =  last_sample_process_experiment.permId
            
            # Automatically select the instrument used in the last sample process task
            for parent in last_sample_process_object.parents:
                parent_object = self.openbis_session.get_object(parent)
                if parent_object.type == "INSTRUMENT":
                    self.instruments_dropdown.value = parent_object.permId
        
        self.sample_details_textbox.value = sample_metadata_string
        if self.method_type.upper() in self.process_sample_types:
            sample_name = sample_object.props['$name']
            self.sample_out_name_textbox.value = f"{sample_name}_{self.method_name_textbox.value}" if self.method_name_textbox.value else sample_name
            
    def load_list(self, type, dropdown, sorting_checkboxes, label):
        if type == "EXPERIMENT":
            items = self.openbis_session.get_collections(type = type)
            items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
        elif type == "PROJECT":
            items = self.openbis_session.get_projects()
            items_names_permids = [(item.attrs.identifier, item.permId) for item in items]
        else:
            items = self.openbis_session.get_objects(type = type)
            if type == "SAMPLE":
                items = self.filter_samples(items)
                items_names_permids = [(f"{item.props['$name']}", item.permId) for item in items]
            elif type == "SUBSTANCE":
                items_names_permids = [(f"{item.props['empa_number']}{item.props['batch']} ({item.permId})", item.permId) for item in items]
            else:
                items_names_permids = [(f"{item.props['$name']} ({item.permId})", item.permId) for item in items]
        
        items_names_permids.insert(0, (f'Select {label}...', -1))
        dropdown.options = items_names_permids
        dropdown.value = -1
        self.sort_dropdown(sorting_checkboxes, dropdown)
        for checkbox in sorting_checkboxes.children[1:3]:
            checkbox.observe(lambda change: self.sort_dropdown(sorting_checkboxes, dropdown), names='value')
    
    def load_dropdown_lists(self):
        # Populate dropdown lists
        self.load_list("SAMPLE", self.samples_dropdown, self.sample_sorting_checkboxes, "sample")
        self.load_list("INSTRUMENT", self.instruments_dropdown, self.instrument_sorting_checkboxes, "instrument")
        if self.method_type.upper() in self.process_sample_types:
            self.load_list("EXPERIMENT", self.experiments_dropdown, self.experiment_sorting_checkboxes, "experiment")
            if self.method_type.upper() == "DEPOSITION":
                self.load_list("SUBSTANCE", self.molecules_dropdown, self.molecule_sorting_checkboxes, "substance")

    def get_parents_recursive(self, object, object_parents_metadata):
        object_parents_metadata.append([object.attrs.type, object.attrs.permId, object.attrs.registrationDate,object.props['$name']])
        
        for parent in object.parents:
            self.get_parents_recursive(self.openbis_session.get_object(parent), object_parents_metadata)
        return object_parents_metadata

    def update_text(self, change):
        if self.samples_dropdown.value != -1:
            selected_sample_name = next(label for label, val in self.samples_dropdown.options if val == self.samples_dropdown.value)
            if len(self.method_name_textbox.value) > 0:
                self.sample_out_name_textbox.value = f"{selected_sample_name}_{self.method_name_textbox.value}"
            else:
                self.sample_out_name_textbox.value = selected_sample_name
    
    def upload_measurements_to_openbis(self, b):
        sample_object = self.openbis_session.get_object(self.samples_dropdown.value)
    
        # Get measurements collection if available
        measurements_collection = next(
            (self.openbis_session.get_object(child_id).collection 
            for child_id in sample_object.children 
            if self.openbis_session.get_object(child_id).type in ["1D_MEASUREMENT", "2D_MEASUREMENT"]), 
            None
        )
        
        # Get sample project from parent
        sample_project = next(
            (self.openbis_session.get_object(parent_id).project 
            for parent_id in sample_object.parents 
            if self.openbis_session.get_object(parent_id).type in self.process_sample_types), 
            None
        )

        # Check folder for valid files
        data_folder = self.folder_selector.selected_path
        data_files = os.listdir(data_folder)
        correct_folder = all(f".{filename.split('.')[-1]}" in self.measurement_file_extensions for filename in data_files)
        
        if correct_folder:
            if measurements_collection:
                measurements_collection = f"{sample_project.identifier}/{measurements_collection}"
                measurements_collection = self.openbis_session.get_collection(measurements_collection)
            else:
                if not sample_project:
                    print("Sample does not belong to any project yet.")
                    return
                collection_props = {"$name": f"Measurements from Sample {sample_object.props['$name']}", "$default_collection_view": "IMAGING_GALLERY_VIEW"}
                measurements_collection = self.create_openbis_collection(
                    code=f"MEASUREMENTS_COLLECTION_{sample_object.code}", 
                    type="COLLECTION", project=sample_project, 
                    props=collection_props
                )
            
            nanonis_importer.upload_measurements_into_openbis(
                self.session_data['url'], data_folder, 
                measurements_collection.permId, sample_object.permId, 
                self.instruments_dropdown.value
            )
            print("Upload successful!")
        else:
            print(f"Folder contains unrecognised files. Please remove files with extensions other than: {', '.join(self.measurement_file_extensions)}.")
    
    def filter_samples(self, samples):
        """
        Function for hiding intermediate samples which result 
        from the different steps of samples preparation.

        Args:
            samples (_type_): _description_
        """
        selected_samples = []
        for sample in samples:
            sample = self.openbis_session.get_object(sample.permId)
            if sample.props["exists"] == "true":
                selected_samples.append(sample)
        return selected_samples