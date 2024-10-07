import json
from pathlib import Path
import os
import ipywidgets as widgets
from IPython.display import display, Javascript, HTML, clear_output
from pybis import Openbis
from ipyfilechooser import FileChooser
import sys
sys.path.append('/home/jovyan/aiida-openbis/Notebooks/importer')
import nanonis_importer
import pandas as pd
from datetime import datetime
import numpy as np

def read_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)

def create_json(json_content, filename):
    with open(filename, 'w') as file:
        json.dump(json_content, file, indent=4)

def close_notebook(b):
    display(Javascript('window.location.replace("home.ipynb")'))

class AppWidgets():
    def __init__(self, method_type, raw_materials_types = [], process_sample_types = [], samples_collection_openbis_path = ""):
        self.openbis_session = None
        self.method_type = method_type
        self.raw_materials_types = raw_materials_types
        self.process_sample_types = process_sample_types
        self.samples_collection_openbis_path = samples_collection_openbis_path
        self.measurement_file_extensions = [".sxm", ".dat"]
        
        self.support_files_uploader = widgets.FileUpload(multiple = True)
        
        self.materials_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '350px'))
        self.material_details_textbox = self.get_textarea_box(description = "", disabled = True, layout = widgets.Layout(width = '425px', height = '200px'))
        self.material_image_box = widgets.Image(value = open("images/white_screen.jpg", "rb").read(), format = 'jpg', width = '200px', height = '300px', layout=widgets.Layout(border='solid 1px #cccccc'))
        self.material_metadata_boxes = widgets.HBox([self.material_details_textbox, self.material_image_box])

        self.material_selection_radio_button = self.get_radio_button(description = '', options=['No material', 'Crystal', 'Wafer substrate', '2D-layer material'], disabled = False, layout = widgets.Layout(width = '300px'), description_width = "100px")
        material_sorting_checkboxes_list = [
            widgets.Label(value = "Sort by:", layout = widgets.Layout(width = '60px')),
            self.get_check_box(description = 'Name', value = False, disabled = False, layout = widgets.Layout(width = '60px')),
            self.get_check_box(description = 'Registration date', value = False, disabled = False, layout = widgets.Layout(width = '200px'))
        ]
        self.material_sorting_checkboxes = widgets.HBox([e for e in material_sorting_checkboxes_list])
        
        self.samples_dropdown = self.get_dropdown_box(description='Sample', disabled=False, layout = widgets.Layout(width = '400px'), description_width = '110px')
        self.sample_details_textbox = self.get_textarea_box(description = "", disabled = True, layout = widgets.Layout(width = '589px', height = '300px'))
        
        sample_sorting_checkboxes_list = [
            widgets.Label(value = "Sort by:", layout = widgets.Layout(width = '130px', justify_content='flex-end')),
            self.get_check_box(description = 'Name', value = False, disabled = False, layout = widgets.Layout(width = '60px')),
            self.get_check_box(description = 'Registration date', value = False, disabled = False, layout = widgets.Layout(width = '200px'))
        ]
        self.sample_sorting_checkboxes = widgets.HBox([e for e in sample_sorting_checkboxes_list])
        
        self.sample_metadata_boxes = widgets.HBox([widgets.VBox([self.samples_dropdown, self.sample_sorting_checkboxes]), self.sample_details_textbox])

        self.instruments_dropdown = self.get_dropdown_box(description='Instrument', disabled=False, layout = widgets.Layout(width = '993px'), description_width = '110px')
        self.experiments_dropdown = self.get_dropdown_box(description='Experiment', disabled=False, layout = widgets.Layout(width = '993px'), description_width = '110px')

        self.duration_value_floatbox = self.get_floattext_box(description='Duration', disabled=False, layout = widgets.Layout(width = '200px'), value = 0, description_width = "110px")
        self.duration_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["sec", "min", "hrs"], value = "sec")
        self.duration_hbox = widgets.HBox([self.duration_value_floatbox, self.duration_unit_dropdown])

        self.pressure_value_floatbox = self.get_floattext_box(description='Pressure', disabled=False, layout = widgets.Layout(width = '200px'), value = 0, description_width = "110px")
        self.pressure_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["mBar", "Bar", "Pa", "kPa"], value = "mBar")
        self.pressure_hbox = widgets.HBox([self.pressure_value_floatbox, self.pressure_unit_dropdown])

        self.discharge_voltage_value_floatbox = self.get_floattext_box(description='Discharge voltage', disabled=False, layout = widgets.Layout(width = '200px'), value = 0, description_width = "110px")
        self.discharge_voltage_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["V", "kV"], value = "kV")
        self.discharge_voltage_hbox = widgets.HBox([self.discharge_voltage_value_floatbox, self.discharge_voltage_unit_dropdown])
        
        self.voltage_value_floatbox = self.get_floattext_box(description='Voltage', disabled=False, layout = widgets.Layout(width = '200px'), value = 0, description_width = "110px")
        self.voltage_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["uV", "mV", "V"], value = "V")
        self.voltage_hbox = widgets.HBox([self.voltage_value_floatbox, self.voltage_unit_dropdown])
        
        self.temperature_value_floatbox = self.get_floattext_box(description='Temperature', disabled=False, layout = widgets.Layout(width = '200px'), value = 0, description_width = "110px")
        self.temperature_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["oF", "oC", "K"], value = "K")
        self.temperature_hbox = widgets.HBox([self.temperature_value_floatbox, self.temperature_unit_dropdown])

        self.current_value_floatbox = self.get_floattext_box(description='Current', disabled=False, layout = widgets.Layout(width = '200px'), value = 0, description_width = "110px")
        self.current_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["uA", "mA", "A"], value = "A")
        self.current_hbox = widgets.HBox([self.current_value_floatbox, self.current_unit_dropdown])
        
        self.angle_value_floatbox = self.get_floattext_box(description='Angle', disabled=False, layout = widgets.Layout(width = '200px'), value = 0, description_width = "110px")
        self.angle_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["deg", "rad"], value = "deg")
        self.angle_hbox = widgets.HBox([self.angle_value_floatbox, self.angle_unit_dropdown])
        
        self.molecules_dropdown = self.get_dropdown_box(description='Molecule', disabled=False, layout = widgets.Layout(width = '350px'), description_width = "110px")
        self.molecule_details_textbox = self.get_textarea_box(description = "", disabled = True, layout = widgets.Layout(width = '415px', height = '250px'))
        self.molecule_image_box = widgets.Image(value = open("images/white_screen.jpg", "rb").read(), format = 'jpg', width = '220px', height = '250px', layout=widgets.Layout(border='solid 1px #cccccc'))
        self.molecule_metadata_boxes = widgets.HBox([self.molecules_dropdown, self.molecule_details_textbox, self.molecule_image_box])

        self.stabilisation_time_value_floatbox = self.get_floattext_box(description='Stabilisation time', disabled=False, layout = widgets.Layout(width = '250px'), value = 0, description_width = "110px")
        self.stabilisation_time_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["sec", "min", "hrs"], value = "sec")
        self.stabilisation_time_hbox = widgets.HBox([self.stabilisation_time_value_floatbox, self.stabilisation_time_unit_dropdown])

        self.deposition_time_value_floatbox = self.get_floattext_box(description='Deposition time', disabled=False, layout = widgets.Layout(width = '250px'), value = 0, description_width = "110px")
        self.deposition_time_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["sec", "min", "hrs"], value = "sec")
        self.deposition_time_hbox = widgets.HBox([self.deposition_time_value_floatbox, self.deposition_time_unit_dropdown])
        
        self.substrate_temperature_value_floatbox = self.get_floattext_box(description='Substrate temperature', disabled=False, layout = widgets.Layout(width = '250px'), value = 0, description_width = "150px")
        self.substrate_temperature_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["oF", "oC", "K"], value = "K")
        self.substrate_temperature_hbox = widgets.HBox([self.substrate_temperature_value_floatbox, self.substrate_temperature_unit_dropdown])

        self.molecule_temperature_value_floatbox = self.get_floattext_box(description='Molecule temperature', disabled=False, layout = widgets.Layout(width = '250px'), value = 0, description_width = "150px")
        self.molecule_temperature_unit_dropdown = self.get_dropdown_box(description='', disabled=False, layout = widgets.Layout(width = '100px'), options = ["oF", "oC", "K"], value = "K")
        self.molecule_temperature_hbox = widgets.HBox([self.molecule_temperature_value_floatbox, self.molecule_temperature_unit_dropdown])

        self.evaporation_slot_value_intslider = widgets.IntSlider(value = 1, description = 'Evaporation slot', min = 1, max = 6, disabled = False, layout = widgets.Layout(width = '300px'))
        self.evaporation_slot_value_intslider.style = {'description_width': '150px'}
        self.evaporation_slot_details_textbox = widgets.Text(value = '', description = '', placeholder= "Write evaporator slot details...", disabled = False, layout = widgets.Layout(width = '300px'))
        self.evaporation_slot_hbox = widgets.HBox([self.evaporation_slot_value_intslider, self.evaporation_slot_details_textbox])
        
        self.molecule_formula_textbox = self.get_text_box(description = "Substance", disabled = False, layout = widgets.Layout(width = '350px'), placeholder = f"Write substance formula here...", description_width = "110px")
        
        if self.method_type.upper() in self.process_sample_types:
            if self.method_type == "sputtering":
                self.properties_on_left = widgets.VBox([self.duration_hbox, self.pressure_hbox, self.discharge_voltage_hbox, self.voltage_hbox])
                self.properties_on_right = widgets.VBox([self.temperature_hbox, self.angle_hbox, self.current_hbox])
            elif self.method_type == "annealing":
                self.properties_on_left = widgets.VBox([self.duration_hbox, self.pressure_hbox, self.voltage_hbox])
                self.properties_on_right = widgets.VBox([self.temperature_hbox, self.current_hbox])
            elif self.method_type == "deposition":
                self.properties_on_left = widgets.VBox([self.stabilisation_time_hbox, self.deposition_time_hbox, self.pressure_hbox])
                self.properties_on_right = widgets.VBox([self.substrate_temperature_hbox, self.molecule_temperature_hbox, self.evaporation_slot_hbox])
            elif self.method_type == "dosing":
                self.properties_on_left = widgets.VBox([self.molecule_formula_textbox, self.pressure_hbox])
                self.properties_on_right = widgets.VBox([self.pressure_hbox, self.temperature_hbox])
        
            self.method_properties = widgets.HBox([self.properties_on_left, self.properties_on_right])

        self.comments_textbox = self.get_textarea_box(description = "Comments", disabled = False, layout = widgets.Layout(width = '993px', height = '200px'), placeholder = "Write comments here...", description_width = "110px")
        self.method_name_textbox = self.get_text_box(description = "Name", disabled = False, layout = widgets.Layout(width = '400px'), placeholder = f"Write {method_type} task name here...", description_width = "110px")
        self.sample_out_name_textbox = self.get_text_box(description = "Name", disabled = False, layout = widgets.Layout(width = '400px'), placeholder = f"Write sample name here...", description_width = "110px")
        
        self.measurements_files_import = widgets.FileUpload(accept = '.sxm, .dat', multiple = True)
        
        self.folder_selector = FileChooser('.', select_default=True, use_dir_icons=True)
        self.folder_selector.show_only_dirs = True  # Set to True to display only folders

        self.create_button = widgets.Button(description = '', disabled = False, button_style = '', tooltip = 'Save', icon = 'save', layout = widgets.Layout(width = '100px', height = '50px'))
        self.quit_button = widgets.Button(description = '', disabled = False, button_style = '', tooltip = 'Main menu', icon = 'home', layout = widgets.Layout(width = '100px', height = '50px'))

        self.bottom_buttons_hbox = widgets.HBox([self.create_button, self.quit_button])

        self.increase_buttons_size = HTML("""
        <style>
            .fa-save {
                font-size: 2em !important; /* Increase icon size */
            }
            .fa-home {
                font-size: 2em !important; /* Increase icon size */
            }
        </style>
        """)
        
        self.openbis_connection_status_htmlbox = widgets.HTML(value = '')
        self.open_notebooks_html_disable_code = '''
        <div id="box" style="display: flex; pointer-events: none;">
            <div style="flex: 1; padding: 10px;">
                <ul>
                    <li><a href="create_sample.ipynb" style="color: gray;" target="_blank">Sample creation</a></li>
                    <li><a href="sputtering.ipynb" style="color: gray;" target="_blank">Sputtering</a></li>
                    <li><a href="annealing.ipynb" style="color: gray;" target="_blank">Annealing</a></li>
                    <li><a href="deposition.ipynb" style="color: gray;" target="_blank">Molecule deposition</a></li>
                    <li><a href="measurement.ipynb" style="color: gray;" target="_blank">Sample measurement</a></li>
                    <li><a href="" style="color: gray;" target="_blank">Other sample preparation task</a></li>
                    <li><a href="dosing.ipynb" style="color: gray;" target="_blank">Dosing</a></li>
                    <li><a href="" style="color: gray;" target="_blank">Sublimation rate determination</a></li>
                    <li><a href="" style="color: gray;" target="_blank">Comments/Notes</a></li>
                    <li><a href="" style="color: gray;" target="_blank">Calibration/Optimisation</a></li>
                </ul>
            </div>
            <div style="flex: 1; padding: 10px;">
                <ul>
                    <li><a href="" style="color: gray;" target="_blank">Edit openBIS objects</a></li>
                </ul>
            </div>
        </div>
        '''

        self.open_notebooks_html_enable_code = '''
            <div style="display: flex;">
                <div style="flex: 1; padding: 10px;">
                    <ul>
                        <li><a href="create_sample.ipynb" target="_blank">Sample creation</a></li>
                        <li><a href="sputtering.ipynb" target="_blank">Sputtering</a></li>
                        <li><a href="annealing.ipynb" target="_blank">Annealing</a></li>
                        <li><a href="deposition.ipynb" target="_blank">Molecule deposition</a></li>
                        <li><a href="measurement.ipynb" target="_blank">Sample measurement</a></li>
                        <li><a href="" target="_blank">Other sample preparation task</a></li>
                        <li><a href="dosing.ipynb" target="_blank">Dosing</a></li>
                        <li><a href="" target="_blank">Sublimation rate determination</a></li>
                        <li><a href="" target="_blank">Comments/Notes</a></li>
                        <li><a href="" target="_blank">Calibration/Optimisation</a></li>
                    </ul>
                </div>
                <div style="flex: 1; padding: 10px;">
                    <ul>
                        <li><a href="" target="_blank">Edit openBIS objects</a></li>
                    </ul>
                </div>
            </div>
            '''
              
        self.open_notebooks_htmlbox = widgets.HTML(self.open_notebooks_html_disable_code)
    
    @staticmethod
    def get_check_box(description, disabled, layout, indent = False, value = ''):
        check_box = widgets.Checkbox(value = value, description = description, disabled = disabled, layout = layout, indent = indent)
        return check_box
    
    @staticmethod
    def get_textarea_box(description, disabled, layout, placeholder = '', description_width = '', value = ''):
        textarea_box = widgets.Textarea(value = value, description = description, disabled = disabled, layout = layout, placeholder = placeholder)
        if description_width != '':
            textarea_box.style = {'description_width': description_width}
        
        return textarea_box
    
    @staticmethod
    def get_text_box(description, disabled, layout, placeholder = '', description_width = '', value = ''):
        text_box = widgets.Text(value = value, description = description, disabled = disabled, layout = layout, placeholder = placeholder)
        if description_width != '':
            text_box.style = {'description_width': description_width}
        
        return text_box
        
    @staticmethod
    def get_radio_button(description, disabled, layout, description_width = '', options = None):
        if options is None:
            radio_button = widgets.RadioButtons(description = description, disabled = disabled, layout = layout)
        else:
            radio_button = widgets.RadioButtons(description = description, disabled = disabled, layout = layout, options = options)
        if description_width != '':
            radio_button.style = {'description_width': description_width}
        return radio_button
    
    @staticmethod
    def get_dropdown_box(description, disabled, layout, description_width = '', options = None, value = ''):
        if options is None:
            dropdown_box = widgets.Dropdown(description=description, disabled=disabled, layout = layout)
        else:
            if value == '':
                dropdown_box = widgets.Dropdown(description=description, disabled=disabled, layout = layout, options = options)
            else:
                dropdown_box = widgets.Dropdown(description=description, disabled=disabled, layout = layout, options = options, value = value)
            
        if description_width == '':
            dropdown_box.style = {'description_width': description_width}
            
        return dropdown_box
    
    @staticmethod
    def get_floattext_box(description, disabled, layout, description_width = '', value = 0):
        floattext_box = widgets.FloatText(description=description, disabled=disabled, layout = layout, value = value)
        if description_width != '':
            floattext_box.style = {'description_width': description_width}
        return floattext_box

    @staticmethod
    def read_file(filename):
        file = open(filename, "rb")
        return file.read()
    
    @staticmethod
    def sort_dropdown(df, columns, ascending_columns):
        df = df.sort_values(by = columns, ascending = ascending_columns)
        return list(df.itertuples(index = False, name = None))
    
    def sort_materials_dropdown(self):
        dropdown_list = self.materials_dropdown.options[1:] # Default -1 message should not be sorted.
        df = pd.DataFrame(dropdown_list, columns = ['Name', 'PermID'])
        
        if self.material_sorting_checkboxes.children[1].value == True and self.material_sorting_checkboxes.children[2].value == True:
            dropdown_list = self.sort_dropdown(df, ['Name', 'PermID'], [True, False])
            
        elif self.material_sorting_checkboxes.children[1].value == True:
            dropdown_list = self.sort_dropdown(df, ['Name'], [True])
            
        elif self.material_sorting_checkboxes.children[2].value == True:
            dropdown_list = self.sort_dropdown(df, ['PermID'], [False])
            
        dropdown_list.insert(0, self.materials_dropdown.options[0])
        self.materials_dropdown.options = dropdown_list
    
    def sort_materials_dropdown_on_change(self, change):
        self.sort_materials_dropdown()
        
    def sort_samples_dropdown(self):
        dropdown_list = list(self.samples_dropdown.options[1:]) # Default -1 message should not be sorted.
        df = pd.DataFrame(dropdown_list, columns = ['Name', 'PermID'])
        
        if self.sample_sorting_checkboxes.children[1].value == True and self.sample_sorting_checkboxes.children[2].value == True:
            dropdown_list = self.sort_dropdown(df, ['Name', 'PermID'], [True, False])
            
        elif self.sample_sorting_checkboxes.children[1].value == True:
            dropdown_list = self.sort_dropdown(df, ['Name'], [True])
            
        elif self.sample_sorting_checkboxes.children[2].value == True:
            dropdown_list = self.sort_dropdown(df, ['PermID'], [False])
        
        dropdown_list.insert(0, self.samples_dropdown.options[0])
        self.samples_dropdown.options = dropdown_list

    def sort_samples_dropdown_on_change(self, change):
        self.sort_samples_dropdown()
    
    def connect_openbis(self):
        eln_config = Path.home() / ".aiidalab" / "aiidalab-eln-config.json"
        eln_config.parent.mkdir(
            parents=True, exist_ok=True
        )  # making sure that the folder exists.

        config = read_json(eln_config)
        eln_url = "https://local.openbis.ch"
        self.session_data = {"url": eln_url, "token": config[eln_url]["token"]}
        self.openbis_session = Openbis(eln_url, verify_certificates = False)
        self.openbis_session.set_token(self.session_data["token"])

    def create_openbis_object(self, object_type, collection_id, object_props, object_parents):
        openbis_object = self.openbis_session.new_object(
            type = object_type, 
            collection = collection_id,
            props = object_props,
            parents = object_parents
        )
        openbis_object.save()
        return openbis_object

    def create_openbis_collection(self, collection_code, collection_type, collection_project, collection_props):
        measurements_collection = self.openbis_session.new_collection(
            code = collection_code,
            type = collection_type,
            project = collection_project,
            props = collection_props
        )
        measurements_collection.save()
        return measurements_collection
    
    # Function to create sample object inside openBIS using information selected in the app
    def create_sample_action(self, b):
        samples = self.openbis_session.get_objects(type = "SAMPLE")
        samples_names = [sample.props["$name"] for sample in samples]
        
        if self.sample_out_name_textbox.value in samples_names:
            display(Javascript(f"alert('{'Sample name already exists!'}')"))
        else:
            if self.materials_dropdown.value == -1 or self.material_selection_radio_button.value == "No material":
                sample_parents = []
            else:
                sample_parents = [self.materials_dropdown.value]
            
            sample_props = {"$name": self.sample_out_name_textbox.value}
            _ = self.create_openbis_object("SAMPLE", self.samples_collection_openbis_path, sample_props, sample_parents)
            print("Upload successful!")
    
    # Function to handle changes in the materials dropdown
    def load_material_metadata(self, change):
        if self.materials_dropdown.value == -1:
            self.material_details_textbox.value = ''
            self.material_image_box.value = self.read_file("images/white_screen.jpg")
        else:
            if self.material_selection_radio_button.value == "Crystal":
                property_list = [("Name", "$name"), ("Material", "material"), ("Face", "face"), ("Sample Plate", "sample_plate"), ("Diameter", "diameter"), ("Height", "height"), ("Specifications", "specifications")]
            elif self.material_selection_radio_button.value == "Wafer substrate":
                property_list = [("Name", "$name"), ("Substrates", "substrates"), ("Sample Plate", "sample_plate"), ("Diameter", "diameter"), ("Height", "height"), ("Thickness", "thickness")]
            elif self.material_selection_radio_button.value == "2D-layer material":
                property_list = [("Name", "$name"), ("Sample Plate", "sample_plate"), ("2D Layers", "layers_2d"), ("Substrates", "substrates"), ("Width", "width"), ("Length", "length"), ("Thickness", "thickness")]
            
            material_object = self.openbis_session.get_object(self.materials_dropdown.value)
            material_dataset = material_object.get_datasets()[0]
            
            if material_dataset is None:
                self.material_image_box.value = self.read_file("images/white_screen.jpg")
            else:
                material_dataset.download(destination = "images")
                material_dataset_filenames = material_dataset.file_list
                material_image_filepath = material_dataset_filenames[0]
                self.material_image_box.value = self.read_file(f"images/{material_dataset.permId}/{material_image_filepath}")
            
            material_metadata = material_object.props.all()
            material_metadata_string = ""
            for property in property_list:
                if property[1] in ["diameter", "height", "width", "thickness"]:
                    if material_metadata[property[1]] is None:
                        material_metadata_string = f"{material_metadata_string} {property[0]}: {material_metadata[property[1]]}\n"
                    else:
                        property_dict = json.loads(material_metadata[property[1]])
                        material_metadata_string = f"{material_metadata_string} {property[0]}: {property_dict['value']} {property_dict['unit']}\n"
                else:
                    material_metadata_string = f"{material_metadata_string} {property[0]}: {material_metadata[property[1]]}\n"
            self.material_details_textbox.value = material_metadata_string

    # Function to handle changes in the radio button
    def select_material_radio_change(self, change):
        self.material_details_textbox.value = ''
        
        if self.material_selection_radio_button.value == "No material":
            clear_output()
            display(self.material_selection_radio_button)
        else: 
            if self.material_selection_radio_button.value == "Crystal":
                materials = self.openbis_session.get_objects(type = "CRYSTAL")
                materials_placeholder = "Select crystal..."

            elif self.material_selection_radio_button.value == "Wafer substrate":
                materials = self.openbis_session.get_objects(type = "WAFER_SUBSTRATE")
                materials_placeholder = "Select wafer substrate..."
                
            elif self.material_selection_radio_button.value == "2D-layer material":
                materials = self.openbis_session.get_objects(type = "2D_LAYERED_MATERIAL")
                materials_placeholder = "Select 2D-layer material..."
            
            materials_names_permids = [(f"{material.props['$name']} ({material.permId})", material.permId) for material in materials]
            materials_names_permids.insert(0, (materials_placeholder, -1))
            self.materials_dropdown.options = materials_names_permids
            self.materials_dropdown.value = -1
            
            # Dropdown list must be sorted when changing material type
            self.sort_materials_dropdown()
            
            self.materials_dropdown.observe(self.load_material_metadata, names = 'value')
            
            self.material_sorting_checkboxes.children[1].observe(self.sort_materials_dropdown_on_change, names = 'value')
            self.material_sorting_checkboxes.children[2].observe(self.sort_materials_dropdown_on_change, names = 'value')
            
            clear_output()
            display(self.material_selection_radio_button)
            display(widgets.HBox([self.materials_dropdown, self.material_sorting_checkboxes]))
            display(self.material_metadata_boxes)
            
        display(self.increase_buttons_size)
    
    def upload_datasets(self, method_object):
        for filename, file_info in self.support_files_uploader.value.items():
            with open(filename, 'wb') as f:  # Save file content
                f.write(file_info['content'])
                
            ds_new = self.openbis_session.new_dataset(
                type       = 'RAW_DATA',
                sample     = method_object,
                files      = [filename]
            )
            ds_new.save()
            
            os.remove(filename)
    
    def create_process_action(self,b):
        samples = self.openbis_session.get_objects(type = "SAMPLE")
        samples_names = [sample.props["$name"] for sample in samples]
        
        if self.sample_out_name_textbox.value in samples_names:
            display(Javascript(f"alert('{'Sample name already exists!'}')"))
        else:
            if self.experiments_dropdown.value == -1:
                print("Select an experiment.")
            elif self.samples_dropdown.value == -1:
                print("Select a sample.")
            elif self.instruments_dropdown.value == -1:
                print("Select an instrument.")
            elif self.molecules_dropdown.value == -1 and self.method_type == "deposition":
                print("Select a molecule.")
            else:
                if self.method_type == "deposition":
                    sample_parents = [self.samples_dropdown.value, self.instruments_dropdown.value, self.molecules_dropdown.value]
                else:
                    sample_parents = [self.samples_dropdown.value, self.instruments_dropdown.value]
                
                object_properties = {}
                object_properties["$name"] = self.method_name_textbox.value
                object_properties["comments"] = self.comments_textbox.value
                object_properties["pressure"] = json.dumps({"value": self.pressure_value_floatbox.value, "unit": self.pressure_unit_dropdown.value})
                
                if self.method_type in ["sputtering", "annealing", "dosing"]:
                    object_properties["duration"] = json.dumps({"value": self.duration_value_floatbox.value, "unit": self.duration_unit_dropdown.value})
                    object_properties["temperature"] = json.dumps({"value": self.temperature_value_floatbox.value, "unit": self.temperature_unit_dropdown.value})
                    
                    if self.method_type in ["sputtering", "annealing"]:
                        object_properties["voltage"] = json.dumps({"value": self.voltage_value_floatbox.value, "unit": self.voltage_unit_dropdown.value})
                        object_properties["current"] = json.dumps({"value": self.current_value_floatbox.value, "unit": self.current_unit_dropdown.value})

                        if self.method_type == "sputtering":
                            object_properties["discharge_voltage"] = json.dumps({"value": self.discharge_voltage_value_floatbox.value, "unit": self.discharge_voltage_unit_dropdown.value})
                            object_properties["angle"] = json.dumps({"value": self.angle_value_floatbox.value, "unit": self.angle_unit_dropdown.value})

                elif self.method_type in ["deposition"]:
                    object_properties["stabilisation_time"] = json.dumps({"value": self.stabilisation_time_value_floatbox.value, "unit": self.stabilisation_time_unit_dropdown.value})
                    object_properties["deposition_time"] = json.dumps({"value": self.deposition_time_value_floatbox.value, "unit": self.deposition_time_unit_dropdown.value})
                    object_properties["substrate_temperature"] = json.dumps({"value": self.substrate_temperature_value_floatbox.value, "unit": self.substrate_temperature_unit_dropdown.value})
                    object_properties["molecule_temperature"] = json.dumps({"value": self.molecule_temperature_value_floatbox.value, "unit": self.molecule_temperature_unit_dropdown.value})
                    object_properties["evaporator_slot"] = json.dumps({"evaporator_number": self.evaporation_slot_value_intslider.value, "details": self.evaporation_slot_details_textbox.value})
                
                method_object = self.create_openbis_object(self.method_type.upper(), self.experiments_dropdown.value, object_properties, sample_parents)
                self.upload_datasets(method_object)
                sample_props = {"$name": self.sample_out_name_textbox.value}
                _ = self.create_openbis_object("SAMPLE", self.samples_collection_openbis_path, sample_props, [method_object])
                print("Upload successful!")
    
    # Function to handle changes in the materials dropdown
    def load_molecule_metadata(self, change):
        
        if self.molecules_dropdown.value == -1:
            self.molecule_details_textbox.value = ''
            self.molecule_image_box.value = self.read_file("images/white_screen.jpg")
        else:
            property_list = [("Name", "$name"), ("IUPAC name", "iupac_name"), ("Sum formula", "sum_formula"), ("SMILES", "smiles"), ("Empa number", "empa_number"), ("Batch", "batch"), ("Amount", "amount"), ("Comments", "comments")]
            
            material_object = self.openbis_session.get_object(self.molecules_dropdown.value)
            material_dataset = material_object.get_datasets(type = "ELN_PREVIEW")[0]
            
            if material_dataset is None:
                self.molecule_image_box.value = self.read_file("images/white_screen.jpg")
            else:
                material_dataset.download(destination = "images")
                material_dataset_filenames = material_dataset.file_list
                material_image_filepath = material_dataset_filenames[0]
                self.molecule_image_box.value = self.read_file(f"images/{material_dataset.permId}/{material_image_filepath}")
            
            material_metadata = material_object.props.all()
            material_metadata_string = ""
            for property in property_list:
                if property[1] in ["amount"]:
                    if material_metadata[property[1]] is None:
                        material_metadata_string = f"{material_metadata_string} {property[0]}: {material_metadata[property[1]]}\n"
                    else:
                        property_dict = json.loads(material_metadata[property[1]])
                        material_metadata_string = f"{material_metadata_string} {property[0]}: {property_dict['value']} {property_dict['unit']}\n"
                else:
                    material_metadata_string = f"{material_metadata_string} {property[0]}: {material_metadata[property[1]]}\n"
            self.molecule_details_textbox.value = material_metadata_string

    def load_sample_metadata(self, change):
        if self.samples_dropdown.value == -1:
            self.experiments_dropdown.value = -1
            self.instruments_dropdown.value = -1
            self.sample_details_textbox.value = ''
            
            if self.method_type.upper() in self.process_sample_types:
                self.sample_out_name_textbox.value = ''
        else:
            sample_object = self.openbis_session.get_object(self.samples_dropdown.value)
            sample_parents_metadata = []
            sample_parents_metadata = self.get_parents_recursive(sample_object, sample_parents_metadata)
            
            last_sample_process = None
            sample_processes_strings = []
            sample_materials_strings = []
            number_parents = len(sample_parents_metadata)
            parent_idx = 0
            while parent_idx < number_parents:
                parent_metadata = sample_parents_metadata[parent_idx]
                if parent_metadata[0] == "DEPOSITION":
                    next_parent_metadata = sample_parents_metadata[parent_idx + 1]
                    sample_metadata_string = f"> {parent_metadata[0]} ({parent_metadata[3]}, {parent_metadata[1]}, {parent_metadata[2]}) [{next_parent_metadata[0]} ({next_parent_metadata[3]}, {next_parent_metadata[1]}, {next_parent_metadata[4]})]"
                    parent_idx += 1
                else:
                    sample_metadata_string = f"> {parent_metadata[0]} ({parent_metadata[3]}, {parent_metadata[1]}, {parent_metadata[2]})"
                
                if parent_metadata[0] in self.process_sample_types:
                    sample_processes_strings.append(sample_metadata_string)
                    # Get the last sample preparation method performed on the sample in order to search the correct experiment where the sample is being used.
                    if last_sample_process is None:
                        last_sample_process = parent_metadata[1]
                else:
                    sample_materials_strings.append(sample_metadata_string)
                    
                parent_idx += 1
            
            sample_metadata_string = f"PermId: {sample_object.attrs.permId}\nMaterial:"
            
            for material_string in sample_materials_strings:
                sample_metadata_string = f"{sample_metadata_string}\n{material_string}"
                
            sample_metadata_string = f"{sample_metadata_string}\nProcesses:"
            for process_string in sample_processes_strings:
                sample_metadata_string = f"{sample_metadata_string}\n{process_string}"

            if last_sample_process is not None:
                last_sample_process_object = self.openbis_session.get_object(last_sample_process)
                
                if self.method_type.upper() in self.process_sample_types:
                    # Automatically select the experiment where the last sample process task was saved
                    last_sample_process_experiment_id = last_sample_process_object.attrs.experiment
                    last_sample_process_experiment = self.openbis_session.get_experiment(last_sample_process_experiment_id)
                    self.experiments_dropdown.value =  last_sample_process_experiment.permId
                
                # Automatically select the instrument used in the last sample process task
                last_sample_process_object_parents = last_sample_process_object.parents
                for parent in last_sample_process_object_parents:
                    parent_object = self.openbis_session.get_object(parent)
                    if parent_object.type == "INSTRUMENT":
                        self.instruments_dropdown.value = parent_object.permId
            
            self.sample_details_textbox.value = sample_metadata_string
            
            if self.method_type.upper() in self.process_sample_types:
                self.sample_out_name_textbox.value = f"{sample_object.props['$name']}_{self.method_name_textbox.value}"
    
    def load_experiment_list(self):
        experiments = self.openbis_session.get_collections(type = "EXPERIMENT")
        experiments_names_permids = [(f"{experiment.props['$name']} ({experiment.attrs.identifier})", experiment.permId) for experiment in experiments]
        experiments_names_permids.insert(0, ('Select an experiment...', -1))
        self.experiments_dropdown.options = experiments_names_permids
        self.experiments_dropdown.value = -1
    
    def load_sample_list(self):
        samples = self.openbis_session.get_objects(type = "SAMPLE")
        samples_names_permids = [(sample.props['$name'], sample.permId) for sample in samples]
        samples_names_permids.insert(0, ('Select an input sample...', -1))
        self.samples_dropdown.options = samples_names_permids
        self.samples_dropdown.value = -1
        self.sort_samples_dropdown()
        
        self.sample_sorting_checkboxes.children[1].observe(self.sort_samples_dropdown_on_change, names = 'value')
        self.sample_sorting_checkboxes.children[2].observe(self.sort_samples_dropdown_on_change, names = 'value')
    
    def load_instrument_list(self):
        instruments = self.openbis_session.get_objects(type = "INSTRUMENT")
        instruments_names_permids = [(f"{instrument.props['$name']} ({instrument.attrs.permId})", instrument.permId) for instrument in instruments]
        instruments_names_permids.insert(0, ('Select an instrument...', -1))
        self.instruments_dropdown.options = instruments_names_permids
        self.instruments_dropdown.value = -1
    
    def load_molecule_list(self):
        molecules = self.openbis_session.get_objects(type = "MOLECULE")
        molecules_names_permids = [(f"{molecule.props['$name']} ({molecule.attrs.permId})", molecule.permId) for molecule in molecules]
        # molecules_empa_numbers = [np.nan if molecule.props['empa_number'] is None else int(molecule.props['empa_number']) for molecule in molecules]
        # molecules_batches = [np.nan if molecule.props['batch'] is None else molecule.props['batch'] for molecule in molecules]
        # molecules_names_permids = [x for _, _, x in sorted(zip(molecules_empa_numbers, molecules_batches, molecules_names_permids), key = lambda item: (item[0], item[1]), reverse=True)]
        molecules_names_permids.insert(0, ('Select a molecule...', -1))
        self.molecules_dropdown.options = molecules_names_permids
        self.molecules_dropdown.value = -1
    
    def load_dropdown_lists(self):
        # Populate dropdown lists
        self.load_sample_list()
        self.load_instrument_list()
        
        if self.method_type.upper() in self.process_sample_types:
            self.load_experiment_list()
            
            if self.method_type.upper() == "DEPOSITION":
                self.load_molecule_list()

    def get_parents_recursive(self, object, object_parents_metadata):
        if object.attrs.type in self.process_sample_types or object.attrs.type in self.raw_materials_types:
            object_parents_metadata.append([object.attrs.type, object.attrs.permId, object.attrs.registrationDate, 
                                            object.props['$name'], object.props['sum_formula']])
        object_parents = object.parents
        for parent in object_parents:
            parent_object = self.openbis_session.get_object(parent)
            self.get_parents_recursive(parent_object, object_parents_metadata)
        return object_parents_metadata

    def update_text(self, change):
        if self.samples_dropdown.value != -1:
            selected_sample_name = next(label for label, val in self.samples_dropdown.options if val == self.samples_dropdown.value)
            self.sample_out_name_textbox.value = f"{selected_sample_name}_{self.method_name_textbox.value}"
    
    def upload_measurements_to_openbis(self, b):
        measurements_collection = None
        sample_object = self.openbis_session.get_object(self.samples_dropdown.value)
        for child_id in sample_object.children:
            child_object = self.openbis_session.get_object(child_id)
            if child_object.type == "1D_MEASUREMENT" or child_object.type == "2D_MEASUREMENT":
                measurements_collection = child_object.collection
                break
        
        sample_project = None
        for parent_id in sample_object.parents:
            parent_object = self.openbis_session.get_object(parent_id)
            if parent_object.type in self.process_sample_types:
                sample_project = parent_object.project
                break
        
        data_folder = self.folder_selector.selected_path
        
        correct_folder = False
        data_files = os.listdir(data_folder)
        for filename in data_files:
            filename_extension = filename.split('.')[-1]
            if f".{filename_extension}" in self.measurement_file_extensions:
                correct_folder = True
            else:
                correct_folder = False
                break
            
        if correct_folder:
            if measurements_collection:
                measurements_collection = f"{sample_project.identifier}/{measurements_collection}"
                nanonis_importer.upload_measurements_into_openbis(self.session_data['url'], data_folder, measurements_collection, sample_object.permId, self.instruments_dropdown.value)
                print("Upload successful!")
                
            elif sample_project:
                collection_props = {"$name": f"Measurements from Sample {sample_object.props['$name']}", "$default_collection_view": "IMAGING_GALLERY_VIEW"}
                measurements_collection = self.create_openbis_collection(f"MEASUREMENTS_COLLECTION_{sample_object.code}", "COLLECTION", sample_project, collection_props)
                nanonis_importer.upload_measurements_into_openbis(self.session_data['url'], data_folder, measurements_collection.permId, sample_object.permId, self.instruments_dropdown.value)
                print("Upload successful!")
                
            else:
                print("Sample does not belong to any project yet.")
                
        else:
            print(f"Folder contains unrecognised files. Please remove files whose extension does not belong to the following list: .sxm, .dat.")