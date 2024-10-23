import json
from pathlib import Path
import os
import ipywidgets as ipw
from IPython.display import display, clear_output
from pybis import Openbis
import ipyfilechooser
import sys
sys.path.append('/home/jovyan/aiida-openbis/Notebooks/importer')
import nanonis_importer
import pandas as pd
import numpy as np
import yaml
import shutil
import utils

class AppWidgets():
    def __init__(self, method_type, config_filename):
        self.method_type = method_type
        self.openbis_session, self.session_data = utils.connect_openbis()
        self.config = utils.read_json(config_filename)
        
        # Necessary for refreshing the widgets needed to create new experiments
        self.select_experiment_output = ipw.Output()
        
        # Home page configuration
        if self.method_type == "Home":
            self.openbis_connection_status_htmlbox = utils.HTMLbox(value = '')
            if self.openbis_session:
                self.openbis_connection_status_htmlbox.value = self.config["home_page"]["enable_status"]
                self.open_notebooks_html_enable_code = ''.join(self.config["home_page"]["enable_links"])
                self.open_notebooks_htmlbox = utils.HTMLbox(value = self.open_notebooks_html_enable_code)
            else:
                self.openbis_connection_status_htmlbox.value = self.config["home_page"]["disable_status"]
                self.open_notebooks_html_disable_code = ''.join(self.config["home_page"]["disable_links"])
                self.open_notebooks_htmlbox = utils.HTMLbox(value = self.open_notebooks_html_disable_code)
        else:
            self.raw_materials_types = [object_info["openbis_object_type"] for _, object_info in self.config["objects"].items() if object_info["object_type"] == "raw_material"]
            self.process_sample_types = [object_info["openbis_object_type"] for _, object_info in self.config["objects"].items() if object_info["object_type"] == "process"]
            self.samples_collection_openbis_path = self.config["samples_collection_openbis_path"]
            self.measurement_file_extensions = self.config["measurement_file_extensions"]
            
            self.support_files_uploader = ipw.FileUpload(multiple = True)
            
            self.materials_dropdown = utils.Dropdown(description='', disabled=False, layout = ipw.Layout(width = '350px'))
            self.material_details_textbox = utils.Textarea(description = "", disabled = True, layout = ipw.Layout(width = '425px', height = '200px'))
            self.material_image_box = utils.Image(value = open("images/white_screen.jpg", "rb").read(), format = 'jpg', width = '200px', height = '300px', layout=ipw.Layout(border='solid 1px #cccccc'))
            self.material_metadata_boxes = ipw.HBox([self.material_details_textbox, self.material_image_box])
            raw_materials_options = [object_key for object_key, object_info in self.config["objects"].items() if object_info["object_type"] == "raw_material"]
            raw_materials_options.insert(0, "No material")
            self.material_selection_radio_button = utils.Radiobuttons(description = '', options = raw_materials_options, disabled = False, layout = ipw.Layout(width = '300px'), style = {'description_width': "100px"})
            self.material_sorting_checkboxes = utils.SortingCheckboxes("50px", "60px", "200px")
            
            self.samples_dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget('Sample', ipw.Layout(width = '385px'), {'description_width': "110px"}, [-1])
            self.sample_details_textbox = utils.Textarea(description = "", disabled = True, layout = ipw.Layout(width = '589px', height = '300px'))
            self.sample_metadata_boxes = ipw.HBox([self.samples_dropdown_boxes, self.sample_details_textbox])

            self.instruments_dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget('Instrument', ipw.Layout(width = '982px'), {'description_width': "110px"}, [-1])
            self.projects_dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget('Project', ipw.Layout(width = '982px'), {'description_width': "110px"}, [-1])
            
            self.create_new_experiment_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Add', icon = 'plus', layout = ipw.Layout(width = '50px', height = '25px'))
            self.save_new_experiment_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Save', icon = 'save', layout = ipw.Layout(width = '50px', height = '35px', margin = '0 0 0 90px'))
            self.cancel_new_experiment_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Cancel', icon = 'times', layout = ipw.Layout(width = '50px', height = '35px', margin = '0 0 0 5px'))
            self.new_experiment_name_textbox = utils.Text(description = "Name", disabled = False, layout = ipw.Layout(width = '400px'), placeholder = f"Write experiment name here...", style = {'description_width': "110px"})
            
            self.experiments_dropdown = utils.Dropdown(description='Experiment', disabled=False, layout = ipw.Layout(width = '993px'), style = {'description_width': "110px"}, options = [-1])
            self.experiment_sorting_checkboxes = utils.SortingCheckboxes("130px", "60px", "200px")
            self.experiments_dropdown_details = ipw.HBox([self.experiments_dropdown, self.create_new_experiment_button])
            self.experiments_dropdown_boxes = ipw.VBox([self.experiments_dropdown_details, self.experiment_sorting_checkboxes])
            
            self.molecules_dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget('Substance', ipw.Layout(width = '335px'), {'description_width': "110px"}, [-1])
            self.molecule_details_textbox = utils.Textarea(description = "", disabled = True, layout = ipw.Layout(width = '415px', height = '250px'))
            self.molecule_image_box = utils.Image(value = open("images/white_screen.jpg", "rb").read(), format = 'jpg', width = '220px', height = '250px', layout=ipw.Layout(border='solid 1px #cccccc'))
            self.molecule_metadata_boxes = ipw.HBox([self.molecules_dropdown_boxes, self.molecule_details_textbox, self.molecule_image_box])

            self.method_name_textbox = utils.Text(description = "Name", disabled = False, layout = ipw.Layout(width = '400px'), placeholder = f"Write task name here...", style = {'description_width': "150px"})
            self.comments_textbox = utils.Textarea(description = "Comments", disabled = False, layout = ipw.Layout(width = '993px', height = '200px'), placeholder = "Write comments here...", style = {'description_width': "150px"})
            
            # Quantity properties widgets
            self.property_widgets = {}
            for prop in self.config["properties"].keys():
                if self.config["properties"][prop]["property_type"] == "quantity_value":
                    self.property_widgets[prop] = utils.FloatTextwithDropdownWidget(
                        self.config["properties"][prop]["title"], ipw.Layout(width = self.config["properties"][prop]["box_layout"]["width"]), 
                        0, {'description_width': self.config["properties"][prop]["box_layout"]["description_width"]}, 
                        ipw.Layout(width = self.config["properties"][prop]["dropdown_layout"]["width"]), self.config["properties"][prop]["units"], 
                        self.config["properties"][prop]["units"][0]
                    )
            
            self.property_widgets["sum_formula"] = utils.Text(
                description = self.config["properties"]["sum_formula"]["title"], 
                disabled = False, 
                layout = ipw.Layout(
                    width = self.config["properties"]["sum_formula"]["box_layout"]["width"]), 
                    placeholder = self.config["properties"]["sum_formula"]["placeholder"],
                    style = {'description_width': self.config["properties"]["sum_formula"]["box_layout"]["description_width"]}
                )
            
            self.property_widgets["evaporator_slot"] = utils.IntSliderwithTextWidget(
                1, self.config["properties"]["evaporator_slot"]["title"], [1,6], 
                ipw.Layout(width = self.config["properties"]["evaporator_slot"]["slider_layout"]["width"]), 
                {'description_width': self.config["properties"]["evaporator_slot"]["slider_layout"]["description_width"]}, 
                self.config["properties"]["evaporator_slot"]["placeholder"], 
                ipw.Layout(width = self.config["properties"]["evaporator_slot"]["box_layout"]["width"])
            )
            
            if self.method_type in self.config["objects"].keys():
                items = [self.property_widgets[prop] for prop in self.config["objects"][self.method_type]["properties"]]
                self.method_properties = ipw.VBox(items)

            self.sample_out_name_textbox = utils.Text(description = "Name", disabled = False, layout = ipw.Layout(width = '400px'), placeholder = f"Write sample name here...", style = {'description_width': "110px"})
            self.create_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Save', icon = 'save', layout = ipw.Layout(width = '100px', height = '50px'))
            self.quit_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Main menu', icon = 'home', layout = ipw.Layout(width = '100px', height = '50px'))
            self.bottom_buttons_hbox = ipw.HBox([self.create_button, self.quit_button])
            
            # Widget to select folder with measurement files
            self.folder_selector = utils.FileChooser(path = '.', select_default=True, use_dir_icons=True, show_only_dirs = True)
            
            # Assign functions to widgets
            self.create_new_experiment_button.on_click(self.create_new_experiment_button_on_click)
            self.cancel_new_experiment_button.on_click(self.cancel_new_experiment_button_on_click)
            self.save_new_experiment_button.on_click(self.save_new_experiment_button_on_click)
            
            # Increase buttons icons' size
            self.increase_buttons_size = utils.HTML(data = ''.join(self.config["save_home_buttons_settings"]))
            display(self.increase_buttons_size)
    
    def create_new_experiment_button_on_click(self, b):
        with self.select_experiment_output:
            clear_output()
            utils.load_openbis_elements_list(self.openbis_session, "PROJECT", self.projects_dropdown_boxes.children[0], self.projects_dropdown_boxes.children[1], "project")
            display_list = [self.experiments_dropdown_boxes, self.new_experiment_name_textbox, self.projects_dropdown_boxes,
                            ipw.HBox([self.save_new_experiment_button, self.cancel_new_experiment_button]), 
                            self.sample_metadata_boxes, self.instruments_dropdown_boxes]
            if self.method_type == "Deposition":
                display_list.append(self.molecule_metadata_boxes)
            display(ipw.VBox(display_list))
    
    def cancel_new_experiment_button_on_click(self, b):
        with self.select_experiment_output:
            clear_output()
            display_list = [self.experiments_dropdown_boxes, self.sample_metadata_boxes, self.instruments_dropdown_boxes]
            if self.method_type == "Deposition":
                display_list.append(self.molecule_metadata_boxes)
            display(ipw.VBox(display_list))
    
    def save_new_experiment_button_on_click(self, b):
        utils.create_experiment_in_openbis(self.projects_dropdown_boxes.children[0].value, self.new_experiment_name_textbox.value)
        with self.select_experiment_output:
            clear_output()
            utils.load_openbis_elements_list(self.openbis_session, "EXPERIMENT", self.experiments_dropdown, self.experiment_sorting_checkboxes, "experiment")
            display_list = [self.experiments_dropdown_boxes, self.sample_metadata_boxes, self.instruments_dropdown_boxes]
            if self.method_type == "Deposition":
                display_list.append(self.molecule_metadata_boxes)
            display(ipw.VBox(display_list))
    
    # Function to create sample object inside openBIS using information selected in the app
    def create_sample_action(self, b):
        samples_names = [sample.props["$name"] for sample in self.openbis_session.get_objects(type = "SAMPLE")]
        if self.sample_out_name_textbox.value in samples_names:
            display(utils.Javascript(data = f"alert('{'Sample name already exists!'}')"))
            return
        else:
            sample_parents = [] if self.materials_dropdown.value == -1 or self.material_selection_radio_button.value == "No material" else [self.materials_dropdown.value]
            sample_props = {"$name": self.sample_out_name_textbox.value, "exists": True}
            utils.create_openbis_object(type="SAMPLE", collection=self.samples_collection_openbis_path, props=sample_props, parents=sample_parents)
            print("Upload successful!")
    
    # Function to handle changes in the materials dropdown
    def load_material_metadata(self, change):
        if self.materials_dropdown.value == -1:
            self.material_details_textbox.value = ''
            self.material_image_box.value = utils.read_file("images/white_screen.jpg")
            return
        
        # Get selected object properties information from config file
        selected_object = self.material_selection_radio_button.value
        selected_object_properties = self.config["objects"][selected_object]["properties"]
        
        # Get material object information and dataset
        material_object = self.openbis_session.get_object(self.materials_dropdown.value)
        material_dataset = material_object.get_datasets()[0]
        
        # Get the object image preview
        if material_dataset:
            material_dataset.download(destination="images")
            self.material_image_box.value = utils.read_file(f"images/{material_dataset.permId}/{material_dataset.file_list[0]}")
            shutil.rmtree(f"images/{material_dataset.permId}")
        else:
            self.material_image_box.value = utils.read_file("images/white_screen.jpg")

        # Make a string with the property values of the object
        material_metadata = material_object.props.all()
        material_metadata_string = ""
        for prop_key in selected_object_properties:
            prop_title = self.config["properties"][prop_key]["title"]
            if self.config["properties"][prop_key]["property_type"] == "quantity_value":
                value = material_metadata.get(prop_key)
                if value:
                    prop_dict = json.loads(value)
                    material_metadata_string += f"{prop_title}: {prop_dict['value']} {prop_dict['unit']}\n"
                else:
                    material_metadata_string += f"{prop_title}: {value}\n"
            else:
                material_metadata_string += f"{prop_title}: {material_metadata.get(prop_key)}\n"

        self.material_details_textbox.value = material_metadata_string
    
    def select_material_radio_change(self, change):
        self.material_details_textbox.value = ''
        clear_output()
        display(self.material_selection_radio_button)
        
        material_types = {}
        for object_key, object_info in self.config["objects"].items():
            if object_info["object_type"] == "raw_material":
                material_types[object_key] = (object_info["openbis_object_type"], object_info["placeholder"])

        material_type = self.material_selection_radio_button.value
        if material_type == "No material":
            return

        material_class, placeholder = material_types.get(material_type, (None, None))
        if material_class:
            materials = self.openbis_session.get_objects(type = material_class)
            materials_names_permids = [(f"{mat.props['$name']} ({mat.permId})", mat.permId) for mat in materials]
            self.materials_dropdown.options = [(placeholder, -1)] + materials_names_permids
            self.materials_dropdown.value = -1

            utils.sort_dropdown(self.material_sorting_checkboxes, self.materials_dropdown)

            self.materials_dropdown.observe(self.load_material_metadata, names='value')
            for checkbox in self.material_sorting_checkboxes.children[1:3]:
                checkbox.observe(lambda change: utils.sort_dropdown(self.material_sorting_checkboxes, self.materials_dropdown), names='value')

            display(ipw.HBox([self.materials_dropdown, self.material_sorting_checkboxes]))
            display(self.material_metadata_boxes)
    
    def upload_datasets(self, method_object):
        for file_info in self.support_files_uploader.value:
            filename = file_info['name']
            utils.save_file(file_info['content'], filename)
            self.openbis_session.new_dataset(type = 'RAW_DATA', sample = method_object, files = [filename]).save()
            os.remove(filename)
    
    def create_process_action(self, b):
        samples_names = [sample.props["$name"] for sample in self.openbis_session.get_objects(type="SAMPLE")]
        
        if self.sample_out_name_textbox.value in samples_names:
            display(utils.Javascript(data = "alert('Sample name already exists!')"))
            return
        
        if self.experiments_dropdown.value == -1:
            print("Select an experiment.")
            return
        
        if self.samples_dropdown_boxes.children[0].value == -1:
            print("Select a sample.")
            return
        
        if self.instruments_dropdown_boxes.children[0].value == -1:
            print("Select an instrument.")
            return
        
        if self.molecules_dropdown_boxes.children[0].value == -1 and self.method_type == "Deposition":
            print("Select a substance.")
            return

        # Prepare sample parents based on method type
        sample_parents = [self.samples_dropdown_boxes.children[0].value, self.instruments_dropdown_boxes.children[0].value]
        if self.method_type == "Deposition":
            sample_parents.append(self.molecules_dropdown_boxes.children[0].value)

        object_properties = {"$name": self.method_name_textbox.value, "comments": self.comments_textbox.value}
        
        for prop in self.config["objects"][self.method_type]["properties"]:
            if prop == "evaporator_slot":
                object_properties[prop] = json.dumps({"evaporator_number": self.property_widgets[prop].children[0].value, "details": self.property_widgets[prop].children[1].value})
            elif prop == "sum_formula":
                object_properties[prop] = self.property_widgets[prop]
            else:
                object_properties[prop] = json.dumps({"has_value": self.property_widgets[prop].children[0].value, "has_unit": self.property_widgets[prop].children[1].value})

        method_object = utils.create_openbis_object(type = self.config["objects"][self.method_type]["openbis_object_type"], collection = self.experiments_dropdown.value, props = object_properties, parents = sample_parents)
        self.upload_datasets(method_object)

        # Turn off sample visibility
        parent_sample = self.openbis_session.get_object(self.samples_dropdown_boxes.children[0].value)
        parent_sample.props["exists"] = False
        parent_sample.save()

        sample_props = {"$name": self.sample_out_name_textbox.value, "exists": True}
        utils.create_openbis_object(type = "SAMPLE", collection = self.samples_collection_openbis_path, props = sample_props, parents = [method_object])
        print("Upload successful!")
    
    # Function to handle changes in the substances dropdown
    def load_molecule_metadata(self, change):
        if self.molecules_dropdown_boxes.children[0].value == -1:
            self.molecule_details_textbox.value = ''
            self.molecule_image_box.value = utils.read_file("images/white_screen.jpg")
            return
        
        # Get substance metadata
        property_list = self.config["objects"]["Substance"]["properties"]
        molecule_property_list = self.config["objects"]["Molecule"]["properties"]
        material_object = self.openbis_session.get_object(self.molecules_dropdown_boxes.children[0].value)
        
        # Get molecule image
        molecule_object = self.openbis_session.get_object(material_object.props.all()['has_molecule'])
        molecule_metadata = molecule_object.props.all()
        molecule_dataset = molecule_object.get_datasets(type="ELN_PREVIEW")[0]

        if molecule_dataset:
            molecule_dataset.download(destination="images")
            material_image_filepath = molecule_dataset.file_list[0]
            self.molecule_image_box.value = utils.read_file(f"images/{molecule_dataset.permId}/{material_image_filepath}")
            shutil.rmtree(f"images/{molecule_dataset.permId}/")
        else:
            self.molecule_image_box.value = utils.read_file("images/white_screen.jpg")

        material_metadata = material_object.props.all()
        material_metadata_string = ""
        for prop_key in property_list:
            prop_title = self.config["properties"][prop_key]["title"]
            value = material_metadata.get(prop_key)
            if self.config["properties"][prop_key]["property_type"] == "quantity_value" and value is not None:
                value = json.loads(value)
                material_metadata_string += f"{prop_title}: {value['value']} {value['unit']}\n"
            elif prop_key == "has_molecule":
                material_metadata_string += f"Molecule:\n"
                for prop_key in molecule_property_list:
                    prop_title = self.config["properties"][prop_key]["title"]
                    prop_value = molecule_metadata.get(prop_key)
                    material_metadata_string += f"- {prop_title}: {prop_value}\n"
            else:
                material_metadata_string += f"{prop_title}: {value}\n"

        self.molecule_details_textbox.value = material_metadata_string

    def load_sample_metadata(self, change):
        if self.samples_dropdown_boxes.children[0].value == -1:
            self.experiments_dropdown.value = -1
            self.instruments_dropdown_boxes.children[0].value = -1
            self.sample_details_textbox.value = ''
            
            if self.method_type in self.config["objects"].keys():
                if self.config["objects"][self.method_type]["openbis_object_type"] in self.process_sample_types:
                    self.sample_out_name_textbox.value = ''
            return
        
        sample_object = self.openbis_session.get_object(self.samples_dropdown_boxes.children[0].value)
        sample_parents_metadata = utils.get_openbis_parents_recursive(self.openbis_session, sample_object, [])
        
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
            
            if self.method_type in self.config["objects"].keys():
                if self.config["objects"][self.method_type]["openbis_object_type"] in self.process_sample_types:
                    # Automatically select the experiment where the last sample process task was saved
                    last_sample_process_experiment = self.openbis_session.get_experiment(last_sample_process_object.attrs.experiment)
                    # Notify the user in case the experiment changed
                    if self.experiments_dropdown.value != last_sample_process_experiment.permId and self.experiments_dropdown.value != -1:
                        dropdown_options_dict = {key: val for val, key in self.experiments_dropdown.options}
                        previous_experiment_name = dropdown_options_dict.get(self.experiments_dropdown.value)
                        new_experiment_name = dropdown_options_dict.get(last_sample_process_experiment.permId)
                        display(utils.Javascript(data = f"alert('{f'Experiment was changed from {previous_experiment_name} to {new_experiment_name}!'}')"))
                    self.experiments_dropdown.value =  last_sample_process_experiment.permId
            
            # Automatically select the instrument used in the last sample process task
            for parent in last_sample_process_object.parents:
                parent_object = self.openbis_session.get_object(parent)
                if parent_object.type == "INSTRUMENT":
                    self.instruments_dropdown_boxes.children[0].value = parent_object.permId
        
        self.sample_details_textbox.value = sample_metadata_string
        if self.method_type in self.config["objects"].keys():
            if self.config["objects"][self.method_type]["openbis_object_type"] in self.process_sample_types:
                sample_name = sample_object.props['$name']
                self.sample_out_name_textbox.value = f"{sample_name}_{self.method_name_textbox.value}" if self.method_name_textbox.value else sample_name
    
    def load_dropdown_lists(self):
        # Populate dropdown lists
        utils.load_openbis_elements_list(self.openbis_session, "SAMPLE", self.samples_dropdown_boxes.children[0], self.samples_dropdown_boxes.children[1], "sample")
        utils.load_openbis_elements_list(self.openbis_session, "INSTRUMENT", self.instruments_dropdown_boxes.children[0], self.instruments_dropdown_boxes.children[1], "instrument")
        if self.method_type in self.config["objects"].keys():
            if self.config["objects"][self.method_type]["openbis_object_type"] in self.process_sample_types:
                utils.load_openbis_elements_list(self.openbis_session, "EXPERIMENT", self.experiments_dropdown, self.experiment_sorting_checkboxes, "experiment")
                if self.config["objects"][self.method_type]["openbis_object_type"] == "DEPOSITION":
                    utils.load_openbis_elements_list(self.openbis_session, "SUBSTANCE", self.molecules_dropdown_boxes.children[0], self.molecules_dropdown_boxes.children[1], "substance")

    def update_text(self, change):
        if self.samples_dropdown_boxes.children[0].value != -1:
            selected_sample_name = next(label for label, val in self.samples_dropdown_boxes.children[0].options if val == self.samples_dropdown_boxes.children[0].value)
            if len(self.method_name_textbox.value) > 0:
                self.sample_out_name_textbox.value = f"{selected_sample_name}_{self.method_name_textbox.value}"
            else:
                self.sample_out_name_textbox.value = selected_sample_name
    
    def upload_measurements_to_openbis(self, b):
        sample_object = self.openbis_session.get_object(self.samples_dropdown_boxes.children[0].value)
    
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
                measurements_collection = utils.create_openbis_collection(
                    code=f"MEASUREMENTS_COLLECTION_{sample_object.code}", 
                    type="COLLECTION", project=sample_project, 
                    props=collection_props
                )
            
            nanonis_importer.upload_measurements_into_openbis(
                self.session_data['url'], data_folder, 
                measurements_collection.permId, sample_object.permId, 
                self.instruments_dropdown_boxes.children[0].value
            )
            print("Upload successful!")
        else:
            print(f"Folder contains unrecognised files. Please remove files with extensions other than: {', '.join(self.measurement_file_extensions)}.")

    def close_notebook(self, b):
        display(utils.Javascript(data = 'window.location.replace("home.ipynb")'))