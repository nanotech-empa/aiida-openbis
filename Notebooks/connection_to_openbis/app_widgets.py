import json
import os
import ipywidgets as ipw
from IPython.display import display, clear_output
import sys
sys.path.append('/home/jovyan/aiida-openbis/Notebooks/importer')
import nanonis_importer
import shutil
import utils

class AppWidgets():
    def __init__(self, method_type, config_filename):
        self.method_type = method_type
        
        # Get configuration data
        self.config = utils.read_json(config_filename)
        self.config_eln = utils.read_json("eln_config.json")
        self.samples_collection_openbis_path = self.config["samples_collection_openbis_path"]
        self.measurement_file_extensions = self.config["measurement_file_extensions"]
        
        # Connect to openBIS
        self.openbis_session, self.session_data = utils.connect_openbis(self.config_eln["url"], self.config_eln["token"])
        
        # Necessary for refreshing the widgets needed to create new experiments in sample preparation page
        self.select_experiment_output = ipw.Output()
        
        # Home page configuration
        self.openbis_connection_status_htmlbox = utils.HTMLbox(value = '')
        if self.openbis_session:
            self.openbis_connection_status_htmlbox.value = self.config["home_page"]["enable_status"]
            self.open_notebooks_html_enable_code = ''.join(self.config["home_page"]["enable_links"])
            self.open_notebooks_htmlbox = utils.HTMLbox(value = self.open_notebooks_html_enable_code)
        else:
            self.openbis_connection_status_htmlbox.value = self.config["home_page"]["disable_status"]
            self.open_notebooks_html_disable_code = ''.join(self.config["home_page"]["disable_links"])
            self.open_notebooks_htmlbox = utils.HTMLbox(value = self.open_notebooks_html_disable_code)
         
        # Get slabs types
        self.slabs_types = [object_info["openbis_object_type"] for _, object_info in self.config["objects"].items() if object_info["object_type"] == "slab"]
        
        # Get sample preparation types
        self.sample_preparation_sample_types = [object_info["openbis_object_type"] for _, object_info in self.config["objects"].items() if object_info["object_type"] == "sample_preparation"]
        self.sample_preparation_options = [object_key for object_key, object_info in self.config["objects"].items() if object_info["object_type"] == "sample_preparation"]
        
        # Get widgets for new experiments
        self.create_new_experiment_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Add', icon = 'plus', layout = ipw.Layout(width = '50px', height = '25px'))
        self.save_new_experiment_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Save', icon = 'save', layout = ipw.Layout(width = '50px', height = '35px', margin = '0 0 0 90px'))
        self.cancel_new_experiment_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Cancel', icon = 'times', layout = ipw.Layout(width = '50px', height = '35px', margin = '0 0 0 5px'))
        self.new_experiment_name_textbox = utils.Text(description = "Name", disabled = False, layout = ipw.Layout(width = '400px'), placeholder = f"Write experiment name here...", style = {'description_width': "110px"})
        
        # Get widgets for slab dropdown
        slabs_options = [object_key for object_key, object_info in self.config["objects"].items() if object_info["object_type"] == "slab"]
        slabs_options.insert(0, "No material")
        self.material_selection_radio_button = utils.Radiobuttons(description = '', options = slabs_options, disabled = False, layout = ipw.Layout(width = '300px'), style = {'description_width': "100px"})
        self.materials_dropdown = utils.Dropdown(description='', disabled=False, layout = ipw.Layout(width = '350px'))
        self.material_details_textbox = utils.Textarea(description = "", disabled = True, layout = ipw.Layout(width = '425px', height = '200px'))
        self.material_image_box = utils.Image(value = utils.read_file(self.config["default_image_filepath"]), format = 'jpg', width = '200px', height = '300px', layout=ipw.Layout(border='solid 1px #cccccc'))
        self.material_metadata_boxes = ipw.HBox([self.material_details_textbox, self.material_image_box])
        self.material_sorting_checkboxes = utils.SortingCheckboxes("50px", "60px", "200px")
        
        # Get widgets for sample dropdown
        self.samples_dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget('Sample', ipw.Layout(width = '385px'), {'description_width': "110px"}, [-1])
        self.sample_details_textbox = utils.Textarea(description = "", disabled = True, layout = ipw.Layout(width = '589px', height = '300px'))
        self.sample_metadata_boxes = ipw.HBox([self.samples_dropdown_boxes, self.sample_details_textbox])

        # Get widgets for instrument dropdown
        self.instruments_dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget('Instrument', ipw.Layout(width = '982px'), {'description_width': "110px"}, [-1])
        
        # Get widgets for project dropdown
        self.projects_dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget('Project', ipw.Layout(width = '982px'), {'description_width': "110px"}, [-1])
        
        # Get widgets for experiment dropdown
        self.experiments_dropdown = utils.Dropdown(description='Experiment', disabled=False, layout = ipw.Layout(width = '993px'), style = {'description_width': "110px"}, options = [-1])
        self.experiment_sorting_checkboxes = utils.SortingCheckboxes("130px", "60px", "200px")
        self.experiments_dropdown_details = ipw.HBox([self.experiments_dropdown, self.create_new_experiment_button])
        self.experiments_dropdown_boxes = ipw.VBox([self.experiments_dropdown_details, self.experiment_sorting_checkboxes])
        
        # Get widgets for molecule dropdown
        self.molecules_dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget('Substance', ipw.Layout(width = '335px'), {'description_width': "110px"}, [-1])
        self.molecule_details_textbox = utils.Textarea(description = "", disabled = True, layout = ipw.Layout(width = '415px', height = '250px'))
        self.molecule_image_box = utils.Image(value = utils.read_file(self.config["default_image_filepath"]), format = 'jpg', width = '220px', height = '250px', layout=ipw.Layout(border='solid 1px #cccccc'))
        self.molecule_metadata_boxes = ipw.HBox([self.molecules_dropdown_boxes, self.molecule_details_textbox, self.molecule_image_box])

        # Get widgets for properties
        self.property_widgets = {}
        self.property_widgets["$name"] = utils.Text(description = "Name", disabled = False, layout = ipw.Layout(width = '400px'), placeholder = f"Write task name here...", style = {'description_width': "150px"})
        self.property_widgets["description"] = utils.Textarea(description = "Description", disabled = False, layout = ipw.Layout(width = '993px', height = '100px'), placeholder = f"Write task description here...", style = {'description_width': "150px"})
        self.property_widgets["comments"] = utils.Textarea(description = "Comments", disabled = False, layout = ipw.Layout(width = '993px', height = '200px'), placeholder = "Write comments here...", style = {'description_width': "150px"})
        
        for prop in self.config["properties"].keys():
            if self.config["properties"][prop]["property_type"] == "quantity_value":
                self.property_widgets[prop] = utils.FloatTextwithDropdownWidget(
                    self.config["properties"][prop]["title"], ipw.Layout(width = self.config["properties"][prop]["box_layout"]["width"]), 
                    self.config["properties"][prop]["default_value"],
                    {'description_width': self.config["properties"][prop]["box_layout"]["description_width"]}, 
                    ipw.Layout(width = self.config["properties"][prop]["dropdown_layout"]["width"]), self.config["properties"][prop]["units"], 
                    self.config["properties"][prop]["default_unit"]
                )
        
        self.property_widgets["pid_controller"] = utils.Text(
            description = self.config["properties"]["pid_controller"]["title"], 
            disabled = False, 
            layout = ipw.Layout(
                width = self.config["properties"]["pid_controller"]["box_layout"]["width"]), 
                placeholder = self.config["properties"]["pid_controller"]["placeholder"],
                style = {'description_width': self.config["properties"]["pid_controller"]["box_layout"]["description_width"]}
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
            self.config["properties"]["evaporator_slot"]["default_value"], 
            self.config["properties"]["evaporator_slot"]["title"], 
            self.config["properties"]["evaporator_slot"]["values_list"], 
            ipw.Layout(width = self.config["properties"]["evaporator_slot"]["slider_layout"]["width"]), 
            {'description_width': self.config["properties"]["evaporator_slot"]["slider_layout"]["description_width"]}, 
            self.config["properties"]["evaporator_slot"]["placeholder"], 
            ipw.Layout(width = self.config["properties"]["evaporator_slot"]["box_layout"]["width"])
        )
        
        # Select the properties for the specific object one is interested
        self.object_widgets = {}
        for key, object in self.config["objects"].items():
            if object["object_type"] in ["sample_preparation", "other"]:
                items = [self.property_widgets[prop] for prop in object["properties"]]
                self.object_widgets[key] = ipw.VBox(items)

        # Get widget for sample name
        self.sample_out_name_textbox = utils.Text(description = "Name", disabled = False, layout = ipw.Layout(width = '400px'), placeholder = f"Write sample name here...", style = {'description_width': "110px"})
        
        # Get widgets for saving/closing buttons
        self.create_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Save', icon = 'save', layout = ipw.Layout(width = '100px', height = '50px'))
        self.quit_button = utils.Button(description = '', disabled = False, button_style = '', tooltip = 'Main menu', icon = 'home', layout = ipw.Layout(width = '100px', height = '50px'))
        self.bottom_buttons_hbox = ipw.HBox([self.create_button, self.quit_button])
        
        self.sample_prep_selector = utils.SelectMultiple(description = 'Processes', disabled = False, 
                                                         layout = ipw.Layout(width = '300px', height = '110px'), 
                                                         style = {'description_width': "65px"},
                                                         options = self.sample_preparation_options)
        self.sample_prep_accordion = ipw.Accordion(children = [])
        
        # Widget to select folder with measurement files
        self.folder_selector = utils.FileChooser(path = '.', select_default=True, use_dir_icons=True, show_only_dirs = True)
        
        # Measurements selector widget
        self.measurements_selector = utils.SelectMultiple(description = 'Measurements', disabled = False, layout = ipw.Layout(width = '800px'), style = {'description_width': "110px"})
        
        # Results selector widget
        self.results_selector = utils.SelectMultiple(description = 'Results', disabled = False, layout = ipw.Layout(width = '800px'), style = {'description_width': "110px"})
        
        # Drafts selector widget
        self.drafts_selector = utils.SelectMultiple(description = 'Drafts', disabled = False, layout = ipw.Layout(width = '800px'), style = {'description_width': "110px"})
        
        # Authors selector widget
        self.authors_selector = utils.SelectMultiple(description = 'Authors', disabled = False, layout = ipw.Layout(width = '800px'), style = {'description_width': "110px"})
        self.load_authors()
        
        # Grants selector widget
        self.grants_selector = utils.SelectMultiple(description = 'Grants', disabled = False, layout = ipw.Layout(width = '800px'), style = {'description_width': "110px"})
        self.load_grants()
        
        # Draft dropdown box widget
        self.draft_type_dropdown = utils.Dropdown(description = 'Draft type', disabled = False, layout = ipw.Layout(width = '300px'), style = {'description_width': "150px"},
                                                    options = [("Select draft type...", -1), ("Preprint", "PREPRINT"), ("Postprint", "POSTPRINT")], value = -1)
        
        # Publication properties
        self.publication_abstract_textarea = utils.Textarea(description = "Abstract", disabled = False, layout = ipw.Layout(width = '993px', height = '100px'), 
                                                            placeholder = f"Write publication abstract here...", style = {'description_width': "150px"})
        self.publication_doi_textbox = utils.Text(description = "DOI", disabled = False, layout = ipw.Layout(width = '993px'), 
                                                    placeholder = f"Write publication DOI here...", style = {'description_width': "150px"})
        self.publication_year_intbox = utils.IntText(description = "Year", disabled = False, layout = ipw.Layout(width = '250px'), 
                                                        placeholder = f"Write publication year here...", style = {'description_width': "150px"})
        self.publication_url_textbox = utils.Text(description = "URL", disabled = False, layout = ipw.Layout(width = '993px'), 
                                                    placeholder = f"Write publication URL here...", style = {'description_width': "150px"})
        self.publication_dataset_url_textbox = utils.Text(description = "URL", disabled = False, layout = ipw.Layout(width = '993px'), 
                                                            placeholder = f"Write publication dataset URL here...", style = {'description_width': "150px"})
        
        self.support_files_uploader = ipw.FileUpload(multiple = True)
        
        # Assign functions to widgets
        self.create_new_experiment_button.on_click(self.create_new_experiment_button_on_click)
        self.cancel_new_experiment_button.on_click(self.cancel_new_experiment_button_on_click)
        self.save_new_experiment_button.on_click(self.save_new_experiment_button_on_click)
        self.sample_prep_selector.observe(self.sample_preparation_widgets, names='value')
        
        # Increase buttons icons' size
        self.increase_buttons_size = utils.HTML(data = ''.join(self.config["save_home_buttons_settings"]))
        display(self.increase_buttons_size)

    def sample_preparation_widgets(self, change):
    
        accordion_options = []
        
        for i, task in enumerate(self.sample_prep_selector.value):
            
            task_widgets = []
            
            if self.config["objects"][task]["uses_molecule"]:
                task_widgets.append(self.molecule_metadata_boxes)
                
            task_widgets.extend([self.object_widgets[task]])
            
            task_widgets.append(self.support_files_uploader)
            
            self.sample_prep_accordion.set_title(i, task)
            accordion_options.append(ipw.VBox(task_widgets))
        
        self.sample_prep_accordion.children = accordion_options
    
    def create_publication_openbis(self, b):
        publication_props = {"$name": self.property_widgets["$name"].value, "comments": self.property_widgets["comments"].value,
                             "abstract": self.publication_abstract_textarea.value, "doi": self.publication_doi_textbox.value,
                             "url": self.publication_url_textbox.value, "year": self.publication_year_intbox.value,
                             "dataset_url": self.publication_dataset_url_textbox.value}
        publication_parents = []
        publication_parents.extend(list(self.authors_selector.value))
        publication_parents.extend(list(self.drafts_selector.value))
        publication_parents.extend(list(self.grants_selector.value))
        publication_object = utils.create_openbis_object(self.openbis_session, type="PUBLICATION_CUSTOM", 
                                                         collection="/PUBLICATIONS/PUBLIC_REPOSITORIES/PUBLICATIONS_COLLECTION", 
                                                         props = publication_props, parents = publication_parents)
        self.upload_datasets(publication_object)
        print("Upload successful!")
    
    def create_draft_openbis(self, b):
        draft_props = {"$name": self.property_widgets["$name"].value, "comments": self.property_widgets["comments"].value}
        project_permid = self.projects_dropdown_boxes.children[0].value
        project_drafts_collection = self.openbis_session.get_collections(project = project_permid, code = "DRAFTS_COLLECTION")
        drafts_parents = list(self.results_selector.value)
        
        if len(project_drafts_collection) == 0:
            drafts_collection = utils.create_openbis_collection(self.openbis_session, project = self.projects_dropdown_boxes.children[0].value, 
                                                                 code = "DRAFTS_COLLECTION", type = "COLLECTION", props = {"$name": "Drafts"})
        else:
            drafts_collection = project_drafts_collection[0]

        draft_object = utils.create_openbis_object(self.openbis_session, type="DRAFT", collection=drafts_collection, 
                                                   props = draft_props, parents = drafts_parents)
        self.upload_datasets(draft_object)
        
        print("Upload successful!")
    
    def create_results_openbis(self, b):
        results_props = {"$name": self.property_widgets["$name"].value, "description": self.property_widgets["description"].value,
                         "comments": self.property_widgets["comments"].value}
        project_permid = self.projects_dropdown_boxes.children[0].value
        project_results_collection = self.openbis_session.get_collections(project = project_permid, code = "RESULTS_COLLECTION")
        results_parents = list(self.measurements_selector.value)
        
        if len(project_results_collection) == 0:
            results_collection = utils.create_openbis_collection(self.openbis_session, project = self.projects_dropdown_boxes.children[0].value, 
                                                                 code = "RESULTS_COLLECTION", type = "COLLECTION", props = {"$name": "Results"})
        else:
            results_collection = project_results_collection[0]

        utils.create_openbis_object(self.openbis_session, type="RESULTS", collection=results_collection, 
                                    props = results_props, parents = results_parents)
        
        print("Upload successful!")
    
    def load_measurements(self, change):
        if self.projects_dropdown_boxes.children[0].value != -1:
            one_d_measurements = self.openbis_session.get_objects(type = "1D_MEASUREMENT", project = self.projects_dropdown_boxes.children[0].value)
            two_d_measurements = self.openbis_session.get_objects(type = "2D_MEASUREMENT", project = self.projects_dropdown_boxes.children[0].value)
            measurements_list = []
            measurements_list.extend([(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in one_d_measurements])
            measurements_list.extend([(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in two_d_measurements])
            self.measurements_selector.options = measurements_list
    
    def load_results(self, change):
        if self.projects_dropdown_boxes.children[0].value != -1:
            results = self.openbis_session.get_objects(type = "RESULTS", project = self.projects_dropdown_boxes.children[0].value)
            results_list = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in results]
            self.results_selector.options = results_list
    
    def load_drafts(self, change):
        if self.projects_dropdown_boxes.children[0].value != -1:
            drafts = self.openbis_session.get_objects(type = "DRAFT", project = self.projects_dropdown_boxes.children[0].value)
            drafts_list = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in drafts]
            self.drafts_selector.options = drafts_list
    
    def load_authors(self):
        authors = self.openbis_session.get_objects(type = "AUTHOR")
        authors_list = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in authors]
        self.authors_selector.options = authors_list
    
    def load_grants(self):
        grants = self.openbis_session.get_objects(type = "GRANT")
        grants_list = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in grants]
        self.grants_selector.options = grants_list
    
    def create_new_experiment_button_on_click(self, b):
        with self.select_experiment_output:
            clear_output()
            display_list = [self.experiments_dropdown_boxes, self.new_experiment_name_textbox, self.projects_dropdown_boxes,
                            ipw.HBox([self.save_new_experiment_button, self.cancel_new_experiment_button]), 
                            self.sample_metadata_boxes, self.instruments_dropdown_boxes]
            if self.config["objects"][self.method_type]["uses_molecule"]:
                display_list.append(self.molecule_metadata_boxes)
            display(ipw.VBox(display_list))
    
    def cancel_new_experiment_button_on_click(self, b):
        with self.select_experiment_output:
            clear_output()
            display_list = [self.experiments_dropdown_boxes, self.sample_metadata_boxes, self.instruments_dropdown_boxes]
            # This flag is needed to check whether molecule interface should be displayed (dropdown and metadata details boxes)
            if self.config["objects"][self.method_type]["uses_molecule"]: 
                display_list.append(self.molecule_metadata_boxes)
            display(ipw.VBox(display_list))
    
    def save_new_experiment_button_on_click(self, b):
        utils.create_experiment_in_openbis(self.openbis_session, self.projects_dropdown_boxes.children[0].value, self.new_experiment_name_textbox.value)
        with self.select_experiment_output:
            clear_output()
            utils.load_openbis_elements_list(self.openbis_session, "EXPERIMENT", self.experiments_dropdown, self.experiment_sorting_checkboxes, "experiment")
            display_list = [self.experiments_dropdown_boxes, self.sample_metadata_boxes, self.instruments_dropdown_boxes]
            # This flag is needed to check whether molecule interface should be displayed (dropdown and metadata details boxes)
            if self.config["objects"][self.method_type]["uses_molecule"]:
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
            utils.create_openbis_object(self.openbis_session, type="SAMPLE", collection=self.samples_collection_openbis_path, props=sample_props, parents=sample_parents)
            print("Upload successful!")
    
    # Function to handle changes in the materials dropdown
    def load_material_metadata(self, change):
        if self.materials_dropdown.value == -1:
            self.material_details_textbox.value = ''
            self.material_image_box.value = utils.read_file(self.config["default_image_filepath"])
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
            # Erase file after downloading it
            shutil.rmtree(f"images/{material_dataset.permId}")
        else:
            self.material_image_box.value = utils.read_file(self.config["default_image_filepath"])

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
        self.sample_out_name_textbox.value = utils.convert_datetime_to_string(utils.get_current_datetime()) + f"_{material_metadata['$name']}"
    
    def select_material_radio_change(self, change):
        self.material_details_textbox.value = ''
        clear_output()
        display(self.material_selection_radio_button)
        
        material_types = {}
        for object_key, object_info in self.config["objects"].items():
            if object_info["object_type"] == "slab":
                material_types[object_key] = (object_info["openbis_object_type"], object_info["placeholder"])

        material_type = self.material_selection_radio_button.value
        if material_type == "No material":
            return

        material_class, placeholder = material_types.get(material_type)
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
        for filename in self.support_files_uploader.value:
            file_info = self.support_files_uploader.value[filename]
            utils.save_file(file_info['content'], filename)
            self.openbis_session.new_dataset(type = 'RAW_DATA', sample = method_object, files = [filename]).save()
            os.remove(filename)
    
    def create_miscellaneous_action(self, b):
        if self.experiments_dropdown.value == -1:
            print("Select an experiment.")
            return

        object_properties = {"$name": self.property_widgets["$name"].value, "comments": self.property_widgets["comments"].value}
        
        method_object = utils.create_openbis_object(self.openbis_session, type = self.config["objects"][self.method_type]["openbis_object_type"],
                                                    collection = self.experiments_dropdown.value, props = object_properties)
        self.upload_datasets(method_object)
        print("Upload successful!")
    
    def create_calibration_optimisation_action(self, b):
        if self.experiments_dropdown.value == -1:
            print("Select an experiment.")
            return
        
        if self.instruments_dropdown_boxes.children[0].value == -1:
            print("Select an instrument.")
            return
        
        if self.molecules_dropdown_boxes.children[0].value == -1 and self.config["objects"][self.method_type]["uses_molecule"]:
            print("Select a substance.")
            return

        # Prepare sample parents based on method type
        object_parents = [self.instruments_dropdown_boxes.children[0].value]
        if self.config["objects"][self.method_type]["uses_molecule"]:
            object_parents.append(self.molecules_dropdown_boxes.children[0].value)

        object_properties = {}
        for prop in self.config["objects"][self.method_type]["properties"]:
            if self.config["properties"][prop]["property_type"] == "string":
                object_properties[prop] = self.property_widgets[prop].value
            elif prop == "evaporator_slot":
                object_properties[prop] = json.dumps({"evaporator_number": self.property_widgets[prop].children[0].value, "details": self.property_widgets[prop].children[1].value})
            elif self.config["properties"][prop]["property_type"] == "quantity_value":
                object_properties[prop] = json.dumps({"has_value": self.property_widgets[prop].children[0].value, "has_unit": self.property_widgets[prop].children[1].value})
        
        method_object = utils.create_openbis_object(self.openbis_session, type = self.config["objects"][self.method_type]["openbis_object_type"], 
                                                    collection = self.experiments_dropdown.value, props = object_properties, parents = object_parents)
        self.upload_datasets(method_object)
        print("Upload successful!")
    
    def create_sample_preparation_action(self, b):
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

        # Prepare sample parents based on method type
        methods_objects = []
        for i, task in enumerate(self.sample_prep_selector.value):
            
            method_parents = [self.samples_dropdown_boxes.children[0].value, self.instruments_dropdown_boxes.children[0].value]
            
            # ---------
            # TODO: I should create a widget per property per object. Otherwise, the duration in sputtering and in annealing is shared (is the same object).
            if self.molecules_dropdown_boxes.children[0].value == -1 and self.config["objects"][task]["uses_molecule"]:
                print("Select a substance.")
                return
            # ---------
        
            if self.config["objects"][task]["uses_molecule"]:
                method_parents.append(self.molecules_dropdown_boxes.children[0].value)

            object_properties = {}
            
            for prop in self.config["objects"][task]["properties"]:
                if prop == "evaporator_slot":
                    object_properties[prop] = json.dumps({"evaporator_number": self.property_widgets[prop].children[0].value, "details": self.property_widgets[prop].children[1].value})
                elif self.config["properties"][prop]["property_type"] == "string":
                    object_properties[prop] = self.property_widgets[prop].value
                elif self.config["properties"][prop]["property_type"] == "quantity_value":
                    object_properties[prop] = json.dumps({"has_value": self.property_widgets[prop].children[0].value, "has_unit": self.property_widgets[prop].children[1].value})

            method_object = utils.create_openbis_object(self.openbis_session, type = self.config["objects"][task]["openbis_object_type"], 
                                                        collection = self.experiments_dropdown.value, props = object_properties, parents = method_parents)
            self.upload_datasets(method_object)
            
            methods_objects.append(method_object)

        # Turn off sample visibility
        parent_sample = self.openbis_session.get_object(self.samples_dropdown_boxes.children[0].value)
        parent_sample.props["exists"] = False
        parent_sample.save()

        sample_props = {"$name": self.sample_out_name_textbox.value, "exists": True}
        utils.create_openbis_object(self.openbis_session, type = "SAMPLE", collection = self.samples_collection_openbis_path, 
                                    props = sample_props, parents = methods_objects)
        print("Upload successful!")
    
    # Function to handle changes in the substances dropdown
    def load_molecule_metadata(self, change):
        if self.molecules_dropdown_boxes.children[0].value == -1:
            self.molecule_details_textbox.value = ''
            self.molecule_image_box.value = utils.read_file(self.config["default_image_filepath"])
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
            self.molecule_image_box.value = utils.read_file(self.config["default_image_filepath"])

        material_metadata = material_object.props.all()
        material_metadata_string = ""
        for prop_key in property_list:
            prop_title = self.config["properties"][prop_key]["title"]
            prop_value = material_metadata.get(prop_key)
            if self.config["properties"][prop_key]["property_type"] == "quantity_value" and prop_value is not None:
                prop_value = json.loads(prop_value)
                material_metadata_string += f"{prop_title}: {prop_value['value']} {prop_value['unit']}\n"
            elif prop_key == "has_molecule":
                material_metadata_string += f"{prop_title}:\n"
                for mol_prop_key in molecule_property_list:
                    mol_prop_title = self.config["properties"][mol_prop_key]["title"]
                    mol_prop_value = molecule_metadata.get(mol_prop_key)
                    material_metadata_string += f"- {mol_prop_title}: {mol_prop_value}\n"
            else:
                material_metadata_string += f"{prop_title}: {prop_value}\n"

        self.molecule_details_textbox.value = material_metadata_string

    def load_sample_metadata(self, change):
        if self.samples_dropdown_boxes.children[0].value == -1:
            self.experiments_dropdown.value = -1
            self.instruments_dropdown_boxes.children[0].value = -1
            self.sample_details_textbox.value = ''
            
            if self.method_type in self.config["objects"].keys():
                if self.config["objects"][self.method_type]["object_type"] in ["sample_preparation", "comments_notes"]:
                    self.sample_out_name_textbox.value = ''
            return
        
        sample_object = self.openbis_session.get_object(self.samples_dropdown_boxes.children[0].value)
        sample_parents_metadata = utils.get_openbis_parents_recursive(self.openbis_session, sample_object, [])
        
        last_sample_preparation = None
        sample_strings = {"sample_preparations": [], "materials": []}
        number_parents = len(sample_parents_metadata)
        parent_idx = 0
        while parent_idx < number_parents:
            parent_metadata = sample_parents_metadata[parent_idx]
            if parent_metadata[0] in ["DEPOSITION"]:
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
            
            if parent_metadata[0] in self.sample_preparation_sample_types:
                sample_strings["sample_preparations"].append(sample_metadata_string)
                # Get the last sample preparation method performed on the sample in order to search the correct experiment where the sample is being used.
                if last_sample_preparation is None:
                    last_sample_preparation = parent_metadata[1]
            elif parent_metadata[0] in self.slabs_types:
                sample_strings["materials"].append(sample_metadata_string)
                
            parent_idx += 1
        
        sample_metadata_string = (f"PermId: {sample_object.attrs.permId}\nMaterial:\n" +
                              "\n".join(sample_strings["materials"]) + 
                              "\nProcesses:\n" + 
                              "\n".join(sample_strings["sample_preparations"]))

        if last_sample_preparation:
            last_sample_preparation_object = self.openbis_session.get_object(last_sample_preparation)
            
            if self.method_type in self.config["objects"].keys():
                if self.config["objects"][self.method_type]["object_type"] in ["sample_preparation", "comments_notes"]:
                # if self.config["objects"][self.method_type]["openbis_object_type"] in self.sample_preparation_sample_types:
                    # Automatically select the experiment where the last sample preparation task was saved
                    last_sample_preparation_experiment = self.openbis_session.get_experiment(last_sample_preparation_object.attrs.experiment)
                    # Notify the user in case the experiment changed
                    if self.experiments_dropdown.value != last_sample_preparation_experiment.permId and self.experiments_dropdown.value != -1:
                        dropdown_options_dict = {key: val for val, key in self.experiments_dropdown.options}
                        previous_experiment_name = dropdown_options_dict.get(self.experiments_dropdown.value)
                        new_experiment_name = dropdown_options_dict.get(last_sample_preparation_experiment.permId)
                        display(utils.Javascript(data = f"alert('{f'Experiment was changed from {previous_experiment_name} to {new_experiment_name}!'}')"))
                    self.experiments_dropdown.value =  last_sample_preparation_experiment.permId
            
            # Automatically select the instrument used in the last sample preparation task
            for parent in last_sample_preparation_object.parents:
                parent_object = self.openbis_session.get_object(parent)
                if parent_object.type == "INSTRUMENT":
                    self.instruments_dropdown_boxes.children[0].value = parent_object.permId
        
        self.sample_details_textbox.value = sample_metadata_string
        if self.method_type in self.config["objects"].keys():
            if self.config["objects"][self.method_type]["object_type"] in ["sample_preparation", "comments_notes"]:
            # if self.config["objects"][self.method_type]["openbis_object_type"] in self.sample_preparation_sample_types:
                sample_name = sample_object.props['$name']
                self.sample_out_name_textbox.value = f"{sample_name}_{self.property_widgets['$name'].value}" if self.property_widgets["$name"].value else sample_name
    
    def load_dropdown_lists(self):
        # Populate dropdown lists
        utils.load_openbis_elements_list(self.openbis_session, "SAMPLE", self.samples_dropdown_boxes.children[0], self.samples_dropdown_boxes.children[1], "sample")
        utils.load_openbis_elements_list(self.openbis_session, "INSTRUMENT", self.instruments_dropdown_boxes.children[0], self.instruments_dropdown_boxes.children[1], "instrument")
        utils.load_openbis_elements_list(self.openbis_session, "EXPERIMENT", self.experiments_dropdown, self.experiment_sorting_checkboxes, "experiment")
        utils.load_openbis_elements_list(self.openbis_session, "PROJECT", self.projects_dropdown_boxes.children[0], self.projects_dropdown_boxes.children[1], "project")
        utils.load_openbis_elements_list(self.openbis_session, "SUBSTANCE", self.molecules_dropdown_boxes.children[0], self.molecules_dropdown_boxes.children[1], "substance")

    def update_text(self, change):
        if self.samples_dropdown_boxes.children[0].value != -1:
            selected_sample_name = next(label for label, val in self.samples_dropdown_boxes.children[0].options if val == self.samples_dropdown_boxes.children[0].value)
            if len(self.property_widgets["$name"].value) > 0:
                self.sample_out_name_textbox.value = f"{selected_sample_name}_{self.property_widgets['$name'].value}"
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
            if self.openbis_session.get_object(parent_id).type in self.sample_preparation_sample_types), 
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
                    self.openbis_session,
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