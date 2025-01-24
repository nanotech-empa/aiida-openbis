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
CONFIG_ELN = utils.get_aiidalab_eln_config()
OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis(CONFIG_ELN["url"], CONFIG_ELN["token"])

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
        utils.create_experiment_in_openbis(OPENBIS_SESSION, self.projects_dropdown_boxes.dropdown.value, self.new_experiment_name_textbox.value)
        self.load_dropdown_box()
        with self.add_experiment_output:
            clear_output()
  
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
            value=utils.read_file(CONFIG["default_image_filepath"]), 
            format='jpg',
            layout=ipw.Layout(border='solid 1px #cccccc', width = '220px', height = '250px')
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

class MultipleSelectorWidget(ipw.VBox):
    def __init__(self, selector_type):
        # Initialize the parent VBox
        super().__init__()
        self.selectors = []
        self.selector_type = selector_type
        
    def add_selector(self):
        if self.selector_type == "molecule":
            selector_config = {
                "dropdown": {"width": "315px"},
                "details": {"width": "415px", "height": "250px"},
                "image": {"width": "220px", "height": "250px"}
            }
            selector = ObjectSelectionWidget("Molecule", selector_config)
            selector.load_dropdown_box("MOLECULE", "molecule")
            selector.dropdown.observe(selector.load_metadata, names = "value")
            
        elif self.selector_type == "molecule_concept":
            selector_config = {
                "dropdown": {"width": "315px"},
                "details": {"width": "415px", "height": "250px"},
                "image": {"width": "220px", "height": "250px"}
            }
            selector = ObjectSelectionWidget("Molecule Concept", selector_config)
            selector.load_dropdown_box("MOLECULE_CONCEPT", "molecule concept")
            selector.dropdown.observe(selector.load_metadata, names = "value")
            
        elif self.selector_type == "product":
            selector_config = {
                "dropdown": {"width": "315px"},
                "details": {"width": "415px", "height": "250px"},
                "image": {"width": "220px", "height": "250px"}
            }
            selector = ObjectSelectionWidget("Reaction Product", selector_config)
            selector.load_dropdown_box("REACTION_PRODUCT", "reaction product")
            selector.dropdown.observe(selector.load_metadata, names = "value")
        
        elif self.selector_type == "product_concept":
            selector_config = {
                "dropdown": {"width": "315px"},
                "details": {"width": "415px", "height": "250px"},
                "image": {"width": "220px", "height": "250px"}
            }
            selector = ObjectSelectionWidget("Reaction Product Concept", selector_config)
            selector.load_dropdown_box("REACTION_PRODUCT_CONCEPT", "reaction product concept")
            selector.dropdown.observe(selector.load_metadata, names = "value")
            
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
            
            elif property["property_widget"] == "MULTIPLE_CHECKBOXES":
                checkboxes_list = []
                label_widget = ipw.Label(property["title"])
                checkboxes_list.append(label_widget)
                for i in range(property["num_elements"]):
                    checkbox_widget = utils.Checkbox(
                        disabled = property["disabled"], 
                        layout = ipw.Layout(width = property["box_layout"]["width"]),
                        style = {'description_width': property["box_layout"]["description_width"]},
                        indent = False, value = False
                    )
                    checkboxes_list.append(checkbox_widget)
                
                prop_widget = ipw.HBox(checkboxes_list)
            
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
            elif property["property_widget"] == "INTTEXT":
                prop_widget = utils.IntText(
                    description = property["title"], disabled = property["disabled"], 
                    layout = ipw.Layout(width = property["box_layout"]["width"]), 
                    placeholder = property["placeholder"], 
                    style = {'description_width': property["box_layout"]["description_width"]}
                )
            
            elif property["property_widget"] == "FLOATTEXT":
                prop_widget = utils.FloatText(
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

class ObjectMultipleSelectionWidget(ipw.HBox):
    def __init__(self, description):
        # Initialize the parent HBox
        super().__init__()
        self.description = description
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = description, disabled = False, 
                                             layout = ipw.Layout(width = '800px'), 
                                             style = {'description_width': "110px"})
        
        self.children = [self.selector]
    
    def load_selector(self, openbis_objects_types, project_id = None):
        
        if isinstance(openbis_objects_types, str):
            openbis_objects_types = [openbis_objects_types]
        
        selector_options = []
        for type in openbis_objects_types:
            if project_id:
                items = utils.get_openbis_objects(OPENBIS_SESSION, type = type, project = project_id)
            else:
                items = utils.get_openbis_objects(OPENBIS_SESSION, type = type)
            
            if type == "2D_MEASUREMENT" and self.description == "Simulations":
                items = [item for item in items if item.props["wfms_uuid"]]
            
            options = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
            selector_options += options
            
        self.selector.options = selector_options

class ObjectSelectionWidget(ipw.HBox):
    def __init__(self, description, widgets_types = {"dropdown": {"width": "982px"}}, checkboxes_descriptions = ["Name", "Registration date"]):
        # Initialize the parent HBox
        super().__init__()
        
        widgets_list = []
        self.dropdown_boxes = None
        self.details_textbox = None
        self.image_box = None
        
        self.description = description
        self.checkboxes_descriptions = checkboxes_descriptions
        
        if "dropdown" in widgets_types:
            self.dropdown = utils.Dropdown(
                description = description, 
                disabled = False, 
                layout = ipw.Layout(width = widgets_types["dropdown"]["width"]), 
                style = {'description_width': "110px"}, 
                options = [-1]
            )
        
            self.sorting_checkboxes_list = ipw.HBox(
                [
                    ipw.Label(value = "Sort by:", layout = ipw.Layout(width = "130px", display = "flex", justify_content='flex-end')),
                    utils.Checkbox(description = self.checkboxes_descriptions[0], value = False, disabled = False, layout = ipw.Layout(width = "60px"), indent = False),
                    utils.Checkbox(description = self.checkboxes_descriptions[1], value = False, disabled = False, layout = ipw.Layout(width = "200px"), indent = False)
                ]
            )
        
            self.dropdown_boxes = ipw.VBox([self.dropdown, self.sorting_checkboxes_list])
            
            widgets_list.append(self.dropdown_boxes)
        
        if "details" in widgets_types:
            # Details textbox (disabled for display purposes)
            self.details_textbox = ipw.Textarea(
                description="", 
                disabled=True, 
                layout=ipw.Layout(width=widgets_types["details"]["width"], height=widgets_types["details"]["height"])
            )
            
            widgets_list.append(self.details_textbox)
        
        if "image" in widgets_types:
            # Image box (displaying a default image)
            self.image_box = ipw.Image(
                value=utils.read_file(CONFIG["default_image_filepath"]), 
                format='jpg',
                layout=ipw.Layout(border='solid 1px #cccccc', width = widgets_types["image"]["width"], height = widgets_types["image"]["height"])
            )
            
            widgets_list.append(self.image_box)
        
        self.children = widgets_list
    
    def load_dropdown_box(self, type = None, placeholder = "", source = "openbis"):
        self.type = type
        
        if source == "openbis":
            items = utils.get_openbis_objects(
                OPENBIS_SESSION,
                type = type
            )
            options = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
            options.insert(0, (f'Select {placeholder}...', -1))
            
        elif source == "aiidalab":
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
            options = []
            for result in results:
                if result[3]:
                    name_pk_string = f"{result[3][:20]} - {result[2]} (PK: {result[0]})"
                else:
                    name_pk_string = f"{result[2]} (PK: {result[0]})"
                    
                name_pk_tuple = (name_pk_string, result[0])
                options.append(name_pk_tuple)
            
            options.insert(0, (f'Select {placeholder}...', -1))
        
        self.dropdown.options = options
        self.dropdown.value = -1
            
        utils.sort_dropdown(
            self.sorting_checkboxes_list,
            self.dropdown,
            self.checkboxes_descriptions,
            [True, False]
        )
        
        for checkbox in self.sorting_checkboxes_list.children[1:]:
            checkbox.observe(
                lambda change: utils.sort_dropdown(
                    self.sorting_checkboxes_list, 
                    self.dropdown,
                    self.checkboxes_descriptions,
                    [True, False]
                ), 
                names='value'
            )
    
    # Function to handle changes in the molecules dropdown
    def load_metadata(self, change):
        if self.dropdown.value == -1:
            if self.details_textbox:
                self.details_textbox.value = ''
            if self.image_box:
                self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])
            return
        
        # Get metadata
        object = OPENBIS_SESSION.get_object(self.dropdown.value)
        
        if self.image_box:
            if self.type in ["CRYSTAL", "MOLECULE"]:
                parents = object.get_parents()
                for parent in parents:
                    if parent.type == f"{self.type}_CONCEPT":
                        concept_object = parent
                        
                # Get image
                object_dataset = concept_object.get_datasets(type="ELN_PREVIEW")[0]
            
            else:
                object_dataset = object.get_datasets(type="ELN_PREVIEW")[0]

            if object_dataset:
                object_dataset.download(destination="images")
                object_image_filepath = object_dataset.file_list[0]
                self.image_box.value = utils.read_file(f"images/{object_dataset.permId}/{object_image_filepath}")
                shutil.rmtree(f"images/{object_dataset.permId}/")
            else:
                self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])

        if self.details_textbox:
            metadata_string = utils.get_metadata_string(OPENBIS_SESSION, object, "", CONFIG)
            self.details_textbox.value = metadata_string

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
        items_names_permids = [(f"{item.code} ({item.identifier})", item.permId) for item in items]
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

class SamplePreparationAccordionWidget(ipw.Accordion):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        self.tasks_properties_widgets = []

class SamplePreparationMultipleSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        self.selector = ipw.SelectMultiple()
        self.sample_preparation_options = [object_key for object_key, object_info in CONFIG["objects"].items() if object_info["object_type"] == "sample_preparation"]
        self.selector.description = "Processes"
        self.selector.disabled = False
        self.selector.style = {'description_width': "65px"}
        self.selector.options = self.sample_preparation_options
        selector_height = len(self.sample_preparation_options) * 18
        selector_height = f"{selector_height}px"
        self.selector.layout = ipw.Layout(width = '220px', height = selector_height)
        self.children = [self.selector]

