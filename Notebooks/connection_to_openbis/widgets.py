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
import datetime
import copy
from aiida import orm

DATA_MODEL = utils.read_yaml("/home/jovyan/aiida-openbis/Notebooks/Metadata_Schemas_LinkML/materialMLinfo.yaml")
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
        items_names_permids = [(f"{item.props['$name']} (Project: {item.project.code})", item.permId) for item in items]
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
    
    def load_dropdown_box(self, type):
        self.schema_object_type = type
        openbis_type = DATA_MODEL["classes"][type]["annotations"]["openbis_label"].replace(" ", "_").upper()
        placeholder = DATA_MODEL["classes"][type]["title"]
        
        items = utils.get_openbis_objects(
            OPENBIS_SESSION,
            type = openbis_type
        )
        items_names_permids = [(f"{item.props['$name']}", item.permId) for item in items]
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

    # Function to handle changes in the materials objects dropdown
    def load_metadata(self, change):
        if self.dropdown.value == -1:
            self.details_textbox.value = ''
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])
            return
        
        # Get material object information and dataset
        object = OPENBIS_SESSION.get_object(self.dropdown.value)
        object_dataset = object.get_datasets()[0]
        
        # Get the object image preview
        if object_dataset:
            object_dataset.download(destination="images")
            object_image_filepath = object_dataset.file_list[0]
            self.image_box.value = utils.read_file(f"images/{object_dataset.permId}/{object_image_filepath}")
            # Erase file after downloading it
            shutil.rmtree(f"images/{object_dataset.permId}")
        else:
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])

        # Make a string with the property values of the object
        metadata_string = utils.get_metadata_string(OPENBIS_SESSION, object, self.schema_object_type, "", DATA_MODEL)
        self.details_textbox.value = metadata_string
    
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
            selector.load_dropdown_box()
            selector.dropdown.observe(selector.load_metadata, names = "value")
            
        elif self.selector_type == "molecule_concept":
            selector_config = {
                "dropdown": {"width": "315px"},
                "details": {"width": "415px", "height": "250px"},
                "image": {"width": "220px", "height": "250px"}
            }
            selector = ObjectSelectionWidget("MoleculeConcept", selector_config)
            selector.load_dropdown_box()
            selector.dropdown.observe(selector.load_metadata, names = "value")
            
        elif self.selector_type == "product":
            selector_config = {
                "dropdown": {"width": "315px"},
                "details": {"width": "415px", "height": "250px"},
                "image": {"width": "220px", "height": "250px"}
            }
            selector = ObjectSelectionWidget("ReactionProduct", selector_config)
            selector.load_dropdown_box()
            selector.dropdown.observe(selector.load_metadata, names = "value")
        
        elif self.selector_type == "product_concept":
            selector_config = {
                "dropdown": {"width": "315px"},
                "details": {"width": "415px", "height": "250px"},
                "image": {"width": "220px", "height": "250px"}
            }
            selector = ObjectSelectionWidget("ReactionProductConcept", selector_config)
            selector.load_dropdown_box()
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
        self.properties_widgets_dict = {}
        self.properties_widgets_detailed_dict = {}
        self.task = task
    
    def get_properties_widgets(self):
        self.all_schema_classes = DATA_MODEL["classes"]
        self.all_schema_slots = DATA_MODEL["slots"]
        object = self.all_schema_classes[self.task]
        object = self.get_properties_recursive(object)
        for property in object["slots"]:
            if self.all_schema_slots[property]["annotations"]["openbis_type"] != "Not used":
                property_widget, property_dict = self.get_property_widget(property)
                if property_widget:
                    if property == "name":
                        new_property = "$name"
                    else:
                        new_property = property
                        
                    self.properties_widgets_dict[new_property] = property_widget
                    self.properties_widgets_detailed_dict[new_property] = property_dict[property]
                
        self.children = list(self.properties_widgets_dict.values())
    
    def get_property_widget(self, property, is_group = False):
        property_widget = None
        property_dict = None
        property_settings = self.all_schema_slots[property]
        property_multivalued = property_settings["multivalued"]
        property_description = property_settings["description"]
        property_range = property_settings["range"]
        property_openbis_type = property_settings["annotations"]["openbis_type"]
        
        # Multivalued properties are a special case of properties
        if property_multivalued:
            return None, None
        
        if property_openbis_type == "VARCHAR":
            if is_group:
                label_widget = ipw.HTML(value = property_description)
            else:
                label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                
            text_widget = utils.Text(
                layout = ipw.Layout(width = "200px"), 
                placeholder = "",
            )
            property_widget = ipw.VBox([label_widget, text_widget])
            property_dict = {property: text_widget}
        
        elif property_openbis_type == "MULTILINE_VARCHAR":
            if is_group:
                label_widget = ipw.HTML(value = property_description)
            else:
                label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                
            textarea_widget = utils.Textarea(
                layout = ipw.Layout(width = "200px", height = "100px"), 
                placeholder = "",
            )
            property_widget = ipw.VBox([label_widget, textarea_widget])
            property_dict = {property: textarea_widget}
        
        elif property_openbis_type == "BOOLEAN":
            if is_group:
                label_widget = ipw.HTML(value = property_description)
            else:
                label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                
            boolean_widget = utils.Checkbox(
                layout = ipw.Layout(width = "200px"),
                value = False,
                indent = False
            )
            property_widget = ipw.VBox([label_widget, boolean_widget])
            property_dict = {property: boolean_widget}
        
        elif property_openbis_type == "DATE":
            if is_group:
                label_widget = ipw.HTML(value = property_description)
            else:
                label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                
            datepicker_widget = ipw.DatePicker(
                layout = ipw.Layout(width = "200px"), 
                value = datetime.date.today()
            )
            property_widget = ipw.VBox([label_widget, datepicker_widget])
            property_dict = {property: datepicker_widget}
        
        elif property_openbis_type == "TIMESTAMP":
            if is_group:
                label_widget = ipw.HTML(value = property_description)
            else:
                label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                
            text_widget = utils.Text(
                layout = ipw.Layout(width = "200px"), 
                placeholder = "",
                value = datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S")
            )
                    
            property_widget = ipw.VBox([label_widget, text_widget])
            property_dict = {property: text_widget}
        
        elif property_openbis_type == "INTEGER":
            if is_group:
                label_widget = ipw.HTML(value = property_description)
            else:
                label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                
            int_widget = utils.Text(
                layout = ipw.Layout(width = "200px")
            )
            
            def validate_input(change):
                """Remove non-numeric characters from input."""
                new_value = change['new']
                if not new_value.isdigit():  # Allow only digits
                    int_widget.value = ''.join(filter(str.isdigit, new_value))
            
            int_widget.observe(validate_input, names = 'value')
            
            property_widget = ipw.VBox([label_widget, int_widget])
            property_dict = {property: int_widget}
        
        elif property_openbis_type == "REAL":
            if is_group:
                label_widget = ipw.HTML(value = property_description)
            else:
                label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                
            float_widget = utils.Text(
                layout = ipw.Layout(width = "200px")
            )
            
            def validate_float_input(change):
                """Ensure input contains only a valid float format."""
                new_value = change['new']
                
                # Allow only numbers and one decimal point
                if new_value and not new_value.replace('.', '', 1).isdigit():
                    float_widget.value = ''.join(filter(lambda c: c.isdigit() or c == '.', new_value))
                    
                # Ensure only one decimal point
                if float_widget.value.count('.') > 1:
                    float_widget.value = float_widget.value[:-1]
            
            float_widget.observe(validate_float_input, names = 'value')
            
            property_widget = ipw.VBox([label_widget, float_widget])
            property_dict = {property: float_widget}

        elif property_openbis_type == "CONTROLLEDVOCABULARY":
            if is_group:
                label_widget = ipw.HTML(value = property_description)
            else:
                label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                
            property_vocabulary = DATA_MODEL["enums"][property_range]["permissible_values"].keys()
            property_vocabulary = list(property_vocabulary)
            dropdown_widget = utils.Dropdown(
                layout = ipw.Layout(width = "200px"), 
                options = property_vocabulary,
                value = property_vocabulary[0]
            )
            property_widget = ipw.VBox([label_widget, dropdown_widget])
            property_dict = {property: dropdown_widget}
        
        elif property_openbis_type == "JSON":
            if "slots" in self.all_schema_classes[property_range]:
                property_widget, property_dict = self.get_property_widgets_recursive(property, property_range)
            else:
                if is_group:
                    label_widget = ipw.HTML(value = property_description)
                else:
                    label_widget = ipw.HTML(value = f"<b>{property_description}</b>")
                    
                text_widget = utils.Text(
                    layout = ipw.Layout(width = "200px"), 
                    placeholder = "",
                )
                property_widget = ipw.VBox([label_widget, text_widget])
                property_dict = {property: text_widget}
        
        return property_widget, property_dict

    def get_property_widgets_recursive(self, property, property_range):
        widget_accordion_children = []
        property_dict = {property: {}}
        for slot in self.all_schema_classes[property_range]["slots"]:
            if self.all_schema_slots[slot]["annotations"]["openbis_type"] != "Not used":
                sub_property_widget, sub_property_dict = self.get_property_widget(slot, True)
                property_dict[property][slot] = sub_property_dict[slot]
                widget_accordion_children.append(sub_property_widget)
        property_widget = ipw.HBox(children = widget_accordion_children)
        property_description = self.all_schema_slots[property]["description"]
        
        return ipw.VBox([ipw.HTML(value = f"<b>{property_description}</b>"), property_widget]), property_dict
    
    def get_properties_recursive(self, object):
        # Get all the property types until the last parent class
        object = copy.deepcopy(object)
        object_copy = object
        while "is_a" in object_copy:
            parent_object = self.all_schema_classes[object_copy["is_a"]]
            # Used when the object inherits all the properties from another object
            if object["slots"]:
                object["slots"] = parent_object["slots"] + object["slots"]
            else:
                object["slots"] = parent_object["slots"]
                
            object_copy = parent_object
        
        return object

    def reset_properties_widgets(self):
        self.get_properties_widgets()

class ObjectMultipleSelectionWidget(ipw.HBox):
    def __init__(self, description):
        # Initialize the parent HBox
        super().__init__()
        self.description = description
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = description, disabled = False, 
                                             layout = ipw.Layout(width = '800px', height = '300px'), 
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
            
            options = [(f"{item.props['$name']} ({item.permId})", item.permId) for item in items]
            selector_options += options
        
        selector_options.sort(key=lambda x: (x[1], x[0]), reverse = True)
        self.selector.options = selector_options

class ObjectSelectionWidget(ipw.HBox):
    def __init__(self, type, widgets_types = {"dropdown": {"width": "982px"}}, checkboxes_descriptions = ["Name", "Registration date"]):
        # Initialize the parent HBox
        super().__init__()
        
        widgets_list = []
        self.dropdown_boxes = None
        self.details_textbox = None
        self.image_box = None
        self.schema_object_type = type
        
        self.description = DATA_MODEL["classes"][type]["title"]
        self.checkboxes_descriptions = checkboxes_descriptions
        
        if "dropdown" in widgets_types:
            self.dropdown = utils.Dropdown(
                description = self.description, 
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
    
    def load_dropdown_box(self, source = "openbis"):
        self.type = DATA_MODEL["classes"][self.schema_object_type]["annotations"]["openbis_label"].replace(" ", "_").upper()
        placeholder = DATA_MODEL["classes"][self.schema_object_type]["title"]
        
        if source == "openbis":
            items = utils.get_openbis_objects(
                OPENBIS_SESSION,
                type = self.type
            )
            if type == "SAMPLE":
                options = [(f"{item.props['$name']}", item.permId) for item in items if item.props["exists"] == "true"]
            else:
                options = [(f"{item.props['$name']}", item.permId) for item in items]
            
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
    
    # Function to handle changes in the objects dropdown
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
            if self.type in ["CRYSTAL", "MOLECULE", "CHEMICAL"]:
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
                # Erase file after downloading it
                shutil.rmtree(f"images/{object_dataset.permId}/")
            else:
                self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])

        if self.details_textbox:
            metadata_string = utils.get_metadata_string(OPENBIS_SESSION, object, self.schema_object_type, "", DATA_MODEL)
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
        items_names_permids = [(f"{item.code}", item.permId) for item in items]
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
    def __init__(self, sample_preparation_tasks):
        # Initialize the parent HBox
        super().__init__()
        self.selector = ipw.SelectMultiple()
        self.sample_preparation_options = [(DATA_MODEL["classes"][sample_preparation]["description"], sample_preparation) for sample_preparation in sample_preparation_tasks]
        self.selector.description = "Processes"
        self.selector.disabled = False
        self.selector.style = {'description_width': "65px"}
        self.selector.options = self.sample_preparation_options
        selector_height = len(self.sample_preparation_options) * 18
        selector_height = f"{selector_height}px"
        self.selector.layout = ipw.Layout(width = '220px', height = selector_height)
        self.children = [self.selector]

