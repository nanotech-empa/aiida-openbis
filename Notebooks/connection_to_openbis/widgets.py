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

CONFIG = utils.read_json("config.json")
CONFIG_ELN = utils.read_json("eln_config.json")
OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis(CONFIG_ELN["url"], CONFIG_ELN["token"])

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
        self.selector.options = utils.load_openbis_elements_list(OPENBIS_SESSION, "ANALYSIS", project_id)

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
        self.selector.options = utils.load_openbis_elements_list(OPENBIS_SESSION, "CODE")

class PublicationSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Dropdown with sorting and checkboxes
        self.dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget(
            'Publication', ipw.Layout(width='982px'), 
            {'description_width': "110px"}, [-1], "vertical"
        )
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        utils.load_dropdown_elements(
            OPENBIS_SESSION, "PUBLICATION_CUSTOM", self.dropdown_boxes.children[0], 
            self.dropdown_boxes.children[1], "publication"
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
        measurements_list.extend(
            utils.load_openbis_elements_list(OPENBIS_SESSION, "1D_MEASUREMENT", project_id)
        )
        measurements_list.extend(
            utils.load_openbis_elements_list(OPENBIS_SESSION, "2D_MEASUREMENT", project_id)
        )
        self.selector.options = measurements_list
        
class DraftSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Select multiple drafts
        self.selector = utils.SelectMultiple(description = 'Drafts', disabled = False, 
                                             layout = ipw.Layout(width = '800px'), 
                                             style = {'description_width': "110px"})
        
        self.children = [self.selector]
    
    def load_selector(self, project_id):
        self.selector.options = utils.load_openbis_elements_list(OPENBIS_SESSION, "DRAFT", project_id)
        
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
        self.selector.options = utils.load_openbis_elements_list(OPENBIS_SESSION, "RESULTS", project_id)
        
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
        self.selector.options = utils.load_openbis_elements_list(OPENBIS_SESSION, "GRANT")

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
        self.selector.options = utils.load_openbis_elements_list(OPENBIS_SESSION, "AUTHOR")
        
class MoleculeSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Dropdown with sorting and checkboxes
        self.dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget(
            'Molecule', ipw.Layout(width='315px'), 
            {'description_width': "110px"}, [-1], "vertical"
        )
        
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
        utils.load_dropdown_elements(
            OPENBIS_SESSION, "MOLECULE", self.dropdown_boxes.children[0], 
            self.dropdown_boxes.children[1], "molecule"
        )

class SubstanceSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Dropdown with sorting and checkboxes
        self.dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget(
            'Substance', ipw.Layout(width='315px'), 
            {'description_width': "110px"}, [-1], "vertical"
        )
        
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
        utils.load_dropdown_elements(
            OPENBIS_SESSION, "SUBSTANCE", self.dropdown_boxes.children[0], 
            self.dropdown_boxes.children[1], "substance"
        )
    
    # Function to handle changes in the substances dropdown
    def load_substance_metadata(self, change):
        if self.dropdown_boxes.children[0].value == -1:
            self.details_textbox.value = ''
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])
            return
        
        # Get substance metadata
        property_list = CONFIG["objects"]["Substance"]["properties"]
        substance_property_list = CONFIG["objects"]["Molecule"]["properties"]
        material_object = OPENBIS_SESSION.get_object(self.dropdown_boxes.children[0].value)
        
        # Get substance image
        substance_object = OPENBIS_SESSION.get_object(material_object.props.all()['has_molecule'])
        substance_metadata = substance_object.props.all()
        substance_dataset = substance_object.get_datasets(type="ELN_PREVIEW")[0]

        if substance_dataset:
            substance_dataset.download(destination="images")
            material_image_filepath = substance_dataset.file_list[0]
            self.image_box.value = utils.read_file(f"images/{substance_dataset.permId}/{material_image_filepath}")
            shutil.rmtree(f"images/{substance_dataset.permId}/")
        else:
            self.image_box.value = utils.read_file(CONFIG["default_image_filepath"])

        material_metadata = material_object.props.all()
        material_metadata_string = ""
        for prop_key in property_list:
            prop_title = CONFIG["properties"][prop_key]["title"]
            prop_value = material_metadata.get(prop_key)
            if CONFIG["properties"][prop_key]["property_type"] == "QUANTITY_VALUE" and prop_value is not None:
                prop_value = json.loads(prop_value)
                material_metadata_string += f"{prop_title}: {prop_value['value']} {prop_value['unit']}\n"
            elif prop_key == "has_molecule":
                material_metadata_string += f"{prop_title}:\n"
                for mol_prop_key in substance_property_list:
                    mol_prop_title = CONFIG["properties"][mol_prop_key]["title"]
                    mol_prop_value = substance_metadata.get(mol_prop_key)
                    material_metadata_string += f"- {mol_prop_title}: {mol_prop_value}\n"
            else:
                material_metadata_string += f"{prop_title}: {prop_value}\n"

        self.details_textbox.value = material_metadata_string
    
class MaterialSelectionWidget(ipw.Output):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Dropdown with sorting and checkboxes
        self.dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget(
            '', ipw.Layout(width='315px'), {'description_width': "110px"}, 
            [-1], "horizontal"
        )
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
        utils.load_dropdown_elements(
            OPENBIS_SESSION, object_type, self.dropdown_boxes.children[0], 
            self.dropdown_boxes.children[1], placeholder
        )

class ProjectSelectionWidget(ipw.VBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        self.dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget(
            'Project', ipw.Layout(width = '982px'),
            {'description_width': "110px"}, [-1], "vertical"
        )
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        utils.load_dropdown_elements(
            OPENBIS_SESSION, "PROJECT", self.dropdown_boxes.children[0], 
            self.dropdown_boxes.children[1], "project"
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
            description='Experiment', disabled=False, layout = ipw.Layout(width = '993px'), 
            style = {'description_width': "110px"}, options = [-1]
        )
        self.sorting_checkboxes = utils.SortingCheckboxes("130px", "60px", "200px")
        self.dropdown_details = ipw.HBox([self.dropdown, self.create_new_experiment_button])
        self.add_experiment_output = ipw.Output()
        self.children = [self.dropdown_details, self.sorting_checkboxes, self.add_experiment_output]
        
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
        utils.load_dropdown_elements(
            OPENBIS_SESSION, "EXPERIMENT", self.dropdown, self.sorting_checkboxes, "experiment"
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
        utils.load_dropdown_elements(OPENBIS_SESSION, "EXPERIMENT", self.dropdown, self.sorting_checkboxes, "experiment")
        with self.add_experiment_output:
            clear_output()

class SampleSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Dropdown with sorting and checkboxes
        self.dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget(
            'Sample', ipw.Layout(width='385px'), 
            {'description_width': "110px"}, [-1], "vertical"
        )
        
        # Details textbox (disabled for display purposes)
        self.details_textbox = ipw.Textarea(
            description="", disabled=True, layout=ipw.Layout(width='589px', height='250px')
        )
        
        self.children = [self.dropdown_boxes, self.details_textbox]
    
    def load_dropdown_box(self):
        utils.load_dropdown_elements(
            OPENBIS_SESSION, "SAMPLE", self.dropdown_boxes.children[0], 
            self.dropdown_boxes.children[1], "sample"
        )

class InstrumentSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Dropdown with sorting and checkboxes
        self.dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget(
            'Instrument', ipw.Layout(width='982px'), 
            {'description_width': "110px"}, [-1], "vertical"
        )
        
        self.children = [self.dropdown_boxes]
    
    def load_dropdown_box(self):
        utils.load_dropdown_elements(
            OPENBIS_SESSION, "INSTRUMENT", self.dropdown_boxes.children[0], 
            self.dropdown_boxes.children[1], "instrument"
        )
            
class ProductSelectionWidget(ipw.HBox):
    def __init__(self):
        # Initialize the parent HBox
        super().__init__()
        
        # Dropdown with sorting and checkboxes
        self.dropdown_boxes = utils.DropdownwithSortingCheckboxesWidget(
            'Product', ipw.Layout(width='315px'), 
            {'description_width': "110px"}, [-1], "vertical"
        )
        
        # Details textbox (disabled for display purposes)
        self.details_textbox = ipw.Textarea(
            description="", disabled=True, layout=ipw.Layout(width='415px', height='250px')
        )
        
        # Image box (displaying a default image)
        self.image_box = ipw.Image(
            value=utils.read_file(self.custom_config["default_image_filepath"]), format='jpg', width='220px', 
            height='250px', layout=ipw.Layout(border='solid 1px #cccccc')
        )
        
        self.children = [self.dropdown_boxes, self.details_textbox, self.image_box]
    
    def load_dropdown_box(self):
        utils.load_dropdown_elements(
            OPENBIS_SESSION, "PRODUCT", self.dropdown_boxes.children[0], 
            self.dropdown_boxes.children[1], "substance"
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

class ObjectPropertiesWidgets(ipw.VBox):
    def __init__(self, task):
        super().__init__()
        
        self.properties_widgets = {}
        properties = CONFIG["objects"][task]["properties"]
        for prop_key in properties:
            property = CONFIG["properties"][prop_key]
            if property["property_widget"] == "TEXT":
                prop_widget = utils.Text(
                    description = property["title"], disabled = False, 
                    layout = ipw.Layout(width = property["box_layout"]["width"]), 
                    placeholder = property["placeholder"], 
                    style = {'description_width': property["box_layout"]["description_width"]}
                )
                
            elif property["property_widget"] == "TEXTAREA":
                prop_widget = utils.Textarea(
                    description = property["title"], disabled = False, 
                    layout = ipw.Layout(width = property["box_layout"]["width"]), 
                    placeholder = property["placeholder"], 
                    style = {'description_width': property["box_layout"]["description_width"]}
                )
            
            elif property["property_widget"] == "INTTEXT":
                prop_widget = utils.IntText(
                    description = property["title"], disabled = False, 
                    layout = ipw.Layout(width = property["box_layout"]["width"]), 
                    placeholder = property["placeholder"], 
                    style = {'description_width': property["box_layout"]["description_width"]}
                )
            
            elif property["property_widget"] == "DROPDOWN":
                prop_widget = utils.Dropdown(
                    description = property["title"], disabled = False, 
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
                
            self.properties_widgets[prop_key] = prop_widget
            
        self.children = list(self.properties_widgets.values())
        