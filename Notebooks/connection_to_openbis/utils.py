import json
from pathlib import Path
import IPython.display
import ipywidgets as ipw
from pybis import Openbis
import ipyfilechooser
import yaml
import pandas as pd
import copy
import datetime
import os
from IPython.display import display

def full_listdir(path):
    return [f"{path}{os.sep}{filepath}" for filepath in os.listdir(path)]

def upload_datasets(openbis_session, method_object, support_files_widget, dataset_type):
    for filename in support_files_widget.value:
        file_info = support_files_widget.value[filename]
        save_file(file_info['content'], filename)
        openbis_session.new_dataset(type = dataset_type, sample = method_object, files = [filename]).save()
        os.remove(filename)

def is_valid_json(string):
    try:
        json.loads(string)
        return True
    except ValueError:
        return False

def get_aiidalab_eln_config():
    eln_config = Path.home() / ".aiidalab" / "aiidalab-eln-config.json"
    eln_config.parent.mkdir(
        parents=True, exist_ok=True
    )  # making sure that the folder exists.

    with open(eln_config) as file:
        config = json.load(file)
        eln_url = config["default"]
        eln_config = config[eln_url]
        eln_token = eln_config["token"]
        return {"url": eln_url, "token": eln_token}

def read_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)

def read_yaml(filepath: str):
    return yaml.safe_load(open(filepath, 'r'))

def create_json(json_content, filename):
    with open(filename, 'w') as file:
        json.dump(json_content, file, indent=4)
        
def DropdownwithSortingCheckboxesWidget(*args):
    dropdown = Dropdown(description=args[0], disabled=False, layout = args[1], style = args[2], options = args[3])
    sorting_checkboxes_list = SortingCheckboxes("130px", "60px", "200px")
    widgets = [dropdown, sorting_checkboxes_list]
    if args[4] == "horizontal":
        return ipw.HBox(widgets)
    else:
        return ipw.VBox(widgets)

def IntSliderwithTextWidget(*args):
    intslider = IntSlider(value = args[0], description = args[1], min = args[2][0], max = args[2][-1], disabled = False, layout = args[3], style = args[4])
    textbox = Text(value = '', description = '', placeholder= args[5], disabled = False, layout = args[6])
    return ipw.HBox([intslider, textbox])

def SortingCheckboxes(*args):
    return ipw.HBox([
        ipw.Label(value = "Sort by:", layout = ipw.Layout(width = args[0], display = "flex", justify_content='flex-end')),
        Checkbox(description = 'Name', value = False, disabled = False, layout = ipw.Layout(width = args[1]), indent = False),
        Checkbox(description = 'Registration date', value = False, disabled = False, layout = ipw.Layout(width = args[2]), indent = False)
    ])

def FloatTextwithDropdownWidget(*args):
    floatbox = FloatText(description=args[0], disabled=False, layout = args[1], value = args[2], style = args[3])
    dropdown = Dropdown(description='', disabled=False, layout = args[4], options = args[5], value = args[6])
    return ipw.HBox([floatbox, dropdown])

def FloatSlider(**kwargs):
    return ipw.FloatSlider(**kwargs)

def IntSlider(**kwargs):
    return ipw.IntSlider(**kwargs)

def HTMLbox(**kwargs):
    return ipw.HTML(**kwargs)

def HTML(**kwargs):
    return IPython.display.HTML(**kwargs)

def Javascript(**kwargs):
    return IPython.display.Javascript(**kwargs)

def Markdown(**kwargs):
    return IPython.display.Markdown(**kwargs)

def FileChooser(**kwargs):
    return ipyfilechooser.FileChooser(**kwargs)

def Image(**kwargs):
    return ipw.Image(**kwargs)

def Button(**kwargs):
    return ipw.Button(**kwargs)

def Checkbox(**kwargs):
    return ipw.Checkbox(**kwargs)

def Textarea(**kwargs):
    return ipw.Textarea(**kwargs)

def Text(**kwargs):
    return ipw.Text(**kwargs)
    
def Radiobuttons(**kwargs):
    return ipw.RadioButtons(**kwargs)

def Dropdown(**kwargs):
    return ipw.Dropdown(**kwargs)

def IntText(**kwargs):
    return ipw.IntText(**kwargs)

def FloatText(**kwargs):
    return ipw.FloatText(**kwargs)

def SelectMultiple(**kwargs):
    return ipw.SelectMultiple(**kwargs)

def read_file(filename):
    return open(filename, "rb").read()

def sort_dataframe(df, columns, ascending_columns):
    df = df.sort_values(by = columns, ascending = ascending_columns)
    return dataframe_to_list_of_tuples(df)

def dataframe_to_list_of_tuples(df):
    return list(df.itertuples(index = False, name = None))

def sort_dropdown(sorting_checkboxes, dropdown, columns_names, column_ascending):
    current_value = dropdown.value
    dropdown_list = list(dropdown.options[1:]) # Default -1 message should not be sorted.
    df = pd.DataFrame(dropdown_list, columns = columns_names)
    
    # Determine sort columns and order based on checkboxes
    columns, ascending = [], []
    for i,e in enumerate(sorting_checkboxes.children[1:]):
        if e.value:
            columns.append(columns_names[i])
            ascending.append(column_ascending[i])
    
    if columns:
        dropdown_list = sort_dataframe(df, columns, ascending)

    dropdown_list.insert(0, dropdown.options[0])
    dropdown.options = dropdown_list
    dropdown.value = current_value

def connect_openbis_aiida():
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

def connect_openbis(eln_url, eln_token):
    try:
        session_data = {"url": eln_url, "token": eln_token}
        openbis_session = Openbis(eln_url, verify_certificates = False)
        openbis_session.set_token(eln_token)
    except ValueError:
        print("Session is no longer valid. Please check if the token is still valid.")
        openbis_session = None
        session_data = {}
    
    return openbis_session, session_data

def save_file(file_content, filename):
    with open(filename, 'wb') as f:  # Save file content
        f.write(file_content)

def get_openbis_parents_recursive(openbis_session, object, object_parents_metadata):
    object_parents_metadata.append([object.attrs.type, object.attrs.permId, object.attrs.registrationDate,object.props['$name']])
    
    for parent in object.parents:
        get_openbis_parents_recursive(openbis_session, openbis_session.get_object(parent), object_parents_metadata)
    return object_parents_metadata

def get_next_experiment_code(openbis_session):
    experiments = openbis_session.get_experiments(type = "EXPERIMENT")
    experiment_number = 1
    if experiments:
        experiment_number = max(int(exp.code.rsplit('_')[-1]) for exp in experiments) + 1
    return f"EXPERIMENT_{experiment_number}"

def create_experiment_in_openbis(openbis_session, project_id, experiment_name):
    experiment_code = get_next_experiment_code(openbis_session)
    experiment = openbis_session.new_experiment(
        code = experiment_code, 
        type = "EXPERIMENT", 
        project = project_id, 
        props = {
            "$name": experiment_name, 
            "$default_collection_view": "IMAGING_GALLERY_VIEW"
        }
    )
    experiment.save()

def create_openbis_dataset(openbis_session, **kwargs):
    openbis_ds = openbis_session.new_dataset(**kwargs)
    openbis_ds.save()

def create_openbis_object(openbis_session, **kwargs):
    openbis_object = openbis_session.new_object(**kwargs)
    openbis_object.save()
    return openbis_object

def create_openbis_collection(openbis_session, **kwargs):
    collection = openbis_session.new_collection(**kwargs)
    collection.save()
    return collection

def get_openbis_collections(openbis_session, **kwargs):
    return openbis_session.get_collections(**kwargs)

def get_openbis_projects(openbis_session, **kwargs):
    return openbis_session.get_projects(**kwargs)

def get_openbis_objects(openbis_session, **kwargs):
    return openbis_session.get_objects(**kwargs)

def get_current_datetime():
    return datetime.datetime.now()

def convert_datetime_to_string(dt):
    return dt.strftime('%Y%m%d%H%M%S')

def remove_digits_from_string(string):
    return ''.join([i for i in string if not i.isdigit()])

def is_nan(var):
    """Function to verify whether it is a NaN."""
    return var != var

def get_parent_child_relationships_nested(openbis_session, selected_object, parent_child_relationships = None, object_history = None):
    """
    Function to get all the parent-child relations together with the information about the objects from openBIS.
    This is a recursive function because the objects inside openBIS are like trees containing multiple
    relations with other objects.

    """
    if parent_child_relationships is None:
        parent_child_relationships = {}
    
    if object_history is None:
        object_history = {}
    
    if selected_object.identifier not in object_history:
        object_history[selected_object.identifier] = selected_object

    # Fetch parents of the current publication
    parents = selected_object.parents

    if parents is not None:
        # Add current parents to the dictionary with child-parent relationships
        parent_child_relationships["has_part"] = []
        for parent_id in parents:
            if parent_id in object_history:
                parent = object_history[parent_id]
            else:
                parent = openbis_session.get_sample(parent_id)
                object_history[parent.identifier] = parent
            parent_props = parent.props.all()
            for prop in parent_props:
                if openbis_session.get_property_type(prop).dataType == "SAMPLE" and parent_props[prop] is not None:
                    
                    if openbis_session.get_property_type(prop).multiValue:
                        prop_objects = openbis_session.get_object(parent_props[prop]) # When multivalued, this function returns a list of openBIS objects
                        parent_props_objects = []
                        for obj in prop_objects:
                            parent_props[prop] = {"perm_id": obj.permId}
                            parent_props[prop]["type"] = str(obj.type)
                            parent_props[prop]["registrationDate"] = obj.registrationDate
                            parent_props[prop].update(obj.props.all())
                            parent_props[prop], _ = get_parent_child_relationships_nested(openbis_session, obj, parent_props[prop], object_history)
                            parent_props_objects.append(parent_props[prop])
                        parent_props[prop] = parent_props_objects
                        
                    else:
                        prop_object = openbis_session.get_object(parent_props[prop])
                        parent_props[prop] = {"perm_id": parent_props[prop]}
                        parent_props[prop]["type"] = str(prop_object.type)
                        parent_props[prop]["registrationDate"] = prop_object.registrationDate
                        parent_props[prop].update(prop_object.props.all())
                        parent_props[prop], _ = get_parent_child_relationships_nested(openbis_session, prop_object, parent_props[prop], object_history)
                    
            parent_props["perm_id"] = parent.permId
            parent_props["type"] = str(parent.type)
            parent_props["registrationDate"] = parent.registrationDate
            parent_props, _ = get_parent_child_relationships_nested(openbis_session, parent, parent_props, object_history)
            parent_child_relationships["has_part"].append(parent_props)

    return parent_child_relationships, object_history

def get_object_type(openbis_object_type, config):
    for object_type in config["objects"]:
        if config["objects"][object_type]["openbis_object_type"] == openbis_object_type:
            return object_type

def get_metadata_string(openbis_session, object, metadata_string, config):
    
    object_type = get_object_type(object.type, config)
    
    for prop_key in config["objects"][object_type]["properties"]:
        prop_title = config["properties"][prop_key]["title"]
        prop_datatype = config["properties"][prop_key]["property_type"]
        prop_string = get_property_string(openbis_session, object, prop_title, prop_key, prop_datatype, config)
        metadata_string += prop_string

    return metadata_string
        
def get_property_string(openbis_session, object, prop_title, prop_key, prop_datatype, config):
    property_string = ""
    object_metadata = object.props.all()
    if prop_datatype == "QUANTITY_VALUE":
        value = object_metadata.get(prop_key)
        if value:
            prop_dict = json.loads(value)
            property_string = f"{prop_title}: {prop_dict['has_value']} {prop_dict['has_unit']}\n"
        else:
            property_string = f"{prop_title}: {value}\n"
    elif prop_datatype == "PARENT_OBJECT":
        parent_object_type = config["properties"][prop_key]["parent_object_type"]
        parents_objects = object.get_parents()
        for parent_object in parents_objects:
            if parent_object.type == parent_object_type:
                # For each property related to an object type there should be no more than one parent
                property_string = f"-------\n{prop_title}:\n" 
                property_string = get_metadata_string(openbis_session, parent_object, property_string, config)
                property_string = f"{property_string}-------\n"
                break
    else:
        property_string = f"{prop_title}: {object_metadata.get(prop_key)}\n"
    
    return property_string

def find_first_object_permid(openbis_session, openbis_object, openbis_type):
    first_object = openbis_object
    object_parents_identifiers = openbis_object.parents
    for parent_identifier in object_parents_identifiers:
        parent_object = openbis_session.get_object(parent_identifier)
        if parent_object.type == openbis_type:
            first_object = find_first_object_permid(openbis_session, parent_object, openbis_type)
    
    return first_object

def get_sample_details(openbis_session, sample_object, sample_preparation_types, slabs_types):
    sample_parents_metadata = get_openbis_parents_recursive(openbis_session, sample_object, [])
    sample_details = {
        "sample_preparations": [], 
        "materials": [], 
        "sample_preparation_datetimes": [],
        "last_sample_preparation_object": None
    }
    number_parents = len(sample_parents_metadata)
    parent_idx = 0
    while parent_idx < number_parents:
        parent_metadata = sample_parents_metadata[parent_idx]
        
        if parent_metadata[0] == "DEPOSITION":
            deposition_object = openbis_session.get_object(parent_metadata[1])
            
            # Get molecule used in the specific deposition
            molecule_metadata = []
            for parent_id in deposition_object.parents:
                deposition_parent_object = openbis_session.get_object(parent_id)
                if deposition_parent_object.type == "MOLECULE":
                    deposition_parent_object_metadata = deposition_parent_object.props.all()
                    molecule_metadata = [deposition_parent_object_metadata["$name"], deposition_parent_object.permId]
                    
            if molecule_metadata: # If deposition does not contain any molecule, the app must not try to display it
                sample_metadata_string = f"> {parent_metadata[0]} ({parent_metadata[3]}, {parent_metadata[1]}, {parent_metadata[2]}) [{molecule_metadata[0]} ({molecule_metadata[1]})]"
            else:
                sample_metadata_string = f"> {parent_metadata[0]} ({parent_metadata[3]}, {parent_metadata[1]}, {parent_metadata[2]})"
            parent_idx += 1
        else:
            sample_metadata_string = f"> {parent_metadata[0]} ({parent_metadata[3]}, {parent_metadata[1]}, {parent_metadata[2]})"
        
        if parent_metadata[0] in sample_preparation_types:
            if sample_metadata_string not in sample_details["sample_preparations"]:
                sample_details["sample_preparations"].append(sample_metadata_string)
                sample_details["sample_preparation_datetimes"].append(parent_metadata[2])
            
            # Get the last sample preparation method performed on the sample in order to search the correct experiment where the sample is being used.
            if sample_details["last_sample_preparation_object"] is None:
                sample_details["last_sample_preparation_object"] = parent_metadata[1]
                
        elif parent_metadata[0] in slabs_types:
            if sample_metadata_string not in sample_details["materials"]:
                sample_details["materials"].append(sample_metadata_string)
            
        parent_idx += 1
    
    # Parse the datetime strings and zip them with sample preparations and materials
    parsed_datetimes = [datetime.datetime.strptime(dt, "%Y-%m-%d %H:%M:%S") for dt in sample_details["sample_preparation_datetimes"]]
    zipped_lists = list(zip(parsed_datetimes, sample_details["sample_preparations"]))

    # Sort by the datetime
    sorted_zipped_lists = sorted(zipped_lists, key=lambda x: x[0], reverse = True)

    # Unpack the sorted values back into sample_strings
    if sample_details["sample_preparations"]:
        _, sample_details["sample_preparations"] = zip(*sorted_zipped_lists)
    
    return sample_details

def get_last_sample_experiment(openbis_session, sample_object, last_sample_preparation_object):
    # Automatically select the experiment where the last sample preparation task was saved
    last_sample_experiment = openbis_session.get_experiment(last_sample_preparation_object.attrs.experiment)
    
    # Check if the sample is being used in a set of measurements that were created in another experiment after the preparation steps
    sample_children = sample_object.get_children()
    for child in sample_children:
        if child.registrationDate > last_sample_preparation_object.registrationDate:
            last_sample_experiment = child.experiment
    
    return last_sample_experiment

def get_last_sample_instrument(openbis_session, last_sample_preparation_object):
    # Automatically select the instrument used in the last sample preparation task
    for parent in last_sample_preparation_object.parents:
        parent_object = openbis_session.get_object(parent)
        if parent_object.type == "INSTRUMENT":
            instrument_object = parent_object
            break
    
    return instrument_object

def get_object_property(property_type, property_widget):
    property_value = None
    if property_type == "QUANTITY_VALUE":
        property_value = json.dumps(
            {
                "has_value": property_widget.children[0].value,
                "has_unit": property_widget.children[1].value
            }
        )
    elif property_type == "JSON":
        property_value = property_widget.value
        property_value = property_value.replace("'", '"')
        if property_value:
            if is_valid_json(property_value) == False:
                display(Javascript(data = "alert('There is a JSON property that is not a valid JSON object.')"))
    elif property_type == "DATE":
        property_value = str(property_widget.value)
    elif property_type == "MULTIVALUE_BOOLEAN":
        property_value = []
        for child in property_widget.children:
            property_value.append(child.value)
    else:
        property_value = property_widget.value
    
    return property_value