import json
from pathlib import Path
import IPython.display
import ipywidgets as ipw
from pybis import Openbis
import ipyfilechooser
import yaml
import pandas as pd
import copy
from datetime import datetime

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

def sort_dropdown(sorting_checkboxes, dropdown_box):
    current_value = dropdown_box.value
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
        dropdown_list = sort_dataframe(df, columns, ascending)

    dropdown_list.insert(0, dropdown_box.options[0])
    dropdown_box.options = dropdown_list
    dropdown_box.value = current_value

def dataframe_to_list_of_tuples(df):
    return list(df.itertuples(index = False, name = None))

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

def filter_samples(openbis_session, samples):
    """
    Function for hiding intermediate samples which result 
    from the different steps of samples preparation.

    Args:
        samples (_type_): _description_
    """
    
    selected_samples = []
    for sample in samples:
        sample = openbis_session.get_object(sample.permId)
        if sample.props["exists"] == "true":
            selected_samples.append(sample)
    return selected_samples

def get_openbis_parents_recursive(openbis_session, object, object_parents_metadata):
    object_parents_metadata.append([object.attrs.type, object.attrs.permId, object.attrs.registrationDate,object.props['$name']])
    
    for parent in object.parents:
        get_openbis_parents_recursive(openbis_session, openbis_session.get_object(parent), object_parents_metadata)
    return object_parents_metadata

def load_openbis_elements_list(openbis_session, type, project = None):
    if type == "EXPERIMENT":
        items = openbis_session.get_collections(type = type, project = project)
        items_names_permids = [(f"{item.props['$name']} ({item.attrs.identifier})", item.permId) for item in items]
    elif type == "PROJECT":
        items = openbis_session.get_projects()
        items_names_permids = [(item.attrs.identifier, item.permId) for item in items]
    else:
        items = openbis_session.get_objects(type = type, project = project)
        if type == "SAMPLE":
            items = filter_samples(openbis_session, items)
            items_names_permids = [(f"{item.props['$name']}", item.permId) for item in items]
        elif type == "SUBSTANCE":
            items_names_permids = [(f"{item.props['empa_number']}{item.props['batch']} ({item.permId})", item.permId) for item in items]
        else:
            items_names_permids = [(f"{item.props['$name']} ({item.permId})", item.permId) for item in items]
    
    return items_names_permids

def load_dropdown_elements(openbis_session, type, dropdown, sorting_checkboxes, label):
    items_names_permids = load_openbis_elements_list(openbis_session, type)
    items_names_permids.insert(0, (f'Select {label}...', -1))
    dropdown.options = items_names_permids
    dropdown.value = -1
    sort_dropdown(sorting_checkboxes, dropdown)
    for checkbox in sorting_checkboxes.children[1:3]:
        checkbox.observe(lambda change: sort_dropdown(sorting_checkboxes, dropdown), names='value')
    
def get_next_experiment_code(openbis_session):
    experiments = openbis_session.get_experiments(type = "EXPERIMENT")
    experiment_number = max(int(exp.code.rsplit('_')[-1]) for exp in experiments) + 1
    return f"{experiments[0].code.rsplit('_')[0]}_{experiment_number}"

def create_experiment_in_openbis(openbis_session, project_id, experiment_name):
    experiment_code = get_next_experiment_code(openbis_session)
    experiment = openbis_session.new_experiment(code = experiment_code, type = "EXPERIMENT", project = project_id, props = {"$name": experiment_name})
    experiment.save()

def create_openbis_object(openbis_session, **kwargs):
    openbis_object = openbis_session.new_object(**kwargs)
    openbis_object.save()
    return openbis_object

def create_openbis_collection(openbis_session, **kwargs):
    collection = openbis_session.new_collection(**kwargs)
    collection.save()
    return collection

def get_current_datetime():
    return datetime.now()

def convert_datetime_to_string(dt):
    return dt.strftime('%Y%m%d%H%M%S')

def remove_digits_from_string(string):
    return ''.join([i for i in string if not i.isdigit()])

def is_nan(var):
    """Function to verify whether it is a NaN."""
    return var != var

# def get_parent_child_relationships_nested(openbis_session, selected_object, parent_child_relationships=None):
#     """
#     Function to get all the parent-child relations together with the information about the objects from openBIS.
#     This is a recursive function because the objects inside openBIS are like trees containing multiple
#     relations with other objects.

#     Parameters
#     ----------
#     selected_object : pybis.sample.Sample
#         Selected openBIS object.
#     parent_child_relationships : dict, optional
#         Dictionary containing the openBIS objects and the relations between them. The default is None.

#     Returns
#     -------
#     parent_child_relationships : TYPE
#         Dictionary containing the openBIS objects and the relations between them.

#     """
#     if parent_child_relationships is None:
#         parent_child_relationships = {}

#     # Fetch parents of the current publication
#     parents = selected_object.parents
    
#     # Get selected object ID
#     selected_object_id = selected_object.permId

#     if parents is not None:
#         # Add current parents to the dictionary with child-parent relationships
#         for parent_id in parents:
#             parent = openbis_session.get_sample(parent_id)
#             parent_id = parent.permId
#             parent_props = parent.props.all()
#             parent_props["eln_object_type"] = str(parent.type)
#             parent_child_relationships[selected_object_id][parent_id] = parent_props
#             get_parent_child_relationships_nested(openbis_session, parent, parent_child_relationships[selected_object_id])

#     return parent_child_relationships

def get_parent_child_relationships_nested(openbis_session, selected_object, parent_child_relationships=None):
    """
    Function to get all the parent-child relations together with the information about the objects from openBIS.
    This is a recursive function because the objects inside openBIS are like trees containing multiple
    relations with other objects.

    Parameters
    ----------
    selected_object : pybis.sample.Sample
        Selected openBIS object.
    parent_child_relationships : dict, optional
        Dictionary containing the openBIS objects and the relations between them. The default is None.

    Returns
    -------
    parent_child_relationships : TYPE
        Dictionary containing the openBIS objects and the relations between them.

    """
    if parent_child_relationships is None:
        parent_child_relationships = {}

    # Fetch parents of the current publication
    parents = selected_object.parents

    if parents is not None:
        # Add current parents to the dictionary with child-parent relationships
        parent_child_relationships["has_part"] = []
        for parent_id in parents:
            parent = openbis_session.get_sample(parent_id)
            parent_props = parent.props.all()
            for prop in parent_props:
                if openbis_session.get_property_type(prop).dataType == "SAMPLE" and parent_props[prop] is not None:
                    prop_object = openbis_session.get_object(parent_props[prop])
                    parent_props[prop] = {"perm_id": parent_props[prop]}
                    parent_props[prop]["type"] = str(prop_object.type)
                    parent_props[prop].update(prop_object.props.all())
                    parent_props[prop] = get_parent_child_relationships_nested(openbis_session, prop_object, parent_props[prop])
                    
            parent_props["perm_id"] = parent.permId
            parent_props["type"] = str(parent.type)
            parent_props = get_parent_child_relationships_nested(openbis_session, parent, parent_props)
            parent_child_relationships["has_part"].append(parent_props)

    return parent_child_relationships