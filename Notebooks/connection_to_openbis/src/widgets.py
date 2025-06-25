import ipywidgets as ipw
import utils
from IPython.display import display
import pandas as pd
import json

class CreateSampleWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_material_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select material</span>"
        )
        
        material_type_options = [
            ("Select material type...", "-1"),
            ("Crystal", "CRYSTAL"),
            ("2D layer material", "2D_LAYER_MATERIAL"),
            ("Wafer substrate", "WAFER_SUBSTRATE")
        ]
        
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
            material_objects = utils.get_openbis_objects(
                self.openbis_session,
                type = material_type
            )
            materials_objects_names_permids = [(obj.props["$name"], obj.permId) for obj in material_objects]
            material_options += materials_objects_names_permids
            material_dropdown.options = material_options
            
            def sort_material_dropdown(change):
                options = material_options[1:]
                
                df = pd.DataFrame(options, columns=["$name", "registration_date"])
                if name_checkbox.value and not registration_date_checkbox.value:
                    df = df.sort_values(by="$name", ascending=True)
                elif not name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by="registration_date", ascending=False)
                elif name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by=["$name", "registration_date"], ascending=[True, False])

                options = list(df.itertuples(index=False, name=None))
                options.insert(0, material_options[0])
                material_dropdown.options = options
            
            def load_material_details(change):
                obj_permid = material_dropdown.value
                if obj_permid == "-1":
                    return
                else:
                    obj = self.openbis_session.get_object(obj_permid)
                    obj_props = obj.props.all()
                    obj_name = obj_props.get("$name", "")
                    obj_details_string = "<div style='border: 1px solid grey; padding: 10px; margin: 10px;'>"
                    for key, value in obj_props.items():
                        if value:
                            prop_type = self.openbis_session.get_property_type(key)
                            prop_label = prop_type.label
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
                sample_props = {
                    "$name": sample_name,
                    "exists": True
                }
                utils.create_openbis_object(
                    self.openbis_session,
                    type = "SAMPLE",
                    collection = "/MATERIALS/SAMPLES/SAMPLE_COLLECTION",
                    props = sample_props,
                    parents = [material_object]
                )
                
                display(utils.Javascript(data = "alert('Sample created successfully!')"))
                
                # Clear interface
                self.material_type_dropdown.value = "-1"
                self.sample_name_textbox.value = ""     

class RegisterProcessWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.select_experiment_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select experiment</span>"
        )
        
        self.select_sample_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Select sample</span>"
        )
        
        self.sample_history_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Sample history</span>"
        )
        
        self.register_new_processes_title = ipw.HTML(
            value = "<span style='font-weight: bold; font-size: 20px;'>Register new steps</span>"
        )
        
        self.select_experiment_dropdown = SelectExperimentWidget(self.openbis_session)
        self.select_sample_dropdown = SelectSampleWidget(self.openbis_session)
        self.sample_history_vbox = SampleHistoryWidget(self.openbis_session)
        self.register_new_processes_vbox = RegisterProcessStepWidget(self.openbis_session)
        
        self.add_process_template_button = ipw.Button(
            description = 'Add process template', 
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
                self.add_process_template_button,
                self.add_process_step_button
            ]
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
            self.select_experiment_title,
            self.select_experiment_dropdown,
            self.select_sample_title,
            self.select_sample_dropdown,
            self.sample_history_title,
            self.sample_history_vbox,
            self.register_new_processes_title,
            self.register_new_processes_vbox,
            self.process_buttons_hbox,
            self.save_button
        ]
        
        self.select_sample_dropdown.sample_dropdown.observe(
            self.load_sample_data,
            names = "value"
        )
    
    def load_sample_data(self, change):
        if self.select_sample_dropdown.sample_dropdown.value == "-1":
            self.sample_history_vbox.sample_history.children = []
            return
        
        sample_identifier = self.select_sample_dropdown.sample_dropdown.value
        sample_object = utils.get_openbis_object(self.openbis_session, sample_ident = sample_identifier)
        sample_object_parents = sample_object.parents
        most_recent_parent = None
        for parent_id in sample_object_parents:
            parent_object = utils.get_openbis_object(self.openbis_session, sample_ident = parent_id)
            parent_type = parent_object.type
            if parent_type == "PROCESS_STEP":
                if most_recent_parent:
                    if parent_object.registrationDate > most_recent_parent.registrationDate:
                        most_recent_parent = parent_object
                else:
                    most_recent_parent = parent_object
        
        if most_recent_parent:
            if most_recent_parent.experiment.permId != self.select_experiment_dropdown.experiment_dropdown.value:
                self.select_experiment_dropdown.experiment_dropdown.value = most_recent_parent.experiment.permId
                display(utils.Javascript(data = "alert('Experiment was changed!')"))
                
            # Load sample history
            self.sample_history_vbox.load_sample_history(sample_object)
        else:
            self.sample_history_vbox.sample_history.children = []

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
        
        self.experiment_options = self.load_experiments()
        self.experiment_dropdown = ipw.Dropdown(
            options = self.experiment_options,
            layout=ipw.Layout(width='500px'),
            value = "-1"
        )
        
        
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
        self.sort_name_checkbox.observe(self.sort_experiment_dropdown)
        self.sort_registration_date_checkbox.observe(self.sort_experiment_dropdown)
    
    def load_experiments(self):
        experiments = utils.get_openbis_collections(
            self.openbis_session,
            type = "EXPERIMENT"
        )
        experiment_options = []
        for exp in experiments:
            if "$name" in exp.props.all():
                exp_option = (f"{exp.props['$name']} from Project {exp.project.code} and Space {exp.project.space}", exp.permId)
            else:
                exp_option = (f"{exp.code} from Project {exp.project.code} and Space {exp.project.space}", exp.permId)
            experiment_options.append(exp_option)
        experiment_options.insert(0, ("Select experiment...", "-1"))
        return experiment_options

    def sort_experiment_dropdown(self, change):
        options = self.experiment_options[1:]
        
        df = pd.DataFrame(options, columns=["$name", "registration_date"])
        if self.sort_name_checkbox.value and not self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="$name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="registration_date", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by=["$name", "registration_date"], ascending=[True, False])

        options = list(df.itertuples(index=False, name=None))
        options.insert(0, self.experiment_options[0])
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
                "$name": new_experiment_name_textbox.value,
                "default_collection_view": "IMAGING_GALLERY_VIEW"
            }
            
            project_id = project_dropdown.value
            
            if project_id == "-1":
                display(utils.Javascript(data = "alert('Select a project.')"))
                return
            else:
                try:
                    utils.create_experiment_in_openbis(
                        self.openbis_session, 
                        project_dropdown.value, 
                        experiment_props
                    )
                    self.create_new_experiment_widgets.children = []
                    self.experiment_dropdown.options = self.load_experiments()
                    display(utils.Javascript(data = "alert('Experiment successfully created!')"))
                except ValueError as e:
                    display(utils.Javascript(data = "alert('Error! Check if experiment already exists (either in ELN or Trash).')"))
        
        def cancel_new_experiment(b):
            self.create_new_experiment_widgets.children = []
            
        def sort_project_dropdown(change):
            options = project_dropdown_options[1:]
            
            df = pd.DataFrame(options, columns=["$name", "registration_date"])
            if sort_project_name_checkbox.value and not sort_project_registration_date_checkbox.value:
                df = df.sort_values(by="$name", ascending=True)
            elif not sort_project_name_checkbox.value and sort_project_registration_date_checkbox.value:
                df = df.sort_values(by="registration_date", ascending=False)
            elif sort_project_name_checkbox.value and sort_project_registration_date_checkbox.value:
                df = df.sort_values(by=["$name", "registration_date"], ascending=[True, False])

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
        
        self.children = [
            self.sample_dropdown_hbox,
            self.sort_sample_hbox
        ]
    
    def load_samples(self):
        samples = utils.get_openbis_objects(
            self.openbis_session,
            type = "SAMPLE"
        )
        sample_options = [(f"{obj.props['$name']}", obj.permId) for obj in samples if obj.props["exists"] == "true"]
        sample_options.insert(0, ("Select sample...", "-1"))
        self.sample_dropdown.options =  sample_options
        self.sample_dropdown.value = "-1"

class SampleHistoryWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
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
                if "PRST" in parent_code or "SAMP" in parent_code:
                    parent_object = utils.get_openbis_object(self.openbis_session, sample_ident=parent)
                    next_parents.extend(parent_object.parents)
                    if "PRST" in parent_code:
                        process_steps.append(parent_object)
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
        if openbis_object_props["$name"]:
            self.name_html.value = openbis_object_props["$name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
            
        self.registration_date = self.openbis_object.registrationDate
        
        instrument_id = openbis_object_props["instrument"]
        
        if instrument_id:
            instrument_object = utils.get_openbis_object(self.openbis_session, sample_ident = instrument_id)
            self.instrument_html.value = instrument_object.props["$name"]
        
        self.load_actions()
        self.load_observables()
    
    def load_actions(self):
        actions_ids = self.openbis_object.props["actions"]
        if actions_ids:
            actions_accordion_children = []
            for i, act_id in enumerate(actions_ids):
                act_object = utils.get_openbis_object(self.openbis_session, sample_ident = act_id)
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
                obs_object = utils.get_openbis_object(self.openbis_session, sample_ident = obs_id)
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
        self.object_type = self.openbis_object.type
        
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
        
        if self.object_type == "ANNEALING":
            self.target_temperature_label = ipw.Label(value = "Name:")
            self.target_temperature_html = ipw.HTML()
            self.target_temperature_hbox = ipw.HBox(children = [self.target_temperature_label, self.target_temperature_html])
            widget_children.append(self.target_temperature_hbox)
            
        elif self.object_type == "COOLDOWN":
            self.target_temperature_label = ipw.Label(value = "Target temperature:")
            self.target_temperature_html = ipw.HTML()
            self.target_temperature_hbox = ipw.HBox(children = [self.target_temperature_label, self.target_temperature_html])
            
            self.cryogen_label = ipw.Label(value = "Cryogen:")
            self.cryogen_html = ipw.HTML()
            self.cryogen_hbox = ipw.HBox(children = [self.cryogen_label, self.cryogen_html])
            
            widget_children.append(self.cryogen_hbox)
            widget_children.append(self.target_temperature_hbox)
        
        elif self.object_type == "DEPOSITION":
            self.substrate_temperature_label = ipw.Label(value = "Substrate temperature:")
            self.substrate_temperature_html = ipw.HTML()
            self.substrate_temperature_hbox = ipw.HBox(children = [self.substrate_temperature_label, self.substrate_temperature_html])
            
            self.substance_label = ipw.Label(value = "Substance:")
            self.substance_html = ipw.HTML()
            self.substance_hbox = ipw.HBox(children = [self.substance_label, self.substance_html])
            
            widget_children.append(self.substance_hbox)
            widget_children.append(self.substrate_temperature_hbox)
        
        elif self.object_type == "DOSING":
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

        elif self.object_type == "SPUTTERING":
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
        if openbis_object_props["$name"]:
            self.name_html.value = openbis_object_props["$name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["duration"]:
            self.duration_html.value = openbis_object_props["duration"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
        
        if "target_temperature" in openbis_object_props:
            if openbis_object_props["target_temperature"]:
                self.target_temperature_html.value = openbis_object_props["target_temperature"]
        
        if "cryogen" in openbis_object_props:
            if openbis_object_props["cryogen"]:
                self.cryogen_html.value = openbis_object_props["cryogen"]
        
        if "substrate_temperature" in openbis_object_props:
            if openbis_object_props["substrate_temperature"]:
                self.substrate_temperature_html.value = openbis_object_props["substrate_temperature"]
        
        if "dosing_gas" in openbis_object_props:
            if openbis_object_props["dosing_gas"]:
                self.dosing_gas_html.value = openbis_object_props["dosing_gas"]
        
        if "pressure" in openbis_object_props:
            if openbis_object_props["pressure"]:
                self.pressure_html.value = openbis_object_props["pressure"]
        
        if "current" in openbis_object_props:
            if openbis_object_props["current"]:
                self.current_html.value = openbis_object_props["current"]
        
        if "angle" in openbis_object_props:
            if openbis_object_props["angle"]:
                self.angle_html.value = openbis_object_props["angle"]

        if "substance" in openbis_object_props:
            substance_id = openbis_object_props["substance"]
            if substance_id:
                substance_object = utils.get_openbis_object(self.openbis_session, sample_ident = substance_id)
                self.substance_html.value = substance_object.props["$name"]
        
        component_id = openbis_object_props["component"]
        if component_id:
            component_object = utils.get_openbis_object(self.openbis_session, sample_ident = component_id)
            self.component_html.value = component_object.props["$name"]
            
        if openbis_object_props["component_settings"]:
            component_settings = json.loads(openbis_object_props["component_settings"])
            component_settings_string = ""
            for prop_key, prop_value in component_settings.items():
                prop_type = self.openbis_session.get_property_type(prop_key)
                prop_label = prop_type.label
                component_settings_string += f"<p>&bull; {prop_label}: {prop_value}</p>"
            
            self.component_settings_html.value = component_settings_string
        
class ObservableHistoryWidget(ipw.VBox):
    def __init__(self, openbis_session, openbis_object):
        super().__init__()
        self.openbis_session = openbis_session
        self.openbis_object = openbis_object
        
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
        if openbis_object_props["$name"]:
            self.name_html.value = openbis_object_props["$name"]
        
        if openbis_object_props["description"]:
            self.description_html.value = openbis_object_props["description"]
        
        if openbis_object_props["channel_name"]:
            self.ch_name_html.value = openbis_object_props["channel_name"]
        
        if openbis_object_props["comments"]:
            self.comments_html.value = openbis_object_props["comments"]
            
        component_id = openbis_object_props["component"]
            
        if component_id:
            component_object = utils.get_openbis_object(self.openbis_session, sample_ident = component_id)
            self.component_html.value = component_object.props["$name"]
            
        if openbis_object_props["component_settings"]:
            component_settings = json.loads(openbis_object_props["component_settings"])
            component_settings_string = ""
            for prop_key, prop_value in component_settings.items():
                prop_type = self.openbis_session.get_property_type(prop_key)
                prop_label = prop_type.label
                component_settings_string += f"<p>&bull; {prop_label}: {prop_value}</p>"
            
            self.component_settings_html.value = component_settings_string

class RegisterProcessStepWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session