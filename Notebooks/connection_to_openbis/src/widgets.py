import ipywidgets as ipw
import utils
from IPython.display import display
import pandas as pd

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
            materials_objects_names_permids = [(obj.props["name"], obj.permId) for obj in material_objects]
            material_options += materials_objects_names_permids
            material_dropdown.options = material_options
            
            def sort_material_dropdown(change):
                options = material_options[1:]
                
                df = pd.DataFrame(options, columns=["name", "registration_date"])
                if name_checkbox.value and not registration_date_checkbox.value:
                    df = df.sort_values(by="name", ascending=True)
                elif not name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by="registration_date", ascending=False)
                elif name_checkbox.value and registration_date_checkbox.value:
                    df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

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
                    obj_name = obj_props.get("name", "")
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
                    "name": sample_name,
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
            return
        
        sample_identifier = self.select_sample_dropdown.sample_dropdown.value
        sample_object = utils.get_openbis_object(self.openbis_session, sample_ident = sample_identifier)
        sample_object_parents = sample_object.parents
        most_recent_parent = None
        for parent_id in sample_object_parents:
            parent_object = utils.get_openbis_object(self.openbis_session, sample_ident = parent_id)
            if most_recent_parent:
                if parent_object.registrationDate > most_recent_parent.registrationDate:
                    most_recent_parent = parent_object
            else:
                most_recent_parent = parent_object
        
        self.select_experiment_dropdown.experiment_dropdown.value = most_recent_parent.experiment.permId
        display(utils.Javascript(data = "alert('Experiment was changed!')"))
        
        # TODO: Get all sample history and populate the widget. Then, get the last experiment where the sample was used
        self.sample_history_vbox.load_sample_history(sample_object)

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
        experiment_options = [(f"{exp.props['name']} from Project {exp.project.code} and Space {exp.project.space}", exp.permId) for exp in experiments]
        experiment_options.insert(0, ("Select experiment...", "-1"))
        return experiment_options

    def sort_experiment_dropdown(self, change):
        options = self.experiment_options[1:]
        
        df = pd.DataFrame(options, columns=["name", "registration_date"])
        if self.sort_name_checkbox.value and not self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="name", ascending=True)
        elif not self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by="registration_date", ascending=False)
        elif self.sort_name_checkbox.value and self.sort_registration_date_checkbox.value:
            df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

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
                "name": new_experiment_name_textbox.value,
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
            
            df = pd.DataFrame(options, columns=["name", "registration_date"])
            if sort_project_name_checkbox.value and not sort_project_registration_date_checkbox.value:
                df = df.sort_values(by="name", ascending=True)
            elif not sort_project_name_checkbox.value and sort_project_registration_date_checkbox.value:
                df = df.sort_values(by="registration_date", ascending=False)
            elif sort_project_name_checkbox.value and sort_project_registration_date_checkbox.value:
                df = df.sort_values(by=["name", "registration_date"], ascending=[True, False])

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
        
        self.sample_options = self.load_samples()
        self.sample_dropdown = ipw.Dropdown(
            options = self.sample_options,
            value = "-1"
        )
        
        
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
        sample_options = [(f"{obj.props['name']}", obj.permId) for obj in samples if obj.props["exists"] == "true"]
        sample_options.insert(0, ("Select sample...", "-1"))
        return sample_options

class SampleHistoryWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session
        
        self.sample_history = ipw.Accordion()
        
        self.children = [
            self.sample_history
        ]
    
    def load_sample_history(self, sample_object):
        pass

class RegisterProcessStepWidget(ipw.VBox):
    def __init__(self, openbis_session):
        super().__init__()
        self.openbis_session = openbis_session