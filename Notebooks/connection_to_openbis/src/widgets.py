import ipywidgets as ipw
from datetime import datetime

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
                layout=ipw.Layout(margin='0px', width='30px')
            )
            
            name_label = ipw.Label(
                value = "Name", 
                layout=ipw.Layout(margin='0px', width='50px'),
                style = {'description_width': 'initial'}
            )
            
            registration_date_checkbox = ipw.Checkbox(
                indent = False,
                layout=ipw.Layout(margin='0px', width='20px')
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
                    name_label,
                    name_checkbox,
                    registration_date_label,
                    registration_date_checkbox
                ]
            )
            
            material_details_html = ipw.HTML()
            
            self.material_details_vbox.children = [
                select_material_box,
                material_details_html
            ]
            
            material_type = self.material_type_dropdown.value
            material_objects = self.openbis_session.get_objects(
                type = material_type
            )
            materials_objects_names_permids = [(obj.props["name"], obj.permId) for obj in material_objects]
            material_options += materials_objects_names_permids
            material_dropdown.options = material_options
            
            def sort_material_dropdown(change):
                options = material_options[1:]
                
                if name_checkbox.value and registration_date_checkbox.value == False:
                    options = sorted(options, key=lambda x: x[0].lower())
                elif name_checkbox.value == False and registration_date_checkbox.value:
                    options = sorted(options, key=lambda x: x[1].lower())
                elif name_checkbox.value and registration_date_checkbox.value:
                    options = sorted(options, key=lambda x: (x[0].lower(), x[1].lower()))

                options.insert(0, material_options[0])
                
                material_dropdown.options = options
            
            def load_material_details(change):
                obj_permid = material_dropdown.value
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
                
                current_datetime = datetime.now()
                current_datetime_str = current_datetime.strftime("%Y%m%d%H%M%S")
                self.sample_name_textbox.value = f"{current_datetime_str}_{obj_name}"
                    
            name_checkbox.observe(sort_material_dropdown, names = "value")
            registration_date_checkbox.observe(sort_material_dropdown, names = "value")
            material_dropdown.observe(load_material_details, names = "value")

class RegisterProcessWidget(ipw.VBox):
    def __init__(self):
        super().__init__()