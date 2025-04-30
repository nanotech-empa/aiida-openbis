import ipywidgets as ipw
from IPython.display import display, clear_output
import utils
import io
import os
import urllib.parse
import rdkit
from rdkit.Chem import AllChem, Draw, rdMolDescriptors

class SummaryObjectWidget(ipw.VBox):
    def __init__(self):
        # Initialize the parent
        super().__init__()
        
        self.names_label = utils.Label(
            value = "Names: "
        )
        
        self.names_textbox = utils.Text(
            layout = ipw.Layout(width = '30%'),
            placeholder = "Write names separated by ;"
        )
        
        self.names_widgets = utils.HBox(
            [
                self.names_label,
                self.names_textbox
            ]
        )
        
        self.structure_file_uploader_label = utils.Label(
            value = "Upload structure:"
        )
        
        self.structure_file_uploader = utils.FileUpload(
            accept = ".cdxml",
            multiple = False
        )
        
        self.structure_file_uploader_widgets = utils.HBox(
            [
                self.structure_file_uploader_label,
                self.structure_file_uploader
            ]
        )
        
        self.structure_sketch_label = utils.Label(
            value = "Structure sketch: "
        )
        
        self.structure_sketch_image_file = "images/white_screen.jpg"
        
        self.structure_sketch_imagebox = utils.Image(
            value = utils.read_file(self.structure_sketch_image_file),
            layout = ipw.Layout(border = 'solid 1px #cccccc', width = '250px', height = '250px')
        )
        
        self.structure_sketch_widgets = ipw.HBox(
            [
                self.structure_sketch_label,
                self.structure_sketch_imagebox
            ]
        )
        
        self.structure_file_label = utils.Label(
            value = "Structure file:"
        )
        
        self.structure_file_hyperlink = utils.HTML(
            value = ""
        )
        
        self.structure_cdxml_file = ""
        
        self.structure_file_widgets = utils.HBox(
            [
                self.structure_file_label,
                self.structure_file_hyperlink
            ]
        )
        
        self.classification_label = utils.Label(
            value = "Classification: "
        )
        
        self.classification_textbox = utils.Text(
            layout = ipw.Layout(width = '30%')
        )
        
        self.classification_widgets = utils.HBox(
            [
                self.classification_label,
                self.classification_textbox
            ]
        )
        
        self.structure_file_uploader.observe(
            self.load_sketch_cdxml
        )
        
        self.children = [
            self.names_widgets,
            self.structure_file_uploader_widgets,
            self.structure_sketch_widgets,
            self.structure_file_widgets,
            self.classification_widgets
        ]
    
    def load_sketch_cdxml(self, change):
        for filename, data in self.structure_file_uploader.value.items():
            content = data["content"]
            cdxml_filepath = f"structures/{filename}"
            self.structure_cdxml_file = cdxml_filepath
            utils.save_file(content, cdxml_filepath)
            encoded_path = urllib.parse.quote(cdxml_filepath)
            self.structure_file_hyperlink.value = f"<a href='{encoded_path}' target='_blank'>Download structure cdxml file</a>"
            
            # Convert CDXML to PNG
            molecules = rdkit.Chem.MolsFromCDXML(content)
            generate_molecule_image = False
            
            if len(molecules) == 1:
                mol = molecules[0] # Get first molecule
                mol_smiles = rdkit.Chem.MolToSmiles(mol) # Canonical Smiles
                chem_mol = rdkit.Chem.MolFromSmiles(mol_smiles)
                
                if chem_mol is not None:
                    generate_molecule_image = True
                    AllChem.Compute2DCoords(chem_mol) # Add coords to the atoms in the molecule
                    image = Draw.MolToImage(chem_mol)
                    buffer = io.BytesIO()
                    image.save(buffer, format='PNG')
                    image_bytes = buffer.getvalue()
                    image_filepath = "images/structure.png"
                    image.save(image_filepath)
                    self.structure_sketch_imagebox.value = image_bytes
                    self.structure_sketch_image_file = image_filepath
                    
            if generate_molecule_image == False:
                self.structure_sketch_imagebox.value = utils.read_file("images/white_screen.jpg")
                print(f"Cannot generate molecule image.")

class SimulationObjectWidget(ipw.VBox):
    def __init__(self, accordion_parent_widget, accordion_parent_item_index):
        # Initialize the parent
        super().__init__()
        
        # Make the connection with the accordion where the Simulations are stored
        self.accordion_parent_item_index = accordion_parent_item_index
        self.accordion_parent_widget = accordion_parent_widget
        
        self.substrate_label = utils.Label(
            value = "Substrate: "
        )
        
        self.substrate_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.substrate_widgets = utils.HBox(
            [
                self.substrate_label,
                self.substrate_textbox
            ]
        )
        
        self.uks_multi_label = utils.Label(
            value = "UKS Multi.: "
        )
        
        self.uks_multi_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.uks_multi_widgets = utils.HBox(
            [
                self.uks_multi_label,
                self.uks_multi_textbox
            ]
        )
        
        self.charge_label = utils.Label(
            value = "Charge: "
        )
        
        self.charge_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.charge_widgets = utils.HBox(
            [
                self.charge_label,
                self.charge_textbox
            ]
        )
        
        self.num_alpha_label = utils.Label(
            value = "Number of alpha: "
        )
        
        self.num_alpha_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.num_alpha_widgets = utils.HBox(
            [
                self.num_alpha_label,
                self.num_alpha_textbox
            ]
        )
        
        self.num_beta_label = utils.Label(
            value = "Number of beta: "
        )
        
        self.num_beta_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.num_beta_widgets = utils.HBox(
            [
                self.num_beta_label,
                self.num_beta_textbox
            ]
        )
        
        self.energy_au_label = utils.Label(
            value = "Energy [au]: "
        )
        
        self.energy_au_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.energy_au_widgets = utils.HBox(
            [
                self.energy_au_label,
                self.energy_au_textbox
            ]
        )
        
        self.energy_ev_label = utils.Label(
            value = "Energy [eV]: "
        )
        
        self.energy_ev_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.energy_ev_widgets = utils.HBox(
            [
                self.energy_ev_label,
                self.energy_ev_textbox
            ]
        )
        
        self.gap_ev_label = utils.Label(
            value = "Gap [eV]: "
        )
        
        self.gap_ev_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.gap_ev_widgets = utils.HBox(
            [
                self.gap_ev_label,
                self.gap_ev_textbox
            ]
        )
        
        self.band_energies_label = utils.Label(
            value = "Band energies (wrt E_F)"
        )
        
        self.band_energies_items = [
            utils.Label(
                value = "Index"
            ),
            utils.Label(
                value = "Label"
            ),
            utils.Label(
                value = "E - E_F (eV)"
            ),
            utils.Label(
                value = "Effective mass (m_e)"
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            )
        ]
        
        self.band_energies_gridbox = utils.GridBox(
            self.band_energies_items,
            layout = ipw.Layout(grid_template_columns = "repeat(4, 15%)")
        )
        
        self.add_row_band_energies_button = utils.Button(
            button_style = "success",
            icon = "plus",
            layout = ipw.Layout(width = "5%")
        )
        
        self.remove_row_band_energies_button = utils.Button(
            button_style = "danger",
            icon = "minus",
            layout = ipw.Layout(width = "5%")
        )
        
        self.change_rows_band_energies_buttons = utils.HBox(
            [
                self.add_row_band_energies_button,
                self.remove_row_band_energies_button
            ]
        )
        
        self.ldos_orbital_maps_label = utils.Label(
            value = "LDOS/Orbital Maps"
        )
        
        self.ldos_orbital_maps_imagebox = utils.Image(
            value = utils.read_file("images/white_screen.jpg"),
            layout = ipw.Layout(border = 'solid 1px #cccccc', width = '250px', height = '250px')
        )
        
        self.ldos_orbital_maps_image_uploader = utils.FileUpload(
            multiple = True
        )
        
        self.ldos_orbital_maps_image_output = utils.Output()
        
        self.ldos_orbital_maps_widgets = utils.VBox([])
        
        self.ldos_orbital_maps_files = []
        
        self.add_row_band_energies_button.on_click(
            self.add_band_energies_row
        )
        
        self.remove_row_band_energies_button.on_click(
            self.remove_band_energies_row
        )
        
        self.remove_simulation_button = utils.Button(
            button_style = "danger",
            description = "Remove simulation"
        )
        
        self.remove_simulation_button.on_click(
            self.remove_simulation
        )
        
        self.ldos_orbital_maps_image_uploader.observe(
            self.load_ldos_orbital_maps_image, 
            names = '_counter'
        )
        
        self.children = [
            self.substrate_widgets,
            self.uks_multi_widgets,
            self.charge_widgets,
            self.num_alpha_widgets,
            self.num_beta_widgets,
            self.energy_au_widgets,
            self.energy_ev_widgets,
            self.gap_ev_widgets,
            self.band_energies_label,
            self.band_energies_gridbox,
            self.change_rows_band_energies_buttons,
            self.ldos_orbital_maps_label,
            self.ldos_orbital_maps_image_uploader,
            self.ldos_orbital_maps_image_output,
            self.remove_simulation_button
        ]

    def load_ldos_orbital_maps_image(self, change):
        with self.ldos_orbital_maps_image_output:
            clear_output()
            images_content = []
            for image_filename, image_data in self.ldos_orbital_maps_image_uploader.value.items():
                image_binary = image_data["content"]
                image_widget = utils.Image(
                    value = image_binary,
                    layout = ipw.Layout(border = 'solid 1px #cccccc', width = '250px', height = '250px')
                )
                
                image_caption = utils.Text(
                    description = "Caption: "
                )
                
                image_widgets = utils.VBox(
                    [
                        image_widget,
                        image_caption
                    ]
                )
                
                image_filepath = f"images/{image_filename}"
                utils.save_file(image_binary, image_filepath)
                images_content.append(image_widgets)
                self.ldos_orbital_maps_files.append(image_filepath)
            
            self.ldos_orbital_maps_widgets.children = images_content
            display(self.ldos_orbital_maps_widgets)
        
    def add_band_energies_row(self, b):
        band_energies_grid_box_items = list(self.band_energies_gridbox.children)
        band_energies_grid_box_items.extend(
            [
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                ),
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                ),
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                ),
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                )
            ]
        )
        self.band_energies_gridbox.children = band_energies_grid_box_items
    
    def remove_band_energies_row(self, b):
        band_energies_grid_box_items = list(self.band_energies_gridbox.children)
        if len(band_energies_grid_box_items) > 8:
            self.band_energies_gridbox.children = band_energies_grid_box_items[:-4]
       
    def remove_simulation(self, b):
        simulations_accordion_titles = self.accordion_parent_widget._titles
        simulations_accordion_children = self.accordion_parent_widget.children
        simulations_accordion_children = list(simulations_accordion_children)
        simulations_accordion_children.pop(self.accordion_parent_item_index)
        titles_simulation_key = str(self.accordion_parent_item_index)
        del simulations_accordion_titles[titles_simulation_key]
        self.accordion_parent_widget.children = simulations_accordion_children
        
        # Update parent accordion indexes
        simulations_accordion_titles = list(simulations_accordion_titles.values())
        for idx, simulation_widget in enumerate(simulations_accordion_children):
            if simulation_widget.accordion_parent_item_index > self.accordion_parent_item_index:
                simulation_widget.accordion_parent_item_index -= 1
            self.accordion_parent_widget.set_title(idx, simulations_accordion_titles[idx])
            
        # Remove the last element title
        new_num_simulations = len(simulations_accordion_children)
        self.accordion_parent_widget.set_title(new_num_simulations, "")

class ExperimentObjectWidget(ipw.VBox):
    def __init__(self, accordion_parent_widget, accordion_parent_item_index):
        # Initialize the parent
        super().__init__()
        
        # Make the connection with the accordion where the Experiments are stored
        self.accordion_parent_item_index = accordion_parent_item_index
        self.accordion_parent_widget = accordion_parent_widget
        
        self.doi_label = utils.Label(
            value = "DOI: "
        )
        
        self.doi_textbox = utils.Text(
            layout = ipw.Layout(width = '40%')
        )
        
        self.doi_widgets = utils.HBox(
            [
                self.doi_label,
                self.doi_textbox
            ]
        )
        
        self.year_label = utils.Label(
            value = "Year: "
        )
        
        self.year_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.year_widgets = utils.HBox(
            [
                self.year_label,
                self.year_textbox
            ]
        )
        
        self.comments_label = utils.Label(
            value = "Comments/flags: "
        )
        
        self.comments_textareabox = utils.Textarea(
            layout = ipw.Layout(width = '20%')
        )
        
        self.comments_widgets = utils.HBox(
            [
                self.comments_label,
                self.comments_textareabox
            ]
        )
        
        self.substrate_label = utils.Label(
            value = "Substrate: "
        )
        
        self.substrate_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.substrate_widgets = utils.HBox(
            [
                self.substrate_label,
                self.substrate_textbox
            ]
        )
        
        self.spin_excitation_energy_ev_label = utils.Label(
            value = "Spin excitation energy [eV]: "
        )
        
        self.spin_excitation_energy_ev_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.spin_excitation_energy_ev_widgets = utils.HBox(
            [
                self.spin_excitation_energy_ev_label,
                self.spin_excitation_energy_ev_textbox
            ]
        )
        
        self.gap_ev_label = utils.Label(
            value = "Gap [eV]: "
        )
        
        self.gap_ev_textbox = utils.Text(
            layout = ipw.Layout(width = '20%')
        )
        
        self.gap_ev_widgets = utils.HBox(
            [
                self.gap_ev_label,
                self.gap_ev_textbox
            ]
        )
        
        self.band_energies_label = utils.Label(
            value = "Band energies"
        )
        
        self.bands_energies_items = [
            utils.Label(
                value = "Index"
            ),
            utils.Label(
                value = "Label"
            ),
            utils.Label(
                value = "E - E_F (eV)"
            ),
            utils.Label(
                value = "Error (eV)"
            ),
            utils.Label(
                value = "Effective mass (m_e)"
            ),
            utils.Label(
                value = "Error (m_e)"
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            ),
            utils.Text(
                layout = ipw.Layout(width = "80%")
            )
        ]
        
        self.band_energies_gridbox = utils.GridBox(
            self.bands_energies_items,
            layout = ipw.Layout(grid_template_columns = "repeat(6, 15%)")
        )
        
        self.add_row_band_energies_button = utils.Button(
            button_style = "success",
            icon = "plus",
            layout = ipw.Layout(width = "5%")
        )
        
        self.remove_row_band_energies_button = utils.Button(
            button_style = "danger",
            icon = "minus",
            layout = ipw.Layout(width = "5%")
        )
        
        self.change_rows_band_energies_buttons = utils.HBox(
            [
                self.add_row_band_energies_button,
                self.remove_row_band_energies_button
            ]
        )
        
        self.add_row_band_energies_button.on_click(
            self.add_band_energies_row
        )
        
        self.remove_row_band_energies_button.on_click(
            self.remove_band_energies_row
        )
        
        self.relevant_stm_label = utils.Label(
            value = "Relevant STMs:"
        )
        
        self.relevant_stm_output = utils.Output()
        
        self.relevant_stm_image_widgets = utils.VBox([])
        
        self.relevant_stm_image_files = []
        
        self.relevant_stm_image_uploader = utils.FileUpload(
            multiple = True
        )
        
        self.relevant_sts_label = utils.Label(
            value = "Relevant STSs:"
        )
        
        self.relevant_sts_output = utils.Output()
        
        self.relevant_sts_image_widgets = utils.VBox([])
        
        self.relevant_sts_image_files = []
        
        self.relevant_sts_image_uploader = utils.FileUpload(
            multiple = True
        )
        
        self.precursor_label = utils.Label(
            value = "Precursor"
        )
        
        self.precursor_image_output = utils.Output()
        
        self.precursor_imagebox = utils.Image(
            value = utils.read_file("images/white_screen.jpg"),
            layout = ipw.Layout(border = 'solid 1px #cccccc', width = '250px', height = '250px')
        )
        
        self.precursor_cdxml_uploader_label = utils.Label(
            value = "Upload precursor CDXML file:"
        )
        
        self.precursor_cdxml_uploader = utils.FileUpload(
            accept = ".cdxml",
            multiple = False
        )
        
        self.precursor_cdxml_uploader_widgets = utils.HBox(
            [
                self.precursor_cdxml_uploader_label,
                self.precursor_cdxml_uploader
            ]
        )
        
        self.precursor_file_hyperlink_label = utils.Label(
            value = "Precursor file:"
        )
        
        self.precursor_file_hyperlink = utils.HTML(
            value = ""
        )
        
        self.precursor_file_hyperlink_widgets = utils.HBox(
            [
                self.precursor_file_hyperlink_label,
                self.precursor_file_hyperlink
            ]
        )
        
        self.precursor_cdxml_file = ""
        self.precursor_image_file = ""
        
        self.precursor_caption_label = utils.Label(
            value = "Precursor caption: "
        )
        
        self.precursor_caption_textbox = utils.Text(
            layout = ipw.Layout(width = "20%")
        )
        
        self.precursor_caption_widgets = utils.HBox(
            [
                self.precursor_caption_label,
                self.precursor_caption_textbox
            ]
        )
        
        self.remove_experiment_button = utils.Button(
            button_style = "danger",
            description = "Remove experiment"
        )
        
        self.remove_experiment_button.on_click(
            self.remove_experiment
        )
        
        self.relevant_stm_image_uploader.observe(
            self.load_relevant_stm,
            names = "_counter"
        )
        
        self.relevant_sts_image_uploader.observe(
            self.load_relevant_sts,
            names = "_counter"
        )
        
        self.precursor_cdxml_uploader.observe(
            self.load_precursor_cdxml,
            names = "_counter"
        )
        
        self.children = [
            self.doi_widgets,
            self.year_widgets,
            self.comments_widgets,
            self.substrate_widgets,
            self.spin_excitation_energy_ev_widgets,
            self.gap_ev_widgets,
            self.band_energies_gridbox,
            self.change_rows_band_energies_buttons,
            self.relevant_stm_label,
            self.relevant_stm_image_uploader,
            self.relevant_stm_image_widgets,
            self.relevant_sts_label,
            self.relevant_sts_image_uploader,
            self.relevant_sts_image_widgets,
            self.precursor_label,
            self.precursor_cdxml_uploader_widgets,
            self.precursor_file_hyperlink_widgets,
            self.precursor_imagebox,
            self.precursor_caption_widgets,
            self.remove_experiment_button
        ]
    
    def load_relevant_stm(self, change):
        images_content = []
        for image_filename, image_data in self.relevant_stm_image_uploader.value.items():
            image_binary = image_data["content"]
            image_widget = utils.Image(
                value = image_binary,
                layout = ipw.Layout(border = 'solid 1px #cccccc', width = '250px', height = '250px')
            )
            
            image_caption = utils.Text(
                description = "Caption: "
            )
            
            image_widgets = utils.VBox(
                [
                    image_widget,
                    image_caption
                ]
            )
            
            image_filepath = f"images/{image_filename}"
            utils.save_file(image_binary, image_filepath)
            images_content.append(image_widgets)
            self.relevant_stm_image_files.append(image_filepath)
        
        self.relevant_stm_image_widgets.children = images_content
    
    def load_relevant_sts(self, change):
        images_content = []
        for image_filename, image_data in self.relevant_sts_image_uploader.value.items():
            image_binary = image_data["content"]
            image_widget = utils.Image(
                value = image_binary,
                layout = ipw.Layout(border = 'solid 1px #cccccc', width = '250px', height = '250px')
            )
            
            image_caption = utils.Text(
                description = "Caption: "
            )
            
            image_widgets = utils.VBox(
                [
                    image_widget,
                    image_caption
                ]
            )
            
            image_filepath = f"images/{image_filename}"
            utils.save_file(image_binary, image_filepath)
            images_content.append(image_widgets)
            self.relevant_sts_image_files.append(image_filepath)
        
        self.relevant_sts_image_widgets.children = images_content
    
    def load_precursor_cdxml(self, change):
        for filename, data in self.precursor_cdxml_uploader.value.items():
            content = data["content"]
            cdxml_filepath = f"structures/{filename}"
            self.precursor_cdxml_file = cdxml_filepath
            utils.save_file(content, cdxml_filepath)
            encoded_path = urllib.parse.quote(cdxml_filepath)
            self.precursor_file_hyperlink.value = f"<a href='{encoded_path}' target='_blank'>Download precursor cdxml file</a>"
            
            # Convert CDXML to PNG
            molecules = rdkit.Chem.MolsFromCDXML(content)
            generate_molecule_image = False
            
            if len(molecules) == 1:
                mol = molecules[0] # Get first molecule
                mol_name = os.path.splitext(filename)[0]
                mol_smiles = rdkit.Chem.MolToSmiles(mol) # Canonical Smiles
                chem_mol = rdkit.Chem.MolFromSmiles(mol_smiles)
                
                if chem_mol is not None:
                    generate_molecule_image = True
                    AllChem.Compute2DCoords(chem_mol) # Add coords to the atoms in the molecule
                    image = Draw.MolToImage(chem_mol)
                    buffer = io.BytesIO()
                    image.save(buffer, format='PNG')
                    image_bytes = buffer.getvalue()
                    image_filepath = f"images/{mol_name}.png"
                    image.save(image_filepath)
                    self.precursor_imagebox.value = image_bytes
                    self.precursor_image_file = image_filepath
                    
            if generate_molecule_image == False:
                self.precursor_imagebox.value = utils.read_file("images/white_screen.jpg")
                print(f"Cannot generate molecule image.")

      
    def add_band_energies_row(self, b):
        band_energies_grid_box_items = list(self.band_energies_gridbox.children)
        band_energies_grid_box_items.extend(
            [
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                ),
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                ),
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                ),
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                ),
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                ),
                utils.Text(
                    layout = ipw.Layout(width = "80%")
                )
            ]
        )
        self.band_energies_gridbox.children = band_energies_grid_box_items
    
    def remove_band_energies_row(self, b):
        band_energies_grid_box_items = list(self.band_energies_gridbox.children)
        if len(band_energies_grid_box_items) > 12:
            self.band_energies_gridbox.children = band_energies_grid_box_items[:-6]

    def remove_experiment(self, b):
        experiments_accordion_titles = self.accordion_parent_widget._titles
        experiments_accordion_children = self.accordion_parent_widget.children
        experiments_accordion_children = list(experiments_accordion_children)
        experiments_accordion_children.pop(self.accordion_parent_item_index)
        titles_experiment_key = str(self.accordion_parent_item_index)
        del experiments_accordion_titles[titles_experiment_key]
        self.accordion_parent_widget.children = experiments_accordion_children
        
        # Update parent accordion indexes
        experiments_accordion_titles = list(experiments_accordion_titles.values())
        for idx, experiment_widget in enumerate(experiments_accordion_children):
            if experiment_widget.accordion_parent_item_index > self.accordion_parent_item_index:
                experiment_widget.accordion_parent_item_index -= 1
            self.accordion_parent_widget.set_title(idx, experiments_accordion_titles[idx])
            
        # Remove the last element title
        new_num_experiments = len(experiments_accordion_children)
        self.accordion_parent_widget.set_title(new_num_experiments, "")

class PublicationObjectWidget(ipw.VBox):
    def __init__(self):
        # Initialize the parent
        super().__init__()
        
        self.summary_info_title_label = ipw.HTML(
            value = '<span style="font-size:20px; font-weight:bold;">Summary of information</span>'
        )
        
        self.summary_info_widgets = SummaryObjectWidget()
        
        self.simulations_title_label = ipw.HTML(
            value = '<span style="font-size:20px; font-weight:bold;">Simulations</span>'
        )

        self.add_simulation_button = utils.Button(
            button_style = "success",
            icon = "plus",
            layout = ipw.Layout(width = "5%")
        )
        
        self.simulations_accordion = utils.Accordion()
        
        self.experiments_title_label = ipw.HTML(
            value = '<span style="font-size:20px; font-weight:bold;">Experiments</span>'
        )
        
        self.add_experiment_button = utils.Button(
            button_style = "success",
            icon = "plus",
            layout = ipw.Layout(width = "5%")
        )
        
        self.experiments_accordion = utils.Accordion()
        
        self.add_simulation_button.on_click(
            self.add_simulation
        )
        
        self.add_experiment_button.on_click(
            self.add_experiment
        )
        
        self.children = [
            self.summary_info_title_label,
            self.summary_info_widgets,
            self.simulations_title_label,
            self.simulations_accordion,
            self.add_simulation_button,
            self.experiments_title_label,
            self.experiments_accordion,
            self.add_experiment_button
        ]
    
    def get_widgets_data(self):
        widgets_data = {}
        
        if self.summary_info_widgets.names_textbox.value:
            names = self.summary_info_widgets.names_textbox.value
            widgets_data["names"] = names.split(";")
        
        if self.summary_info_widgets.structure_sketch_image_file:
            widgets_data["structure_sketch_file"] = self.summary_info_widgets.structure_sketch_image_file
        
        if self.summary_info_widgets.structure_cdxml_file:
            widgets_data["structure_file"] = self.summary_info_widgets.structure_cdxml_file
        
        if self.summary_info_widgets.classification_textbox.value:
            widgets_data["classification"] = self.summary_info_widgets.classification_textbox.value
        
        simulations_data = []
        simulations_objects = self.simulations_accordion.children
        for simulation in simulations_objects:
            simulation_dict = {}
            
            if simulation.substrate_textbox.value:
                simulation_dict["substrate"] = simulation.substrate_textbox.value
            
            if simulation.uks_multi_textbox.value:
                uks_multiplicity = simulation.uks_multi_textbox.value
                if utils.is_integer(uks_multiplicity):
                    simulation_dict["uks_multiplicity"] = int(uks_multiplicity)
            
            if simulation.charge_textbox.value:
                charge = simulation.charge_textbox.value
                if utils.is_integer(charge):
                    simulation_dict["charge"] = int(charge)
            
            if simulation.num_alpha_textbox.value:
                num_alpha = simulation.num_alpha_textbox.value
                if utils.is_integer(num_alpha):
                    simulation_dict["number_alpha_electrons"] = int(num_alpha)
            
            if simulation.num_beta_textbox.value:
                num_beta = simulation.num_beta_textbox.value
                if utils.is_integer(num_beta):
                    simulation_dict["number_beta_electrons"] = int(num_beta)
            
            if simulation.energy_au_textbox.value:
                energy_au = simulation.energy_au_textbox.value
                if utils.is_number(energy_au):
                    simulation_dict["energy_au"] = float(energy_au)
            
            if simulation.energy_ev_textbox.value:
                energy_ev = simulation.energy_ev_textbox.value
                if utils.is_number(energy_ev):
                    simulation_dict["energy_ev"] = float(energy_ev)
            
            if simulation.gap_ev_textbox.value:
                gap_ev = simulation.gap_ev_textbox.value
                if utils.is_number(gap_ev):
                    simulation_dict["gap_ev"] = float(gap_ev)
            
            band_energies_data = []
            band_energies_objects = simulation.band_energies_gridbox.children
            for idx in range(4, len(band_energies_objects), 4):
                band_energy_dict = {}
                if band_energies_objects[idx].value:
                    index = band_energies_objects[idx].value
                    if utils.is_integer(index):
                        band_energy_dict["index"] = int(index)
                
                if band_energies_objects[idx + 1].value:
                    band_energy_dict["label"] = band_energies_objects[idx + 1].value
                
                if band_energies_objects[idx + 2].value:
                    e_ef_ev = band_energies_objects[idx + 2].value
                    if utils.is_number(e_ef_ev):
                        band_energy_dict["e_ef_ev"] = float(e_ef_ev)
                
                if band_energies_objects[idx + 3].value:
                    effective_mass_m_e = band_energies_objects[idx + 3].value
                    if utils.is_number(effective_mass_m_e):
                        band_energy_dict["effective_mass_m_e"] = float(effective_mass_m_e)
                
                band_energies_data.append(band_energy_dict)

            if band_energies_data:
                simulation_dict["band_energies"] = band_energies_data
            
            ldos_orbital_maps = []
            ldos_orbital_maps_objects = simulation.ldos_orbital_maps_widgets.children
            for idx, ldos_orbital_map in enumerate(ldos_orbital_maps_objects):
                ldos_orbital_map_dict = {}
                ldos_orbital_map_dict["file"] = simulation.ldos_orbital_maps_files[idx]
                
                if ldos_orbital_map.children[1].value:
                    ldos_orbital_map_dict["caption"] = ldos_orbital_map.children[1].value
                
                ldos_orbital_maps.append(ldos_orbital_map_dict)
            
            if ldos_orbital_maps:
                simulation_dict["ldos_orbital_maps"] = ldos_orbital_maps
            
            if simulation_dict:
                simulations_data.append(simulation_dict)
        
        if simulations_data:
            widgets_data["simulations"] = simulations_data
        
        # Experiment data
        experiments_data = []
        experiments_objects = self.experiments_accordion.children
        
        for experiment in experiments_objects:
            experiment_dict = {}
            
            if experiment.doi_textbox.value:
                experiment_dict["doi"] = experiment.doi_textbox.value
                
            if experiment.year_textbox.value:
                year = experiment.year_textbox.value
                if utils.is_integer(year):
                    experiment_dict["year"] = int(year)
            
            if experiment.comments_textareabox.value:
                experiment_dict["comments_flags"] = experiment.comments_textareabox.value
            
            if experiment.substrate_textbox.value:
                experiment_dict["substrate"] = experiment.substrate_textbox.value
            
            if experiment.spin_excitation_energy_ev_textbox.value:
                spin_excitation_energy_ev = experiment.spin_excitation_energy_ev_textbox.value
                if utils.is_number(spin_excitation_energy_ev):
                    experiment_dict["spin_excitation_energy_ev"] = float(spin_excitation_energy_ev)
            
            if experiment.gap_ev_textbox.value:
                gap_ev = experiment.gap_ev_textbox.value
                if utils.is_number(gap_ev):
                    experiment_dict["gap_ev"] = float(gap_ev)
            
            band_energies_data = []
            band_energies_objects = experiment.band_energies_gridbox.children
            for idx in range(6, len(band_energies_objects), 6):
                band_energy_dict = {}
                if band_energies_objects[idx].value:
                    index = band_energies_objects[idx].value
                    if utils.is_integer(index):
                        band_energy_dict["index"] = int(index)
                
                if band_energies_objects[idx + 1].value:
                    band_energy_dict["label"] = band_energies_objects[idx + 1].value
                
                if band_energies_objects[idx + 2].value:
                    e_ef_ev = band_energies_objects[idx + 2].value
                    if utils.is_number(e_ef_ev):
                        band_energy_dict["e_ef_ev"] = float(e_ef_ev)
                
                if band_energies_objects[idx + 3].value:
                    error_ev = band_energies_objects[idx + 3].value
                    if utils.is_number(error_ev):
                        band_energy_dict["error_ev"] = float(error_ev)
                
                if band_energies_objects[idx + 4].value:
                    effective_mass_m_e = band_energies_objects[idx + 4].value
                    if utils.is_number(effective_mass_m_e):
                        band_energy_dict["effective_mass_m_e"] = float(effective_mass_m_e)
                
                if band_energies_objects[idx + 5].value:
                    error_m_e = band_energies_objects[idx + 5].value
                    if utils.is_number(error_m_e):
                        band_energy_dict["error_m_e"] = float(error_m_e)
                
                band_energies_data.append(band_energy_dict)

            if band_energies_data:
                experiment_dict["band_energies"] = band_energies_data
            
            relevant_stm_images = []
            stm_images_objects = experiment.relevant_stm_image_widgets.children
            for idx, stm_image in enumerate(stm_images_objects):
                stm_image_dict = {}
                stm_image_dict["file"] = experiment.relevant_stm_image_files[idx]
                
                if stm_image.children[1].value:
                    stm_image_dict["caption"] = stm_image.children[1].value
                
                relevant_stm_images.append(stm_image_dict)

            relevant_sts_images = []
            sts_images_objects = experiment.relevant_sts_image_widgets.children
            for idx, sts_image in enumerate(sts_images_objects):
                sts_image_dict = {}
                sts_image_dict["file"] = experiment.relevant_sts_image_files[idx]
                
                if sts_image.children[1].value:
                    sts_image_dict["caption"] = sts_image.children[1].value
                
                relevant_sts_images.append(sts_image_dict)
            
            precursor_dict = {}
            if experiment.precursor_image_file:
                precursor_dict["image_file"] = experiment.precursor_image_file
            
            if experiment.precursor_cdxml_file:
                precursor_dict["structure_file"] = experiment.precursor_cdxml_file
            
            if experiment.precursor_caption_textbox.value:
                precursor_dict["caption"] = experiment.precursor_caption_textbox.value
            
            if relevant_stm_images:
                experiment_dict["relevant_stm_images"] = relevant_stm_images
            
            if relevant_sts_images:
                experiment_dict["relevant_sts_images"] = relevant_sts_images
            
            if precursor_dict:
                experiment_dict["precursor"] = precursor_dict
            
            if experiment_dict:
                experiments_data.append(experiment_dict)
        
        if experiments_data:
            widgets_data["experiments"] = experiments_data
        
        return widgets_data
    
    def add_simulation(self, b):
        simulations_accordion_children = list(self.simulations_accordion.children)
        num_simulations = len(simulations_accordion_children)
        simulation_index = num_simulations
        simulation_object_widget = SimulationObjectWidget(self.simulations_accordion, simulation_index)
        simulations_accordion_children.append(simulation_object_widget)
        self.simulations_accordion.set_title(num_simulations, f"Simulation {num_simulations + 1}")
        self.simulations_accordion.children = simulations_accordion_children
    
    def add_experiment(self, b):
        experiments_accordion_children = list(self.experiments_accordion.children)
        num_experiments = len(experiments_accordion_children)
        experiment_index = num_experiments
        experiment_object_widget = ExperimentObjectWidget(self.experiments_accordion, experiment_index)
        experiments_accordion_children.append(experiment_object_widget)
        self.experiments_accordion.set_title(experiment_index, f"Experiment {num_experiments + 1}")
        self.experiments_accordion.children = experiments_accordion_children