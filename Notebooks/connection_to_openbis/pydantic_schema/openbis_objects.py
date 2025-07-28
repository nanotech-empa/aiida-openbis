from pydantic import BaseModel, Field, model_validator, field_validator, conint
from typing import List, Union, Dict, Literal, Optional
from enum import Enum

# Enums

class ComponentMainCategoryEnum(str, Enum):
    Auxiliary = "Auxiliary"
    Microscope_Core_Components = "Microscope Core Components"
    Sample_Preparation_Handling = "Sample Preparation & Handling"
    Vacuum_System = "Vacuum System"

class ComponentSubCategoryEnum(str, Enum):
    Auxiliary = "Auxiliary"
    Control_Data_Acquisition = "Control and Data Acquisition"
    Gauge = "Gauge"
    Mechanical_Component = "Mechanical Component"
    Preparation_Tools = "Preparation Tools"
    Pump = "Pump"
    STM = "STM"
    Temperature_Environment_Control = "Temperature and Environment Control"
    Vacuum_Chamber = "Vacuum Chamber"

class DraftTypeEnum(str, Enum):
    Preprint = "Preprint"
    Postprint = "Postprint"

class ObjectStatusEnum(str, Enum):
    Active = "Active"
    Broken = "Broken"
    Disposed = "Disposed"
    Inactive = "Inactive"

class AngleUnitEnum(str, Enum):
    Degree = "deg"
    Radian = "rad"

class CurrentUnitEnum(str, Enum):
    Ampere = "A"
    Milliampere = "mA"
    Microampere = "uA"

class DensityUnitEnum(str, Enum):
    Grams_per_Cubic_Centimiter = "g/cm3"

class EnergyUnitEnum(str, Enum):
    Joule = "J"
    ElectronVolt = "eV"

class MassUnitEnum(str, Enum):
    Gram = "g"
    Milligram = "mg"
    Microgram = "ug"

class TemperatureUnitEnum(str, Enum):
    Celsius = "C"
    Kelvin = "K"
    Fahrenheit = "F"

class PressureUnitEnum(str, Enum):
    Millibar = "mBar"
    Bar = "Bar"
    Pascal = "Pa"
    KiloPascal = "kPa"

class VoltageUnitEnum(str, Enum):
    Kilovolt = "kV"
    Volt = "V"
    Millivolt = "mV"
    Microvolt = "uV"

class VolumeUnitEnum(str, Enum):
    Liter = "l"
    Milliliter = "ml"
    Microliter = "ul"

class SpecialStorageConditionsEnum(str, Enum):
    Dry = "Dry"
    No_Oxygen = "No Oxygen"
    Fridge = "Fridge"
    Freezer = "Freezer"
    Dark = "Dark"
    Flammable = "Flammable"
    Poisenous = "Poisenous"

class CurrencyUnitEnum(str, Enum):
    Euro = "EUR"
    SwissFranc = "CHF"
    USDollar = "USD"

# Complex properties

class AngleValue(BaseModel):
    value: float
    unit: AngleUnitEnum
    
class CurrentValue(BaseModel):
    value: float
    unit: CurrentUnitEnum

class CurrencyValue(BaseModel):
    value: float
    unit: CurrencyUnitEnum

class DensityValue(BaseModel):
    value: float
    unit: DensityUnitEnum
    
class EnergyValue(BaseModel):
    value: float
    unit: EnergyUnitEnum

class MassValue(BaseModel):
    value: float
    unit: MassUnitEnum

class PressureValue(BaseModel):
    value: float
    unit: PressureUnitEnum

class TemperatureValue(BaseModel):
    value: float
    unit: TemperatureUnitEnum

class VoltageValue(BaseModel):
    value: float
    unit: VoltageUnitEnum

class VolumeValue(BaseModel):
    value: float
    unit: VolumeUnitEnum

# Objects

class OpenBISObject(BaseModel):
    permId: str = Field(default="", description="Permanent identifier for the OpenBIS object")
    name: str = Field(default="", description="Name of the OpenBIS object")
    description: str = Field(default="", description="Description of the OpenBIS object")
    registration_date: str = Field(default="", description="Date when the object was registered")
    comments: str = Field(default="", description="Additional comments about the object")

class Manufacturer(OpenBISObject):
    pass

class Grant(OpenBISObject):
    project_id: str
    acknowledgement_sentence: str
    budget: CurrencyValue
    start_date: str
    end_date: str
    funder_name: str
    
class Organisation(OpenBISObject):
    address: str = Field(default="", description="Address of the organisation")
    email: str = Field(default="", description="Email address of the organisation")
    work_phone: str = Field(default="", description="Work phone number of the organisation")
    organisation: "Organisation" = Field(default=None, description="Parent organisation (if applicable)")

class Person(OpenBISObject):
    username: str = Field(default="", description="Username of the person")
    email: str = Field(default="", description="Email address of the person")
    mobile_phone: str = Field(default="", description="Mobile phone number")
    work_phone: str = Field(default="", description="Work phone number")
    work_status: Literal["Active", "Inactive"] = Field(default="Active", description="Current work status")
    organisations: List[Organisation] = Field(default_factory=list, description="List of organisations the person belongs to")

class Supplier(OpenBISObject):
    pass

class Molecule(OpenBISObject):
    smiles: str = Field(default="", description = "SMILES string for the substance, e.g. CCO")
    sum_formula: str = Field(default="", description = "Molecular sum formula, e.g. CH4")
    cas_number: str = Field(default="", description="Unique numerical identifier assigned to a specific chemical substance by the Chemical Abstracts Service (CAS), e.g. 58-08-2")
    iupac_name: str = Field(default="", description="IUPAC standardised chemical name, e.g. 1,3,7-Trimethyl-3,7-dihydro-1H-purine-2,6-dione for caffeine.")

class Substance(OpenBISObject):
    molecules: List[Molecule] = Field(default_factory=list, description = "List of molecules that substance contains")
    empa_number: int = Field(default=1, ge=1, description = "Integer value given to new substances inside Empa")
    batch: str = Field(default="", description = "Letter given to the batch of substances, e.g., a for the first, b for the second.")
    vial: str = Field(default="", description = "Letter given to the vial of the batch of substances, e.g. i, ii, iii.")
    evaporation_temperatures: Dict = Field(default_factory=dict, description="Evaporation temperatures")
    purity: float = Field(default=0.0, description="Purity of the substance")
    # substance_type: str = Field(default="", description="Type of the substance") What is this?
    amount: Union[MassValue, VolumeValue] = Field(default=None, description="Amount of substance")
    chemist_own_name: str = Field(default="", description="Chemist's own name for the substance")
    location: "Location" = Field(default=None, description="Storage location")
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions required")
    package_opening_date: str = Field(default="", description="Date when package was opened")
    object_status: ObjectStatusEnum = Field(default=ObjectStatusEnum.Active, description="Current status of the object")
    supplier: Supplier = Field(default=None, description="Supplier information")
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the substance")
    supplier_own_name: str = Field(default="", description="Supplier's own name for the substance")
    receive_date: str = Field(default="", description="Date when substance was received")

class CrystalConcept(OpenBISObject):
    face: str = Field(default="", description="Crystal face information")
    material: str = Field(default="", description="Material of the crystal")
    crystal_lattice_parameters: Dict = Field(default_factory=dict, description="Dictionary of crystal lattice parameters")
    crystal_space_group: int = Field(default=0, description="Crystal space group number")

class Crystal(OpenBISObject):
    crystal_concept: CrystalConcept = Field(default=None, description="Crystal concept information")
    dimensions: Dict = Field(default_factory=dict, description="Crystal dimensions")
    reference_number: str = Field(default="", description="Reference number for the crystal")
    location: "Location" = Field(default=None, description="Storage location") #TODO: Correct this
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions required")
    package_opening_date: str = Field(default="", description="Date when package was opened")
    object_status: ObjectStatusEnum = Field(default=ObjectStatusEnum.Active, description="Current status of the object")
    supplier: Supplier = Field(default=None, description="Supplier information")
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the crystal")
    supplier_own_name: str = Field(default="", description="Supplier's own name for the crystal")
    receive_date: str = Field(default="", description="Date when crystal was received")

class TwoDLayerMaterial(OpenBISObject):
    top_layer_material: str = Field(default="", description="Material of the top layer")
    layer_count: str = Field(default="", description="Number of layers as string")
    impurities: str = Field(default="", description="Information about impurities present")
    heretostructure_stack: Dict = Field(default_factory=dict, description="Heterostructure stack")
    substrate: str = Field(default="", description="Substrate material")
    growth_method: str = Field(default="", description="Method used for growth")
    dimensions: Dict = Field(default_factory=dict, description="Dimensions")
    location: "Location" = Field(default=None, description="Storage location")
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions required")
    package_opening_date: str = Field(default="", description="Date when package was opened")
    sample_plate: str = Field(default="", description="Sample plate identifier")
    object_status: ObjectStatusEnum = Field(default=ObjectStatusEnum.Active, description="Current status of the object")
    supplier: Supplier = Field(default=None, description="Supplier information")
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the sample")
    supplier_own_name: str = Field(default="", description="Supplier's own name for the material")
    receive_date: str = Field(default="", description="Date when sample was received")

class Wafer(OpenBISObject):
    material: str
    face: str
    doping: str
    material_coating: Dict
    dimensions: Dict
    location: "Location"
    special_storage_conditions: List[SpecialStorageConditionsEnum]
    package_opening_date: str
    sample_plate: str
    object_status: ObjectStatusEnum
    supplier: Supplier
    synthesised_by: List[Person]
    supplier_own_name: str
    receive_date: str

class WaferSubstrate(OpenBISObject):
    material: str
    face: str
    material_coating: Dict
    dimensions: Dict
    location: "Location"
    special_storage_conditions: List[SpecialStorageConditionsEnum]
    package_opening_date: str
    sample_plate: str
    object_status: ObjectStatusEnum
    supplier: Supplier
    synthesised_by: List[Person]
    supplier_own_name: str
    receive_date: str

class Wire(OpenBISObject):
    material: str
    purity: float
    dimensions: Dict
    location: "Location"
    special_storage_conditions: List[SpecialStorageConditionsEnum]
    package_opening_date: str
    object_status: ObjectStatusEnum
    supplier: Supplier
    synthesised_by: List[Person]
    supplier_own_name: str
    receive_date: str

class Sample(OpenBISObject):
    exists: bool
    crystal: Crystal
    two_d_layer_material: TwoDLayerMaterial
    wafer_substrate: WaferSubstrate
    process_step: "ProcessStep"

class ReactionProductConcept(OpenBISObject):
    sum_formula: str
    molecules: List[Molecule]
    crystal_concepts: List[CrystalConcept]

class ReactionProduct(OpenBISObject):
    temperature: TemperatureValue
    reaction_yield: float
    reaction_product_concept: ReactionProductConcept
    sample: Sample

class DeviceSubstrate(OpenBISObject):
    location: OpenBISObject # TODO: correct
    special_storage_conditions: List[SpecialStorageConditionsEnum]
    sample_plate: str
    object_status: ObjectStatusEnum
    supplier: Supplier
    synthesied_by: List[Person]
    supplier_own_name: str
    receive_date: str

class WaferSample(OpenBISObject):
    pass

class Component(OpenBISObject):
    component_main_category: ComponentMainCategoryEnum
    component_sub_category: ComponentSubCategoryEnum
    manufacturer: Manufacturer
    model: str
    serial_number: str
    receive_date: str
    location: OpenBISObject #TODO: Merge all the possible locations
    object_status: ObjectStatusEnum
    actions_settings: dict #TODO: Define this object
    observables_settings: dict #TODO: Define this object

class Analyser(Component):
    density: DensityValue

class Chamber(Component):
    pass

class Cryostat(Component):
    pass

class Electronics(Component):
    pass

class EvaporatorSlot(Component):
    target_temperature: TemperatureValue
    slot_number: int = Field(..., strict = True, ge = 1, le = 6)
    p_value: float
    i_value: float
    ep_percentage: float = Field(..., strict = True, ge = 0, le = 100)

class Evaporator(Component):
    evaporator_slots: List[EvaporatorSlot]

class IonGauge(Component):
    filament: str
    filament_current: CurrentValue

class IonPump(Component):
    pass

class PBNStage(Component):
    target_temperature: TemperatureValue

class SputterGun(Component):
    bias_voltage: VoltageValue
    discharge_voltage: VoltageValue
    discharge_current: CurrentValue

class ScrollPump(Component):
    pass

class STM_AFM_TipSensor(Component):
    pass

class Thyracont(Component):
    pass

class TurboPump(Component):
    pass

class Valve(Component):
    pass

class Instrument(OpenBISObject):
    manufacturer: Manufacturer
    model: str
    serial_number: str
    empa_id: str
    receive_date: str
    location: OpenBISObject #TODO: Merge all the possible locations
    object_status: ObjectStatusEnum
    responsibles: List[Person]

class InstrumentSTM(Instrument):
    # Vacuum system
    pumps: List[Union[IonPump, ScrollPump, TurboPump]]
    gauges: List[Union[IonGauge]]
    vacuum_chambers: List[Chamber]
    ports_valves: List[Union[Component, Valve]]
    # Sample Preparation & Handling
    preparation_tools: List[Union[Component, PBNStage, SputterGun, Valve, Evaporator, Cryostat]]
    analysers: List[Analyser]
    mechanical_components: List[Component]
    # Microscope Core Components
    stm_components: List[Union[Component, Electronics, STM_AFM_TipSensor]]
    control_data_acquisition: List[Electronics]
    temperature_environment_control: List[Union[Component, Cryostat]]
    # Auxiliary
    auxiliary_components: List[Component]
    # Consumables
    consumables: List
    # Materials
    substances: List[Substance]
    single_crystals: List[Crystal]
    wafer_samples: List[WaferSample]
    # Accessories
    tip_sensors: List[Component]
    accessories: List[Component]
    # Logbook
    logbook_entries: List[Component]

class Code(OpenBISObject):
    version: str
    filepath_executable: str
    repository_url: str

class MeasurementSession(OpenBISObject):
    pass

class Software(OpenBISObject):
    version: str

class Simulation(OpenBISObject):
    wfms_uuid: str

class AtomisticModel(Simulation):
    cell: Dict
    dimensionality: int
    periodic_boundary_conditions: List[bool] = Field(min_items = 3, max_items = 3)
    volume: float
    crystal_concepts: List[CrystalConcept]
    geometry_optimisation: "GeometryOptimisation"
    molecules: List[Molecule]
    reaction_product_concepts: List[ReactionProductConcept]

class BandStructure(Simulation):
    band_gap: EnergyValue
    level_theory_method: str
    level_theory_parameters: Dict
    input_parameters: Dict
    output_parameters: Dict
    codes: List[Code]

class GeometryOptimisation(Simulation):
    cell_opt_constraints: Dict
    cell_optimised: bool
    driver_code: str
    constrained: bool
    force_convergence_threshold: Dict
    level_theory_method: str
    level_theory_parameters: Dict
    input_parameters: Dict
    output_parameters: Dict
    atomistic_models: List["AtomisticModel"]
    codes: List[Code]

class MeanFieldHubbard(OpenBISObject):
    charge: int
    convergence_coefficient: float
    density_tolerance: float
    down_densities: List[int]
    down_electrons: List[int]
    number_k_points: int
    on_site_hubbard_repulsion_u_value: float
    on_site_hubbard_repulsion_u_unit: str
    up_densities: List[int]
    up_electrons: List[int]
    atomsitic_models: List[AtomisticModel]

class MinimumEnergyPotential(OpenBISObject):
    energies: List[EnergyValue]
    energy_barrier: EnergyValue
    geometry_constraints: List[float]
    geometry_constraints_inscrements: List[float]
    method_type: Literal["NEB", "String method"]
    number_geometries: int
    atomistic_models: List[AtomisticModel]

class PDOS(Simulation):
    level_theory_method: str
    level_theory_parameters: Dict
    input_parameters: Dict
    output_parameters: Dict
    atomistic_models: List["AtomisticModel"]
    codes: List[Code]

class UnclassifiedSimulation(Simulation):
    input_parameters: Dict
    output_parameters: Dict
    atomistic_models: List["AtomisticModel"]
    codes: List[Code]

class VibrationalSpectroscopy(Simulation):
    level_theory_method: str
    level_theory_parameters: Dict
    input_parameters: Dict
    output_parameters: Dict
    atomistic_models: List["AtomisticModel"]
    codes: List[Code]

class PotentialEnergyCalculation(Simulation):
    energy: EnergyValue
    level_theory_method: str
    level_theory_parameters: Dict
    input_parameters: Dict
    output_parameters: Dict
    atomistic_models: List["AtomisticModel"]
    codes: List[Code]

class Action(OpenBISObject):
    duration: str
    component: Component
    component_settings: dict
    
class Annealing(Action):
    target_temperature: TemperatureValue

class Coating(Action):
    pass

class Delamination(Action):
    pass

class Deposition(Action):
    substrate_temperature: TemperatureValue
    substance: Substance

class Dosing(Action):
    dosing_gas: Substance
    pressure: PressureValue
    substrate_temperature: TemperatureValue

class Etching(Action):
    pass

class Fishing(Action):
    pass

class FieldEmission(Action):
    pass

class LightIrradiation(Action):
    pass

class MechanicalPressing(Action):
    pass

class Rinse(Action):
    pass

class Sputtering(Action):
    sputter_ion: str
    pressure: PressureValue
    current: CurrencyValue
    angle: AngleValue
    substrate_temperature: TemperatureValue

class Location(OpenBISObject):
    organisation: Organisation

class Analysis(OpenBISObject):
    codes: List[Code]
    measurements: List[MeasurementSession]
    software: List[Software]

class Result(OpenBISObject):
    analysis: List[Analysis]
    simulations: List[Union[BandStructure, GeometryOptimisation, PDOS, VibrationalSpectroscopy]]
    measurements: List[MeasurementSession]

class Draft(OpenBISObject):
    draft_type: DraftTypeEnum
    results: List[Result]

class AiidaNode(OpenBISObject):
    pass

class Observable(OpenBISObject):
    channel_name: str
    component: Component
    component_settings: Dict

class CurrentObservable(Observable):
    pass

class ElementalCompositionObservable(Observable):
    pass

class FluxObservable(Observable):
    pass

class ForceObservable(Observable):
    pass

class InductanceObservable(Observable):
    pass

class PHValueObservable(Observable):
    pass

class PressureObservable(Observable):
    pass

class ResistanceObservable(Observable):
    pass

class SpeedObservable(Observable):
    pass

class TemperatureObservable(Observable):
    pass

class VoltageObservable(Observable):
    pass

class MeasurementSession(OpenBISObject):
    measurement_folder_path: str
    default_object_view: Literal["IMAGING_GALLERY_VIEW"]
    measurement_session: "MeasurementSession"

class Process(OpenBISObject):
    short_name: str
    instrument: Union[Instrument, InstrumentSTM]
    process_steps_settings: Dict

class ProcessStep(OpenBISObject):
    actions: List[Action]
    observables: List[Observable]
    instrument: Union[Instrument, InstrumentSTM]

class Preparation(OpenBISObject):
    process_steps: List[ProcessStep]

class Publication(OpenBISObject):
    abstract: str
    doi: str
    url: str
    dataset_url: str
    year: int
    authors: List[Person]
    drafts: List[Draft]
    grants: List[Grant]

# atomistic_model = AtomisticModel(name = "aa")
# GeometryOptimisation(atomistic_models = [atomistic_model])