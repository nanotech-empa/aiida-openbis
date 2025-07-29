from pydantic import BaseModel, Field, model_validator, field_validator, conint
from typing import List, Union, Dict, Literal, Optional
from enum import Enum
import re
from datetime import datetime

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

class LengthUnitEnum(str, Enum):
    Meter = "m"
    Millimeter = "mm"
    Micrometer = "um"
    Nanometer = "nm"
    Angstorm = "Ang"

class EnergyDensityUnitEnum(str, Enum):
    ElectronVolt_by_CubicBohr = "eV/Bohr^3"

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

class EnergyDensityValue(BaseModel):
    value: float
    unit: EnergyDensityUnitEnum

class LengthValue(BaseModel):
    value: float
    unit: LengthUnitEnum

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

class ComponentActionSettings(BaseModel):
    action_type: str = Field(default="", description="Type of action")
    component_properties: List[str] = Field(default_factory=list, description="Properties of the component that are used by this action type")

class ObservableReading(BaseModel):
    reading_source: str = Field(default="", description="Reading source name")
    reading_unit: str = Field(default="", description="Reading unit")
    reading_variable: str = Field(default="", description="Reading variable") # TODO: Ask Bruno what is this

class ComponentObservableSettings(BaseModel):
    observable_type: str = Field(default="", description="Type of observable")
    component_properties: List[str] = Field(default_factory=list, description="List of properties of the component that are used by this observable type")
    observable_readings: List[ObservableReading] = Field(default_factory=list, description="List of readings (channels) details")

class CrystalLatticeParameters(BaseModel):
    values: List[float] = Field(default_factory=list, description="List of lattice values")
    units: List[Union[LengthValue, AngleValue]] = Field(default_factory=list, description="List of lattice values units")

class Dimensions(BaseModel):
    shape: Literal["rectangular", "round"] = Field(default="rectangular", description="Shape of the object")
    diameter: LengthValue = Field(default=None, description="Diameter for round shapes")
    width: LengthValue = Field(default=None, description="Width of the object")
    height: LengthValue = Field(default=None, description="Height of the object")
    thickness: LengthValue = Field(default=None, description="Thickness of the object")
    
    @model_validator(mode='after')
    def validate_shape_dimensions(self):
        if self.shape == "round" and self.diameter is None:
            raise ValueError("Diameter is required for round shapes")
        return self

class MaterialCoating(BaseModel):
    material: str
    face: str
    thickness: LengthValue

class HeterostructureStack(BaseModel):
    impurities: str
    material: str
    thickness: LengthValue
    source: str

class ActionSettings(BaseModel):
    action_type: str
    action: Union["Annealing", "Coating", "Cooldown", "Deposition", "Delamination", "Dosing", "Etching", "Fishing", "FieldEmission", "LightIrradiation", "MechanicalPressing", "Rinse", "Sputtering"]

class ObservableSettings(BaseModel):
    observable_type: str
    observable: Union["CurrentObservable", "ElementalCompositionObservable", "FluxObservable", "ForceObservable", "InductanceObservable", "PHValueObservable", "PressureObservable", "ResistanceObservable", "SpeedObservable", "TemperatureObservable", "VoltageObservable"]

class ProcessStepSettings(BaseModel):
    name: str
    actions_settings: List[ActionSettings]
    observables_settings: List[ObservableSettings]

# Objects

class OpenBISObject(BaseModel):
    permId: str = Field(default="", description="Permanent identifier for the OpenBIS object", metadata={"openbis_type": "No type"})
    name: str = Field(default="", description="Name of the OpenBIS object", metadata={"openbis_type": "VARCHAR"})
    description: str = Field(default="", description="Description of the OpenBIS object", metadata={"openbis_type": "VARCHAR"})
    registration_date: str = Field(default="", description="Date when the object was registered", metadata={"openbis_type": "No type"})
    comments: str = Field(default="", description="Additional comments about the object", metadata={"openbis_type": "MULTILINE_VARCHAR"})

class Manufacturer(OpenBISObject):
    pass

class Grant(OpenBISObject):
    project_id: str = Field(default=None, description="Identifier of the project", metadata={"openbis_type": "VARCHAR"})
    acknowledgement_sentence: str = Field(default=None, description="Acknowledgement sentence", metadata={"openbis_type": "MULTILINE_VARCHAR"})
    budget: CurrencyValue = Field(default=None, description="Budget allocated", metadata={"openbis_type": "JSON"})
    start_date: str = Field(default=None, description="Start date in 'YYYY-MM-DD' format", metadata={"openbis_type": "DATE"})
    end_date: str = Field(default=None, description="End date in 'YYYY-MM-DD' format", metadata={"openbis_type": "DATE"})
    funder_name: str = Field(default=None, description="Name of the funder", metadata={"openbis_type": "VARCHAR"})
    
    @field_validator('start_date', 'end_date')
    @classmethod
    def validate_date_format(cls, v: Optional[str]) -> Optional[str]:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("Date must be in 'YYYY-MM-DD' format.")
        return v

    @model_validator(mode='after')
    def check_date_order(cls, model):
        if model.start_date and model.end_date:
            start = datetime.strptime(model.start_date, "%Y-%m-%d")
            end = datetime.strptime(model.end_date, "%Y-%m-%d")
            if end < start:
                raise ValueError("end_date must be greater than or equal to start_date.")
        return model
    
class Organisation(OpenBISObject):
    address: str = Field(default="", description="Address of the organisation", metadata={"openbis_type": "VARCHAR"})
    email: str = Field(default="", description="Email address of the organisation", metadata={"openbis_type": "VARCHAR"})
    work_phone: str = Field(default="", description="Work phone number of the organisation", metadata={"openbis_type": "VARCHAR"})
    organisation: "Organisation" = Field(default=None, description="Parent organisation (if applicable)", metadata={"openbis_type": "SAMPLE"})

class Person(OpenBISObject):
    username: str = Field(default="", description="Username of the person", metadata={"openbis_type": "VARCHAR"})
    email: str = Field(default="", description="Email address of the person", metadata={"openbis_type": "VARCHAR"})
    mobile_phone: str = Field(default="", description="Mobile phone number", metadata={"openbis_type": "VARCHAR"})
    work_phone: str = Field(default="", description="Work phone number", metadata={"openbis_type": "VARCHAR"})
    work_status: Literal["Active", "Inactive"] = Field(default="Active", description="Current work status", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    organisations: List[Organisation] = Field(default_factory=list, description="List of organisations the person belongs to", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class Supplier(OpenBISObject):
    pass

class Molecule(OpenBISObject):
    smiles: str = Field(default="", description = "SMILES string for the substance, e.g. CCO", metadata={"openbis_type": "VARCHAR"})
    sum_formula: str = Field(default="", description = "Molecular sum formula, e.g. CH4", metadata={"openbis_type": "VARCHAR"})
    cas_number: str = Field(default="", description="Unique numerical identifier assigned to a specific chemical substance by the Chemical Abstracts Service (CAS), e.g. 58-08-2", metadata={"openbis_type": "VARCHAR"})
    iupac_name: str = Field(default="", description="IUPAC standardised chemical name, e.g. 1,3,7-Trimethyl-3,7-dihydro-1H-purine-2,6-dione for caffeine.", metadata={"openbis_type": "VARCHAR"})

class Substance(OpenBISObject):
    molecules: List[Molecule] = Field(default_factory=list, description = "List of molecules that substance contains", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    empa_number: int = Field(default=0, ge = 1, description = "Integer value given to new substances inside Empa (0 = unassigned)", metadata={"openbis_type": "INTEGER"})
    batch: str = Field(default="", description = "Letter given to the batch of substances, e.g., a for the first, b for the second.", metadata={"openbis_type": "VARCHAR"})
    vial: str = Field(default="", description = "Letter given to the vial of the batch of substances, e.g. i, ii, iii.", metadata={"openbis_type": "VARCHAR"})
    evaporation_temperatures: Dict = Field(default_factory=dict, description="Evaporation temperatures", metadata={"openbis_type": "XML", "custom_widget": "Spreadsheet"})
    purity: float = Field(default=0.0, description="Purity of the substance", metadata={"openbis_type": "REAL"})
    # substance_type: str = Field(default="", description="Type of the substance") What is this?
    amount: Union[MassValue, VolumeValue] = Field(default=None, description="Amount of substance", metadata={"openbis_type": "JSON"})
    chemist_own_name: str = Field(default="", description="Chemist's own name for the substance", metadata={"openbis_type": "VARCHAR"})
    location: Union["Instrument", "InstrumentSTM", "Location"] = Field(default=None, description="Storage location", metadata={"openbis_type": "SAMPLE"})
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions required", metadata={"openbis_type": "CONTROLLEDVOCABULARY", "openbis_multivalue": True})
    package_opening_date: str = Field(default="", description="Date when package was opened", metadata={"openbis_type": "DATE"})
    object_status: ObjectStatusEnum = Field(default=ObjectStatusEnum.Active, description="Current status of the object", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    supplier: Supplier = Field(default=None, description="Supplier information", metadata={"openbis_type": "SAMPLE"})
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the substance", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    supplier_own_name: str = Field(default="", description="Supplier's own name for the substance", metadata={"openbis_type": "VARCHAR"})
    receive_date: str = Field(default=None, description="Date when the substance was received", metadata={"openbis_type": "DATE"})

    @field_validator("receive_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v

class CrystalConcept(OpenBISObject):
    face: str = Field(default="", description="Crystal face information", metadata={"openbis_type": "VARCHAR"})
    material: str = Field(default="", description="Material of the crystal", metadata={"openbis_type": "VARCHAR"})
    crystal_lattice_parameters: CrystalLatticeParameters = Field(default=None, description="Crystal lattice parameters", metadata={"openbis_type": "JSON"})
    crystal_space_group: int = Field(default=0, description="Crystal space group number", metadata={"openbis_type": "INTEGER"})

class Crystal(OpenBISObject):
    crystal_concept: CrystalConcept = Field(default=None, description="Crystal concept information", metadata={"openbis_type": "SAMPLE"})
    dimensions: Dimensions = Field(default=None, description="Crystal dimensions", metadata={"openbis_type": "JSON"})
    reference_number: str = Field(default="", description="Reference number for the crystal", metadata={"openbis_type": "VARCHAR"})
    location: Union["Location", "Instrument", "InstrumentSTM"] = Field(default=None, description="Storage location", metadata={"openbis_type": "SAMPLE"})
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions required", metadata={"openbis_type": "CONTROLLEDVOCABULARY", "openbis_multivalue": True})
    package_opening_date: str = Field(default="", description="Date when package was opened", metadata={"openbis_type": "DATE"})
    object_status: ObjectStatusEnum = Field(default=ObjectStatusEnum.Active, description="Current status of the object", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    supplier: Supplier = Field(default=None, description="Supplier information", metadata={"openbis_type": "SAMPLE"})
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the crystal", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    supplier_own_name: str = Field(default="", description="Supplier's own name for the crystal", metadata={"openbis_type": "VARCHAR"})
    receive_date: str = Field(default="", description="Date when crystal was received", metadata={"openbis_type": "DATE"})

    @field_validator("receive_date", "package_opening_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v

class TwoDLayerMaterial(OpenBISObject):
    top_layer_material: str = Field(default="", description="Material of the top layer", metadata={"openbis_type": "VARCHAR"})
    layer_count: str = Field(default="", description="Number of layers as string", metadata={"openbis_type": "VARCHAR"})
    impurities: str = Field(default="", description="Information about impurities present", metadata={"openbis_type": "VARCHAR"})
    heterostructure_stack: List[HeterostructureStack] = Field(default_factory=list, description="Heterostructure stack", metadata={"openbis_type": "XML", "custom_widget": "Spreadsheet"})
    substrate: str = Field(default="", description="Substrate material", metadata={"openbis_type": "VARCHAR"})
    growth_method: str = Field(default="", description="Method used for growth", metadata={"openbis_type": "VARCHAR"})
    dimensions: Dimensions = Field(default_factory=dict, description="Dimensions", metadata={"openbis_type": "JSON"})
    location: Union["Instrument", "InstrumentSTM", "Location"] = Field(default=None, description="Storage location", metadata={"openbis_type": "SAMPLE"})
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions required", metadata={"openbis_type": "CONTROLLEDVOCABULARY", "openbis_multivalue": True})
    package_opening_date: str = Field(default="", description="Date when package was opened", metadata={"openbis_type": "DATE"})
    sample_plate: str = Field(default="", description="Sample plate identifier", metadata={"openbis_type": "VARCHAR"})
    object_status: ObjectStatusEnum = Field(default=ObjectStatusEnum.Active, description="Current status of the object", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    supplier: Supplier = Field(default=None, description="Supplier information", metadata={"openbis_type": "SAMPLE"})
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the sample", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    supplier_own_name: str = Field(default="", description="Supplier's own name for the material", metadata={"openbis_type": "VARCHAR"})
    receive_date: str = Field(default="", description="Date when sample was received", metadata={"openbis_type": "DATE"})

    @field_validator("receive_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v

class Wafer(OpenBISObject):
    material: str = Field(default=None, description="Material of the wafer", metadata={"openbis_type": "VARCHAR"})
    face: str = Field(default=None, description="Face orientation of the wafer", metadata={"openbis_type": "VARCHAR"})
    doping: str = Field(default=None, description="Doping type or level", metadata={"openbis_type": "VARCHAR"})
    material_coating: List[MaterialCoating] = Field(default_factory=list, description="Material coating details", metadata={"openbis_type": "JSON"})
    dimensions: Dimensions = Field(default_factory=Dimensions, description="Dimensions of the wafer", metadata={"openbis_type": "JSON"})
    location: Union["Instrument", "InstrumentSTM", "Location"] = Field(default=None, description="Storage location", metadata={"openbis_type": "SAMPLE"})
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions", metadata={"openbis_type": "CONTROLLEDVOCABULARY", "openbis_multivalue": True})
    package_opening_date: str = Field(default=None, description="Date when the package was opened", metadata={"openbis_type": "DATE"})
    sample_plate: str = Field(default=None, description="Sample plate identifier", metadata={"openbis_type": "VARCHAR"})
    object_status: ObjectStatusEnum = Field(default=None, description="Status of the wafer", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    supplier: Supplier = Field(default=None, description="Supplier information", metadata={"openbis_type": "SAMPLE"})
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the wafer", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    supplier_own_name: str = Field(default=None, description="Supplier's own name for the wafer", metadata={"openbis_type": "VARCHAR"})
    receive_date: str = Field(default="", description="Date when wafer was received", metadata={"openbis_type": "VARCHAR"})

    @field_validator("receive_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v

class WaferSubstrate(OpenBISObject):
    material: str = Field(default=None, description="Material of the wafer", metadata={"openbis_type": "VARCHAR"})
    face: str = Field(default=None, description="Face orientation of the wafer", metadata={"openbis_type": "VARCHAR"})
    material_coating: List[MaterialCoating] = Field(default_factory=list, description="Material coating details", metadata={"openbis_type": "JSON"})
    dimensions: Dimensions = Field(default_factory=Dimensions, description="Dimensions of the wafer", metadata={"openbis_type": "JSON"})
    location: Union["Instrument", "InstrumentSTM", "Location"] = Field(default=None, description="Storage location", metadata={"openbis_type": "SAMPLE"})
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions", metadata={"openbis_type": "CONTROLLEDVOCABULARY", "openbis_multivalue": True})
    package_opening_date: str = Field(default=None, description="Date when the package was opened", metadata={"openbis_type": "VARCHAR"})
    sample_plate: str = Field(default=None, description="Sample plate identifier", metadata={"openbis_type": "VARCHAR"})
    object_status: ObjectStatusEnum = Field(default=None, description="Status of the wafer", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    supplier: Supplier = Field(default=None, description="Supplier information", metadata={"openbis_type": "SAMPLE"})
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the wafer", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    supplier_own_name: str = Field(default=None, description="Supplier's own name for the wafer", metadata={"openbis_type": "VARCHAR"})
    receive_date: str = Field(default="", description="Date when wafer substrate was received", metadata={"openbis_type": "DATE"})

    @field_validator("receive_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v

class Wire(OpenBISObject):
    material: str = Field(default=None, description="Material of the wire", metadata={"openbis_type": "VARCHAR"})
    purity: float = Field(default=None, description="Purity level of the wire", metadata={"openbis_type": "REAL"})
    dimensions: Dimensions = Field(default=None, description="Dimensions of the wafer", metadata={"openbis_type": "JSON"})
    location: Union["Instrument", "InstrumentSTM", "Location"] = Field(default=None, description="Storage location", metadata={"openbis_type": "SAMPLE"})
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="Special storage conditions", metadata={"openbis_type": "CONTROLLEDVOCABULARY", "openbis_multivalue": True})
    package_opening_date: str = Field(default=None, description="Date when the package was opened", metadata={"openbis_type": "DATE"})
    object_status: ObjectStatusEnum = Field(default=None, description="Status of the wafer", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    supplier: Supplier = Field(default=None, description="Supplier information", metadata={"openbis_type": "SAMPLE"})
    synthesised_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the wafer", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    supplier_own_name: str = Field(default=None, description="Supplier's own name for the wafer", metadata={"openbis_type": "VARCHAR"})
    receive_date: str = Field(default="", description="Date when wire was received", metadata={"openbis_type": "DATE"})

    @field_validator("receive_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v

class Sample(OpenBISObject):
    exists: bool = Field(default=True, description="Whether the sample physically exists", metadata={"openbis_type": "BOOLEAN"})
    crystal: Crystal = Field(default=None, description="Crystal structure of the sample", metadata={"openbis_type": "PARENT"})
    two_d_layer_material: TwoDLayerMaterial = Field(default=None, description="2D layer material associated with the sample", metadata={"openbis_type": "PARENT"})
    wafer_substrate: WaferSubstrate = Field(default=None, description="Wafer substrate used for the sample", metadata={"openbis_type": "PARENT"})
    process_step: "ProcessStep" = Field(default=None, description="Process step used to prepare the sample", metadata={"openbis_type": "PARENT"})

class ReactionProductConcept(OpenBISObject):
    sum_formula: str = Field(default=None, description="Sum formula of the reaction product concept", metadata={"openbis_type": "VARCHAR"})
    molecules: List[Molecule] = Field(default_factory=list, description="List of molecules involved", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    crystal_concepts: List[CrystalConcept] = Field(default_factory=list, description="Associated crystal concepts", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class ReactionProduct(OpenBISObject):
    temperature: TemperatureValue = Field(default=None, description="Temperature at which the product was obtained", metadata={"openbis_type": "JSON"})
    reaction_yield: float = Field(default=None, description="Yield of the reaction product", metadata={"openbis_type": "REAL"})
    reaction_product_concept: ReactionProductConcept = Field(default=None, description="Concept associated with the reaction product", metadata={"openbis_type": "SAMPLE"})
    sample: Sample = Field(default=None, description="Sample associated with the reaction product", metadata={"openbis_type": "PARENT"})

class DeviceSubstrate(OpenBISObject):
    location: Union["Location", "Instrument", "InstrumentSTM"] = Field(default=None, description="Location of the substrate", metadata={"openbis_type": "SAMPLE"})
    special_storage_conditions: List[SpecialStorageConditionsEnum] = Field(default_factory=list, description="List of special storage requirements", metadata={"openbis_type": "CONTROLLEDVOCABULARY", "openbis_multivalue": True})
    sample_plate: str = Field(default=None, description="Identifier for the sample plate", metadata={"openbis_type": "VARCHAR"})
    object_status: ObjectStatusEnum = Field(default="Active", description="Current status of the object", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    supplier: Supplier = Field(default=None, description="Supplier information", metadata={"openbis_type": "SAMPLE"})
    synthesied_by: List[Person] = Field(default_factory=list, description="List of people who synthesized the device substrate", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    supplier_own_name: str = Field(default=None, description="Supplier-specific name for the substrate", metadata={"openbis_type": "VARCHAR"})
    receive_date: str = Field(default=None, description="Date when the device substrate was received", metadata={"openbis_type": "DATE"})

    @field_validator("receive_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v
    
class WaferSample(OpenBISObject):
    pass

class Component(OpenBISObject):
    component_main_category: ComponentMainCategoryEnum = Field(default="", description="Main category of the component", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    component_sub_category: ComponentSubCategoryEnum = Field(default="", description="Sub category of the component", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    manufacturer: Manufacturer = Field(default=None, description="Manufacturer of the component", metadata={"openbis_type": "SAMPLE"})
    model: str = Field(default="", description="Model of the component", metadata={"openbis_type": "VARCHAR"})
    serial_number: str = Field(default="", description="Serial number of the component", metadata={"openbis_type": "VARCHAR"})
    location: Union["Location", "Instrument", "InstrumentSTM"] = Field(default=None, description="Current location of the component", metadata={"openbis_type": "SAMPLE"})
    object_status: ObjectStatusEnum = Field(default=ObjectStatusEnum.Active, description="Current status of the component", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    component_actions_settings: List[ComponentActionSettings] = Field(default_factory=list, description="Action types that use this component together with the names of the component properties that the action types use.", metadata={"openbis_type": "JSON"})
    component_observables_settings: List[ComponentObservableSettings]  = Field(default_factory=list, description="Observable types that use this component together with the names of the component properties that the observable types use and the readings that the observable type comprise.", metadata={"openbis_type": "JSON"})
    receive_date: str = Field(default="", description="Date when component was received", metadata={"openbis_type": "DATE"})

    @field_validator("receive_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v

class Analyser(Component):
    density: DensityValue = Field(default=None, description="Density value", metadata={"openbis_type": "JSON"})

class Chamber(Component):
    pass

class Cryostat(Component):
    pass

class Electronics(Component):
    pass

class EvaporatorSlot(Component):
    target_temperature: TemperatureValue = Field(default = None, description = "Target temperature for evaporation of the substance in it.", metadata={"openbis_type": "JSON"})
    slot_number: int = Field(default = 1, strict = True, ge = 1, le = 6, description = "Number of the evaporator slot.", metadata={"openbis_type": "INTEGER"})
    p_value: float = Field(default = 0, strict = True, ge = 0, description = "P-value", metadata={"openbis_type": "REAL"})
    i_value: float = Field(default = 0, strict = True, ge = 0, description = "I-value", metadata={"openbis_type": "REAL"})
    ep_percentage: float = Field(default = 0, strict = True, ge = 0, le = 100, description = "EP percentage", metadata={"openbis_type": "REAL"})

class Evaporator(Component):
    evaporator_slots: List[EvaporatorSlot] = Field(default_factory = list, description = "List of evaporator slots attached to the evaporator.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class IonGauge(Component):
    filament: str = Field(default = "", description = "Filament material", metadata={"openbis_type": "VARCHAR"})
    filament_current: CurrentValue = Field(default = None, description = "Filament current", metadata={"openbis_type": "JSON"})

class IonPump(Component):
    pass

class PBNStage(Component):
    target_temperature: TemperatureValue = Field(default = None, description = "Target temperature for evaporation of the substance in it.", metadata={"openbis_type": "JSON"})

class SputterGun(Component):
    bias_voltage: VoltageValue = Field(default = None, description = "Bias voltage set in the sputter gun.", metadata={"openbis_type": "JSON"})
    discharge_voltage: VoltageValue = Field(default = None, description = "Discharge voltage set in the sputter gun.", metadata={"openbis_type": "JSON"})
    discharge_current: CurrentValue = Field(default = None, description = "Discharge current set in the sputter gun.", metadata={"openbis_type": "JSON"})

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
    manufacturer: Manufacturer = Field(default=None, description="Manufacturer of the instrument", metadata={"openbis_type": "SAMPLE"})
    model: str = Field(default=None, description="Model number", metadata={"openbis_type": "VARCHAR"})
    serial_number: str = Field(default=None, description="Serial number", metadata={"openbis_type": "VARCHAR"})
    empa_id: str = Field(default=None, description="Empa identifier", metadata={"openbis_type": "VARCHAR"})
    location: "Location" = Field(default=None, description="Physical location", metadata={"openbis_type": "SAMPLE"})
    object_status: ObjectStatusEnum = Field(default=None, description="Status of the instrument", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    responsibles: List[Person] = Field(default_factory=list, description="List of responsible persons", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    receive_date: str = Field(default=None, description="Date when the instrument was received", metadata={"openbis_type": "DATE"})

    @field_validator("receive_date")
    @classmethod
    def validate_receive_date(cls, v: str) -> str:
        if v is None:
            return v
        try:
            datetime.strptime(v, "%Y-%m-%d")
        except ValueError:
            raise ValueError("receive_date must be in 'YYYY-MM-DD' format.")
        return v

class InstrumentSTM(Instrument):
    # Vacuum system
    pumps: List[Union[IonPump, ScrollPump, TurboPump]] = Field(default_factory = list, description = "Pumps attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    gauges: List[Union[IonGauge]] = Field(default_factory = list, description = "Gauges attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    vacuum_chambers: List[Chamber] = Field(default_factory = list, description = "Vacuum chambers attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    ports_valves: List[Union[Component, Valve]] = Field(default_factory = list, description = "Ports and valves attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    # Sample Preparation & Handling
    preparation_tools: List[Union[Component, PBNStage, SputterGun, Valve, Evaporator, Cryostat]]  = Field(default_factory = list, description = "Preparation tools attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    analysers: List[Analyser] = Field(default_factory = list, description = "Analysers attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    mechanical_components: List[Component] = Field(default_factory = list, description = "Mechanical components attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    # Microscope Core Components
    stm_components: List[Union[Component, Electronics, STM_AFM_TipSensor]] = Field(default_factory = list, description = "STM components attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    control_data_acquisition: List[Electronics] = Field(default_factory = list, description = "Control and data acquisition components attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    temperature_environment_control: List[Union[Component, Cryostat]] = Field(default_factory = list, description = "Temperature and environment control components attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    # Auxiliary
    auxiliary_components: List[Component] = Field(default_factory = list, description = "Auxiliary components attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    # Consumables
    consumables: List = Field(default_factory = list, description = "Consumables inside the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    # Materials
    substances: List[Substance] = Field(default_factory = list, description = "Substances inside the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    single_crystals: List[Crystal] = Field(default_factory = list, description = "Crystals inside the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    wafer_samples: List[WaferSample] = Field(default_factory = list, description = "Wafer samples inside the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    # Accessories
    tip_sensors: List[Component] = Field(default_factory = list, description = "Tip sensors attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    accessories: List[Component] = Field(default_factory = list, description = "Other accessories attached to the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    # Logbook
    logbook_entries: List[Component] = Field(default_factory = list, description = "Logbook entries of the instrument.", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class Code(OpenBISObject):
    version: str = Field(default="", description="Version of the code", metadata={"openbis_type": "VARCHAR"})
    filepath_executable: str = Field(default="", description="File path to the executable", metadata={"openbis_type": "VARCHAR"})
    repository_url: str = Field(default="", description="URL of the code repository", metadata={"openbis_type": "VARCHAR"})

class Software(OpenBISObject):
    version: str = Field(default=None, description="Software version", metadata={"openbis_type": "VARCHAR"})

class Simulation(OpenBISObject):
    wfms_uuid: str = Field(default=None, description="Workflow management system UUID", metadata={"openbis_type": "VARCHAR"})

class AtomisticModel(Simulation):
    cell: Dict = Field(default_factory=dict, description="Dictionary describing the simulation cell", metadata={"openbis_type": "JSON"})
    dimensionality: int = Field(default=0, description="Dimensionality of the model (e.g., 1D, 2D, 3D)", metadata={"openbis_type": "INTEGER"})
    periodic_boundary_conditions: List[bool] = Field(default_factory=lambda: [False, False, False], min_length=3, max_length=3, description="Periodic boundary conditions for x, y, z directions", metadata={"openbis_type": "BOOLEAN", "openbis_multivalue": True})
    volume: float = Field(default=0.0, description="Volume of the simulation cell", metadata={"openbis_type": "REAL"})
    crystal_concepts: List[CrystalConcept] = Field(default_factory=list, description="List of crystal concepts in the model", metadata={"openbis_type": "PARENT"})
    geometry_optimisation: "GeometryOptimisation" = Field(default=None, description="Geometry optimization settings", metadata={"openbis_type": "PARENT"})
    molecules: List[Molecule] = Field(default_factory=list, description="List of molecules in the model", metadata={"openbis_type": "PARENT"})
    reaction_product_concepts: List[ReactionProductConcept] = Field(default_factory=list, description="List of reaction product concepts", metadata={"openbis_type": "PARENT"})

class BandStructure(Simulation):
    band_gap: EnergyValue = Field(default=None, description="Band gap energy value", metadata={"openbis_type": "JSON"})
    level_theory_method: str = Field(default="", description="Method used for level of theory calculations", metadata={"openbis_type": "VARCHAR"})
    level_theory_parameters: Dict = Field(default_factory=dict, description="Parameters for the level of theory method", metadata={"openbis_type": "JSON"})
    input_parameters: Dict = Field(default_factory=dict, description="Input parameters for the simulation", metadata={"openbis_type": "JSON"})
    output_parameters: Dict = Field(default_factory=dict, description="Output parameters from the simulation", metadata={"openbis_type": "JSON"})
    atomistic_models: List["AtomisticModel"] = Field(default_factory=list, description="List of atomistic models", metadata={"openbis_type": "PARENT"})
    codes: List[Code] = Field(default_factory=list, description="List of codes used in the band structure calculation", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class GeometryOptimisation(Simulation):
    cell_opt_constraints: str = Field(default=None, description="Constraints applied to the cell optimisation", metadata={"openbis_type": "VARCHAR"})
    cell_optimised: bool = Field(default=None, description="Whether the cell was optimised", metadata={"openbis_type": "BOOLEAN"})
    driver_code: str = Field(default=None, description="Driver code used for the simulation", metadata={"openbis_type": "VARCHAR"})
    constrained: bool = Field(default=None, description="Whether the simulation is constrained", metadata={"openbis_type": "BOOLEAN"})
    force_convergence_threshold: EnergyDensityValue = Field(default=None, description="Force convergence threshold value", metadata={"openbis_type": "JSON"})
    level_theory_method: str = Field(default=None, description="Method for the level of theory", metadata={"openbis_type": "VARCHAR"})
    level_theory_parameters: Dict = Field(default_factory=dict, description="Parameters for the level of theory method", metadata={"openbis_type": "JSON"})
    input_parameters: Dict = Field(default_factory=dict, description="Input parameters for the simulation", metadata={"openbis_type": "JSON"})
    output_parameters: Dict = Field(default_factory=dict, description="Output parameters from the simulation", metadata={"openbis_type": "JSON"})
    atomistic_models: List["AtomisticModel"] = Field(default_factory=list, description="List of atomistic models", metadata={"openbis_type": "PARENT"})
    codes: List[Code] = Field(default_factory=list, description="List of codes used", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class MeanFieldHubbard(OpenBISObject):
    charge: int = Field(default=None, description="Charge", metadata={"openbis_type": "INTEGER"})
    convergence_coefficient: float = Field(default=None, description="Convergence coefficient", metadata={"openbis_type": "REAL"})
    density_tolerance: float = Field(default=None, description="Density tolerance", metadata={"openbis_type": "REAL"})
    down_densities: List[int] = Field(default_factory=list, description="Down densities", metadata={"openbis_type": "INTEGER", "openbis_multivalue": True})
    down_electrons: List[int] = Field(default_factory=list, description="Down electrons", metadata={"openbis_type": "INTEGER", "openbis_multivalue": True})
    number_k_points: int = Field(default=None, description="Number of k-points", metadata={"openbis_type": "INTEGER"})
    on_site_hubbard_repulsion_u_value: float = Field(default=None, description="Hubbard U value", metadata={"openbis_type": "REAL"})
    on_site_hubbard_repulsion_u_unit: str = Field(default=None, description="Unit for Hubbard U", metadata={"openbis_type": "VARCHAR"})
    up_densities: List[int] = Field(default_factory=list, description="Up densities", metadata={"openbis_type": "INTEGER", "openbis_multivalue": True})
    up_electrons: List[int] = Field(default_factory=list, description="Up electrons", metadata={"openbis_type": "INTEGER", "openbis_multivalue": True})
    atomsitic_models: List[AtomisticModel] = Field(default_factory=list, description="Atomistic models", metadata={"openbis_type": "PARENT"})

class MinimumEnergyPotential(OpenBISObject):
    energies: List[EnergyValue] = Field(default_factory=list, description="List of energy values", metadata={"openbis_type": "JSON"})
    energy_barrier: EnergyValue = Field(default=None, description="Energy barrier value", metadata={"openbis_type": "JSON"})
    geometry_constraints: List[float] = Field(default_factory=list, description="Geometry constraints", metadata={"openbis_type": "REAL", "openbis_multivalue": True})
    geometry_constraints_inscrements: List[float] = Field(default_factory=list, description="Geometry constraints increments", metadata={"openbis_type": "REAL", "openbis_multivalue": True})
    method_type: Literal["NEB", "String method"] = Field(default=None, description="Method type", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    number_geometries: int = Field(default=None, description="Number of geometries", metadata={"openbis_type": "INTEGER"})
    atomistic_models: List[AtomisticModel] = Field(default_factory=list, description="List of atomistic models", metadata={"openbis_type": "PARENT"})

class PDOS(Simulation):
    level_theory_method: str = Field(default="", description="Method used for level of theory calculations", metadata={"openbis_type": "VARCHAR"})
    level_theory_parameters: Dict = Field(default_factory=dict, description="Parameters for the level of theory method", metadata={"openbis_type": "JSON"})
    input_parameters: Dict = Field(default_factory=dict, description="Input parameters for the simulation", metadata={"openbis_type": "JSON"})
    output_parameters: Dict = Field(default_factory=dict, description="Output parameters from the simulation", metadata={"openbis_type": "JSON"})
    atomistic_models: List["AtomisticModel"] = Field(default_factory=list, description="List of atomistic models", metadata={"openbis_type": "PARENT"})
    codes: List[Code] = Field(default_factory=list, description="List of codes used in the PDOS calculation", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class UnclassifiedSimulation(Simulation):
    input_parameters: Dict = Field(default_factory=dict, description="Input parameters for the simulation", metadata={"openbis_type": "JSON"})
    output_parameters: Dict = Field(default_factory=dict, description="Output parameters from the simulation", metadata={"openbis_type": "JSON"})
    atomistic_models: List["AtomisticModel"] = Field(default_factory=list, description="List of atomistic models", metadata={"openbis_type": "PARENT"})
    codes: List[Code] = Field(default_factory=list, description="List of codes used in the simulation", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class VibrationalSpectroscopy(Simulation):
    level_theory_method: str = Field(default="", description="Method used for level of theory calculations", metadata={"openbis_type": "VARCHAR"})
    level_theory_parameters: Dict = Field(default_factory=dict, description="Parameters for the level of theory method", metadata={"openbis_type": "JSON"})
    input_parameters: Dict = Field(default_factory=dict, description="Input parameters for the simulation", metadata={"openbis_type": "JSON"})
    output_parameters: Dict = Field(default_factory=dict, description="Output parameters from the simulation", metadata={"openbis_type": "JSON"})
    atomistic_models: List["AtomisticModel"] = Field(default_factory=list, description="List of atomistic models", metadata={"openbis_type": "PARENT"})
    codes: List[Code] = Field(default_factory=list, description="List of codes used in the vibrational spectroscopy calculation", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class PotentialEnergyCalculation(Simulation):
    energy: EnergyValue = Field(default=None, description="Energy value", metadata={"openbis_type": "JSON"})
    level_theory_method: str = Field(default="", description="Method used for level of theory calculations", metadata={"openbis_type": "VARCHAR"})
    level_theory_parameters: Dict = Field(default_factory=dict, description="Parameters for the level of theory method", metadata={"openbis_type": "JSON"})
    input_parameters: Dict = Field(default_factory=dict, description="Input parameters for the simulation", metadata={"openbis_type": "JSON"})
    output_parameters: Dict = Field(default_factory=dict, description="Output parameters from the simulation", metadata={"openbis_type": "JSON"})
    atomistic_models: List["AtomisticModel"] = Field(default_factory=list, description="List of atomistic models", metadata={"openbis_type": "PARENT"})
    codes: List[Code] = Field(default_factory=list, description="List of codes used in the potential energy calculation", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class Action(OpenBISObject):
    duration: str = Field(default="0 00:00:00", description="Duration of the action in format 'days hh:mm:ss'", metadata={"openbis_type": "VARCHAR"})
    component: Component = Field(default=None, description="Component used in the action", metadata={"openbis_type": "SAMPLE"})
    component_settings: Dict = Field(default_factory=dict, description="Settings for the component. The keys are the names of the component properties that are used in the specific observable type and the values are the properties values.", metadata={"openbis_type": "JSON"})
    
    @field_validator('duration')
    @classmethod
    def validate_duration_format(cls, v):
        # Pattern for "days hh:mm:ss" format
        pattern = r'^\d+\s+([01]?\d|2[0-3]):([0-5]?\d):([0-5]?\d)$'
        if not re.match(pattern, v):
            raise ValueError('Duration must be in format "days hh:mm:ss" (e.g., "1 12:30:45")')
        return v
    
class Annealing(Action):
    target_temperature: TemperatureValue = Field(default=None, description="Target temperature for the annealing process", metadata={"openbis_type": "JSON"})

class Coating(Action):
    pass

class Cooldown(Action):
    target_temperature: TemperatureValue = Field(default=None, description="Target temperature for the cooldown process", metadata={"openbis_type": "JSON"})
    cryogen: str = Field(default="", description="Name of the used cryogen", metadata={"openbis_type": "VARCHAR"})

class Delamination(Action):
    pass

class Deposition(Action):
    substrate_temperature: TemperatureValue = Field(default=None, description="Temperature of the substrate during deposition", metadata={"openbis_type": "JSON"})
    substance: Substance = Field(default=None, description="Substance being deposited", metadata={"openbis_type": "SAMPLE"})

class Dosing(Action):
    dosing_gas: Substance = Field(default=None, description="Gas used for dosing", metadata={"openbis_type": "SAMPLE"})
    pressure: PressureValue = Field(default=None, description="Pressure during dosing", metadata={"openbis_type": "JSON"})
    substrate_temperature: TemperatureValue = Field(default=None, description="Substrate temperature during dosing", metadata={"openbis_type": "JSON"})

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
    sputter_ion: str = Field(default=None, description="Ion used for sputtering", metadata={"openbis_type": "VARCHAR"})
    pressure: PressureValue = Field(default=None, description="Pressure during sputtering", metadata={"openbis_type": "JSON"})
    current: CurrencyValue = Field(default=None, description="Current applied during sputtering", metadata={"openbis_type": "JSON"})
    angle: AngleValue = Field(default=None, description="Angle of sputtering", metadata={"openbis_type": "JSON"})
    substrate_temperature: TemperatureValue = Field(default=None, description="Substrate temperature during sputtering", metadata={"openbis_type": "JSON"})

class Location(OpenBISObject):
    organisation: Organisation = Field(default=None, description="Organisation associated with the location", metadata={"openbis_type": "SAMPLE"})

class Analysis(OpenBISObject):
    codes: List[Code] = Field(default_factory=list, description="List of codes used in the analysis", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    measurements: List["MeasurementSession"] = Field(default_factory=list, description="List of measurement sessions used in the analysis", metadata={"openbis_type": "PARENT"})
    software: List[Software] = Field(default_factory=list, description="List of software used in the analysis", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

class Result(OpenBISObject):
    analysis: List[Analysis] = Field(default_factory=list, description="List of analysis", metadata={"openbis_type": "PARENT"})
    simulations: List[Union[BandStructure, GeometryOptimisation, PDOS, VibrationalSpectroscopy]] = Field(default_factory=list, description="List of simulations", metadata={"openbis_type": "PARENT"})
    measurements: List["MeasurementSession"] = Field(default_factory=list, description="List of associated measurement sessions", metadata={"openbis_type": "PARENT"})

class Draft(OpenBISObject):
    draft_type: DraftTypeEnum = Field(default=None, description="Type of the draft", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    results: List[Result] = Field(default_factory=list, description="List of associated results", metadata={"openbis_type": "PARENT"})

class AiidaNode(OpenBISObject):
    pass

class Observable(OpenBISObject):
    channel_name: str = Field(default=None, description="Name of the reading channel", metadata={"openbis_type": "VARCHAR"})
    component: Component = Field(default=None, description="Used component", metadata={"openbis_type": "SAMPLE"})
    component_settings: Dict = Field(default_factory=dict, description="Settings for the component. The keys are the names of the component properties that are used in the specific observable type and the values are the properties values.", metadata={"openbis_type": "JSON"})

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
    measurement_folder_path: str = Field(default=None, description="Path to the measurement folder", metadata={"openbis_type": "VARCHAR"})
    default_object_view: Literal["IMAGING_GALLERY_VIEW"] = Field(default=None, description="Default object view", metadata={"openbis_type": "CONTROLLEDVOCABULARY"})
    measurement_session: "MeasurementSession" = Field(default=None, description="Reference to another measurement session", metadata={"openbis_type": "PARENT"})

class Process(OpenBISObject):
    short_name: str = Field(default=None, description="Short name for the process", metadata={"openbis_type": "VARCHAR"})
    instrument: Union[Instrument, InstrumentSTM] = Field(default=None, description="Instrument used in the process", metadata={"openbis_type": "SAMPLE"})
    process_steps_settings: List[ProcessStepSettings] = Field(default_factory=list, description="Settings for each process step", metadata={"openbis_type": "JSON"})

class ProcessStep(OpenBISObject):
    actions: List[Action] = Field(default_factory=list, description="List of actions performed in this step", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    observables: List[Observable] = Field(default_factory=list, description="List of observables recorded in this step", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    instrument: Union[Instrument, InstrumentSTM] = Field(default=None, description="Instrument used for this step", metadata={"openbis_type": "SAMPLE"})

class Preparation(OpenBISObject):
    process_steps: List[ProcessStep] = Field(default_factory=list, description="List of process steps", metadata={"openbis_type": "CHILDREN"})

class Publication(OpenBISObject):
    abstract: str = Field(default=None, description="Abstract of the publication", metadata={"openbis_type": "VARCHAR"})
    doi: str = Field(default=None, description="DOI of the publication", metadata={"openbis_type": "VARCHAR"})
    url: str = Field(default=None, description="URL of the publication", metadata={"openbis_type": "VARCHAR"})
    dataset_url: str = Field(default=None, description="URL to the associated dataset", metadata={"openbis_type": "VARCHAR"})
    year: int = Field(default=None, description="Year of publication", metadata={"openbis_type": "INTEGER"})
    authors: List[Person] = Field(default_factory=list, description="List of authors", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})
    drafts: List[Draft] = Field(default_factory=list, description="Related drafts", metadata={"openbis_type": "PARENT"})
    grants: List[Grant] = Field(default_factory=list, description="Grants supporting the publication", metadata={"openbis_type": "SAMPLE", "openbis_multivalue": True})

import json
atomistic_model = AtomisticModel(name = "aa")
geoopt = GeometryOptimisation(atomistic_models = [atomistic_model])
schema = GeometryOptimisation.model_json_schema()
with open("simulation_schema.json", "w") as f:
    json.dump(schema, f, indent=4)