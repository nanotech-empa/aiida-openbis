# Auto generated from materialMLinfo.yaml by pythongen.py version: 0.0.1
# Generation date: 2024-07-19T16:01:21
# Schema: materialsML
#
# id: https://w3id.org/linkml/examples/materialsML
# description:
# license: https://creativecommons.org/publicdomain/zero/1.0/

import dataclasses
import re
from jsonasobj2 import JsonObj, as_dict
from typing import Optional, List, Union, Dict, ClassVar, Any
from dataclasses import dataclass
from datetime import date, datetime
from linkml_runtime.linkml_model.meta import EnumDefinition, PermissibleValue, PvFormulaOptions

from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.metamodelcore import empty_list, empty_dict, bnode
from linkml_runtime.utils.yamlutils import YAMLRoot, extended_str, extended_float, extended_int
from linkml_runtime.utils.dataclass_extensions_376 import dataclasses_init_fn_with_kwargs
from linkml_runtime.utils.formatutils import camelcase, underscore, sfx
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from rdflib import Namespace, URIRef
from linkml_runtime.utils.curienamespace import CurieNamespace
from linkml_runtime.linkml_model.types import Boolean, Date, Double, Integer, String
from linkml_runtime.utils.metamodelcore import Bool, XSDDate

metamodel_version = "1.7.0"
version = None

# Overwrite dataclasses _init_fn to add **kwargs in __init__
dataclasses._init_fn = dataclasses_init_fn_with_kwargs

# Namespaces
AFE = CurieNamespace('afe', 'http://purl.allotrope.org/ontologies/equipment#')
APOLLO = CurieNamespace('apollo', 'http://purl.obolibrary.org/obo/apollo_sv.owl#')
CHEBI = CurieNamespace('chebi', 'http://purl.obolibrary.org/obo/chebi.owl#')
CHMO = CurieNamespace('chmo', 'http://purl.obolibrary.org/obo/chmo.owl#')
DCAT = CurieNamespace('dcat', 'http://www.w3.org/ns/dcat#')
EMMO = CurieNamespace('emmo', 'https://w3id.org/emmo#')
EMMO_CHAMEO = CurieNamespace('emmo_chameo', 'https://w3id.org/emmo-chameo/chameo#')
EMMO_DA = CurieNamespace('emmo_da', 'http://emmo.info/domain-atomistic#')
EMMO_ELECTROCHEMISTRY = CurieNamespace('emmo_electrochemistry', 'https://w3id.org/emmo/domain/electrochemistry#')
EMMO_NANO = CurieNamespace('emmo_nano', 'http://emmo.info/emmo/nanoind#')
FAIRMAT_STS = CurieNamespace('fairmat_sts', 'https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXsts.html#')
FOAF = CurieNamespace('foaf', 'http://xmlns.com/foaf/0.1/')
IAO = CurieNamespace('iao', 'http://purl.obolibrary.org/obo/iao.owl#')
LINKML = CurieNamespace('linkml', 'https://w3id.org/linkml/')
MATERIALSML = CurieNamespace('materialsML', 'https://w3id.org/linkml/examples/materialsML/')
NCIT = CurieNamespace('ncit', 'http://purl.obolibrary.org/obo/ncit.owl#')
OBI = CurieNamespace('obi', 'http://purl.obolibrary.org/obo/obi.owl#')
OCCO = CurieNamespace('occo', 'http://purl.obolibrary.org/obo/occo.owl#')
OME_ONTO = CurieNamespace('ome_onto', 'http://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2016-06/ome_xsd.html#')
OPMI = CurieNamespace('opmi', 'http://purl.obolibrary.org/obo/opmi.owl#')
QUDT = CurieNamespace('qudt', 'http://qudt.org/schema/qudt/')
SCHEMA = CurieNamespace('schema', 'http://schema.org/')
SIO = CurieNamespace('sio', 'http://semanticscience.org/ontology/sio.owl#')
UNIT = CurieNamespace('unit', 'http://qudt.org/vocab/unit/')
XSD = CurieNamespace('xsd', 'http://www.w3.org/2001/XMLSchema#')
DEFAULT_ = MATERIALSML


# Types
class CustomDatetime(str):
    """ custom datetime. It replaces linkML datetime because it was not possible to use it to validate a data file and to generate other schema file formats. For validating data it required datetime format while for converting it had be a string. """
    type_class_uri = XSD["dateTime"]
    type_class_curie = "xsd:dateTime"
    type_name = "custom_datetime"
    type_model_uri = MATERIALSML.CustomDatetime


# Class references



@dataclass
class AFM(YAMLRoot):
    """
    Atomic Force Microscopy
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = CHMO["CHMO_0000113"]
    class_class_curie: ClassVar[str] = "chmo:CHMO_0000113"
    class_name: ClassVar[str] = "AFM"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.AFM

    perm_id: Optional[str] = None
    name: Optional[str] = None
    start_time: Optional[str] = None
    duration: Optional[Union[dict, "QuantityValue"]] = None
    bias_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    bias_calibration_factor: Optional[Union[dict, "QuantityValue"]] = None
    bias_calibration_offset: Optional[Union[dict, "QuantityValue"]] = None
    current_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    current_calibration_factor: Optional[Union[dict, "QuantityValue"]] = None
    current_calibration_offset: Optional[Union[dict, "QuantityValue"]] = None
    current_gain: Optional[Union[dict, "QuantityValue"]] = None
    z_position: Optional[Union[dict, "QuantityValue"]] = None
    feedback_active: Optional[Union[bool, Bool]] = None
    feedback_type: Optional[Union[str, "FeedbackTypeEnum"]] = None
    z_controller_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_p_gain: Optional[float] = None
    z_controller_i_gain: Optional[float] = None
    z_controller_time_constant: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_tip_lift: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_switch_off_delay: Optional[Union[dict, "QuantityValue"]] = None
    piezo_configuration_settings: Optional[Union[dict, "PiezoConfigurationSettings"]] = None
    scan_settings: Optional[Union[dict, "ScanSettings"]] = None
    dwell_time: Optional[Union[dict, "QuantityValue"]] = None
    sample_temperature: Optional[Union[dict, "QuantityValue"]] = None
    recording_temperature: Optional[Union[dict, "QuantityValue"]] = None
    oscillation_control_settings: Optional[Union[dict, "OscillationControlSettings"]] = None
    comments: Optional[str] = None
    scan_dx: Optional[Union[dict, "QuantityValue"]] = None
    scan_z_min: Optional[Union[dict, "QuantityValue"]] = None
    scan_z_max: Optional[Union[dict, "QuantityValue"]] = None
    amplitude: Optional[Union[dict, "QuantityValue"]] = None
    resonance_frequency: Optional[Union[dict, "QuantityValue"]] = None
    p_model: Optional[Union[str, "PModelEnum"]] = None
    wfms_uuid: Optional[str] = None
    annealing: Optional[Union[dict, "Annealing"]] = None
    deposition: Optional[Union[dict, "Deposition"]] = None
    dft: Optional[Union[dict, "DFT"]] = None
    tb: Optional[Union[dict, "TightBinding"]] = None
    instrument: Optional[Union[dict, "Instrument"]] = None
    sample: Optional[Union[dict, "Sample"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.start_time is not None and not isinstance(self.start_time, str):
            self.start_time = str(self.start_time)

        if self.duration is not None and not isinstance(self.duration, QuantityValue):
            self.duration = QuantityValue(**as_dict(self.duration))

        if self.bias_setpoint is not None and not isinstance(self.bias_setpoint, QuantityValue):
            self.bias_setpoint = QuantityValue(**as_dict(self.bias_setpoint))

        if self.bias_calibration_factor is not None and not isinstance(self.bias_calibration_factor, QuantityValue):
            self.bias_calibration_factor = QuantityValue(**as_dict(self.bias_calibration_factor))

        if self.bias_calibration_offset is not None and not isinstance(self.bias_calibration_offset, QuantityValue):
            self.bias_calibration_offset = QuantityValue(**as_dict(self.bias_calibration_offset))

        if self.current_setpoint is not None and not isinstance(self.current_setpoint, QuantityValue):
            self.current_setpoint = QuantityValue(**as_dict(self.current_setpoint))

        if self.current_calibration_factor is not None and not isinstance(self.current_calibration_factor, QuantityValue):
            self.current_calibration_factor = QuantityValue(**as_dict(self.current_calibration_factor))

        if self.current_calibration_offset is not None and not isinstance(self.current_calibration_offset, QuantityValue):
            self.current_calibration_offset = QuantityValue(**as_dict(self.current_calibration_offset))

        if self.current_gain is not None and not isinstance(self.current_gain, QuantityValue):
            self.current_gain = QuantityValue(**as_dict(self.current_gain))

        if self.z_position is not None and not isinstance(self.z_position, QuantityValue):
            self.z_position = QuantityValue(**as_dict(self.z_position))

        if self.feedback_active is not None and not isinstance(self.feedback_active, Bool):
            self.feedback_active = Bool(self.feedback_active)

        if self.feedback_type is not None and not isinstance(self.feedback_type, FeedbackTypeEnum):
            self.feedback_type = FeedbackTypeEnum(self.feedback_type)

        if self.z_controller_setpoint is not None and not isinstance(self.z_controller_setpoint, QuantityValue):
            self.z_controller_setpoint = QuantityValue(**as_dict(self.z_controller_setpoint))

        if self.z_controller_p_gain is not None and not isinstance(self.z_controller_p_gain, float):
            self.z_controller_p_gain = float(self.z_controller_p_gain)

        if self.z_controller_i_gain is not None and not isinstance(self.z_controller_i_gain, float):
            self.z_controller_i_gain = float(self.z_controller_i_gain)

        if self.z_controller_time_constant is not None and not isinstance(self.z_controller_time_constant, QuantityValue):
            self.z_controller_time_constant = QuantityValue(**as_dict(self.z_controller_time_constant))

        if self.z_controller_tip_lift is not None and not isinstance(self.z_controller_tip_lift, QuantityValue):
            self.z_controller_tip_lift = QuantityValue(**as_dict(self.z_controller_tip_lift))

        if self.z_controller_switch_off_delay is not None and not isinstance(self.z_controller_switch_off_delay, QuantityValue):
            self.z_controller_switch_off_delay = QuantityValue(**as_dict(self.z_controller_switch_off_delay))

        if self.piezo_configuration_settings is not None and not isinstance(self.piezo_configuration_settings, PiezoConfigurationSettings):
            self.piezo_configuration_settings = PiezoConfigurationSettings(**as_dict(self.piezo_configuration_settings))

        if self.scan_settings is not None and not isinstance(self.scan_settings, ScanSettings):
            self.scan_settings = ScanSettings(**as_dict(self.scan_settings))

        if self.dwell_time is not None and not isinstance(self.dwell_time, QuantityValue):
            self.dwell_time = QuantityValue(**as_dict(self.dwell_time))

        if self.sample_temperature is not None and not isinstance(self.sample_temperature, QuantityValue):
            self.sample_temperature = QuantityValue(**as_dict(self.sample_temperature))

        if self.recording_temperature is not None and not isinstance(self.recording_temperature, QuantityValue):
            self.recording_temperature = QuantityValue(**as_dict(self.recording_temperature))

        if self.oscillation_control_settings is not None and not isinstance(self.oscillation_control_settings, OscillationControlSettings):
            self.oscillation_control_settings = OscillationControlSettings(**as_dict(self.oscillation_control_settings))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.scan_dx is not None and not isinstance(self.scan_dx, QuantityValue):
            self.scan_dx = QuantityValue(**as_dict(self.scan_dx))

        if self.scan_z_min is not None and not isinstance(self.scan_z_min, QuantityValue):
            self.scan_z_min = QuantityValue(**as_dict(self.scan_z_min))

        if self.scan_z_max is not None and not isinstance(self.scan_z_max, QuantityValue):
            self.scan_z_max = QuantityValue(**as_dict(self.scan_z_max))

        if self.amplitude is not None and not isinstance(self.amplitude, QuantityValue):
            self.amplitude = QuantityValue(**as_dict(self.amplitude))

        if self.resonance_frequency is not None and not isinstance(self.resonance_frequency, QuantityValue):
            self.resonance_frequency = QuantityValue(**as_dict(self.resonance_frequency))

        if self.p_model is not None and not isinstance(self.p_model, PModelEnum):
            self.p_model = PModelEnum(self.p_model)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.annealing is not None and not isinstance(self.annealing, Annealing):
            self.annealing = Annealing(**as_dict(self.annealing))

        if self.deposition is not None and not isinstance(self.deposition, Deposition):
            self.deposition = Deposition(**as_dict(self.deposition))

        if self.dft is not None and not isinstance(self.dft, DFT):
            self.dft = DFT(**as_dict(self.dft))

        if self.tb is not None and not isinstance(self.tb, TightBinding):
            self.tb = TightBinding(**as_dict(self.tb))

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        if self.sample is not None and not isinstance(self.sample, Sample):
            self.sample = Sample(**as_dict(self.sample))

        super().__post_init__(**kwargs)


@dataclass
class Annealing(YAMLRoot):
    """
    Annealing
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = CHMO["CHMO_0001465"]
    class_class_curie: ClassVar[str] = "chmo:CHMO_0001465"
    class_name: ClassVar[str] = "Annealing"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Annealing

    perm_id: Optional[str] = None
    name: Optional[str] = None
    duration: Optional[Union[dict, "QuantityValue"]] = None
    pressure: Optional[Union[dict, "QuantityValue"]] = None
    voltage: Optional[Union[dict, "QuantityValue"]] = None
    temperature: Optional[Union[dict, "QuantityValue"]] = None
    current: Optional[Union[dict, "QuantityValue"]] = None
    comments: Optional[str] = None
    sputtering: Optional[Union[dict, "Sputtering"]] = None
    deposition: Optional[Union[dict, "Deposition"]] = None
    uhv_components: Optional[Union[Union[dict, "UHVComponent"], List[Union[dict, "UHVComponent"]]]] = empty_list()
    instrument: Optional[Union[dict, "Instrument"]] = None
    sample: Optional[Union[dict, "Sample"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.duration is not None and not isinstance(self.duration, QuantityValue):
            self.duration = QuantityValue(**as_dict(self.duration))

        if self.pressure is not None and not isinstance(self.pressure, QuantityValue):
            self.pressure = QuantityValue(**as_dict(self.pressure))

        if self.voltage is not None and not isinstance(self.voltage, QuantityValue):
            self.voltage = QuantityValue(**as_dict(self.voltage))

        if self.temperature is not None and not isinstance(self.temperature, QuantityValue):
            self.temperature = QuantityValue(**as_dict(self.temperature))

        if self.current is not None and not isinstance(self.current, QuantityValue):
            self.current = QuantityValue(**as_dict(self.current))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.sputtering is not None and not isinstance(self.sputtering, Sputtering):
            self.sputtering = Sputtering(**as_dict(self.sputtering))

        if self.deposition is not None and not isinstance(self.deposition, Deposition):
            self.deposition = Deposition(**as_dict(self.deposition))

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        if self.sample is not None and not isinstance(self.sample, Sample):
            self.sample = Sample(**as_dict(self.sample))

        super().__post_init__(**kwargs)


@dataclass
class AtomisticModel(YAMLRoot):
    """
    Atomistic Model
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = EMMO["EMMO_84cadc45_6758_46f2_ba2a_5ead65c70213"]
    class_class_curie: ClassVar[str] = "emmo:EMMO_84cadc45_6758_46f2_ba2a_5ead65c70213"
    class_name: ClassVar[str] = "AtomisticModel"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.AtomisticModel

    perm_id: Optional[str] = None
    name: Optional[str] = None
    atoms_positions: Optional[Union[Union[dict, "AtomsPositions"], List[Union[dict, "AtomsPositions"]]]] = empty_list()
    comments: Optional[str] = None
    cell_vectors: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    periodic_boundary_conditions: Optional[Union[dict, "PeriodicBoundaryConditions"]] = None
    optimised: Optional[Union[bool, Bool]] = None
    wfms_uuid: Optional[str] = None
    crystal: Optional[Union[dict, "Crystal"]] = None
    molecule: Optional[Union[dict, "Molecule"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if not isinstance(self.atoms_positions, list):
            self.atoms_positions = [self.atoms_positions] if self.atoms_positions is not None else []
        self.atoms_positions = [v if isinstance(v, AtomsPositions) else AtomsPositions(**as_dict(v)) for v in self.atoms_positions]

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.cell_vectors, list):
            self.cell_vectors = [self.cell_vectors] if self.cell_vectors is not None else []
        self.cell_vectors = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.cell_vectors]

        if self.periodic_boundary_conditions is not None and not isinstance(self.periodic_boundary_conditions, PeriodicBoundaryConditions):
            self.periodic_boundary_conditions = PeriodicBoundaryConditions(**as_dict(self.periodic_boundary_conditions))

        if self.optimised is not None and not isinstance(self.optimised, Bool):
            self.optimised = Bool(self.optimised)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.crystal is not None and not isinstance(self.crystal, Crystal):
            self.crystal = Crystal(**as_dict(self.crystal))

        if self.molecule is not None and not isinstance(self.molecule, Molecule):
            self.molecule = Molecule(**as_dict(self.molecule))

        super().__post_init__(**kwargs)


@dataclass
class Author(YAMLRoot):
    """
    Author
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C42781"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C42781"
    class_name: ClassVar[str] = "Author"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Author

    perm_id: Optional[str] = None
    name: Optional[str] = None
    person: Optional[Union[dict, "Person"]] = None
    institution: Optional[Union[dict, "Institution"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.person is not None and not isinstance(self.person, Person):
            self.person = Person(**as_dict(self.person))

        if self.institution is not None and not isinstance(self.institution, Institution):
            self.institution = Institution(**as_dict(self.institution))

        super().__post_init__(**kwargs)


@dataclass
class BandStructure(YAMLRoot):
    """
    Band Structure
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["BandStructure"]
    class_class_curie: ClassVar[str] = "materialsML:BandStructure"
    class_name: ClassVar[str] = "BandStructure"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.BandStructure

    perm_id: Optional[str] = None
    name: Optional[str] = None
    k_points_conditions: Optional[Union[dict, "KPointsConditions"]] = None
    bs_energies: Optional[Union[dict, "BSEnergy"]] = None
    wfms_uuid: Optional[str] = None
    wfms_url: Optional[str] = None
    comments: Optional[str] = None
    atomistic_models: Optional[Union[Union[dict, AtomisticModel], List[Union[dict, AtomisticModel]]]] = empty_list()
    dft: Optional[Union[dict, "DFT"]] = None
    tb: Optional[Union[dict, "TightBinding"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.k_points_conditions is not None and not isinstance(self.k_points_conditions, KPointsConditions):
            self.k_points_conditions = KPointsConditions(**as_dict(self.k_points_conditions))

        if self.bs_energies is not None and not isinstance(self.bs_energies, BSEnergy):
            self.bs_energies = BSEnergy(**as_dict(self.bs_energies))

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.wfms_url is not None and not isinstance(self.wfms_url, str):
            self.wfms_url = str(self.wfms_url)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.atomistic_models, list):
            self.atomistic_models = [self.atomistic_models] if self.atomistic_models is not None else []
        self.atomistic_models = [v if isinstance(v, AtomisticModel) else AtomisticModel(**as_dict(v)) for v in self.atomistic_models]

        if self.dft is not None and not isinstance(self.dft, DFT):
            self.dft = DFT(**as_dict(self.dft))

        if self.tb is not None and not isinstance(self.tb, TightBinding):
            self.tb = TightBinding(**as_dict(self.tb))

        super().__post_init__(**kwargs)


@dataclass
class CanGasBottle(YAMLRoot):
    """
    Bottle/Can of Gas
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["CanGasBottle"]
    class_class_curie: ClassVar[str] = "materialsML:CanGasBottle"
    class_name: ClassVar[str] = "CanGasBottle"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.CanGasBottle

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    pressure: Optional[Union[dict, "QuantityValue"]] = None
    state: Optional[str] = None
    supplier: Optional[Union[dict, "Supplier"]] = None
    comments: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.pressure is not None and not isinstance(self.pressure, QuantityValue):
            self.pressure = QuantityValue(**as_dict(self.pressure))

        if self.state is not None and not isinstance(self.state, str):
            self.state = str(self.state)

        if self.supplier is not None and not isinstance(self.supplier, Supplier):
            self.supplier = Supplier(**as_dict(self.supplier))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        super().__post_init__(**kwargs)


@dataclass
class Chemist(YAMLRoot):
    """
    Chemist/Synthesiser
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = OCCO["OCCO_19203100"]
    class_class_curie: ClassVar[str] = "occo:OCCO_19203100"
    class_name: ClassVar[str] = "Chemist"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Chemist

    perm_id: Optional[str] = None
    name: Optional[str] = None
    person: Optional[Union[dict, "Person"]] = None
    supplier: Optional[Union[dict, "Supplier"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.person is not None and not isinstance(self.person, Person):
            self.person = Person(**as_dict(self.person))

        if self.supplier is not None and not isinstance(self.supplier, Supplier):
            self.supplier = Supplier(**as_dict(self.supplier))

        super().__post_init__(**kwargs)


@dataclass
class Crystal(YAMLRoot):
    """
    Crystal
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = EMMO["EMMO_0bb3b434_73aa_428f_b4e8_2a2468648e19"]
    class_class_curie: ClassVar[str] = "emmo:EMMO_0bb3b434_73aa_428f_b4e8_2a2468648e19"
    class_name: ClassVar[str] = "Crystal"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Crystal

    perm_id: Optional[str] = None
    name: Optional[str] = None
    material: Optional[str] = None
    face: Optional[str] = None
    sample_plate: Optional[str] = None
    diameter: Optional[Union[dict, "QuantityValue"]] = None
    height: Optional[Union[dict, "QuantityValue"]] = None
    hazardous: Optional[Union[bool, Bool]] = None
    hazardous_specification: Optional[str] = None
    fridge: Optional[Union[bool, Bool]] = None
    no_light: Optional[Union[bool, Bool]] = None
    dry: Optional[Union[bool, Bool]] = None
    no_oxygen: Optional[Union[bool, Bool]] = None
    other_storage_condition: Optional[Union[bool, Bool]] = None
    other_storage_condition_specification: Optional[str] = None
    specifications: Optional[str] = None
    receive_date: Optional[Union[str, XSDDate]] = None
    reference_number: Optional[str] = None
    comments: Optional[str] = None
    storage: Optional[Union[dict, "Storage"]] = None
    supplier: Optional[Union[dict, "Supplier"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.material is not None and not isinstance(self.material, str):
            self.material = str(self.material)

        if self.face is not None and not isinstance(self.face, str):
            self.face = str(self.face)

        if self.sample_plate is not None and not isinstance(self.sample_plate, str):
            self.sample_plate = str(self.sample_plate)

        if self.diameter is not None and not isinstance(self.diameter, QuantityValue):
            self.diameter = QuantityValue(**as_dict(self.diameter))

        if self.height is not None and not isinstance(self.height, QuantityValue):
            self.height = QuantityValue(**as_dict(self.height))

        if self.hazardous is not None and not isinstance(self.hazardous, Bool):
            self.hazardous = Bool(self.hazardous)

        if self.hazardous_specification is not None and not isinstance(self.hazardous_specification, str):
            self.hazardous_specification = str(self.hazardous_specification)

        if self.fridge is not None and not isinstance(self.fridge, Bool):
            self.fridge = Bool(self.fridge)

        if self.no_light is not None and not isinstance(self.no_light, Bool):
            self.no_light = Bool(self.no_light)

        if self.dry is not None and not isinstance(self.dry, Bool):
            self.dry = Bool(self.dry)

        if self.no_oxygen is not None and not isinstance(self.no_oxygen, Bool):
            self.no_oxygen = Bool(self.no_oxygen)

        if self.other_storage_condition is not None and not isinstance(self.other_storage_condition, Bool):
            self.other_storage_condition = Bool(self.other_storage_condition)

        if self.other_storage_condition_specification is not None and not isinstance(self.other_storage_condition_specification, str):
            self.other_storage_condition_specification = str(self.other_storage_condition_specification)

        if self.specifications is not None and not isinstance(self.specifications, str):
            self.specifications = str(self.specifications)

        if self.receive_date is not None and not isinstance(self.receive_date, XSDDate):
            self.receive_date = XSDDate(self.receive_date)

        if self.reference_number is not None and not isinstance(self.reference_number, str):
            self.reference_number = str(self.reference_number)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.storage is not None and not isinstance(self.storage, Storage):
            self.storage = Storage(**as_dict(self.storage))

        if self.supplier is not None and not isinstance(self.supplier, Supplier):
            self.supplier = Supplier(**as_dict(self.supplier))

        super().__post_init__(**kwargs)


@dataclass
class Deposition(YAMLRoot):
    """
    Deposition
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = CHMO["CHMO_0001310"]
    class_class_curie: ClassVar[str] = "chmo:CHMO_0001310"
    class_name: ClassVar[str] = "Deposition"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Deposition

    perm_id: Optional[str] = None
    name: Optional[str] = None
    stabilisation_time: Optional[Union[dict, "QuantityValue"]] = None
    deposition_time: Optional[Union[dict, "QuantityValue"]] = None
    pressure: Optional[Union[dict, "QuantityValue"]] = None
    substrate_temperature: Optional[Union[dict, "QuantityValue"]] = None
    molecule_temperature: Optional[Union[dict, "QuantityValue"]] = None
    evaporator_slot: Optional[Union[dict, "EvaporatorSlot"]] = None
    comments: Optional[str] = None
    annealing: Optional[Union[dict, Annealing]] = None
    molecule: Optional[Union[dict, "Molecule"]] = None
    uhv_components: Optional[Union[Union[dict, "UHVComponent"], List[Union[dict, "UHVComponent"]]]] = empty_list()
    instrument: Optional[Union[dict, "Instrument"]] = None
    sample: Optional[Union[dict, "Sample"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.stabilisation_time is not None and not isinstance(self.stabilisation_time, QuantityValue):
            self.stabilisation_time = QuantityValue(**as_dict(self.stabilisation_time))

        if self.deposition_time is not None and not isinstance(self.deposition_time, QuantityValue):
            self.deposition_time = QuantityValue(**as_dict(self.deposition_time))

        if self.pressure is not None and not isinstance(self.pressure, QuantityValue):
            self.pressure = QuantityValue(**as_dict(self.pressure))

        if self.substrate_temperature is not None and not isinstance(self.substrate_temperature, QuantityValue):
            self.substrate_temperature = QuantityValue(**as_dict(self.substrate_temperature))

        if self.molecule_temperature is not None and not isinstance(self.molecule_temperature, QuantityValue):
            self.molecule_temperature = QuantityValue(**as_dict(self.molecule_temperature))

        if self.evaporator_slot is not None and not isinstance(self.evaporator_slot, EvaporatorSlot):
            self.evaporator_slot = EvaporatorSlot(**as_dict(self.evaporator_slot))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.annealing is not None and not isinstance(self.annealing, Annealing):
            self.annealing = Annealing(**as_dict(self.annealing))

        if self.molecule is not None and not isinstance(self.molecule, Molecule):
            self.molecule = Molecule(**as_dict(self.molecule))

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        if self.sample is not None and not isinstance(self.sample, Sample):
            self.sample = Sample(**as_dict(self.sample))

        super().__post_init__(**kwargs)


@dataclass
class DFT(YAMLRoot):
    """
    Density Functional Theory
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = EMMO_DA["83a69cf2-00c1-58b8-915b-76cb3549890a"]
    class_class_curie: ClassVar[str] = "emmo_da:83a69cf2-00c1-58b8-915b-76cb3549890a"
    class_name: ClassVar[str] = "DFT"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.DFT

    perm_id: Optional[str] = None
    name: Optional[str] = None
    wfms_uuid: Optional[str] = None
    wfms_url: Optional[str] = None
    fermi_energy: Optional[float] = None
    scf_convergence_threshold: Optional[Union[dict, "QuantityValue"]] = None
    xc_functional: Optional[Union[str, "XCFunctionalEnum"]] = None
    spin_polarised_calculation: Optional[Union[bool, Bool]] = None
    initial_spin_guesses: Optional[Union[float, List[float]]] = empty_list()
    spin_multiplicity: Optional[Union[bool, Bool]] = None
    absolute_magnetisation: Optional[Union[dict, "QuantityValue"]] = None
    total_magnetisation: Optional[Union[dict, "QuantityValue"]] = None
    total_spin_squared: Optional[float] = None
    net_charge: Optional[float] = None
    smearing: Optional[Union[dict, "QuantityValue"]] = None
    fermi_dirac_temperature: Optional[Union[dict, "QuantityValue"]] = None
    force_multiplicity: Optional[Union[bool, Bool]] = None
    output_mulliken_population_analysis: Optional[Union[dict, "OutputPopulationAnalysis"]] = None
    output_hirshfeld_population_analysis: Optional[Union[dict, "OutputPopulationAnalysis"]] = None
    comments: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.wfms_url is not None and not isinstance(self.wfms_url, str):
            self.wfms_url = str(self.wfms_url)

        if self.fermi_energy is not None and not isinstance(self.fermi_energy, float):
            self.fermi_energy = float(self.fermi_energy)

        if self.scf_convergence_threshold is not None and not isinstance(self.scf_convergence_threshold, QuantityValue):
            self.scf_convergence_threshold = QuantityValue(**as_dict(self.scf_convergence_threshold))

        if self.xc_functional is not None and not isinstance(self.xc_functional, XCFunctionalEnum):
            self.xc_functional = XCFunctionalEnum(self.xc_functional)

        if self.spin_polarised_calculation is not None and not isinstance(self.spin_polarised_calculation, Bool):
            self.spin_polarised_calculation = Bool(self.spin_polarised_calculation)

        if not isinstance(self.initial_spin_guesses, list):
            self.initial_spin_guesses = [self.initial_spin_guesses] if self.initial_spin_guesses is not None else []
        self.initial_spin_guesses = [v if isinstance(v, float) else float(v) for v in self.initial_spin_guesses]

        if self.spin_multiplicity is not None and not isinstance(self.spin_multiplicity, Bool):
            self.spin_multiplicity = Bool(self.spin_multiplicity)

        if self.absolute_magnetisation is not None and not isinstance(self.absolute_magnetisation, QuantityValue):
            self.absolute_magnetisation = QuantityValue(**as_dict(self.absolute_magnetisation))

        if self.total_magnetisation is not None and not isinstance(self.total_magnetisation, QuantityValue):
            self.total_magnetisation = QuantityValue(**as_dict(self.total_magnetisation))

        if self.total_spin_squared is not None and not isinstance(self.total_spin_squared, float):
            self.total_spin_squared = float(self.total_spin_squared)

        if self.net_charge is not None and not isinstance(self.net_charge, float):
            self.net_charge = float(self.net_charge)

        if self.smearing is not None and not isinstance(self.smearing, QuantityValue):
            self.smearing = QuantityValue(**as_dict(self.smearing))

        if self.fermi_dirac_temperature is not None and not isinstance(self.fermi_dirac_temperature, QuantityValue):
            self.fermi_dirac_temperature = QuantityValue(**as_dict(self.fermi_dirac_temperature))

        if self.force_multiplicity is not None and not isinstance(self.force_multiplicity, Bool):
            self.force_multiplicity = Bool(self.force_multiplicity)

        if self.output_mulliken_population_analysis is not None and not isinstance(self.output_mulliken_population_analysis, OutputPopulationAnalysis):
            self.output_mulliken_population_analysis = OutputPopulationAnalysis(**as_dict(self.output_mulliken_population_analysis))

        if self.output_hirshfeld_population_analysis is not None and not isinstance(self.output_hirshfeld_population_analysis, OutputPopulationAnalysis):
            self.output_hirshfeld_population_analysis = OutputPopulationAnalysis(**as_dict(self.output_hirshfeld_population_analysis))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        super().__post_init__(**kwargs)


@dataclass
class Dosing(YAMLRoot):
    """
    Dosing a molecule
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = CHMO["CHMO_0001383"]
    class_class_curie: ClassVar[str] = "chmo:CHMO_0001383"
    class_name: ClassVar[str] = "Dosing"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Dosing

    perm_id: Optional[str] = None
    name: Optional[str] = None
    sum_formula: Optional[str] = None
    duration: Optional[Union[dict, "QuantityValue"]] = None
    pressure: Optional[Union[dict, "QuantityValue"]] = None
    temperature: Optional[Union[dict, "QuantityValue"]] = None
    comments: Optional[str] = None
    instrument: Optional[Union[dict, "Instrument"]] = None
    sample: Optional[Union[dict, "Sample"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.sum_formula is not None and not isinstance(self.sum_formula, str):
            self.sum_formula = str(self.sum_formula)

        if self.duration is not None and not isinstance(self.duration, QuantityValue):
            self.duration = QuantityValue(**as_dict(self.duration))

        if self.pressure is not None and not isinstance(self.pressure, QuantityValue):
            self.pressure = QuantityValue(**as_dict(self.pressure))

        if self.temperature is not None and not isinstance(self.temperature, QuantityValue):
            self.temperature = QuantityValue(**as_dict(self.temperature))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        if self.sample is not None and not isinstance(self.sample, Sample):
            self.sample = Sample(**as_dict(self.sample))

        super().__post_init__(**kwargs)


@dataclass
class Draft(YAMLRoot):
    """
    Publication draft/Manuscript
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C85255"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C85255"
    class_name: ClassVar[str] = "Draft"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Draft

    perm_id: Optional[str] = None
    name: Optional[str] = None
    draft_type: Optional[Union[str, "DraftTypeEnum"]] = None
    comments: Optional[str] = None
    results: Optional[Union[Union[dict, "Result"], List[Union[dict, "Result"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.draft_type is not None and not isinstance(self.draft_type, DraftTypeEnum):
            self.draft_type = DraftTypeEnum(self.draft_type)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.results, list):
            self.results = [self.results] if self.results is not None else []
        self.results = [v if isinstance(v, Result) else Result(**as_dict(v)) for v in self.results]

        super().__post_init__(**kwargs)


@dataclass
class EmpiricalModelling(YAMLRoot):
    """
    Empirical Modelling
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["EmpiricalModelling"]
    class_class_curie: ClassVar[str] = "materialsML:EmpiricalModelling"
    class_name: ClassVar[str] = "EmpiricalModelling"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.EmpiricalModelling

    perm_id: Optional[str] = None
    name: Optional[str] = None
    wfms_uuid: Optional[str] = None
    wfms_url: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.wfms_url is not None and not isinstance(self.wfms_url, str):
            self.wfms_url = str(self.wfms_url)

        super().__post_init__(**kwargs)


@dataclass
class Errors(YAMLRoot):
    """
    Errors and problems
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C43369"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C43369"
    class_name: ClassVar[str] = "Errors"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Errors

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    instrument: Optional[Union[dict, "Instrument"]] = None
    uhv_components: Optional[Union[Union[dict, "UHVComponent"], List[Union[dict, "UHVComponent"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        super().__post_init__(**kwargs)


@dataclass
class FillCryostat(YAMLRoot):
    """
    Fill the cryostat
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["FillCryostat"]
    class_class_curie: ClassVar[str] = "materialsML:FillCryostat"
    class_name: ClassVar[str] = "FillCryostat"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.FillCryostat

    perm_id: Optional[str] = None
    name: Optional[str] = None
    weight_before: Optional[Union[dict, "QuantityValue"]] = None
    weight_after: Optional[Union[dict, "QuantityValue"]] = None
    substance: Optional[str] = None
    dewar: Optional[Union[str, "DewarEnum"]] = None
    comments: Optional[str] = None
    uhv_components: Optional[Union[Union[dict, "UHVComponent"], List[Union[dict, "UHVComponent"]]]] = empty_list()
    instrument: Optional[Union[dict, "Instrument"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.weight_before is not None and not isinstance(self.weight_before, QuantityValue):
            self.weight_before = QuantityValue(**as_dict(self.weight_before))

        if self.weight_after is not None and not isinstance(self.weight_after, QuantityValue):
            self.weight_after = QuantityValue(**as_dict(self.weight_after))

        if self.substance is not None and not isinstance(self.substance, str):
            self.substance = str(self.substance)

        if self.dewar is not None and not isinstance(self.dewar, DewarEnum):
            self.dewar = DewarEnum(self.dewar)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        super().__post_init__(**kwargs)


@dataclass
class GeometryOptimisation(YAMLRoot):
    """
    Geometry optimisation
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["GeometryOptimisation"]
    class_class_curie: ClassVar[str] = "materialsML:GeometryOptimisation"
    class_name: ClassVar[str] = "GeometryOptimisation"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.GeometryOptimisation

    perm_id: Optional[str] = None
    name: Optional[str] = None
    wfms_uuid: Optional[str] = None
    wfms_url: Optional[str] = None
    force_convergence_threshold: Optional[Union[dict, "ForceConvergenceThreshold"]] = None
    geometry_constraints: Optional[Union[float, List[float]]] = empty_list()
    comments: Optional[str] = None
    atomistic_models: Optional[Union[Union[dict, AtomisticModel], List[Union[dict, AtomisticModel]]]] = empty_list()
    dft: Optional[Union[dict, DFT]] = None
    empirical_modelling: Optional[Union[dict, EmpiricalModelling]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.wfms_url is not None and not isinstance(self.wfms_url, str):
            self.wfms_url = str(self.wfms_url)

        if self.force_convergence_threshold is not None and not isinstance(self.force_convergence_threshold, ForceConvergenceThreshold):
            self.force_convergence_threshold = ForceConvergenceThreshold(**as_dict(self.force_convergence_threshold))

        if not isinstance(self.geometry_constraints, list):
            self.geometry_constraints = [self.geometry_constraints] if self.geometry_constraints is not None else []
        self.geometry_constraints = [v if isinstance(v, float) else float(v) for v in self.geometry_constraints]

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.atomistic_models, list):
            self.atomistic_models = [self.atomistic_models] if self.atomistic_models is not None else []
        self.atomistic_models = [v if isinstance(v, AtomisticModel) else AtomisticModel(**as_dict(v)) for v in self.atomistic_models]

        if self.dft is not None and not isinstance(self.dft, DFT):
            self.dft = DFT(**as_dict(self.dft))

        if self.empirical_modelling is not None and not isinstance(self.empirical_modelling, EmpiricalModelling):
            self.empirical_modelling = EmpiricalModelling(**as_dict(self.empirical_modelling))

        super().__post_init__(**kwargs)


@dataclass
class Grant(YAMLRoot):
    """
    Grant
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["Grant"]
    class_class_curie: ClassVar[str] = "schema:Grant"
    class_name: ClassVar[str] = "Grant"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Grant

    perm_id: Optional[str] = None
    name: Optional[str] = None
    funder: Optional[str] = None
    start_date: Optional[Union[str, XSDDate]] = None
    end_date: Optional[Union[str, XSDDate]] = None
    budget: Optional[Union[dict, "QuantityValue"]] = None
    project_id: Optional[str] = None
    acknowledgement_sentence: Optional[str] = None
    comments: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.funder is not None and not isinstance(self.funder, str):
            self.funder = str(self.funder)

        if self.start_date is not None and not isinstance(self.start_date, XSDDate):
            self.start_date = XSDDate(self.start_date)

        if self.end_date is not None and not isinstance(self.end_date, XSDDate):
            self.end_date = XSDDate(self.end_date)

        if self.budget is not None and not isinstance(self.budget, QuantityValue):
            self.budget = QuantityValue(**as_dict(self.budget))

        if self.project_id is not None and not isinstance(self.project_id, str):
            self.project_id = str(self.project_id)

        if self.acknowledgement_sentence is not None and not isinstance(self.acknowledgement_sentence, str):
            self.acknowledgement_sentence = str(self.acknowledgement_sentence)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        super().__post_init__(**kwargs)


@dataclass
class HydrogenCracker(YAMLRoot):
    """
    Splitting gas molecules into elements
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["HydrogenCracker"]
    class_class_curie: ClassVar[str] = "materialsML:HydrogenCracker"
    class_name: ClassVar[str] = "HydrogenCracker"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.HydrogenCracker

    perm_id: Optional[str] = None
    name: Optional[str] = None
    sum_formula: Optional[str] = None
    comments: Optional[str] = None
    uhv_components: Optional[Union[Union[dict, "UHVComponent"], List[Union[dict, "UHVComponent"]]]] = empty_list()
    instrument: Optional[Union[dict, "Instrument"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.sum_formula is not None and not isinstance(self.sum_formula, str):
            self.sum_formula = str(self.sum_formula)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        super().__post_init__(**kwargs)


@dataclass
class Institution(YAMLRoot):
    """
    Institution
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["Organization"]
    class_class_curie: ClassVar[str] = "schema:Organization"
    class_name: ClassVar[str] = "Institution"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Institution

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    address: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.address is not None and not isinstance(self.address, str):
            self.address = str(self.address)

        super().__post_init__(**kwargs)


@dataclass
class Instrument(YAMLRoot):
    """
    Instrument
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C16742"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C16742"
    class_name: ClassVar[str] = "Instrument"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Instrument

    perm_id: Optional[str] = None
    name: Optional[str] = None
    serial_number: Optional[str] = None
    empa_id: Optional[str] = None
    documentation_website: Optional[str] = None
    receive_date: Optional[Union[str, XSDDate]] = None
    comments: Optional[str] = None
    uhv_components: Optional[Union[Union[dict, "UHVComponent"], List[Union[dict, "UHVComponent"]]]] = empty_list()
    room: Optional[Union[dict, "Room"]] = None
    manufacturer: Optional[Union[dict, "Manufacturer"]] = None
    tip_sensor: Optional[Union[dict, "TipSensor"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.serial_number is not None and not isinstance(self.serial_number, str):
            self.serial_number = str(self.serial_number)

        if self.empa_id is not None and not isinstance(self.empa_id, str):
            self.empa_id = str(self.empa_id)

        if self.documentation_website is not None and not isinstance(self.documentation_website, str):
            self.documentation_website = str(self.documentation_website)

        if self.receive_date is not None and not isinstance(self.receive_date, XSDDate):
            self.receive_date = XSDDate(self.receive_date)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        if self.room is not None and not isinstance(self.room, Room):
            self.room = Room(**as_dict(self.room))

        if self.manufacturer is not None and not isinstance(self.manufacturer, Manufacturer):
            self.manufacturer = Manufacturer(**as_dict(self.manufacturer))

        if self.tip_sensor is not None and not isinstance(self.tip_sensor, TipSensor):
            self.tip_sensor = TipSensor(**as_dict(self.tip_sensor))

        super().__post_init__(**kwargs)


@dataclass
class Layered2DMaterial(YAMLRoot):
    """
    Layered 2D material
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["Layered2DMaterial"]
    class_class_curie: ClassVar[str] = "materialsML:Layered2DMaterial"
    class_name: ClassVar[str] = "Layered2DMaterial"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Layered2DMaterial

    perm_id: Optional[str] = None
    name: Optional[str] = None
    empa_number: Optional[int] = None
    batch: Optional[str] = None
    sample_plate: Optional[str] = None
    layers_2d: Optional[Union[Union[dict, "Layers2DDetails"], List[Union[dict, "Layers2DDetails"]]]] = empty_list()
    substrates: Optional[Union[Union[dict, "SubstratesDetails"], List[Union[dict, "SubstratesDetails"]]]] = empty_list()
    width: Optional[Union[dict, "QuantityValue"]] = None
    length: Optional[Union[dict, "QuantityValue"]] = None
    thickness: Optional[Union[dict, "QuantityValue"]] = None
    hazardous: Optional[Union[bool, Bool]] = None
    hazardous_specification: Optional[str] = None
    fridge: Optional[Union[bool, Bool]] = None
    no_light: Optional[Union[bool, Bool]] = None
    dry: Optional[Union[bool, Bool]] = None
    no_oxygen: Optional[Union[bool, Bool]] = None
    other_storage_condition: Optional[Union[bool, Bool]] = None
    other_storage_condition_specification: Optional[str] = None
    receive_date: Optional[Union[str, XSDDate]] = None
    supplier_own_name: Optional[str] = None
    comments: Optional[str] = None
    storage: Optional[Union[dict, "Storage"]] = None
    chemist: Optional[Union[dict, Chemist]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.empa_number is not None and not isinstance(self.empa_number, int):
            self.empa_number = int(self.empa_number)

        if self.batch is not None and not isinstance(self.batch, str):
            self.batch = str(self.batch)

        if self.sample_plate is not None and not isinstance(self.sample_plate, str):
            self.sample_plate = str(self.sample_plate)

        if not isinstance(self.layers_2d, list):
            self.layers_2d = [self.layers_2d] if self.layers_2d is not None else []
        self.layers_2d = [v if isinstance(v, Layers2DDetails) else Layers2DDetails(**as_dict(v)) for v in self.layers_2d]

        if not isinstance(self.substrates, list):
            self.substrates = [self.substrates] if self.substrates is not None else []
        self.substrates = [v if isinstance(v, SubstratesDetails) else SubstratesDetails(**as_dict(v)) for v in self.substrates]

        if self.width is not None and not isinstance(self.width, QuantityValue):
            self.width = QuantityValue(**as_dict(self.width))

        if self.length is not None and not isinstance(self.length, QuantityValue):
            self.length = QuantityValue(**as_dict(self.length))

        if self.thickness is not None and not isinstance(self.thickness, QuantityValue):
            self.thickness = QuantityValue(**as_dict(self.thickness))

        if self.hazardous is not None and not isinstance(self.hazardous, Bool):
            self.hazardous = Bool(self.hazardous)

        if self.hazardous_specification is not None and not isinstance(self.hazardous_specification, str):
            self.hazardous_specification = str(self.hazardous_specification)

        if self.fridge is not None and not isinstance(self.fridge, Bool):
            self.fridge = Bool(self.fridge)

        if self.no_light is not None and not isinstance(self.no_light, Bool):
            self.no_light = Bool(self.no_light)

        if self.dry is not None and not isinstance(self.dry, Bool):
            self.dry = Bool(self.dry)

        if self.no_oxygen is not None and not isinstance(self.no_oxygen, Bool):
            self.no_oxygen = Bool(self.no_oxygen)

        if self.other_storage_condition is not None and not isinstance(self.other_storage_condition, Bool):
            self.other_storage_condition = Bool(self.other_storage_condition)

        if self.other_storage_condition_specification is not None and not isinstance(self.other_storage_condition_specification, str):
            self.other_storage_condition_specification = str(self.other_storage_condition_specification)

        if self.receive_date is not None and not isinstance(self.receive_date, XSDDate):
            self.receive_date = XSDDate(self.receive_date)

        if self.supplier_own_name is not None and not isinstance(self.supplier_own_name, str):
            self.supplier_own_name = str(self.supplier_own_name)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.storage is not None and not isinstance(self.storage, Storage):
            self.storage = Storage(**as_dict(self.storage))

        if self.chemist is not None and not isinstance(self.chemist, Chemist):
            self.chemist = Chemist(**as_dict(self.chemist))

        super().__post_init__(**kwargs)


@dataclass
class Maintenance(YAMLRoot):
    """
    Maintenance performed in the instruments
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C53297"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C53297"
    class_name: ClassVar[str] = "Maintenance"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Maintenance

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    comments: Optional[str] = None
    uhv_components: Optional[Union[Union[dict, "UHVComponent"], List[Union[dict, "UHVComponent"]]]] = empty_list()
    instrument: Optional[Union[dict, Instrument]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        super().__post_init__(**kwargs)


@dataclass
class Manufacturer(YAMLRoot):
    """
    Manufacturer
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C25392"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C25392"
    class_name: ClassVar[str] = "Manufacturer"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Manufacturer

    perm_id: Optional[str] = None
    name: Optional[str] = None
    email: Optional[str] = None
    work_phone: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.email is not None and not isinstance(self.email, str):
            self.email = str(self.email)

        if self.work_phone is not None and not isinstance(self.work_phone, str):
            self.work_phone = str(self.work_phone)

        super().__post_init__(**kwargs)


@dataclass
class MinimumEnergyPotential(YAMLRoot):
    """
    Minimum Energy Potential
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["MinimumEnergyPotential"]
    class_class_curie: ClassVar[str] = "materialsML:MinimumEnergyPotential"
    class_name: ClassVar[str] = "MinimumEnergyPotential"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.MinimumEnergyPotential

    perm_id: Optional[str] = None
    name: Optional[str] = None
    wfms_uuid: Optional[str] = None
    wfms_url: Optional[str] = None
    method_type: Optional[Union[str, "MethodTypeEnum"]] = None
    number_geometries: Optional[int] = None
    energy_barrier: Optional[Union[dict, "QuantityValue"]] = None
    energies: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    geometry_constraints: Optional[Union[float, List[float]]] = empty_list()
    geometry_constraints_increments: Optional[Union[float, List[float]]] = empty_list()
    comments: Optional[str] = None
    atomistic_models: Optional[Union[Union[dict, AtomisticModel], List[Union[dict, AtomisticModel]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.wfms_url is not None and not isinstance(self.wfms_url, str):
            self.wfms_url = str(self.wfms_url)

        if self.method_type is not None and not isinstance(self.method_type, MethodTypeEnum):
            self.method_type = MethodTypeEnum(self.method_type)

        if self.number_geometries is not None and not isinstance(self.number_geometries, int):
            self.number_geometries = int(self.number_geometries)

        if self.energy_barrier is not None and not isinstance(self.energy_barrier, QuantityValue):
            self.energy_barrier = QuantityValue(**as_dict(self.energy_barrier))

        if not isinstance(self.energies, list):
            self.energies = [self.energies] if self.energies is not None else []
        self.energies = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.energies]

        if not isinstance(self.geometry_constraints, list):
            self.geometry_constraints = [self.geometry_constraints] if self.geometry_constraints is not None else []
        self.geometry_constraints = [v if isinstance(v, float) else float(v) for v in self.geometry_constraints]

        if not isinstance(self.geometry_constraints_increments, list):
            self.geometry_constraints_increments = [self.geometry_constraints_increments] if self.geometry_constraints_increments is not None else []
        self.geometry_constraints_increments = [v if isinstance(v, float) else float(v) for v in self.geometry_constraints_increments]

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.atomistic_models, list):
            self.atomistic_models = [self.atomistic_models] if self.atomistic_models is not None else []
        self.atomistic_models = [v if isinstance(v, AtomisticModel) else AtomisticModel(**as_dict(v)) for v in self.atomistic_models]

        super().__post_init__(**kwargs)


@dataclass
class Molecule(YAMLRoot):
    """
    Molecule
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["MolecularEntity"]
    class_class_curie: ClassVar[str] = "schema:MolecularEntity"
    class_name: ClassVar[str] = "Molecule"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Molecule

    name: str = None
    perm_id: Optional[str] = None
    iupac_name: Optional[str] = None
    sum_formula: Optional[str] = None
    smiles: Optional[str] = None
    cas_number: Optional[str] = None
    empa_number: Optional[int] = None
    batch: Optional[str] = None
    vial: Optional[str] = None
    hazardous: Optional[Union[bool, Bool]] = None
    hazardous_specification: Optional[str] = None
    evaporation_temperatures: Optional[Union[Union[dict, "EvaporationTemperature"], List[Union[dict, "EvaporationTemperature"]]]] = empty_list()
    fridge: Optional[Union[bool, Bool]] = None
    no_light: Optional[Union[bool, Bool]] = None
    dry: Optional[Union[bool, Bool]] = None
    no_oxygen: Optional[Union[bool, Bool]] = None
    other_storage_condition: Optional[Union[bool, Bool]] = None
    other_storage_condition_specification: Optional[str] = None
    chemist_molecule_name: Optional[str] = None
    amount: Optional[Union[dict, "QuantityValue"]] = None
    receive_date: Optional[Union[str, XSDDate]] = None
    comments: Optional[str] = None
    chemist: Optional[Union[dict, Chemist]] = None
    storage: Optional[Union[dict, "Storage"]] = None
    precursor_molecules: Optional[Union[Union[dict, "Molecule"], List[Union[dict, "Molecule"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.name):
            self.MissingRequiredField("name")
        if not isinstance(self.name, str):
            self.name = str(self.name)

        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.iupac_name is not None and not isinstance(self.iupac_name, str):
            self.iupac_name = str(self.iupac_name)

        if self.sum_formula is not None and not isinstance(self.sum_formula, str):
            self.sum_formula = str(self.sum_formula)

        if self.smiles is not None and not isinstance(self.smiles, str):
            self.smiles = str(self.smiles)

        if self.cas_number is not None and not isinstance(self.cas_number, str):
            self.cas_number = str(self.cas_number)

        if self.empa_number is not None and not isinstance(self.empa_number, int):
            self.empa_number = int(self.empa_number)

        if self.batch is not None and not isinstance(self.batch, str):
            self.batch = str(self.batch)

        if self.vial is not None and not isinstance(self.vial, str):
            self.vial = str(self.vial)

        if self.hazardous is not None and not isinstance(self.hazardous, Bool):
            self.hazardous = Bool(self.hazardous)

        if self.hazardous_specification is not None and not isinstance(self.hazardous_specification, str):
            self.hazardous_specification = str(self.hazardous_specification)

        if not isinstance(self.evaporation_temperatures, list):
            self.evaporation_temperatures = [self.evaporation_temperatures] if self.evaporation_temperatures is not None else []
        self.evaporation_temperatures = [v if isinstance(v, EvaporationTemperature) else EvaporationTemperature(**as_dict(v)) for v in self.evaporation_temperatures]

        if self.fridge is not None and not isinstance(self.fridge, Bool):
            self.fridge = Bool(self.fridge)

        if self.no_light is not None and not isinstance(self.no_light, Bool):
            self.no_light = Bool(self.no_light)

        if self.dry is not None and not isinstance(self.dry, Bool):
            self.dry = Bool(self.dry)

        if self.no_oxygen is not None and not isinstance(self.no_oxygen, Bool):
            self.no_oxygen = Bool(self.no_oxygen)

        if self.other_storage_condition is not None and not isinstance(self.other_storage_condition, Bool):
            self.other_storage_condition = Bool(self.other_storage_condition)

        if self.other_storage_condition_specification is not None and not isinstance(self.other_storage_condition_specification, str):
            self.other_storage_condition_specification = str(self.other_storage_condition_specification)

        if self.chemist_molecule_name is not None and not isinstance(self.chemist_molecule_name, str):
            self.chemist_molecule_name = str(self.chemist_molecule_name)

        if self.amount is not None and not isinstance(self.amount, QuantityValue):
            self.amount = QuantityValue(**as_dict(self.amount))

        if self.receive_date is not None and not isinstance(self.receive_date, XSDDate):
            self.receive_date = XSDDate(self.receive_date)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.chemist is not None and not isinstance(self.chemist, Chemist):
            self.chemist = Chemist(**as_dict(self.chemist))

        if self.storage is not None and not isinstance(self.storage, Storage):
            self.storage = Storage(**as_dict(self.storage))

        if not isinstance(self.precursor_molecules, list):
            self.precursor_molecules = [self.precursor_molecules] if self.precursor_molecules is not None else []
        self.precursor_molecules = [v if isinstance(v, Molecule) else Molecule(**as_dict(v)) for v in self.precursor_molecules]

        super().__post_init__(**kwargs)


@dataclass
class Notes(YAMLRoot):
    """
    General Notes
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C42619"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C42619"
    class_name: ClassVar[str] = "Notes"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Notes

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        super().__post_init__(**kwargs)


@dataclass
class PDOS(YAMLRoot):
    """
    Projected Density of States
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["PDOS"]
    class_class_curie: ClassVar[str] = "materialsML:PDOS"
    class_name: ClassVar[str] = "PDOS"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.PDOS

    perm_id: Optional[str] = None
    name: Optional[str] = None
    wfms_uuid: Optional[str] = None
    wfms_url: Optional[str] = None
    atomic_selections: Optional[Union[float, List[float]]] = empty_list()
    energies: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    orbitals_selections: Optional[Union[float, List[float]]] = empty_list()
    amplitude: Optional[Union[dict, "QuantityValue"]] = None
    plot_contributions: Optional[Union[str, "PlotContributionsEnum"]] = None
    plot_group_by: Optional[Union[str, "PlotGroupByEnum"]] = None
    degauss_energy: Optional[Union[dict, "QuantityValue"]] = None
    comments: Optional[str] = None
    atomistic_models: Optional[Union[Union[dict, AtomisticModel], List[Union[dict, AtomisticModel]]]] = empty_list()
    dft: Optional[Union[dict, DFT]] = None
    tb: Optional[Union[dict, "TightBinding"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.wfms_url is not None and not isinstance(self.wfms_url, str):
            self.wfms_url = str(self.wfms_url)

        if not isinstance(self.atomic_selections, list):
            self.atomic_selections = [self.atomic_selections] if self.atomic_selections is not None else []
        self.atomic_selections = [v if isinstance(v, float) else float(v) for v in self.atomic_selections]

        if not isinstance(self.energies, list):
            self.energies = [self.energies] if self.energies is not None else []
        self.energies = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.energies]

        if not isinstance(self.orbitals_selections, list):
            self.orbitals_selections = [self.orbitals_selections] if self.orbitals_selections is not None else []
        self.orbitals_selections = [v if isinstance(v, float) else float(v) for v in self.orbitals_selections]

        if self.amplitude is not None and not isinstance(self.amplitude, QuantityValue):
            self.amplitude = QuantityValue(**as_dict(self.amplitude))

        if self.plot_contributions is not None and not isinstance(self.plot_contributions, PlotContributionsEnum):
            self.plot_contributions = PlotContributionsEnum(self.plot_contributions)

        if self.plot_group_by is not None and not isinstance(self.plot_group_by, PlotGroupByEnum):
            self.plot_group_by = PlotGroupByEnum(self.plot_group_by)

        if self.degauss_energy is not None and not isinstance(self.degauss_energy, QuantityValue):
            self.degauss_energy = QuantityValue(**as_dict(self.degauss_energy))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.atomistic_models, list):
            self.atomistic_models = [self.atomistic_models] if self.atomistic_models is not None else []
        self.atomistic_models = [v if isinstance(v, AtomisticModel) else AtomisticModel(**as_dict(v)) for v in self.atomistic_models]

        if self.dft is not None and not isinstance(self.dft, DFT):
            self.dft = DFT(**as_dict(self.dft))

        if self.tb is not None and not isinstance(self.tb, TightBinding):
            self.tb = TightBinding(**as_dict(self.tb))

        super().__post_init__(**kwargs)


@dataclass
class PotentialEnergyCalculation(YAMLRoot):
    """
    Potential Energy Calculation
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["PotentialEnergyCalculation"]
    class_class_curie: ClassVar[str] = "materialsML:PotentialEnergyCalculation"
    class_name: ClassVar[str] = "PotentialEnergyCalculation"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.PotentialEnergyCalculation

    perm_id: Optional[str] = None
    name: Optional[str] = None
    wfms_uuid: Optional[str] = None
    wfms_url: Optional[str] = None
    energy: Optional[Union[dict, "QuantityValue"]] = None
    dft: Optional[Union[dict, DFT]] = None
    tb: Optional[Union[dict, "TightBinding"]] = None
    empirical_modelling: Optional[Union[dict, EmpiricalModelling]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.wfms_url is not None and not isinstance(self.wfms_url, str):
            self.wfms_url = str(self.wfms_url)

        if self.energy is not None and not isinstance(self.energy, QuantityValue):
            self.energy = QuantityValue(**as_dict(self.energy))

        if self.dft is not None and not isinstance(self.dft, DFT):
            self.dft = DFT(**as_dict(self.dft))

        if self.tb is not None and not isinstance(self.tb, TightBinding):
            self.tb = TightBinding(**as_dict(self.tb))

        if self.empirical_modelling is not None and not isinstance(self.empirical_modelling, EmpiricalModelling):
            self.empirical_modelling = EmpiricalModelling(**as_dict(self.empirical_modelling))

        super().__post_init__(**kwargs)


@dataclass
class Person(YAMLRoot):
    """
    Person
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["Person"]
    class_class_curie: ClassVar[str] = "schema:Person"
    class_name: ClassVar[str] = "Person"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Person

    perm_id: Optional[str] = None
    name: Optional[str] = None
    short_name: Optional[str] = None
    work_status: Optional[Union[bool, Bool]] = None
    email: Optional[str] = None
    work_phone: Optional[str] = None
    mobile_phone: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.short_name is not None and not isinstance(self.short_name, str):
            self.short_name = str(self.short_name)

        if self.work_status is not None and not isinstance(self.work_status, Bool):
            self.work_status = Bool(self.work_status)

        if self.email is not None and not isinstance(self.email, str):
            self.email = str(self.email)

        if self.work_phone is not None and not isinstance(self.work_phone, str):
            self.work_phone = str(self.work_phone)

        if self.mobile_phone is not None and not isinstance(self.mobile_phone, str):
            self.mobile_phone = str(self.mobile_phone)

        super().__post_init__(**kwargs)


@dataclass
class Protocol(YAMLRoot):
    """
    Protocol
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = OBI["OBI_0000272"]
    class_class_curie: ClassVar[str] = "obi:OBI_0000272"
    class_name: ClassVar[str] = "Protocol"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Protocol

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        super().__post_init__(**kwargs)


@dataclass
class Publication(YAMLRoot):
    """
    Publication
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["PublicationIssue"]
    class_class_curie: ClassVar[str] = "schema:PublicationIssue"
    class_name: ClassVar[str] = "Publication"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Publication

    perm_id: Optional[str] = None
    name: Optional[str] = None
    abstract: Optional[str] = None
    doi: Optional[str] = None
    year: Optional[int] = None
    url: Optional[str] = None
    dataset_url: Optional[str] = None
    comments: Optional[str] = None
    grants: Optional[Union[Union[dict, Grant], List[Union[dict, Grant]]]] = empty_list()
    drafts: Optional[Union[Union[dict, Draft], List[Union[dict, Draft]]]] = empty_list()
    authors: Optional[Union[Union[dict, Author], List[Union[dict, Author]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.abstract is not None and not isinstance(self.abstract, str):
            self.abstract = str(self.abstract)

        if self.doi is not None and not isinstance(self.doi, str):
            self.doi = str(self.doi)

        if self.year is not None and not isinstance(self.year, int):
            self.year = int(self.year)

        if self.url is not None and not isinstance(self.url, str):
            self.url = str(self.url)

        if self.dataset_url is not None and not isinstance(self.dataset_url, str):
            self.dataset_url = str(self.dataset_url)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.grants, list):
            self.grants = [self.grants] if self.grants is not None else []
        self.grants = [v if isinstance(v, Grant) else Grant(**as_dict(v)) for v in self.grants]

        if not isinstance(self.drafts, list):
            self.drafts = [self.drafts] if self.drafts is not None else []
        self.drafts = [v if isinstance(v, Draft) else Draft(**as_dict(v)) for v in self.drafts]

        if not isinstance(self.authors, list):
            self.authors = [self.authors] if self.authors is not None else []
        self.authors = [v if isinstance(v, Author) else Author(**as_dict(v)) for v in self.authors]

        super().__post_init__(**kwargs)


@dataclass
class Result(YAMLRoot):
    """
    Results
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = EMMO["EMMO_0f6f0120_c079_4d95_bb11_4ddee05e530e"]
    class_class_curie: ClassVar[str] = "emmo:EMMO_0f6f0120_c079_4d95_bb11_4ddee05e530e"
    class_name: ClassVar[str] = "Result"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Result

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    comments: Optional[str] = None
    stms: Optional[Union[Union[dict, "STM"], List[Union[dict, "STM"]]]] = empty_list()
    afms: Optional[Union[Union[dict, AFM], List[Union[dict, AFM]]]] = empty_list()
    stss: Optional[Union[Union[dict, "STS"], List[Union[dict, "STS"]]]] = empty_list()
    geometry_optimisations: Optional[Union[Union[dict, GeometryOptimisation], List[Union[dict, GeometryOptimisation]]]] = empty_list()
    softwares: Optional[Union[Union[dict, "Software"], List[Union[dict, "Software"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.stms, list):
            self.stms = [self.stms] if self.stms is not None else []
        self.stms = [v if isinstance(v, STM) else STM(**as_dict(v)) for v in self.stms]

        if not isinstance(self.afms, list):
            self.afms = [self.afms] if self.afms is not None else []
        self.afms = [v if isinstance(v, AFM) else AFM(**as_dict(v)) for v in self.afms]

        if not isinstance(self.stss, list):
            self.stss = [self.stss] if self.stss is not None else []
        self.stss = [v if isinstance(v, STS) else STS(**as_dict(v)) for v in self.stss]

        if not isinstance(self.geometry_optimisations, list):
            self.geometry_optimisations = [self.geometry_optimisations] if self.geometry_optimisations is not None else []
        self.geometry_optimisations = [v if isinstance(v, GeometryOptimisation) else GeometryOptimisation(**as_dict(v)) for v in self.geometry_optimisations]

        if not isinstance(self.softwares, list):
            self.softwares = [self.softwares] if self.softwares is not None else []
        self.softwares = [v if isinstance(v, Software) else Software(**as_dict(v)) for v in self.softwares]

        super().__post_init__(**kwargs)


@dataclass
class Room(YAMLRoot):
    """
    Room
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["Room"]
    class_class_curie: ClassVar[str] = "schema:Room"
    class_name: ClassVar[str] = "Room"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Room

    perm_id: Optional[str] = None
    name: Optional[str] = None
    comments: Optional[str] = None
    institution: Optional[Union[dict, Institution]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.institution is not None and not isinstance(self.institution, Institution):
            self.institution = Institution(**as_dict(self.institution))

        super().__post_init__(**kwargs)


@dataclass
class Sample(YAMLRoot):
    """
    Sample being measured
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["Sample"]
    class_class_curie: ClassVar[str] = "materialsML:Sample"
    class_name: ClassVar[str] = "Sample"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Sample

    perm_id: Optional[str] = None
    name: Optional[str] = None
    crystal: Optional[Union[dict, Crystal]] = None
    layered_2d_material: Optional[Union[dict, Layered2DMaterial]] = None
    wafer_substrate: Optional[Union[dict, "WaferSubstrate"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.crystal is not None and not isinstance(self.crystal, Crystal):
            self.crystal = Crystal(**as_dict(self.crystal))

        if self.layered_2d_material is not None and not isinstance(self.layered_2d_material, Layered2DMaterial):
            self.layered_2d_material = Layered2DMaterial(**as_dict(self.layered_2d_material))

        if self.wafer_substrate is not None and not isinstance(self.wafer_substrate, WaferSubstrate):
            self.wafer_substrate = WaferSubstrate(**as_dict(self.wafer_substrate))

        super().__post_init__(**kwargs)


@dataclass
class Software(YAMLRoot):
    """
    Used software script
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["SoftwareSourceCode"]
    class_class_curie: ClassVar[str] = "schema:SoftwareSourceCode"
    class_name: ClassVar[str] = "Software"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Software

    perm_id: Optional[str] = None
    name: Optional[str] = None
    repository_url: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.repository_url is not None and not isinstance(self.repository_url, str):
            self.repository_url = str(self.repository_url)

        super().__post_init__(**kwargs)


@dataclass
class Sputtering(YAMLRoot):
    """
    Sputtering
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = CHMO["CHMO_0001569"]
    class_class_curie: ClassVar[str] = "chmo:CHMO_0001569"
    class_name: ClassVar[str] = "Sputtering"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Sputtering

    perm_id: Optional[str] = None
    name: Optional[str] = None
    duration: Optional[Union[dict, "QuantityValue"]] = None
    pressure: Optional[Union[dict, "QuantityValue"]] = None
    discharge_voltage: Optional[Union[dict, "QuantityValue"]] = None
    voltage: Optional[Union[dict, "QuantityValue"]] = None
    temperature: Optional[Union[dict, "QuantityValue"]] = None
    angle: Optional[Union[dict, "QuantityValue"]] = None
    current: Optional[Union[dict, "QuantityValue"]] = None
    comments: Optional[str] = None
    uhv_components: Optional[Union[Union[dict, "UHVComponent"], List[Union[dict, "UHVComponent"]]]] = empty_list()
    annealing: Optional[Union[dict, Annealing]] = None
    deposition: Optional[Union[dict, Deposition]] = None
    instrument: Optional[Union[dict, Instrument]] = None
    sample: Optional[Union[dict, Sample]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.duration is not None and not isinstance(self.duration, QuantityValue):
            self.duration = QuantityValue(**as_dict(self.duration))

        if self.pressure is not None and not isinstance(self.pressure, QuantityValue):
            self.pressure = QuantityValue(**as_dict(self.pressure))

        if self.discharge_voltage is not None and not isinstance(self.discharge_voltage, QuantityValue):
            self.discharge_voltage = QuantityValue(**as_dict(self.discharge_voltage))

        if self.voltage is not None and not isinstance(self.voltage, QuantityValue):
            self.voltage = QuantityValue(**as_dict(self.voltage))

        if self.temperature is not None and not isinstance(self.temperature, QuantityValue):
            self.temperature = QuantityValue(**as_dict(self.temperature))

        if self.angle is not None and not isinstance(self.angle, QuantityValue):
            self.angle = QuantityValue(**as_dict(self.angle))

        if self.current is not None and not isinstance(self.current, QuantityValue):
            self.current = QuantityValue(**as_dict(self.current))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        if self.annealing is not None and not isinstance(self.annealing, Annealing):
            self.annealing = Annealing(**as_dict(self.annealing))

        if self.deposition is not None and not isinstance(self.deposition, Deposition):
            self.deposition = Deposition(**as_dict(self.deposition))

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        if self.sample is not None and not isinstance(self.sample, Sample):
            self.sample = Sample(**as_dict(self.sample))

        super().__post_init__(**kwargs)


@dataclass
class SRD(YAMLRoot):
    """
    Sublimation Rate Determination
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["SRD"]
    class_class_curie: ClassVar[str] = "materialsML:SRD"
    class_name: ClassVar[str] = "SRD"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.SRD

    perm_id: Optional[str] = None
    name: Optional[str] = None
    comments: Optional[str] = None
    molecule: Optional[Union[dict, Molecule]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.molecule is not None and not isinstance(self.molecule, Molecule):
            self.molecule = Molecule(**as_dict(self.molecule))

        super().__post_init__(**kwargs)


@dataclass
class STM(YAMLRoot):
    """
    Scanning Tunneling Microscopy
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = CHMO["CHMO_0000132"]
    class_class_curie: ClassVar[str] = "chmo:CHMO_0000132"
    class_name: ClassVar[str] = "STM"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.STM

    perm_id: Optional[str] = None
    name: Optional[str] = None
    start_time: Optional[str] = None
    duration: Optional[Union[dict, "QuantityValue"]] = None
    bias_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    bias_calibration_factor: Optional[Union[dict, "QuantityValue"]] = None
    bias_calibration_offset: Optional[Union[dict, "QuantityValue"]] = None
    current_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    current_calibration_factor: Optional[Union[dict, "QuantityValue"]] = None
    current_calibration_offset: Optional[Union[dict, "QuantityValue"]] = None
    current_gain: Optional[Union[dict, "QuantityValue"]] = None
    z_position: Optional[Union[dict, "QuantityValue"]] = None
    feedback_active: Optional[Union[bool, Bool]] = None
    feedback_type: Optional[Union[str, "FeedbackTypeEnum"]] = None
    z_controller_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_p_gain: Optional[float] = None
    z_controller_i_gain: Optional[float] = None
    z_controller_time_constant: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_tip_lift: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_switch_off_delay: Optional[Union[dict, "QuantityValue"]] = None
    piezo_configuration_settings: Optional[Union[dict, "PiezoConfigurationSettings"]] = None
    scan_settings: Optional[Union[dict, "ScanSettings"]] = None
    dwell_time: Optional[Union[dict, "QuantityValue"]] = None
    sample_temperature: Optional[Union[dict, "QuantityValue"]] = None
    recording_temperature: Optional[Union[dict, "QuantityValue"]] = None
    comments: Optional[str] = None
    e_min: Optional[Union[dict, "QuantityValue"]] = None
    e_max: Optional[Union[dict, "QuantityValue"]] = None
    de: Optional[Union[dict, "QuantityValue"]] = None
    fwhm: Optional[Union[dict, "QuantityValue"]] = None
    extrap_plane: Optional[Union[dict, "QuantityValue"]] = None
    constant_height: Optional[Union[dict, "QuantityValue"]] = None
    constant_current: Optional[Union[dict, "QuantityValue"]] = None
    wfms_uuid: Optional[str] = None
    annealing: Optional[Union[dict, Annealing]] = None
    deposition: Optional[Union[dict, Deposition]] = None
    dft: Optional[Union[dict, DFT]] = None
    tb: Optional[Union[dict, "TightBinding"]] = None
    instrument: Optional[Union[dict, Instrument]] = None
    sample: Optional[Union[dict, Sample]] = None
    tip_sensor: Optional[Union[dict, "TipSensor"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.start_time is not None and not isinstance(self.start_time, str):
            self.start_time = str(self.start_time)

        if self.duration is not None and not isinstance(self.duration, QuantityValue):
            self.duration = QuantityValue(**as_dict(self.duration))

        if self.bias_setpoint is not None and not isinstance(self.bias_setpoint, QuantityValue):
            self.bias_setpoint = QuantityValue(**as_dict(self.bias_setpoint))

        if self.bias_calibration_factor is not None and not isinstance(self.bias_calibration_factor, QuantityValue):
            self.bias_calibration_factor = QuantityValue(**as_dict(self.bias_calibration_factor))

        if self.bias_calibration_offset is not None and not isinstance(self.bias_calibration_offset, QuantityValue):
            self.bias_calibration_offset = QuantityValue(**as_dict(self.bias_calibration_offset))

        if self.current_setpoint is not None and not isinstance(self.current_setpoint, QuantityValue):
            self.current_setpoint = QuantityValue(**as_dict(self.current_setpoint))

        if self.current_calibration_factor is not None and not isinstance(self.current_calibration_factor, QuantityValue):
            self.current_calibration_factor = QuantityValue(**as_dict(self.current_calibration_factor))

        if self.current_calibration_offset is not None and not isinstance(self.current_calibration_offset, QuantityValue):
            self.current_calibration_offset = QuantityValue(**as_dict(self.current_calibration_offset))

        if self.current_gain is not None and not isinstance(self.current_gain, QuantityValue):
            self.current_gain = QuantityValue(**as_dict(self.current_gain))

        if self.z_position is not None and not isinstance(self.z_position, QuantityValue):
            self.z_position = QuantityValue(**as_dict(self.z_position))

        if self.feedback_active is not None and not isinstance(self.feedback_active, Bool):
            self.feedback_active = Bool(self.feedback_active)

        if self.feedback_type is not None and not isinstance(self.feedback_type, FeedbackTypeEnum):
            self.feedback_type = FeedbackTypeEnum(self.feedback_type)

        if self.z_controller_setpoint is not None and not isinstance(self.z_controller_setpoint, QuantityValue):
            self.z_controller_setpoint = QuantityValue(**as_dict(self.z_controller_setpoint))

        if self.z_controller_p_gain is not None and not isinstance(self.z_controller_p_gain, float):
            self.z_controller_p_gain = float(self.z_controller_p_gain)

        if self.z_controller_i_gain is not None and not isinstance(self.z_controller_i_gain, float):
            self.z_controller_i_gain = float(self.z_controller_i_gain)

        if self.z_controller_time_constant is not None and not isinstance(self.z_controller_time_constant, QuantityValue):
            self.z_controller_time_constant = QuantityValue(**as_dict(self.z_controller_time_constant))

        if self.z_controller_tip_lift is not None and not isinstance(self.z_controller_tip_lift, QuantityValue):
            self.z_controller_tip_lift = QuantityValue(**as_dict(self.z_controller_tip_lift))

        if self.z_controller_switch_off_delay is not None and not isinstance(self.z_controller_switch_off_delay, QuantityValue):
            self.z_controller_switch_off_delay = QuantityValue(**as_dict(self.z_controller_switch_off_delay))

        if self.piezo_configuration_settings is not None and not isinstance(self.piezo_configuration_settings, PiezoConfigurationSettings):
            self.piezo_configuration_settings = PiezoConfigurationSettings(**as_dict(self.piezo_configuration_settings))

        if self.scan_settings is not None and not isinstance(self.scan_settings, ScanSettings):
            self.scan_settings = ScanSettings(**as_dict(self.scan_settings))

        if self.dwell_time is not None and not isinstance(self.dwell_time, QuantityValue):
            self.dwell_time = QuantityValue(**as_dict(self.dwell_time))

        if self.sample_temperature is not None and not isinstance(self.sample_temperature, QuantityValue):
            self.sample_temperature = QuantityValue(**as_dict(self.sample_temperature))

        if self.recording_temperature is not None and not isinstance(self.recording_temperature, QuantityValue):
            self.recording_temperature = QuantityValue(**as_dict(self.recording_temperature))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.e_min is not None and not isinstance(self.e_min, QuantityValue):
            self.e_min = QuantityValue(**as_dict(self.e_min))

        if self.e_max is not None and not isinstance(self.e_max, QuantityValue):
            self.e_max = QuantityValue(**as_dict(self.e_max))

        if self.de is not None and not isinstance(self.de, QuantityValue):
            self.de = QuantityValue(**as_dict(self.de))

        if self.fwhm is not None and not isinstance(self.fwhm, QuantityValue):
            self.fwhm = QuantityValue(**as_dict(self.fwhm))

        if self.extrap_plane is not None and not isinstance(self.extrap_plane, QuantityValue):
            self.extrap_plane = QuantityValue(**as_dict(self.extrap_plane))

        if self.constant_height is not None and not isinstance(self.constant_height, QuantityValue):
            self.constant_height = QuantityValue(**as_dict(self.constant_height))

        if self.constant_current is not None and not isinstance(self.constant_current, QuantityValue):
            self.constant_current = QuantityValue(**as_dict(self.constant_current))

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.annealing is not None and not isinstance(self.annealing, Annealing):
            self.annealing = Annealing(**as_dict(self.annealing))

        if self.deposition is not None and not isinstance(self.deposition, Deposition):
            self.deposition = Deposition(**as_dict(self.deposition))

        if self.dft is not None and not isinstance(self.dft, DFT):
            self.dft = DFT(**as_dict(self.dft))

        if self.tb is not None and not isinstance(self.tb, TightBinding):
            self.tb = TightBinding(**as_dict(self.tb))

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        if self.sample is not None and not isinstance(self.sample, Sample):
            self.sample = Sample(**as_dict(self.sample))

        if self.tip_sensor is not None and not isinstance(self.tip_sensor, TipSensor):
            self.tip_sensor = TipSensor(**as_dict(self.tip_sensor))

        super().__post_init__(**kwargs)


@dataclass
class Storage(YAMLRoot):
    """
    Storage
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C16143"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C16143"
    class_name: ClassVar[str] = "Storage"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Storage

    perm_id: Optional[str] = None
    name: Optional[str] = None
    comments: Optional[str] = None
    room: Optional[Union[dict, Room]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.room is not None and not isinstance(self.room, Room):
            self.room = Room(**as_dict(self.room))

        super().__post_init__(**kwargs)


@dataclass
class STS(YAMLRoot):
    """
    Scanning Tunneling Spectroscopy
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = FAIRMAT_STS["nxsts"]
    class_class_curie: ClassVar[str] = "fairmat_sts:nxsts"
    class_name: ClassVar[str] = "STS"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.STS

    perm_id: Optional[str] = None
    name: Optional[str] = None
    start_time: Optional[str] = None
    duration: Optional[Union[dict, "QuantityValue"]] = None
    experiment_type: Optional[Union[str, "ExperimentTypeEnum"]] = None
    acquisition_coordinate_x: Optional[Union[dict, "QuantityValue"]] = None
    acquisition_coordinate_y: Optional[Union[dict, "QuantityValue"]] = None
    acquisition_coordinate_z: Optional[Union[dict, "QuantityValue"]] = None
    final_z: Optional[Union[dict, "QuantityValue"]] = None
    filter_type: Optional[Union[str, "FilterTypeEnum"]] = None
    filter_order: Optional[int] = None
    filter_cutoff: Optional[Union[dict, "QuantityValue"]] = None
    bias_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    bias_calibration_factor: Optional[Union[dict, "QuantityValue"]] = None
    bias_calibration_offset: Optional[Union[dict, "QuantityValue"]] = None
    num_pixel: Optional[int] = None
    z_avg_time: Optional[Union[dict, "QuantityValue"]] = None
    first_settling_time: Optional[Union[dict, "QuantityValue"]] = None
    settling_time: Optional[Union[dict, "QuantityValue"]] = None
    integration_time: Optional[Union[dict, "QuantityValue"]] = None
    end_settling_time: Optional[Union[dict, "QuantityValue"]] = None
    z_control_time: Optional[Union[dict, "QuantityValue"]] = None
    max_slew_rate: Optional[Union[dict, "QuantityValue"]] = None
    backward_sweep: Optional[Union[bool, Bool]] = None
    num_sweeps: Optional[int] = None
    channel_names: Optional[Union[str, List[str]]] = empty_list()
    record_final_z: Optional[Union[bool, Bool]] = None
    bias_spectroscopy_reset_bias: Optional[Union[bool, Bool]] = None
    bias_spectroscopy_sweep_start: Optional[Union[dict, "QuantityValue"]] = None
    bias_spectroscopy_sweep_end: Optional[Union[dict, "QuantityValue"]] = None
    bias_spectroscopy_z_offset: Optional[Union[dict, "QuantityValue"]] = None
    bias_spectroscopy_z_controller_hold: Optional[Union[bool, Bool]] = None
    z_spectroscopy_reset_z: Optional[Union[bool, Bool]] = None
    z_spectroscopy_initial_z_offset: Optional[Union[dict, "QuantityValue"]] = None
    z_spectroscopy_sweep_distance: Optional[Union[dict, "QuantityValue"]] = None
    z_spectroscopy_time_between_forward_backward: Optional[Union[dict, "QuantityValue"]] = None
    current_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    current_calibration_factor: Optional[Union[dict, "QuantityValue"]] = None
    current_calibration_offset: Optional[Union[dict, "QuantityValue"]] = None
    current_gain: Optional[Union[dict, "QuantityValue"]] = None
    z_position: Optional[Union[dict, "QuantityValue"]] = None
    feedback_active: Optional[Union[bool, Bool]] = None
    feedback_type: Optional[Union[str, "FeedbackTypeEnum"]] = None
    z_controller_setpoint: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_p_gain: Optional[float] = None
    z_controller_i_gain: Optional[float] = None
    z_controller_time_constant: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_tip_lift: Optional[Union[dict, "QuantityValue"]] = None
    z_controller_switch_off_delay: Optional[Union[dict, "QuantityValue"]] = None
    lock_in_status: Optional[Union[bool, Bool]] = None
    lock_in_modulated_signal: Optional[str] = None
    lock_in_frequency: Optional[Union[dict, "QuantityValue"]] = None
    lock_in_amplitude: Optional[Union[dict, "QuantityValue"]] = None
    lock_in_hp_filter_cutoffs: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    lock_in_hp_filter_orders: Optional[Union[int, List[int]]] = empty_list()
    lock_in_demodulated_signal: Optional[Union[str, "LockInDemodulatedSignalEnum"]] = None
    lock_in_harmonics: Optional[Union[int, List[int]]] = empty_list()
    lock_in_reference_phases: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    lock_in_lp_filter_cutoffs: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    lock_in_lp_filter_orders: Optional[Union[int, List[int]]] = empty_list()
    lock_in_sync_filter: Optional[Union[bool, Bool]] = None
    sample_temperature: Optional[Union[dict, "QuantityValue"]] = None
    wfms_uuid: Optional[str] = None
    comments: Optional[str] = None
    annealing: Optional[Union[dict, Annealing]] = None
    deposition: Optional[Union[dict, Deposition]] = None
    instrument: Optional[Union[dict, Instrument]] = None
    sample: Optional[Union[dict, Sample]] = None
    tip_sensor: Optional[Union[dict, "TipSensor"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.start_time is not None and not isinstance(self.start_time, str):
            self.start_time = str(self.start_time)

        if self.duration is not None and not isinstance(self.duration, QuantityValue):
            self.duration = QuantityValue(**as_dict(self.duration))

        if self.experiment_type is not None and not isinstance(self.experiment_type, ExperimentTypeEnum):
            self.experiment_type = ExperimentTypeEnum(self.experiment_type)

        if self.acquisition_coordinate_x is not None and not isinstance(self.acquisition_coordinate_x, QuantityValue):
            self.acquisition_coordinate_x = QuantityValue(**as_dict(self.acquisition_coordinate_x))

        if self.acquisition_coordinate_y is not None and not isinstance(self.acquisition_coordinate_y, QuantityValue):
            self.acquisition_coordinate_y = QuantityValue(**as_dict(self.acquisition_coordinate_y))

        if self.acquisition_coordinate_z is not None and not isinstance(self.acquisition_coordinate_z, QuantityValue):
            self.acquisition_coordinate_z = QuantityValue(**as_dict(self.acquisition_coordinate_z))

        if self.final_z is not None and not isinstance(self.final_z, QuantityValue):
            self.final_z = QuantityValue(**as_dict(self.final_z))

        if self.filter_type is not None and not isinstance(self.filter_type, FilterTypeEnum):
            self.filter_type = FilterTypeEnum(self.filter_type)

        if self.filter_order is not None and not isinstance(self.filter_order, int):
            self.filter_order = int(self.filter_order)

        if self.filter_cutoff is not None and not isinstance(self.filter_cutoff, QuantityValue):
            self.filter_cutoff = QuantityValue(**as_dict(self.filter_cutoff))

        if self.bias_setpoint is not None and not isinstance(self.bias_setpoint, QuantityValue):
            self.bias_setpoint = QuantityValue(**as_dict(self.bias_setpoint))

        if self.bias_calibration_factor is not None and not isinstance(self.bias_calibration_factor, QuantityValue):
            self.bias_calibration_factor = QuantityValue(**as_dict(self.bias_calibration_factor))

        if self.bias_calibration_offset is not None and not isinstance(self.bias_calibration_offset, QuantityValue):
            self.bias_calibration_offset = QuantityValue(**as_dict(self.bias_calibration_offset))

        if self.num_pixel is not None and not isinstance(self.num_pixel, int):
            self.num_pixel = int(self.num_pixel)

        if self.z_avg_time is not None and not isinstance(self.z_avg_time, QuantityValue):
            self.z_avg_time = QuantityValue(**as_dict(self.z_avg_time))

        if self.first_settling_time is not None and not isinstance(self.first_settling_time, QuantityValue):
            self.first_settling_time = QuantityValue(**as_dict(self.first_settling_time))

        if self.settling_time is not None and not isinstance(self.settling_time, QuantityValue):
            self.settling_time = QuantityValue(**as_dict(self.settling_time))

        if self.integration_time is not None and not isinstance(self.integration_time, QuantityValue):
            self.integration_time = QuantityValue(**as_dict(self.integration_time))

        if self.end_settling_time is not None and not isinstance(self.end_settling_time, QuantityValue):
            self.end_settling_time = QuantityValue(**as_dict(self.end_settling_time))

        if self.z_control_time is not None and not isinstance(self.z_control_time, QuantityValue):
            self.z_control_time = QuantityValue(**as_dict(self.z_control_time))

        if self.max_slew_rate is not None and not isinstance(self.max_slew_rate, QuantityValue):
            self.max_slew_rate = QuantityValue(**as_dict(self.max_slew_rate))

        if self.backward_sweep is not None and not isinstance(self.backward_sweep, Bool):
            self.backward_sweep = Bool(self.backward_sweep)

        if self.num_sweeps is not None and not isinstance(self.num_sweeps, int):
            self.num_sweeps = int(self.num_sweeps)

        if not isinstance(self.channel_names, list):
            self.channel_names = [self.channel_names] if self.channel_names is not None else []
        self.channel_names = [v if isinstance(v, str) else str(v) for v in self.channel_names]

        if self.record_final_z is not None and not isinstance(self.record_final_z, Bool):
            self.record_final_z = Bool(self.record_final_z)

        if self.bias_spectroscopy_reset_bias is not None and not isinstance(self.bias_spectroscopy_reset_bias, Bool):
            self.bias_spectroscopy_reset_bias = Bool(self.bias_spectroscopy_reset_bias)

        if self.bias_spectroscopy_sweep_start is not None and not isinstance(self.bias_spectroscopy_sweep_start, QuantityValue):
            self.bias_spectroscopy_sweep_start = QuantityValue(**as_dict(self.bias_spectroscopy_sweep_start))

        if self.bias_spectroscopy_sweep_end is not None and not isinstance(self.bias_spectroscopy_sweep_end, QuantityValue):
            self.bias_spectroscopy_sweep_end = QuantityValue(**as_dict(self.bias_spectroscopy_sweep_end))

        if self.bias_spectroscopy_z_offset is not None and not isinstance(self.bias_spectroscopy_z_offset, QuantityValue):
            self.bias_spectroscopy_z_offset = QuantityValue(**as_dict(self.bias_spectroscopy_z_offset))

        if self.bias_spectroscopy_z_controller_hold is not None and not isinstance(self.bias_spectroscopy_z_controller_hold, Bool):
            self.bias_spectroscopy_z_controller_hold = Bool(self.bias_spectroscopy_z_controller_hold)

        if self.z_spectroscopy_reset_z is not None and not isinstance(self.z_spectroscopy_reset_z, Bool):
            self.z_spectroscopy_reset_z = Bool(self.z_spectroscopy_reset_z)

        if self.z_spectroscopy_initial_z_offset is not None and not isinstance(self.z_spectroscopy_initial_z_offset, QuantityValue):
            self.z_spectroscopy_initial_z_offset = QuantityValue(**as_dict(self.z_spectroscopy_initial_z_offset))

        if self.z_spectroscopy_sweep_distance is not None and not isinstance(self.z_spectroscopy_sweep_distance, QuantityValue):
            self.z_spectroscopy_sweep_distance = QuantityValue(**as_dict(self.z_spectroscopy_sweep_distance))

        if self.z_spectroscopy_time_between_forward_backward is not None and not isinstance(self.z_spectroscopy_time_between_forward_backward, QuantityValue):
            self.z_spectroscopy_time_between_forward_backward = QuantityValue(**as_dict(self.z_spectroscopy_time_between_forward_backward))

        if self.current_setpoint is not None and not isinstance(self.current_setpoint, QuantityValue):
            self.current_setpoint = QuantityValue(**as_dict(self.current_setpoint))

        if self.current_calibration_factor is not None and not isinstance(self.current_calibration_factor, QuantityValue):
            self.current_calibration_factor = QuantityValue(**as_dict(self.current_calibration_factor))

        if self.current_calibration_offset is not None and not isinstance(self.current_calibration_offset, QuantityValue):
            self.current_calibration_offset = QuantityValue(**as_dict(self.current_calibration_offset))

        if self.current_gain is not None and not isinstance(self.current_gain, QuantityValue):
            self.current_gain = QuantityValue(**as_dict(self.current_gain))

        if self.z_position is not None and not isinstance(self.z_position, QuantityValue):
            self.z_position = QuantityValue(**as_dict(self.z_position))

        if self.feedback_active is not None and not isinstance(self.feedback_active, Bool):
            self.feedback_active = Bool(self.feedback_active)

        if self.feedback_type is not None and not isinstance(self.feedback_type, FeedbackTypeEnum):
            self.feedback_type = FeedbackTypeEnum(self.feedback_type)

        if self.z_controller_setpoint is not None and not isinstance(self.z_controller_setpoint, QuantityValue):
            self.z_controller_setpoint = QuantityValue(**as_dict(self.z_controller_setpoint))

        if self.z_controller_p_gain is not None and not isinstance(self.z_controller_p_gain, float):
            self.z_controller_p_gain = float(self.z_controller_p_gain)

        if self.z_controller_i_gain is not None and not isinstance(self.z_controller_i_gain, float):
            self.z_controller_i_gain = float(self.z_controller_i_gain)

        if self.z_controller_time_constant is not None and not isinstance(self.z_controller_time_constant, QuantityValue):
            self.z_controller_time_constant = QuantityValue(**as_dict(self.z_controller_time_constant))

        if self.z_controller_tip_lift is not None and not isinstance(self.z_controller_tip_lift, QuantityValue):
            self.z_controller_tip_lift = QuantityValue(**as_dict(self.z_controller_tip_lift))

        if self.z_controller_switch_off_delay is not None and not isinstance(self.z_controller_switch_off_delay, QuantityValue):
            self.z_controller_switch_off_delay = QuantityValue(**as_dict(self.z_controller_switch_off_delay))

        if self.lock_in_status is not None and not isinstance(self.lock_in_status, Bool):
            self.lock_in_status = Bool(self.lock_in_status)

        if self.lock_in_modulated_signal is not None and not isinstance(self.lock_in_modulated_signal, str):
            self.lock_in_modulated_signal = str(self.lock_in_modulated_signal)

        if self.lock_in_frequency is not None and not isinstance(self.lock_in_frequency, QuantityValue):
            self.lock_in_frequency = QuantityValue(**as_dict(self.lock_in_frequency))

        if self.lock_in_amplitude is not None and not isinstance(self.lock_in_amplitude, QuantityValue):
            self.lock_in_amplitude = QuantityValue(**as_dict(self.lock_in_amplitude))

        if not isinstance(self.lock_in_hp_filter_cutoffs, list):
            self.lock_in_hp_filter_cutoffs = [self.lock_in_hp_filter_cutoffs] if self.lock_in_hp_filter_cutoffs is not None else []
        self.lock_in_hp_filter_cutoffs = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.lock_in_hp_filter_cutoffs]

        if not isinstance(self.lock_in_hp_filter_orders, list):
            self.lock_in_hp_filter_orders = [self.lock_in_hp_filter_orders] if self.lock_in_hp_filter_orders is not None else []
        self.lock_in_hp_filter_orders = [v if isinstance(v, int) else int(v) for v in self.lock_in_hp_filter_orders]

        if self.lock_in_demodulated_signal is not None and not isinstance(self.lock_in_demodulated_signal, LockInDemodulatedSignalEnum):
            self.lock_in_demodulated_signal = LockInDemodulatedSignalEnum(self.lock_in_demodulated_signal)

        if not isinstance(self.lock_in_harmonics, list):
            self.lock_in_harmonics = [self.lock_in_harmonics] if self.lock_in_harmonics is not None else []
        self.lock_in_harmonics = [v if isinstance(v, int) else int(v) for v in self.lock_in_harmonics]

        if not isinstance(self.lock_in_reference_phases, list):
            self.lock_in_reference_phases = [self.lock_in_reference_phases] if self.lock_in_reference_phases is not None else []
        self.lock_in_reference_phases = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.lock_in_reference_phases]

        if not isinstance(self.lock_in_lp_filter_cutoffs, list):
            self.lock_in_lp_filter_cutoffs = [self.lock_in_lp_filter_cutoffs] if self.lock_in_lp_filter_cutoffs is not None else []
        self.lock_in_lp_filter_cutoffs = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.lock_in_lp_filter_cutoffs]

        if not isinstance(self.lock_in_lp_filter_orders, list):
            self.lock_in_lp_filter_orders = [self.lock_in_lp_filter_orders] if self.lock_in_lp_filter_orders is not None else []
        self.lock_in_lp_filter_orders = [v if isinstance(v, int) else int(v) for v in self.lock_in_lp_filter_orders]

        if self.lock_in_sync_filter is not None and not isinstance(self.lock_in_sync_filter, Bool):
            self.lock_in_sync_filter = Bool(self.lock_in_sync_filter)

        if self.sample_temperature is not None and not isinstance(self.sample_temperature, QuantityValue):
            self.sample_temperature = QuantityValue(**as_dict(self.sample_temperature))

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.annealing is not None and not isinstance(self.annealing, Annealing):
            self.annealing = Annealing(**as_dict(self.annealing))

        if self.deposition is not None and not isinstance(self.deposition, Deposition):
            self.deposition = Deposition(**as_dict(self.deposition))

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        if self.sample is not None and not isinstance(self.sample, Sample):
            self.sample = Sample(**as_dict(self.sample))

        if self.tip_sensor is not None and not isinstance(self.tip_sensor, TipSensor):
            self.tip_sensor = TipSensor(**as_dict(self.tip_sensor))

        super().__post_init__(**kwargs)


@dataclass
class Supplier(YAMLRoot):
    """
    Supplier
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C43530"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C43530"
    class_name: ClassVar[str] = "Supplier"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Supplier

    perm_id: Optional[str] = None
    name: Optional[str] = None
    email: Optional[str] = None
    work_phone: Optional[str] = None
    address: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.email is not None and not isinstance(self.email, str):
            self.email = str(self.email)

        if self.work_phone is not None and not isinstance(self.work_phone, str):
            self.work_phone = str(self.work_phone)

        if self.address is not None and not isinstance(self.address, str):
            self.address = str(self.address)

        super().__post_init__(**kwargs)


@dataclass
class SystemStatus(YAMLRoot):
    """
    System Status
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["SystemStatus"]
    class_class_curie: ClassVar[str] = "materialsML:SystemStatus"
    class_name: ClassVar[str] = "SystemStatus"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.SystemStatus

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    instrument: Optional[Union[dict, Instrument]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.instrument is not None and not isinstance(self.instrument, Instrument):
            self.instrument = Instrument(**as_dict(self.instrument))

        super().__post_init__(**kwargs)


@dataclass
class TightBinding(YAMLRoot):
    """
    Tight Binding
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["TightBinding"]
    class_class_curie: ClassVar[str] = "materialsML:TightBinding"
    class_name: ClassVar[str] = "TightBinding"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.TightBinding

    perm_id: Optional[str] = None
    name: Optional[str] = None
    wfms_url: Optional[str] = None
    wfms_uuid: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.wfms_url is not None and not isinstance(self.wfms_url, str):
            self.wfms_url = str(self.wfms_url)

        if self.wfms_uuid is not None and not isinstance(self.wfms_uuid, str):
            self.wfms_uuid = str(self.wfms_uuid)

        super().__post_init__(**kwargs)


@dataclass
class TipPreparation(YAMLRoot):
    """
    Tip Preparation
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["TipPreparation"]
    class_class_curie: ClassVar[str] = "materialsML:TipPreparation"
    class_name: ClassVar[str] = "TipPreparation"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.TipPreparation

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    tip_preparation_type: Optional[Union[str, "TipPreparationTypeEnum"]] = None
    tip_sensor: Optional[Union[dict, "TipSensor"]] = None
    authors: Optional[Union[Union[dict, Author], List[Union[dict, Author]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.tip_preparation_type is not None and not isinstance(self.tip_preparation_type, TipPreparationTypeEnum):
            self.tip_preparation_type = TipPreparationTypeEnum(self.tip_preparation_type)

        if self.tip_sensor is not None and not isinstance(self.tip_sensor, TipSensor):
            self.tip_sensor = TipSensor(**as_dict(self.tip_sensor))

        if not isinstance(self.authors, list):
            self.authors = [self.authors] if self.authors is not None else []
        self.authors = [v if isinstance(v, Author) else Author(**as_dict(v)) for v in self.authors]

        super().__post_init__(**kwargs)


@dataclass
class TipSensor(YAMLRoot):
    """
    Tip Sensor
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = EMMO_CHAMEO["Probe"]
    class_class_curie: ClassVar[str] = "emmo_chameo:Probe"
    class_name: ClassVar[str] = "TipSensor"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.TipSensor

    perm_id: Optional[str] = None
    name: Optional[str] = None
    tip_sensor_type: Optional[Union[str, "TipSensorTypeEnum"]] = None
    wire_material: Optional[str] = None
    holder_id: Optional[str] = None
    thickness: Optional[Union[dict, "QuantityValue"]] = None
    length: Optional[Union[dict, "QuantityValue"]] = None
    resonance_frequency: Optional[Union[dict, "QuantityValue"]] = None
    q_factor: Optional[float] = None
    reference_number: Optional[str] = None
    device_parameters: Optional[str] = None
    receive_date: Optional[Union[str, XSDDate]] = None
    manufacturer: Optional[Union[dict, Manufacturer]] = None
    comments: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.tip_sensor_type is not None and not isinstance(self.tip_sensor_type, TipSensorTypeEnum):
            self.tip_sensor_type = TipSensorTypeEnum(self.tip_sensor_type)

        if self.wire_material is not None and not isinstance(self.wire_material, str):
            self.wire_material = str(self.wire_material)

        if self.holder_id is not None and not isinstance(self.holder_id, str):
            self.holder_id = str(self.holder_id)

        if self.thickness is not None and not isinstance(self.thickness, QuantityValue):
            self.thickness = QuantityValue(**as_dict(self.thickness))

        if self.length is not None and not isinstance(self.length, QuantityValue):
            self.length = QuantityValue(**as_dict(self.length))

        if self.resonance_frequency is not None and not isinstance(self.resonance_frequency, QuantityValue):
            self.resonance_frequency = QuantityValue(**as_dict(self.resonance_frequency))

        if self.q_factor is not None and not isinstance(self.q_factor, float):
            self.q_factor = float(self.q_factor)

        if self.reference_number is not None and not isinstance(self.reference_number, str):
            self.reference_number = str(self.reference_number)

        if self.device_parameters is not None and not isinstance(self.device_parameters, str):
            self.device_parameters = str(self.device_parameters)

        if self.receive_date is not None and not isinstance(self.receive_date, XSDDate):
            self.receive_date = XSDDate(self.receive_date)

        if self.manufacturer is not None and not isinstance(self.manufacturer, Manufacturer):
            self.manufacturer = Manufacturer(**as_dict(self.manufacturer))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        super().__post_init__(**kwargs)


@dataclass
class Transfer(YAMLRoot):
    """
    Transference
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = NCIT["NCIT_C48167"]
    class_class_curie: ClassVar[str] = "ncit:NCIT_C48167"
    class_name: ClassVar[str] = "Transfer"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Transfer

    perm_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    location_before: Optional[Union[dict, Room]] = None
    location_after: Optional[Union[dict, Room]] = None
    sample: Optional[Union[dict, Sample]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.location_before is not None and not isinstance(self.location_before, Room):
            self.location_before = Room(**as_dict(self.location_before))

        if self.location_after is not None and not isinstance(self.location_after, Room):
            self.location_after = Room(**as_dict(self.location_after))

        if self.sample is not None and not isinstance(self.sample, Sample):
            self.sample = Sample(**as_dict(self.sample))

        super().__post_init__(**kwargs)


@dataclass
class UHVComponent(YAMLRoot):
    """
    Ultra High Vacuum Component
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["UHVComponent"]
    class_class_curie: ClassVar[str] = "materialsML:UHVComponent"
    class_name: ClassVar[str] = "UHVComponent"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.UHVComponent

    perm_id: Optional[str] = None
    name: Optional[str] = None
    uhv_component_type: Optional[Union[str, "UHVComponentTypeEnum"]] = None
    empa_id: Optional[str] = None
    model: Optional[str] = None
    serial_number: Optional[str] = None
    receive_date: Optional[Union[str, XSDDate]] = None
    comments: Optional[str] = None
    manufacturer: Optional[Union[dict, Manufacturer]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.uhv_component_type is not None and not isinstance(self.uhv_component_type, UHVComponentTypeEnum):
            self.uhv_component_type = UHVComponentTypeEnum(self.uhv_component_type)

        if self.empa_id is not None and not isinstance(self.empa_id, str):
            self.empa_id = str(self.empa_id)

        if self.model is not None and not isinstance(self.model, str):
            self.model = str(self.model)

        if self.serial_number is not None and not isinstance(self.serial_number, str):
            self.serial_number = str(self.serial_number)

        if self.receive_date is not None and not isinstance(self.receive_date, XSDDate):
            self.receive_date = XSDDate(self.receive_date)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.manufacturer is not None and not isinstance(self.manufacturer, Manufacturer):
            self.manufacturer = Manufacturer(**as_dict(self.manufacturer))

        super().__post_init__(**kwargs)


@dataclass
class WaferSubstrate(YAMLRoot):
    """
    Wafer Substrate
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["WaferSubstrate"]
    class_class_curie: ClassVar[str] = "materialsML:WaferSubstrate"
    class_name: ClassVar[str] = "WaferSubstrate"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.WaferSubstrate

    perm_id: Optional[str] = None
    name: Optional[str] = None
    batch: Optional[str] = None
    sample_plate: Optional[str] = None
    substrates: Optional[Union[Union[dict, "SubstratesDetails"], List[Union[dict, "SubstratesDetails"]]]] = empty_list()
    diameter: Optional[Union[dict, "QuantityValue"]] = None
    height: Optional[Union[dict, "QuantityValue"]] = None
    thickness: Optional[Union[dict, "QuantityValue"]] = None
    hazardous: Optional[Union[bool, Bool]] = None
    hazardous_specification: Optional[str] = None
    fridge: Optional[Union[bool, Bool]] = None
    no_light: Optional[Union[bool, Bool]] = None
    dry: Optional[Union[bool, Bool]] = None
    no_oxygen: Optional[Union[bool, Bool]] = None
    other_storage_condition: Optional[Union[bool, Bool]] = None
    other_storage_condition_specification: Optional[str] = None
    receive_date: Optional[Union[str, XSDDate]] = None
    comments: Optional[str] = None
    storage: Optional[Union[dict, Storage]] = None
    chemist: Optional[Union[dict, Chemist]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.perm_id is not None and not isinstance(self.perm_id, str):
            self.perm_id = str(self.perm_id)

        if self.name is not None and not isinstance(self.name, str):
            self.name = str(self.name)

        if self.batch is not None and not isinstance(self.batch, str):
            self.batch = str(self.batch)

        if self.sample_plate is not None and not isinstance(self.sample_plate, str):
            self.sample_plate = str(self.sample_plate)

        if not isinstance(self.substrates, list):
            self.substrates = [self.substrates] if self.substrates is not None else []
        self.substrates = [v if isinstance(v, SubstratesDetails) else SubstratesDetails(**as_dict(v)) for v in self.substrates]

        if self.diameter is not None and not isinstance(self.diameter, QuantityValue):
            self.diameter = QuantityValue(**as_dict(self.diameter))

        if self.height is not None and not isinstance(self.height, QuantityValue):
            self.height = QuantityValue(**as_dict(self.height))

        if self.thickness is not None and not isinstance(self.thickness, QuantityValue):
            self.thickness = QuantityValue(**as_dict(self.thickness))

        if self.hazardous is not None and not isinstance(self.hazardous, Bool):
            self.hazardous = Bool(self.hazardous)

        if self.hazardous_specification is not None and not isinstance(self.hazardous_specification, str):
            self.hazardous_specification = str(self.hazardous_specification)

        if self.fridge is not None and not isinstance(self.fridge, Bool):
            self.fridge = Bool(self.fridge)

        if self.no_light is not None and not isinstance(self.no_light, Bool):
            self.no_light = Bool(self.no_light)

        if self.dry is not None and not isinstance(self.dry, Bool):
            self.dry = Bool(self.dry)

        if self.no_oxygen is not None and not isinstance(self.no_oxygen, Bool):
            self.no_oxygen = Bool(self.no_oxygen)

        if self.other_storage_condition is not None and not isinstance(self.other_storage_condition, Bool):
            self.other_storage_condition = Bool(self.other_storage_condition)

        if self.other_storage_condition_specification is not None and not isinstance(self.other_storage_condition_specification, str):
            self.other_storage_condition_specification = str(self.other_storage_condition_specification)

        if self.receive_date is not None and not isinstance(self.receive_date, XSDDate):
            self.receive_date = XSDDate(self.receive_date)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.storage is not None and not isinstance(self.storage, Storage):
            self.storage = Storage(**as_dict(self.storage))

        if self.chemist is not None and not isinstance(self.chemist, Chemist):
            self.chemist = Chemist(**as_dict(self.chemist))

        super().__post_init__(**kwargs)


@dataclass
class PiezoConfigurationSettings(YAMLRoot):
    """
    Piezo configuration settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["PiezoConfigurationSettings"]
    class_class_curie: ClassVar[str] = "materialsML:PiezoConfigurationSettings"
    class_name: ClassVar[str] = "PiezoConfigurationSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.PiezoConfigurationSettings

    active_calib: Optional[Union[str, "ActiveCalibEnum"]] = None
    sensitivity_settings: Optional[Union[dict, "SensitivitySettings"]] = None
    hv_gain_settings: Optional[Union[dict, "HVGainSettings"]] = None
    scan_slope: Optional[Union[dict, "ScanSlope"]] = None
    curvature_radius: Optional[Union[dict, "CurvatureRadius"]] = None
    second_order_correction: Optional[Union[dict, "SecondOrderCorrection"]] = None
    drift_settings: Optional[Union[dict, "DriftSettings"]] = None
    drift_correction_status: Optional[Union[bool, Bool]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.active_calib is not None and not isinstance(self.active_calib, ActiveCalibEnum):
            self.active_calib = ActiveCalibEnum(self.active_calib)

        if self.sensitivity_settings is not None and not isinstance(self.sensitivity_settings, SensitivitySettings):
            self.sensitivity_settings = SensitivitySettings(**as_dict(self.sensitivity_settings))

        if self.hv_gain_settings is not None and not isinstance(self.hv_gain_settings, HVGainSettings):
            self.hv_gain_settings = HVGainSettings(**as_dict(self.hv_gain_settings))

        if self.scan_slope is not None and not isinstance(self.scan_slope, ScanSlope):
            self.scan_slope = ScanSlope(**as_dict(self.scan_slope))

        if self.curvature_radius is not None and not isinstance(self.curvature_radius, CurvatureRadius):
            self.curvature_radius = CurvatureRadius(**as_dict(self.curvature_radius))

        if self.second_order_correction is not None and not isinstance(self.second_order_correction, SecondOrderCorrection):
            self.second_order_correction = SecondOrderCorrection(**as_dict(self.second_order_correction))

        if self.drift_settings is not None and not isinstance(self.drift_settings, DriftSettings):
            self.drift_settings = DriftSettings(**as_dict(self.drift_settings))

        if self.drift_correction_status is not None and not isinstance(self.drift_correction_status, Bool):
            self.drift_correction_status = Bool(self.drift_correction_status)

        super().__post_init__(**kwargs)


@dataclass
class SensitivitySettings(YAMLRoot):
    """
    Sensitivity settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["SensitivitySettings"]
    class_class_curie: ClassVar[str] = "materialsML:SensitivitySettings"
    class_name: ClassVar[str] = "SensitivitySettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.SensitivitySettings

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    z_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        if self.z_coordinate is not None and not isinstance(self.z_coordinate, QuantityValue):
            self.z_coordinate = QuantityValue(**as_dict(self.z_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class HVGainSettings(YAMLRoot):
    """
    HV gain settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["HVGainSettings"]
    class_class_curie: ClassVar[str] = "materialsML:HVGainSettings"
    class_name: ClassVar[str] = "HVGainSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.HVGainSettings

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    z_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        if self.z_coordinate is not None and not isinstance(self.z_coordinate, QuantityValue):
            self.z_coordinate = QuantityValue(**as_dict(self.z_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class ScanSlope(YAMLRoot):
    """
    Scan slope axis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["ScanSlope"]
    class_class_curie: ClassVar[str] = "materialsML:ScanSlope"
    class_name: ClassVar[str] = "ScanSlope"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.ScanSlope

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class CurvatureRadius(YAMLRoot):
    """
    Curvature radius axis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["CurvatureRadius"]
    class_class_curie: ClassVar[str] = "materialsML:CurvatureRadius"
    class_name: ClassVar[str] = "CurvatureRadius"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.CurvatureRadius

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class SecondOrderCorrection(YAMLRoot):
    """
    Second order correction axis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["SecondOrderCorrection"]
    class_class_curie: ClassVar[str] = "materialsML:SecondOrderCorrection"
    class_name: ClassVar[str] = "SecondOrderCorrection"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.SecondOrderCorrection

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class DriftSettings(YAMLRoot):
    """
    Drift settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["DriftSettings"]
    class_class_curie: ClassVar[str] = "materialsML:DriftSettings"
    class_name: ClassVar[str] = "DriftSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.DriftSettings

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    z_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        if self.z_coordinate is not None and not isinstance(self.z_coordinate, QuantityValue):
            self.z_coordinate = QuantityValue(**as_dict(self.z_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class ScanSettings(YAMLRoot):
    """
    Scan settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["ScanSettings"]
    class_class_curie: ClassVar[str] = "materialsML:ScanSettings"
    class_name: ClassVar[str] = "ScanSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.ScanSettings

    field: Optional[Union[float, List[float]]] = empty_list()
    angle: Optional[Union[dict, "QuantityValue"]] = None
    scan_offset: Optional[Union[dict, "ScanOffset"]] = None
    scan_time: Optional[Union[dict, "ScanTime"]] = None
    scan_range: Optional[Union[dict, "ScanRange"]] = None
    scan_pixels: Optional[Union[dict, "ScanPixels"]] = None
    scan_speed: Optional[Union[dict, "ScanSpeed"]] = None
    direction_speed: Optional[Union[str, "DirectionSpeedEnum"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if not isinstance(self.field, list):
            self.field = [self.field] if self.field is not None else []
        self.field = [v if isinstance(v, float) else float(v) for v in self.field]

        if self.angle is not None and not isinstance(self.angle, QuantityValue):
            self.angle = QuantityValue(**as_dict(self.angle))

        if self.scan_offset is not None and not isinstance(self.scan_offset, ScanOffset):
            self.scan_offset = ScanOffset(**as_dict(self.scan_offset))

        if self.scan_time is not None and not isinstance(self.scan_time, ScanTime):
            self.scan_time = ScanTime(**as_dict(self.scan_time))

        if self.scan_range is not None and not isinstance(self.scan_range, ScanRange):
            self.scan_range = ScanRange(**as_dict(self.scan_range))

        if self.scan_pixels is not None and not isinstance(self.scan_pixels, ScanPixels):
            self.scan_pixels = ScanPixels(**as_dict(self.scan_pixels))

        if self.scan_speed is not None and not isinstance(self.scan_speed, ScanSpeed):
            self.scan_speed = ScanSpeed(**as_dict(self.scan_speed))

        if self.direction_speed is not None and not isinstance(self.direction_speed, DirectionSpeedEnum):
            self.direction_speed = DirectionSpeedEnum(self.direction_speed)

        super().__post_init__(**kwargs)


@dataclass
class ScanOffset(YAMLRoot):
    """
    Scan offset axis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["ScanOffset"]
    class_class_curie: ClassVar[str] = "materialsML:ScanOffset"
    class_name: ClassVar[str] = "ScanOffset"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.ScanOffset

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class ScanTime(YAMLRoot):
    """
    Scan time by axis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["ScanTime"]
    class_class_curie: ClassVar[str] = "materialsML:ScanTime"
    class_name: ClassVar[str] = "ScanTime"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.ScanTime

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class ScanRange(YAMLRoot):
    """
    Scan range by axis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["ScanRange"]
    class_class_curie: ClassVar[str] = "materialsML:ScanRange"
    class_name: ClassVar[str] = "ScanRange"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.ScanRange

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class ScanPixels(YAMLRoot):
    """
    Scan pixels by axis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["ScanPixels"]
    class_class_curie: ClassVar[str] = "materialsML:ScanPixels"
    class_name: ClassVar[str] = "ScanPixels"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.ScanPixels

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class ScanSpeed(YAMLRoot):
    """
    Speed speed by axis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["ScanSpeed"]
    class_class_curie: ClassVar[str] = "materialsML:ScanSpeed"
    class_name: ClassVar[str] = "ScanSpeed"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.ScanSpeed

    forward_speed: Optional[Union[dict, "QuantityValue"]] = None
    backward_speed: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.forward_speed is not None and not isinstance(self.forward_speed, QuantityValue):
            self.forward_speed = QuantityValue(**as_dict(self.forward_speed))

        if self.backward_speed is not None and not isinstance(self.backward_speed, QuantityValue):
            self.backward_speed = QuantityValue(**as_dict(self.backward_speed))

        super().__post_init__(**kwargs)


@dataclass
class OscillationControlSettings(YAMLRoot):
    """
    Oscillation control settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["OscillationControlSettings"]
    class_class_curie: ClassVar[str] = "materialsML:OscillationControlSettings"
    class_name: ClassVar[str] = "OscillationControlSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.OscillationControlSettings

    differential_input: Optional[Union[bool, Bool]] = None
    one_ten: Optional[Union[bool, Bool]] = None
    input_calibration: Optional[Union[dict, "QuantityValue"]] = None
    input_range: Optional[Union[dict, "QuantityValue"]] = None
    center_frequency: Optional[Union[dict, "QuantityValue"]] = None
    range: Optional[Union[dict, "QuantityValue"]] = None
    demodulation_settings: Optional[Union[dict, "DemodulationSettings"]] = None
    phase_settings: Optional[Union[dict, "PhaseSettings"]] = None
    frequency_shift: Optional[Union[dict, "QuantityValue"]] = None
    amplitude_settings: Optional[Union[dict, "AmplitudeSettings"]] = None
    excitation: Optional[Union[dict, "QuantityValue"]] = None
    output_settings: Optional[Union[dict, "OutputSettings"]] = None
    pll_setup_settings: Optional[Union[dict, "PLLSetupSettings"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.differential_input is not None and not isinstance(self.differential_input, Bool):
            self.differential_input = Bool(self.differential_input)

        if self.one_ten is not None and not isinstance(self.one_ten, Bool):
            self.one_ten = Bool(self.one_ten)

        if self.input_calibration is not None and not isinstance(self.input_calibration, QuantityValue):
            self.input_calibration = QuantityValue(**as_dict(self.input_calibration))

        if self.input_range is not None and not isinstance(self.input_range, QuantityValue):
            self.input_range = QuantityValue(**as_dict(self.input_range))

        if self.center_frequency is not None and not isinstance(self.center_frequency, QuantityValue):
            self.center_frequency = QuantityValue(**as_dict(self.center_frequency))

        if self.range is not None and not isinstance(self.range, QuantityValue):
            self.range = QuantityValue(**as_dict(self.range))

        if self.demodulation_settings is not None and not isinstance(self.demodulation_settings, DemodulationSettings):
            self.demodulation_settings = DemodulationSettings(**as_dict(self.demodulation_settings))

        if self.phase_settings is not None and not isinstance(self.phase_settings, PhaseSettings):
            self.phase_settings = PhaseSettings(**as_dict(self.phase_settings))

        if self.frequency_shift is not None and not isinstance(self.frequency_shift, QuantityValue):
            self.frequency_shift = QuantityValue(**as_dict(self.frequency_shift))

        if self.amplitude_settings is not None and not isinstance(self.amplitude_settings, AmplitudeSettings):
            self.amplitude_settings = AmplitudeSettings(**as_dict(self.amplitude_settings))

        if self.excitation is not None and not isinstance(self.excitation, QuantityValue):
            self.excitation = QuantityValue(**as_dict(self.excitation))

        if self.output_settings is not None and not isinstance(self.output_settings, OutputSettings):
            self.output_settings = OutputSettings(**as_dict(self.output_settings))

        if self.pll_setup_settings is not None and not isinstance(self.pll_setup_settings, PLLSetupSettings):
            self.pll_setup_settings = PLLSetupSettings(**as_dict(self.pll_setup_settings))

        super().__post_init__(**kwargs)


@dataclass
class DemodulationSettings(YAMLRoot):
    """
    Demodulation settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["DemodulationSettings"]
    class_class_curie: ClassVar[str] = "materialsML:DemodulationSettings"
    class_name: ClassVar[str] = "DemodulationSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.DemodulationSettings

    inputs: Optional[Union[float, List[float]]] = empty_list()
    frequencies: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    reference_phases: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    cutoff_frequencies: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    harmonics: Optional[Union[Union[dict, "QuantityValue"], List[Union[dict, "QuantityValue"]]]] = empty_list()
    filters_orders: Optional[Union[float, List[float]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if not isinstance(self.inputs, list):
            self.inputs = [self.inputs] if self.inputs is not None else []
        self.inputs = [v if isinstance(v, float) else float(v) for v in self.inputs]

        if not isinstance(self.frequencies, list):
            self.frequencies = [self.frequencies] if self.frequencies is not None else []
        self.frequencies = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.frequencies]

        if not isinstance(self.reference_phases, list):
            self.reference_phases = [self.reference_phases] if self.reference_phases is not None else []
        self.reference_phases = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.reference_phases]

        if not isinstance(self.cutoff_frequencies, list):
            self.cutoff_frequencies = [self.cutoff_frequencies] if self.cutoff_frequencies is not None else []
        self.cutoff_frequencies = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.cutoff_frequencies]

        if not isinstance(self.harmonics, list):
            self.harmonics = [self.harmonics] if self.harmonics is not None else []
        self.harmonics = [v if isinstance(v, QuantityValue) else QuantityValue(**as_dict(v)) for v in self.harmonics]

        if not isinstance(self.filters_orders, list):
            self.filters_orders = [self.filters_orders] if self.filters_orders is not None else []
        self.filters_orders = [v if isinstance(v, float) else float(v) for v in self.filters_orders]

        super().__post_init__(**kwargs)


@dataclass
class PhaseSettings(YAMLRoot):
    """
    Phase settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["PhaseSettings"]
    class_class_curie: ClassVar[str] = "materialsML:PhaseSettings"
    class_name: ClassVar[str] = "PhaseSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.PhaseSettings

    p_gain: Optional[Union[dict, "QuantityValue"]] = None
    i_gain: Optional[Union[dict, "QuantityValue"]] = None
    controller_on: Optional[Union[bool, Bool]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.p_gain is not None and not isinstance(self.p_gain, QuantityValue):
            self.p_gain = QuantityValue(**as_dict(self.p_gain))

        if self.i_gain is not None and not isinstance(self.i_gain, QuantityValue):
            self.i_gain = QuantityValue(**as_dict(self.i_gain))

        if self.controller_on is not None and not isinstance(self.controller_on, Bool):
            self.controller_on = Bool(self.controller_on)

        super().__post_init__(**kwargs)


@dataclass
class AmplitudeSettings(YAMLRoot):
    """
    Amplitude settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["AmplitudeSettings"]
    class_class_curie: ClassVar[str] = "materialsML:AmplitudeSettings"
    class_name: ClassVar[str] = "AmplitudeSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.AmplitudeSettings

    setpoint: Optional[Union[dict, "QuantityValue"]] = None
    p_gain: Optional[Union[dict, "QuantityValue"]] = None
    i_gain: Optional[Union[dict, "QuantityValue"]] = None
    controller_on: Optional[Union[bool, Bool]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.setpoint is not None and not isinstance(self.setpoint, QuantityValue):
            self.setpoint = QuantityValue(**as_dict(self.setpoint))

        if self.p_gain is not None and not isinstance(self.p_gain, QuantityValue):
            self.p_gain = QuantityValue(**as_dict(self.p_gain))

        if self.i_gain is not None and not isinstance(self.i_gain, QuantityValue):
            self.i_gain = QuantityValue(**as_dict(self.i_gain))

        if self.controller_on is not None and not isinstance(self.controller_on, Bool):
            self.controller_on = Bool(self.controller_on)

        super().__post_init__(**kwargs)


@dataclass
class OutputSettings(YAMLRoot):
    """
    Output settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["OutputSettings"]
    class_class_curie: ClassVar[str] = "materialsML:OutputSettings"
    class_name: ClassVar[str] = "OutputSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.OutputSettings

    amplitude_range: Optional[Union[dict, "QuantityValue"]] = None
    output_active: Optional[Union[bool, Bool]] = None
    output_add: Optional[Union[bool, Bool]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.amplitude_range is not None and not isinstance(self.amplitude_range, QuantityValue):
            self.amplitude_range = QuantityValue(**as_dict(self.amplitude_range))

        if self.output_active is not None and not isinstance(self.output_active, Bool):
            self.output_active = Bool(self.output_active)

        if self.output_add is not None and not isinstance(self.output_add, Bool):
            self.output_add = Bool(self.output_add)

        super().__post_init__(**kwargs)


@dataclass
class PLLSetupSettings(YAMLRoot):
    """
    PLL setup settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["PLLSetupSettings"]
    class_class_curie: ClassVar[str] = "materialsML:PLLSetupSettings"
    class_name: ClassVar[str] = "PLLSetupSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.PLLSetupSettings

    q_factor: Optional[float] = None
    demodulation_bw_amp: Optional[Union[dict, "QuantityValue"]] = None
    demodulation_bw_pha: Optional[Union[dict, "QuantityValue"]] = None
    amplitude_excitation_ratio: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.q_factor is not None and not isinstance(self.q_factor, float):
            self.q_factor = float(self.q_factor)

        if self.demodulation_bw_amp is not None and not isinstance(self.demodulation_bw_amp, QuantityValue):
            self.demodulation_bw_amp = QuantityValue(**as_dict(self.demodulation_bw_amp))

        if self.demodulation_bw_pha is not None and not isinstance(self.demodulation_bw_pha, QuantityValue):
            self.demodulation_bw_pha = QuantityValue(**as_dict(self.demodulation_bw_pha))

        if self.amplitude_excitation_ratio is not None and not isinstance(self.amplitude_excitation_ratio, QuantityValue):
            self.amplitude_excitation_ratio = QuantityValue(**as_dict(self.amplitude_excitation_ratio))

        super().__post_init__(**kwargs)


@dataclass
class PeriodicBoundaryConditions(YAMLRoot):
    """
    Periodic boundary conditions
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["PeriodicBoundaryConditions"]
    class_class_curie: ClassVar[str] = "materialsML:PeriodicBoundaryConditions"
    class_name: ClassVar[str] = "PeriodicBoundaryConditions"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.PeriodicBoundaryConditions

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    z_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        if self.z_coordinate is not None and not isinstance(self.z_coordinate, QuantityValue):
            self.z_coordinate = QuantityValue(**as_dict(self.z_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class AtomsPositions(YAMLRoot):
    """
    Atoms positions and atoms symbols
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["AtomsPositions"]
    class_class_curie: ClassVar[str] = "materialsML:AtomsPositions"
    class_name: ClassVar[str] = "AtomsPositions"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.AtomsPositions

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    z_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    atom_symbol: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        if self.z_coordinate is not None and not isinstance(self.z_coordinate, QuantityValue):
            self.z_coordinate = QuantityValue(**as_dict(self.z_coordinate))

        if self.atom_symbol is not None and not isinstance(self.atom_symbol, str):
            self.atom_symbol = str(self.atom_symbol)

        super().__post_init__(**kwargs)


@dataclass
class KPointsConditions(YAMLRoot):
    """
    K-points conditions
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["KPointsConditions"]
    class_class_curie: ClassVar[str] = "materialsML:KPointsConditions"
    class_name: ClassVar[str] = "KPointsConditions"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.KPointsConditions

    x_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    y_coordinate: Optional[Union[dict, "QuantityValue"]] = None
    z_coordinate: Optional[Union[dict, "QuantityValue"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.x_coordinate is not None and not isinstance(self.x_coordinate, QuantityValue):
            self.x_coordinate = QuantityValue(**as_dict(self.x_coordinate))

        if self.y_coordinate is not None and not isinstance(self.y_coordinate, QuantityValue):
            self.y_coordinate = QuantityValue(**as_dict(self.y_coordinate))

        if self.z_coordinate is not None and not isinstance(self.z_coordinate, QuantityValue):
            self.z_coordinate = QuantityValue(**as_dict(self.z_coordinate))

        super().__post_init__(**kwargs)


@dataclass
class BSEnergy(YAMLRoot):
    """
    Band structure energy
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["BSEnergy"]
    class_class_curie: ClassVar[str] = "materialsML:BSEnergy"
    class_name: ClassVar[str] = "BSEnergy"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.BSEnergy

    charge: Optional[float] = None
    magnetic_momentum: Optional[float] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.charge is not None and not isinstance(self.charge, float):
            self.charge = float(self.charge)

        if self.magnetic_momentum is not None and not isinstance(self.magnetic_momentum, float):
            self.magnetic_momentum = float(self.magnetic_momentum)

        super().__post_init__(**kwargs)


@dataclass
class EvaporatorSlot(YAMLRoot):
    """
    Evaporator slot
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["EvaporatorSlot"]
    class_class_curie: ClassVar[str] = "materialsML:EvaporatorSlot"
    class_name: ClassVar[str] = "EvaporatorSlot"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.EvaporatorSlot

    evaporator_number: Optional[int] = None
    details: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.evaporator_number is not None and not isinstance(self.evaporator_number, int):
            self.evaporator_number = int(self.evaporator_number)

        if self.details is not None and not isinstance(self.details, str):
            self.details = str(self.details)

        super().__post_init__(**kwargs)


@dataclass
class OutputPopulationAnalysis(YAMLRoot):
    """
    Output population analysis details
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["OutputPopulationAnalysis"]
    class_class_curie: ClassVar[str] = "materialsML:OutputPopulationAnalysis"
    class_name: ClassVar[str] = "OutputPopulationAnalysis"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.OutputPopulationAnalysis

    charge: Optional[float] = None
    magnetic_momentum: Optional[float] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.charge is not None and not isinstance(self.charge, float):
            self.charge = float(self.charge)

        if self.magnetic_momentum is not None and not isinstance(self.magnetic_momentum, float):
            self.magnetic_momentum = float(self.magnetic_momentum)

        super().__post_init__(**kwargs)


@dataclass
class ForceConvergenceThreshold(YAMLRoot):
    """
    Force convergence threshold
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["ForceConvergenceThreshold"]
    class_class_curie: ClassVar[str] = "materialsML:ForceConvergenceThreshold"
    class_name: ClassVar[str] = "ForceConvergenceThreshold"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.ForceConvergenceThreshold

    charge: Optional[float] = None
    magnetic_momentum: Optional[float] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.charge is not None and not isinstance(self.charge, float):
            self.charge = float(self.charge)

        if self.magnetic_momentum is not None and not isinstance(self.magnetic_momentum, float):
            self.magnetic_momentum = float(self.magnetic_momentum)

        super().__post_init__(**kwargs)


@dataclass
class Layers2DDetails(YAMLRoot):
    """
    details of 2D layers
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["Layers2DDetails"]
    class_class_curie: ClassVar[str] = "materialsML:Layers2DDetails"
    class_name: ClassVar[str] = "Layers2DDetails"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Layers2DDetails

    material: Optional[str] = None
    number_layers: Optional[int] = None
    dopants: Optional[str] = None
    coverage: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.material is not None and not isinstance(self.material, str):
            self.material = str(self.material)

        if self.number_layers is not None and not isinstance(self.number_layers, int):
            self.number_layers = int(self.number_layers)

        if self.dopants is not None and not isinstance(self.dopants, str):
            self.dopants = str(self.dopants)

        if self.coverage is not None and not isinstance(self.coverage, str):
            self.coverage = str(self.coverage)

        super().__post_init__(**kwargs)


@dataclass
class SubstratesDetails(YAMLRoot):
    """
    details of substrates
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["SubstratesDetails"]
    class_class_curie: ClassVar[str] = "materialsML:SubstratesDetails"
    class_name: ClassVar[str] = "SubstratesDetails"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.SubstratesDetails

    material: Optional[str] = None
    thickness: Optional[Union[dict, "QuantityValue"]] = None
    dopants: Optional[str] = None
    coverage: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.material is not None and not isinstance(self.material, str):
            self.material = str(self.material)

        if self.thickness is not None and not isinstance(self.thickness, QuantityValue):
            self.thickness = QuantityValue(**as_dict(self.thickness))

        if self.dopants is not None and not isinstance(self.dopants, str):
            self.dopants = str(self.dopants)

        if self.coverage is not None and not isinstance(self.coverage, str):
            self.coverage = str(self.coverage)

        super().__post_init__(**kwargs)


@dataclass
class QuantityValue(YAMLRoot):
    """
    Value together with unit
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = QUDT["quantityValue"]
    class_class_curie: ClassVar[str] = "qudt:quantityValue"
    class_name: ClassVar[str] = "QuantityValue"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.QuantityValue

    has_value: Optional[float] = None
    has_unit: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.has_value is not None and not isinstance(self.has_value, float):
            self.has_value = float(self.has_value)

        if self.has_unit is not None and not isinstance(self.has_unit, str):
            self.has_unit = str(self.has_unit)

        super().__post_init__(**kwargs)


@dataclass
class EvaporationTemperature(YAMLRoot):
    """
    evaporation temperature details
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["EvaporationTemperature"]
    class_class_curie: ClassVar[str] = "materialsML:EvaporationTemperature"
    class_name: ClassVar[str] = "EvaporationTemperature"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.EvaporationTemperature

    temperature: Optional[Union[dict, QuantityValue]] = None
    device: Optional[str] = None
    datetime: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.temperature is not None and not isinstance(self.temperature, QuantityValue):
            self.temperature = QuantityValue(**as_dict(self.temperature))

        if self.device is not None and not isinstance(self.device, str):
            self.device = str(self.device)

        if self.datetime is not None and not isinstance(self.datetime, str):
            self.datetime = str(self.datetime)

        super().__post_init__(**kwargs)


@dataclass
class VacuumSystemStatusSettings(YAMLRoot):
    """
    vacuum system settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["VacuumSystemStatusSettings"]
    class_class_curie: ClassVar[str] = "materialsML:VacuumSystemStatusSettings"
    class_name: ClassVar[str] = "VacuumSystemStatusSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.VacuumSystemStatusSettings

    principal_chamber: Optional[Union[dict, "PrincipalChamberSettings"]] = None
    analysis_chamber: Optional[Union[dict, "AnalysisChamberSettings"]] = None
    fel_chamber: Optional[Union[dict, "FELChamberSettings"]] = None
    comments: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.principal_chamber is not None and not isinstance(self.principal_chamber, PrincipalChamberSettings):
            self.principal_chamber = PrincipalChamberSettings(**as_dict(self.principal_chamber))

        if self.analysis_chamber is not None and not isinstance(self.analysis_chamber, AnalysisChamberSettings):
            self.analysis_chamber = AnalysisChamberSettings(**as_dict(self.analysis_chamber))

        if self.fel_chamber is not None and not isinstance(self.fel_chamber, FELChamberSettings):
            self.fel_chamber = FELChamberSettings(**as_dict(self.fel_chamber))

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        super().__post_init__(**kwargs)


@dataclass
class CryoSystemStatusSettings(YAMLRoot):
    """
    cryo system settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["CryoSystemStatusSettings"]
    class_class_curie: ClassVar[str] = "materialsML:CryoSystemStatusSettings"
    class_name: ClassVar[str] = "CryoSystemStatusSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.CryoSystemStatusSettings

    temperature: Optional[Union[dict, QuantityValue]] = None
    liquid_helium: Optional[Union[dict, QuantityValue]] = None
    liquid_nitrogen_refill_date: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.temperature is not None and not isinstance(self.temperature, QuantityValue):
            self.temperature = QuantityValue(**as_dict(self.temperature))

        if self.liquid_helium is not None and not isinstance(self.liquid_helium, QuantityValue):
            self.liquid_helium = QuantityValue(**as_dict(self.liquid_helium))

        if self.liquid_nitrogen_refill_date is not None and not isinstance(self.liquid_nitrogen_refill_date, str):
            self.liquid_nitrogen_refill_date = str(self.liquid_nitrogen_refill_date)

        super().__post_init__(**kwargs)


@dataclass
class PrincipalChamberSettings(YAMLRoot):
    """
    principal chamber settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["PrincipalChamberSettings"]
    class_class_curie: ClassVar[str] = "materialsML:PrincipalChamberSettings"
    class_name: ClassVar[str] = "PrincipalChamberSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.PrincipalChamberSettings

    pressure: Optional[Union[dict, QuantityValue]] = None
    tsp_filament_used: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.pressure is not None and not isinstance(self.pressure, QuantityValue):
            self.pressure = QuantityValue(**as_dict(self.pressure))

        if self.tsp_filament_used is not None and not isinstance(self.tsp_filament_used, str):
            self.tsp_filament_used = str(self.tsp_filament_used)

        super().__post_init__(**kwargs)


@dataclass
class AnalysisChamberSettings(YAMLRoot):
    """
    analysis chamber settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["AnalysisChamberSettings"]
    class_class_curie: ClassVar[str] = "materialsML:AnalysisChamberSettings"
    class_name: ClassVar[str] = "AnalysisChamberSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.AnalysisChamberSettings

    pressure: Optional[Union[dict, QuantityValue]] = None
    tsp_filament_used: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.pressure is not None and not isinstance(self.pressure, QuantityValue):
            self.pressure = QuantityValue(**as_dict(self.pressure))

        if self.tsp_filament_used is not None and not isinstance(self.tsp_filament_used, str):
            self.tsp_filament_used = str(self.tsp_filament_used)

        super().__post_init__(**kwargs)


@dataclass
class FELChamberSettings(YAMLRoot):
    """
    FEL chamber settings
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["FELChamberSettings"]
    class_class_curie: ClassVar[str] = "materialsML:FELChamberSettings"
    class_name: ClassVar[str] = "FELChamberSettings"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.FELChamberSettings

    pressure: Optional[Union[dict, QuantityValue]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.pressure is not None and not isinstance(self.pressure, QuantityValue):
            self.pressure = QuantityValue(**as_dict(self.pressure))

        super().__post_init__(**kwargs)


@dataclass
class Container(YAMLRoot):
    """
    Contains multiple objects
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MATERIALSML["Container"]
    class_class_curie: ClassVar[str] = "materialsML:Container"
    class_name: ClassVar[str] = "Container"
    class_model_uri: ClassVar[URIRef] = MATERIALSML.Container

    afms: Optional[Union[Union[dict, AFM], List[Union[dict, AFM]]]] = empty_list()
    annealings: Optional[Union[Union[dict, Annealing], List[Union[dict, Annealing]]]] = empty_list()
    atomisticmodels: Optional[Union[Union[dict, AtomisticModel], List[Union[dict, AtomisticModel]]]] = empty_list()
    authors: Optional[Union[Union[dict, Author], List[Union[dict, Author]]]] = empty_list()
    band_structures: Optional[Union[Union[dict, BandStructure], List[Union[dict, BandStructure]]]] = empty_list()
    can_gas_bottles: Optional[Union[Union[dict, CanGasBottle], List[Union[dict, CanGasBottle]]]] = empty_list()
    chemists: Optional[Union[Union[dict, Chemist], List[Union[dict, Chemist]]]] = empty_list()
    crystals: Optional[Union[Union[dict, Crystal], List[Union[dict, Crystal]]]] = empty_list()
    depositions: Optional[Union[Union[dict, Deposition], List[Union[dict, Deposition]]]] = empty_list()
    dfts: Optional[Union[Union[dict, DFT], List[Union[dict, DFT]]]] = empty_list()
    dosings: Optional[Union[Union[dict, Dosing], List[Union[dict, Dosing]]]] = empty_list()
    drafts: Optional[Union[Union[dict, Draft], List[Union[dict, Draft]]]] = empty_list()
    empirical_modellings: Optional[Union[Union[dict, EmpiricalModelling], List[Union[dict, EmpiricalModelling]]]] = empty_list()
    errors: Optional[Union[Union[dict, Errors], List[Union[dict, Errors]]]] = empty_list()
    fillcryostats: Optional[Union[Union[dict, FillCryostat], List[Union[dict, FillCryostat]]]] = empty_list()
    geoopts: Optional[Union[Union[dict, GeometryOptimisation], List[Union[dict, GeometryOptimisation]]]] = empty_list()
    grants: Optional[Union[Union[dict, Grant], List[Union[dict, Grant]]]] = empty_list()
    hydrogen_crackers: Optional[Union[Union[dict, HydrogenCracker], List[Union[dict, HydrogenCracker]]]] = empty_list()
    institutions: Optional[Union[Union[dict, Institution], List[Union[dict, Institution]]]] = empty_list()
    instruments: Optional[Union[Union[dict, Instrument], List[Union[dict, Instrument]]]] = empty_list()
    layered2dmaterials: Optional[Union[Union[dict, Layered2DMaterial], List[Union[dict, Layered2DMaterial]]]] = empty_list()
    maintenances: Optional[Union[Union[dict, Maintenance], List[Union[dict, Maintenance]]]] = empty_list()
    manufacturers: Optional[Union[Union[dict, Manufacturer], List[Union[dict, Manufacturer]]]] = empty_list()
    meps: Optional[Union[Union[dict, MinimumEnergyPotential], List[Union[dict, MinimumEnergyPotential]]]] = empty_list()
    molecules: Optional[Union[Union[dict, Molecule], List[Union[dict, Molecule]]]] = empty_list()
    notes: Optional[Union[Union[dict, Notes], List[Union[dict, Notes]]]] = empty_list()
    pdoss: Optional[Union[Union[dict, PDOS], List[Union[dict, PDOS]]]] = empty_list()
    persons: Optional[Union[Union[dict, Person], List[Union[dict, Person]]]] = empty_list()
    protocols: Optional[Union[Union[dict, Protocol], List[Union[dict, Protocol]]]] = empty_list()
    publications: Optional[Union[Union[dict, Publication], List[Union[dict, Publication]]]] = empty_list()
    results: Optional[Union[Union[dict, Result], List[Union[dict, Result]]]] = empty_list()
    rooms: Optional[Union[Union[dict, Room], List[Union[dict, Room]]]] = empty_list()
    samples: Optional[Union[Union[dict, Sample], List[Union[dict, Sample]]]] = empty_list()
    softwares: Optional[Union[Union[dict, Software], List[Union[dict, Software]]]] = empty_list()
    sputterings: Optional[Union[Union[dict, Sputtering], List[Union[dict, Sputtering]]]] = empty_list()
    srds: Optional[Union[Union[dict, SRD], List[Union[dict, SRD]]]] = empty_list()
    stms: Optional[Union[Union[dict, STM], List[Union[dict, STM]]]] = empty_list()
    storages: Optional[Union[Union[dict, Storage], List[Union[dict, Storage]]]] = empty_list()
    stss: Optional[Union[Union[dict, STS], List[Union[dict, STS]]]] = empty_list()
    suppliers: Optional[Union[Union[dict, Supplier], List[Union[dict, Supplier]]]] = empty_list()
    tight_bindings: Optional[Union[Union[dict, TightBinding], List[Union[dict, TightBinding]]]] = empty_list()
    tip_preparations: Optional[Union[Union[dict, TipPreparation], List[Union[dict, TipPreparation]]]] = empty_list()
    tip_sensors: Optional[Union[Union[dict, TipSensor], List[Union[dict, TipSensor]]]] = empty_list()
    transfers: Optional[Union[Union[dict, Transfer], List[Union[dict, Transfer]]]] = empty_list()
    uhv_components: Optional[Union[Union[dict, UHVComponent], List[Union[dict, UHVComponent]]]] = empty_list()
    wafer_substrates: Optional[Union[Union[dict, WaferSubstrate], List[Union[dict, WaferSubstrate]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if not isinstance(self.afms, list):
            self.afms = [self.afms] if self.afms is not None else []
        self.afms = [v if isinstance(v, AFM) else AFM(**as_dict(v)) for v in self.afms]

        if not isinstance(self.annealings, list):
            self.annealings = [self.annealings] if self.annealings is not None else []
        self.annealings = [v if isinstance(v, Annealing) else Annealing(**as_dict(v)) for v in self.annealings]

        if not isinstance(self.atomisticmodels, list):
            self.atomisticmodels = [self.atomisticmodels] if self.atomisticmodels is not None else []
        self.atomisticmodels = [v if isinstance(v, AtomisticModel) else AtomisticModel(**as_dict(v)) for v in self.atomisticmodels]

        if not isinstance(self.authors, list):
            self.authors = [self.authors] if self.authors is not None else []
        self.authors = [v if isinstance(v, Author) else Author(**as_dict(v)) for v in self.authors]

        if not isinstance(self.band_structures, list):
            self.band_structures = [self.band_structures] if self.band_structures is not None else []
        self.band_structures = [v if isinstance(v, BandStructure) else BandStructure(**as_dict(v)) for v in self.band_structures]

        if not isinstance(self.can_gas_bottles, list):
            self.can_gas_bottles = [self.can_gas_bottles] if self.can_gas_bottles is not None else []
        self.can_gas_bottles = [v if isinstance(v, CanGasBottle) else CanGasBottle(**as_dict(v)) for v in self.can_gas_bottles]

        if not isinstance(self.chemists, list):
            self.chemists = [self.chemists] if self.chemists is not None else []
        self.chemists = [v if isinstance(v, Chemist) else Chemist(**as_dict(v)) for v in self.chemists]

        if not isinstance(self.crystals, list):
            self.crystals = [self.crystals] if self.crystals is not None else []
        self.crystals = [v if isinstance(v, Crystal) else Crystal(**as_dict(v)) for v in self.crystals]

        if not isinstance(self.depositions, list):
            self.depositions = [self.depositions] if self.depositions is not None else []
        self.depositions = [v if isinstance(v, Deposition) else Deposition(**as_dict(v)) for v in self.depositions]

        if not isinstance(self.dfts, list):
            self.dfts = [self.dfts] if self.dfts is not None else []
        self.dfts = [v if isinstance(v, DFT) else DFT(**as_dict(v)) for v in self.dfts]

        if not isinstance(self.dosings, list):
            self.dosings = [self.dosings] if self.dosings is not None else []
        self.dosings = [v if isinstance(v, Dosing) else Dosing(**as_dict(v)) for v in self.dosings]

        if not isinstance(self.drafts, list):
            self.drafts = [self.drafts] if self.drafts is not None else []
        self.drafts = [v if isinstance(v, Draft) else Draft(**as_dict(v)) for v in self.drafts]

        if not isinstance(self.empirical_modellings, list):
            self.empirical_modellings = [self.empirical_modellings] if self.empirical_modellings is not None else []
        self.empirical_modellings = [v if isinstance(v, EmpiricalModelling) else EmpiricalModelling(**as_dict(v)) for v in self.empirical_modellings]

        if not isinstance(self.errors, list):
            self.errors = [self.errors] if self.errors is not None else []
        self.errors = [v if isinstance(v, Errors) else Errors(**as_dict(v)) for v in self.errors]

        if not isinstance(self.fillcryostats, list):
            self.fillcryostats = [self.fillcryostats] if self.fillcryostats is not None else []
        self.fillcryostats = [v if isinstance(v, FillCryostat) else FillCryostat(**as_dict(v)) for v in self.fillcryostats]

        if not isinstance(self.geoopts, list):
            self.geoopts = [self.geoopts] if self.geoopts is not None else []
        self.geoopts = [v if isinstance(v, GeometryOptimisation) else GeometryOptimisation(**as_dict(v)) for v in self.geoopts]

        if not isinstance(self.grants, list):
            self.grants = [self.grants] if self.grants is not None else []
        self.grants = [v if isinstance(v, Grant) else Grant(**as_dict(v)) for v in self.grants]

        if not isinstance(self.hydrogen_crackers, list):
            self.hydrogen_crackers = [self.hydrogen_crackers] if self.hydrogen_crackers is not None else []
        self.hydrogen_crackers = [v if isinstance(v, HydrogenCracker) else HydrogenCracker(**as_dict(v)) for v in self.hydrogen_crackers]

        if not isinstance(self.institutions, list):
            self.institutions = [self.institutions] if self.institutions is not None else []
        self.institutions = [v if isinstance(v, Institution) else Institution(**as_dict(v)) for v in self.institutions]

        if not isinstance(self.instruments, list):
            self.instruments = [self.instruments] if self.instruments is not None else []
        self.instruments = [v if isinstance(v, Instrument) else Instrument(**as_dict(v)) for v in self.instruments]

        if not isinstance(self.layered2dmaterials, list):
            self.layered2dmaterials = [self.layered2dmaterials] if self.layered2dmaterials is not None else []
        self.layered2dmaterials = [v if isinstance(v, Layered2DMaterial) else Layered2DMaterial(**as_dict(v)) for v in self.layered2dmaterials]

        if not isinstance(self.maintenances, list):
            self.maintenances = [self.maintenances] if self.maintenances is not None else []
        self.maintenances = [v if isinstance(v, Maintenance) else Maintenance(**as_dict(v)) for v in self.maintenances]

        if not isinstance(self.manufacturers, list):
            self.manufacturers = [self.manufacturers] if self.manufacturers is not None else []
        self.manufacturers = [v if isinstance(v, Manufacturer) else Manufacturer(**as_dict(v)) for v in self.manufacturers]

        if not isinstance(self.meps, list):
            self.meps = [self.meps] if self.meps is not None else []
        self.meps = [v if isinstance(v, MinimumEnergyPotential) else MinimumEnergyPotential(**as_dict(v)) for v in self.meps]

        if not isinstance(self.molecules, list):
            self.molecules = [self.molecules] if self.molecules is not None else []
        self.molecules = [v if isinstance(v, Molecule) else Molecule(**as_dict(v)) for v in self.molecules]

        if not isinstance(self.notes, list):
            self.notes = [self.notes] if self.notes is not None else []
        self.notes = [v if isinstance(v, Notes) else Notes(**as_dict(v)) for v in self.notes]

        if not isinstance(self.pdoss, list):
            self.pdoss = [self.pdoss] if self.pdoss is not None else []
        self.pdoss = [v if isinstance(v, PDOS) else PDOS(**as_dict(v)) for v in self.pdoss]

        if not isinstance(self.persons, list):
            self.persons = [self.persons] if self.persons is not None else []
        self.persons = [v if isinstance(v, Person) else Person(**as_dict(v)) for v in self.persons]

        if not isinstance(self.protocols, list):
            self.protocols = [self.protocols] if self.protocols is not None else []
        self.protocols = [v if isinstance(v, Protocol) else Protocol(**as_dict(v)) for v in self.protocols]

        if not isinstance(self.publications, list):
            self.publications = [self.publications] if self.publications is not None else []
        self.publications = [v if isinstance(v, Publication) else Publication(**as_dict(v)) for v in self.publications]

        if not isinstance(self.results, list):
            self.results = [self.results] if self.results is not None else []
        self.results = [v if isinstance(v, Result) else Result(**as_dict(v)) for v in self.results]

        if not isinstance(self.rooms, list):
            self.rooms = [self.rooms] if self.rooms is not None else []
        self.rooms = [v if isinstance(v, Room) else Room(**as_dict(v)) for v in self.rooms]

        if not isinstance(self.samples, list):
            self.samples = [self.samples] if self.samples is not None else []
        self.samples = [v if isinstance(v, Sample) else Sample(**as_dict(v)) for v in self.samples]

        if not isinstance(self.softwares, list):
            self.softwares = [self.softwares] if self.softwares is not None else []
        self.softwares = [v if isinstance(v, Software) else Software(**as_dict(v)) for v in self.softwares]

        if not isinstance(self.sputterings, list):
            self.sputterings = [self.sputterings] if self.sputterings is not None else []
        self.sputterings = [v if isinstance(v, Sputtering) else Sputtering(**as_dict(v)) for v in self.sputterings]

        if not isinstance(self.srds, list):
            self.srds = [self.srds] if self.srds is not None else []
        self.srds = [v if isinstance(v, SRD) else SRD(**as_dict(v)) for v in self.srds]

        if not isinstance(self.stms, list):
            self.stms = [self.stms] if self.stms is not None else []
        self.stms = [v if isinstance(v, STM) else STM(**as_dict(v)) for v in self.stms]

        if not isinstance(self.storages, list):
            self.storages = [self.storages] if self.storages is not None else []
        self.storages = [v if isinstance(v, Storage) else Storage(**as_dict(v)) for v in self.storages]

        if not isinstance(self.stss, list):
            self.stss = [self.stss] if self.stss is not None else []
        self.stss = [v if isinstance(v, STS) else STS(**as_dict(v)) for v in self.stss]

        if not isinstance(self.suppliers, list):
            self.suppliers = [self.suppliers] if self.suppliers is not None else []
        self.suppliers = [v if isinstance(v, Supplier) else Supplier(**as_dict(v)) for v in self.suppliers]

        if not isinstance(self.tight_bindings, list):
            self.tight_bindings = [self.tight_bindings] if self.tight_bindings is not None else []
        self.tight_bindings = [v if isinstance(v, TightBinding) else TightBinding(**as_dict(v)) for v in self.tight_bindings]

        if not isinstance(self.tip_preparations, list):
            self.tip_preparations = [self.tip_preparations] if self.tip_preparations is not None else []
        self.tip_preparations = [v if isinstance(v, TipPreparation) else TipPreparation(**as_dict(v)) for v in self.tip_preparations]

        if not isinstance(self.tip_sensors, list):
            self.tip_sensors = [self.tip_sensors] if self.tip_sensors is not None else []
        self.tip_sensors = [v if isinstance(v, TipSensor) else TipSensor(**as_dict(v)) for v in self.tip_sensors]

        if not isinstance(self.transfers, list):
            self.transfers = [self.transfers] if self.transfers is not None else []
        self.transfers = [v if isinstance(v, Transfer) else Transfer(**as_dict(v)) for v in self.transfers]

        if not isinstance(self.uhv_components, list):
            self.uhv_components = [self.uhv_components] if self.uhv_components is not None else []
        self.uhv_components = [v if isinstance(v, UHVComponent) else UHVComponent(**as_dict(v)) for v in self.uhv_components]

        if not isinstance(self.wafer_substrates, list):
            self.wafer_substrates = [self.wafer_substrates] if self.wafer_substrates is not None else []
        self.wafer_substrates = [v if isinstance(v, WaferSubstrate) else WaferSubstrate(**as_dict(v)) for v in self.wafer_substrates]

        super().__post_init__(**kwargs)


# Enumerations
class PModelEnum(EnumDefinitionImpl):
    """
    p-model values
    """
    single_probe = PermissibleValue(
        text="single_probe",
        description="single probe")
    two_pp_resp_ptcta = PermissibleValue(
        text="two_pp_resp_ptcta",
        description="2PP RESP PTCTA")
    two_pp_resp_pentacene = PermissibleValue(
        text="two_pp_resp_pentacene",
        description="2PP RESP pentacene")

    _defn = EnumDefinition(
        name="PModelEnum",
        description="p-model values",
    )

class DirectionSpeedEnum(EnumDefinitionImpl):
    """
    speed direction values
    """
    down = PermissibleValue(
        text="down",
        description="down")
    both = PermissibleValue(
        text="both",
        description="both")

    _defn = EnumDefinition(
        name="DirectionSpeedEnum",
        description="speed direction values",
    )

class ActiveCalibEnum(EnumDefinitionImpl):
    """
    active calib values
    """
    l_he = PermissibleValue(
        text="l_he",
        description="LHe")

    _defn = EnumDefinition(
        name="ActiveCalibEnum",
        description="active calib values",
    )

class FeedbackTypeEnum(EnumDefinitionImpl):
    """
    feedback type values
    """
    log_current = PermissibleValue(
        text="log_current",
        description="log current")
    abs_current = PermissibleValue(
        text="abs_current",
        description="abs current")

    _defn = EnumDefinition(
        name="FeedbackTypeEnum",
        description="feedback type values",
    )

class XCFunctionalEnum(EnumDefinitionImpl):
    """
    xc functional values
    """
    exchange_correlation_function = PermissibleValue(
        text="exchange_correlation_function",
        description="exchange correlation function")

    _defn = EnumDefinition(
        name="XCFunctionalEnum",
        description="xc functional values",
    )

class DraftTypeEnum(EnumDefinitionImpl):
    """
    draft type values
    """
    preprint = PermissibleValue(
        text="preprint",
        description="preprint")
    postprint = PermissibleValue(
        text="postprint",
        description="postprint")

    _defn = EnumDefinition(
        name="DraftTypeEnum",
        description="draft type values",
    )

class DewarEnum(EnumDefinitionImpl):
    """
    dewar values
    """
    ln2 = PermissibleValue(
        text="ln2",
        description="dewar liquid molecular nitrogen")

    _defn = EnumDefinition(
        name="DewarEnum",
        description="dewar values",
    )

    @classmethod
    def _addvals(cls):
        setattr(cls, "1",
            PermissibleValue(
                text="1",
                description="dewar 1"))
        setattr(cls, "2",
            PermissibleValue(
                text="2",
                description="dewar 2"))
        setattr(cls, "3",
            PermissibleValue(
                text="3",
                description="dewar 3"))
        setattr(cls, "4",
            PermissibleValue(
                text="4",
                description="dewar 4"))
        setattr(cls, "5",
            PermissibleValue(
                text="5",
                description="dewar 5"))
        setattr(cls, "6",
            PermissibleValue(
                text="6",
                description="dewar 6"))
        setattr(cls, "7",
            PermissibleValue(
                text="7",
                description="dewar 7"))
        setattr(cls, "8",
            PermissibleValue(
                text="8",
                description="dewar 8"))
        setattr(cls, "9",
            PermissibleValue(
                text="9",
                description="dewar 9"))
        setattr(cls, "10",
            PermissibleValue(
                text="10",
                description="dewar 10"))
        setattr(cls, "11",
            PermissibleValue(
                text="11",
                description="dewar 11"))
        setattr(cls, "12",
            PermissibleValue(
                text="12",
                description="dewar 12"))

class MethodTypeEnum(EnumDefinitionImpl):
    """
    method type values
    """
    neb = PermissibleValue(
        text="neb",
        description="NEB")
    string_method = PermissibleValue(
        text="string_method",
        description="string method")

    _defn = EnumDefinition(
        name="MethodTypeEnum",
        description="method type values",
    )

class PlotGroupByEnum(EnumDefinitionImpl):
    """
    group by (plot parameter) values
    """
    index = PermissibleValue(
        text="index",
        description="index")
    element = PermissibleValue(
        text="element",
        description="element")

    _defn = EnumDefinition(
        name="PlotGroupByEnum",
        description="group by (plot parameter) values",
    )

class PlotContributionsEnum(EnumDefinitionImpl):
    """
    contributions (plot parameter) values
    """
    total = PermissibleValue(
        text="total",
        description="total")
    orbital = PermissibleValue(
        text="orbital",
        description="orbital")
    angular_momentum = PermissibleValue(
        text="angular_momentum",
        description="angular momentum")

    _defn = EnumDefinition(
        name="PlotContributionsEnum",
        description="contributions (plot parameter) values",
    )

class ExperimentTypeEnum(EnumDefinitionImpl):
    """
    experiment type values
    """
    bias_spectroscopy = PermissibleValue(
        text="bias_spectroscopy",
        description="bias spectroscopy")
    z_spectroscopy = PermissibleValue(
        text="z_spectroscopy",
        description="z-spectroscopy")

    _defn = EnumDefinition(
        name="ExperimentTypeEnum",
        description="experiment type values",
    )

class FilterTypeEnum(EnumDefinitionImpl):
    """
    filter type values
    """
    butterworth = PermissibleValue(
        text="butterworth",
        description="Butterworth")

    _defn = EnumDefinition(
        name="FilterTypeEnum",
        description="filter type values",
    )

class LockInDemodulatedSignalEnum(EnumDefinitionImpl):
    """
    lock-in demodulated signal values
    """
    current = PermissibleValue(
        text="current",
        description="current")

    _defn = EnumDefinition(
        name="LockInDemodulatedSignalEnum",
        description="lock-in demodulated signal values",
    )

class TipPreparationTypeEnum(EnumDefinitionImpl):
    """
    tip preparation type values
    """
    fib = PermissibleValue(
        text="fib",
        description="FIB")
    etched = PermissibleValue(
        text="etched",
        description="Etched")

    _defn = EnumDefinition(
        name="TipPreparationTypeEnum",
        description="tip preparation type values",
    )

class TipSensorTypeEnum(EnumDefinitionImpl):
    """
    tip sensor type values
    """
    tuning_fork = PermissibleValue(
        text="tuning_fork",
        description="tuning fork")

    _defn = EnumDefinition(
        name="TipSensorTypeEnum",
        description="tip sensor type values",
    )

class UHVComponentTypeEnum(EnumDefinitionImpl):
    """
    ultra high vacuum component type values
    """
    pump = PermissibleValue(
        text="pump",
        description="pump")
    evaporator = PermissibleValue(
        text="evaporator",
        description="evaporator")
    analyser = PermissibleValue(
        text="analyser",
        description="analyser")
    gauge = PermissibleValue(
        text="gauge",
        description="gauge")
    ion_photon_stage = PermissibleValue(
        text="ion_photon_stage",
        description="ion/photon stage")
    storage = PermissibleValue(
        text="storage",
        description="storage")
    gas_panel = PermissibleValue(
        text="gas_panel",
        description="gas panel")
    port_valve = PermissibleValue(
        text="port_valve",
        description="port/valve")
    manipulator_wobble_stick = PermissibleValue(
        text="manipulator_wobble_stick",
        description="manipulator/wobble stick")
    cryostat = PermissibleValue(
        text="cryostat",
        description="cryostat")
    spm = PermissibleValue(
        text="spm",
        description="scanning probe microscope")
    electronics = PermissibleValue(
        text="electronics",
        description="electronics")
    various = PermissibleValue(
        text="various",
        description="various")

    _defn = EnumDefinition(
        name="UHVComponentTypeEnum",
        description="ultra high vacuum component type values",
    )

class MassUnitsEnum(EnumDefinitionImpl):
    """
    mass units
    """
    mg = PermissibleValue(
        text="mg",
        description="milligrams",
        meaning=UNIT["MilliGM"])
    g = PermissibleValue(
        text="g",
        description="grams",
        meaning=UNIT["GM"])
    kg = PermissibleValue(
        text="kg",
        description="kilograms",
        meaning=UNIT["KiloGM"])

    _defn = EnumDefinition(
        name="MassUnitsEnum",
        description="mass units",
    )

# Slots
class slots:
    pass

slots.perm_id = Slot(uri=SCHEMA.identifier, name="perm_id", curie=SCHEMA.curie('identifier'),
                   model_uri=MATERIALSML.perm_id, domain=None, range=Optional[str])

slots.name = Slot(uri=SCHEMA.name, name="name", curie=SCHEMA.curie('name'),
                   model_uri=MATERIALSML.name, domain=None, range=Optional[str])

slots.description = Slot(uri=SCHEMA.description, name="description", curie=SCHEMA.curie('description'),
                   model_uri=MATERIALSML.description, domain=None, range=Optional[str])

slots.work_phone = Slot(uri=SCHEMA.telephone, name="work_phone", curie=SCHEMA.curie('telephone'),
                   model_uri=MATERIALSML.work_phone, domain=None, range=Optional[str])

slots.email = Slot(uri=SCHEMA.email, name="email", curie=SCHEMA.curie('email'),
                   model_uri=MATERIALSML.email, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$'))

slots.address = Slot(uri=SCHEMA.address, name="address", curie=SCHEMA.curie('address'),
                   model_uri=MATERIALSML.address, domain=None, range=Optional[str])

slots.iupac_name = Slot(uri=SCHEMA.iupacName, name="iupac_name", curie=SCHEMA.curie('iupacName'),
                   model_uri=MATERIALSML.iupac_name, domain=None, range=Optional[str])

slots.sum_formula = Slot(uri=SCHEMA.molecularFormula, name="sum_formula", curie=SCHEMA.curie('molecularFormula'),
                   model_uri=MATERIALSML.sum_formula, domain=None, range=Optional[str])

slots.smiles = Slot(uri=SCHEMA.smiles, name="smiles", curie=SCHEMA.curie('smiles'),
                   model_uri=MATERIALSML.smiles, domain=None, range=Optional[str])

slots.cas_number = Slot(uri=EMMO.EMMO_d2a47cd8_662f_438f_855a_b4378eb992ff, name="cas_number", curie=EMMO.curie('EMMO_d2a47cd8_662f_438f_855a_b4378eb992ff'),
                   model_uri=MATERIALSML.cas_number, domain=None, range=Optional[str])

slots.empa_number = Slot(uri=MATERIALSML.empa_number, name="empa_number", curie=MATERIALSML.curie('empa_number'),
                   model_uri=MATERIALSML.empa_number, domain=None, range=Optional[int])

slots.batch = Slot(uri=NCIT.NCIT_C67073, name="batch", curie=NCIT.curie('NCIT_C67073'),
                   model_uri=MATERIALSML.batch, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[a-z]$'))

slots.vial = Slot(uri=AFE.AFE_0000329, name="vial", curie=AFE.curie('AFE_0000329'),
                   model_uri=MATERIALSML.vial, domain=None, range=Optional[str],
                   pattern=re.compile(r'^(i|ii|iii|iv|v|vi|vii|viii|ix|x)$'))

slots.hazardous = Slot(uri=NCIT.NCIT_C73538, name="hazardous", curie=NCIT.curie('NCIT_C73538'),
                   model_uri=MATERIALSML.hazardous, domain=None, range=Optional[Union[bool, Bool]])

slots.hazardous_specification = Slot(uri=MATERIALSML.hazardous_specification, name="hazardous_specification", curie=MATERIALSML.curie('hazardous_specification'),
                   model_uri=MATERIALSML.hazardous_specification, domain=None, range=Optional[str])

slots.evaporation_temperatures = Slot(uri=MATERIALSML.evaporation_temperatures, name="evaporation_temperatures", curie=MATERIALSML.curie('evaporation_temperatures'),
                   model_uri=MATERIALSML.evaporation_temperatures, domain=None, range=Optional[Union[Union[dict, EvaporationTemperature], List[Union[dict, EvaporationTemperature]]]])

slots.fridge = Slot(uri=MATERIALSML.fridge, name="fridge", curie=MATERIALSML.curie('fridge'),
                   model_uri=MATERIALSML.fridge, domain=None, range=Optional[Union[bool, Bool]])

slots.no_light = Slot(uri=MATERIALSML.no_light, name="no_light", curie=MATERIALSML.curie('no_light'),
                   model_uri=MATERIALSML.no_light, domain=None, range=Optional[Union[bool, Bool]])

slots.dry = Slot(uri=MATERIALSML.dry, name="dry", curie=MATERIALSML.curie('dry'),
                   model_uri=MATERIALSML.dry, domain=None, range=Optional[Union[bool, Bool]])

slots.no_oxygen = Slot(uri=MATERIALSML.no_oxygen, name="no_oxygen", curie=MATERIALSML.curie('no_oxygen'),
                   model_uri=MATERIALSML.no_oxygen, domain=None, range=Optional[Union[bool, Bool]])

slots.other_storage_condition = Slot(uri=MATERIALSML.other_storage_condition, name="other_storage_condition", curie=MATERIALSML.curie('other_storage_condition'),
                   model_uri=MATERIALSML.other_storage_condition, domain=None, range=Optional[Union[bool, Bool]])

slots.other_storage_condition_specification = Slot(uri=MATERIALSML.other_storage_condition_specification, name="other_storage_condition_specification", curie=MATERIALSML.curie('other_storage_condition_specification'),
                   model_uri=MATERIALSML.other_storage_condition_specification, domain=None, range=Optional[str])

slots.chemist_molecule_name = Slot(uri=SCHEMA.alternateName, name="chemist_molecule_name", curie=SCHEMA.curie('alternateName'),
                   model_uri=MATERIALSML.chemist_molecule_name, domain=None, range=Optional[str])

slots.amount = Slot(uri=MATERIALSML.amount, name="amount", curie=MATERIALSML.curie('amount'),
                   model_uri=MATERIALSML.amount, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.receive_date = Slot(uri=NCIT.NCIT_C164590, name="receive_date", curie=NCIT.curie('NCIT_C164590'),
                   model_uri=MATERIALSML.receive_date, domain=None, range=Optional[Union[str, XSDDate]])

slots.comments = Slot(uri=SCHEMA.comment, name="comments", curie=SCHEMA.curie('comment'),
                   model_uri=MATERIALSML.comments, domain=None, range=Optional[str])

slots.has_value = Slot(uri=QUDT.hasQuantityKind, name="has_value", curie=QUDT.curie('hasQuantityKind'),
                   model_uri=MATERIALSML.has_value, domain=None, range=Optional[float])

slots.has_unit = Slot(uri=QUDT.hasUnit, name="has_unit", curie=QUDT.curie('hasUnit'),
                   model_uri=MATERIALSML.has_unit, domain=None, range=Optional[str])

slots.device = Slot(uri=MATERIALSML.device, name="device", curie=MATERIALSML.curie('device'),
                   model_uri=MATERIALSML.device, domain=None, range=Optional[str])

slots.temperature = Slot(uri=NCIT.NCIT_C25206, name="temperature", curie=NCIT.curie('NCIT_C25206'),
                   model_uri=MATERIALSML.temperature, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.datetime = Slot(uri=MATERIALSML.datetime, name="datetime", curie=MATERIALSML.curie('datetime'),
                   model_uri=MATERIALSML.datetime, domain=None, range=Optional[str])

slots.duration = Slot(uri=NCIT.NCIT_C25330, name="duration", curie=NCIT.curie('NCIT_C25330'),
                   model_uri=MATERIALSML.duration, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.pressure = Slot(uri=EMMO.EMMO_50a44256_9dc5_434b_bad4_74a4d9a29989, name="pressure", curie=EMMO.curie('EMMO_50a44256_9dc5_434b_bad4_74a4d9a29989'),
                   model_uri=MATERIALSML.pressure, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.voltage = Slot(uri=EMMO.EMMO_4f2d3939_91b1_4001_b8ab_7d19074bf845, name="voltage", curie=EMMO.curie('EMMO_4f2d3939_91b1_4001_b8ab_7d19074bf845'),
                   model_uri=MATERIALSML.voltage, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.current = Slot(uri=EMMO.EMMO_c995ae70_3b84_4ebb_bcfc_69e6a281bb88, name="current", curie=EMMO.curie('EMMO_c995ae70_3b84_4ebb_bcfc_69e6a281bb88'),
                   model_uri=MATERIALSML.current, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.start_time = Slot(uri=NCIT.NCIT_C78441, name="start_time", curie=NCIT.curie('NCIT_C78441'),
                   model_uri=MATERIALSML.start_time, domain=None, range=Optional[str])

slots.bias_setpoint = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-bias-field'], name="bias_setpoint", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-bias-field'),
                   model_uri=MATERIALSML.bias_setpoint, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.bias_calibration_factor = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-bias-calibration-field'], name="bias_calibration_factor", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-bias-calibration-field'),
                   model_uri=MATERIALSML.bias_calibration_factor, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.bias_calibration_offset = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-bias-offset-field'], name="bias_calibration_offset", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-bias-offset-field'),
                   model_uri=MATERIALSML.bias_calibration_offset, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.current_setpoint = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-current-field'], name="current_setpoint", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-current-field'),
                   model_uri=MATERIALSML.current_setpoint, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.current_calibration_factor = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-current-calibration-field'], name="current_calibration_factor", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-current-calibration-field'),
                   model_uri=MATERIALSML.current_calibration_factor, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.current_calibration_offset = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-current-offset-field'], name="current_calibration_offset", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-current-offset-field'),
                   model_uri=MATERIALSML.current_calibration_offset, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.current_gain = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-current-gain-field'], name="current_gain", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-current-gain-field'),
                   model_uri=MATERIALSML.current_gain, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_position = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-position-z-controller-z-field'], name="z_position", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-position-z-controller-z-field'),
                   model_uri=MATERIALSML.z_position, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.feedback_active = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-controller-status-field'], name="feedback_active", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-controller-status-field'),
                   model_uri=MATERIALSML.feedback_active, domain=None, range=Optional[Union[bool, Bool]])

slots.feedback_type = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-controller-name-field'], name="feedback_type", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-controller-name-field'),
                   model_uri=MATERIALSML.feedback_type, domain=None, range=Optional[Union[str, "FeedbackTypeEnum"]])

slots.piezo_configuration_settings = Slot(uri=MATERIALSML.piezo_configuration_settings, name="piezo_configuration_settings", curie=MATERIALSML.curie('piezo_configuration_settings'),
                   model_uri=MATERIALSML.piezo_configuration_settings, domain=None, range=Optional[Union[dict, PiezoConfigurationSettings]])

slots.active_calib = Slot(uri=MATERIALSML.active_calib, name="active_calib", curie=MATERIALSML.curie('active_calib'),
                   model_uri=MATERIALSML.active_calib, domain=None, range=Optional[Union[str, "ActiveCalibEnum"]])

slots.sensitivity_settings = Slot(uri=MATERIALSML.sensitivity_settings, name="sensitivity_settings", curie=MATERIALSML.curie('sensitivity_settings'),
                   model_uri=MATERIALSML.sensitivity_settings, domain=None, range=Optional[Union[dict, SensitivitySettings]])

slots.hv_gain_settings = Slot(uri=MATERIALSML.hv_gain_settings, name="hv_gain_settings", curie=MATERIALSML.curie('hv_gain_settings'),
                   model_uri=MATERIALSML.hv_gain_settings, domain=None, range=Optional[Union[dict, HVGainSettings]])

slots.scan_slope = Slot(uri=MATERIALSML.scan_slope, name="scan_slope", curie=MATERIALSML.curie('scan_slope'),
                   model_uri=MATERIALSML.scan_slope, domain=None, range=Optional[Union[dict, ScanSlope]])

slots.curvature_radius = Slot(uri=MATERIALSML.curvature_radius, name="curvature_radius", curie=MATERIALSML.curie('curvature_radius'),
                   model_uri=MATERIALSML.curvature_radius, domain=None, range=Optional[Union[dict, CurvatureRadius]])

slots.second_order_correction = Slot(uri=MATERIALSML.second_order_correction, name="second_order_correction", curie=MATERIALSML.curie('second_order_correction'),
                   model_uri=MATERIALSML.second_order_correction, domain=None, range=Optional[Union[dict, SecondOrderCorrection]])

slots.drift_settings = Slot(uri=MATERIALSML.drift_settings, name="drift_settings", curie=MATERIALSML.curie('drift_settings'),
                   model_uri=MATERIALSML.drift_settings, domain=None, range=Optional[Union[dict, DriftSettings]])

slots.drift_correction_status = Slot(uri=MATERIALSML.drift_correction_status, name="drift_correction_status", curie=MATERIALSML.curie('drift_correction_status'),
                   model_uri=MATERIALSML.drift_correction_status, domain=None, range=Optional[Union[bool, Bool]])

slots.scan_settings = Slot(uri=MATERIALSML.scan_settings, name="scan_settings", curie=MATERIALSML.curie('scan_settings'),
                   model_uri=MATERIALSML.scan_settings, domain=None, range=Optional[Union[dict, ScanSettings]])

slots.field = Slot(uri=MATERIALSML.field, name="field", curie=MATERIALSML.curie('field'),
                   model_uri=MATERIALSML.field, domain=None, range=Optional[Union[float, List[float]]])

slots.angle = Slot(uri=EMMO.EMMO_f3dd74c0_f480_49e8_9764_33b78638c235, name="angle", curie=EMMO.curie('EMMO_f3dd74c0_f480_49e8_9764_33b78638c235'),
                   model_uri=MATERIALSML.angle, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.scan_offset = Slot(uri=MATERIALSML.scan_offset, name="scan_offset", curie=MATERIALSML.curie('scan_offset'),
                   model_uri=MATERIALSML.scan_offset, domain=None, range=Optional[Union[dict, ScanOffset]])

slots.x_coordinate = Slot(uri=NCIT.NCIT_C44477, name="x_coordinate", curie=NCIT.curie('NCIT_C44477'),
                   model_uri=MATERIALSML.x_coordinate, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.y_coordinate = Slot(uri=NCIT.NCIT_C44478, name="y_coordinate", curie=NCIT.curie('NCIT_C44478'),
                   model_uri=MATERIALSML.y_coordinate, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_coordinate = Slot(uri=NCIT.NCIT_C44479, name="z_coordinate", curie=NCIT.curie('NCIT_C44479'),
                   model_uri=MATERIALSML.z_coordinate, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.scan_time = Slot(uri=MATERIALSML.scan_time, name="scan_time", curie=MATERIALSML.curie('scan_time'),
                   model_uri=MATERIALSML.scan_time, domain=None, range=Optional[Union[dict, ScanTime]])

slots.scan_range = Slot(uri=MATERIALSML.scan_range, name="scan_range", curie=MATERIALSML.curie('scan_range'),
                   model_uri=MATERIALSML.scan_range, domain=None, range=Optional[Union[dict, ScanRange]])

slots.scan_pixels = Slot(uri=MATERIALSML.scan_pixels, name="scan_pixels", curie=MATERIALSML.curie('scan_pixels'),
                   model_uri=MATERIALSML.scan_pixels, domain=None, range=Optional[Union[dict, ScanPixels]])

slots.scan_speed = Slot(uri=MATERIALSML.scan_speed, name="scan_speed", curie=MATERIALSML.curie('scan_speed'),
                   model_uri=MATERIALSML.scan_speed, domain=None, range=Optional[Union[dict, ScanSpeed]])

slots.forward_speed = Slot(uri=MATERIALSML.forward_speed, name="forward_speed", curie=MATERIALSML.curie('forward_speed'),
                   model_uri=MATERIALSML.forward_speed, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.backward_speed = Slot(uri=MATERIALSML.backward_speed, name="backward_speed", curie=MATERIALSML.curie('backward_speed'),
                   model_uri=MATERIALSML.backward_speed, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.direction_speed = Slot(uri=MATERIALSML.direction_speed, name="direction_speed", curie=MATERIALSML.curie('direction_speed'),
                   model_uri=MATERIALSML.direction_speed, domain=None, range=Optional[Union[str, "DirectionSpeedEnum"]])

slots.oscillation_control_settings = Slot(uri=MATERIALSML.oscillation_control_settings, name="oscillation_control_settings", curie=MATERIALSML.curie('oscillation_control_settings'),
                   model_uri=MATERIALSML.oscillation_control_settings, domain=None, range=Optional[Union[dict, OscillationControlSettings]])

slots.differential_input = Slot(uri=MATERIALSML.differential_input, name="differential_input", curie=MATERIALSML.curie('differential_input'),
                   model_uri=MATERIALSML.differential_input, domain=None, range=Optional[Union[bool, Bool]])

slots.one_ten = Slot(uri=MATERIALSML.one_ten, name="one_ten", curie=MATERIALSML.curie('one_ten'),
                   model_uri=MATERIALSML.one_ten, domain=None, range=Optional[Union[bool, Bool]])

slots.input_calibration = Slot(uri=MATERIALSML.input_calibration, name="input_calibration", curie=MATERIALSML.curie('input_calibration'),
                   model_uri=MATERIALSML.input_calibration, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.input_range = Slot(uri=MATERIALSML.input_range, name="input_range", curie=MATERIALSML.curie('input_range'),
                   model_uri=MATERIALSML.input_range, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.center_frequency = Slot(uri=MATERIALSML.center_frequency, name="center_frequency", curie=MATERIALSML.curie('center_frequency'),
                   model_uri=MATERIALSML.center_frequency, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.range = Slot(uri=MATERIALSML.range, name="range", curie=MATERIALSML.curie('range'),
                   model_uri=MATERIALSML.range, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.demodulation_settings = Slot(uri=MATERIALSML.demodulation_settings, name="demodulation_settings", curie=MATERIALSML.curie('demodulation_settings'),
                   model_uri=MATERIALSML.demodulation_settings, domain=None, range=Optional[Union[dict, DemodulationSettings]])

slots.inputs = Slot(uri=MATERIALSML.inputs, name="inputs", curie=MATERIALSML.curie('inputs'),
                   model_uri=MATERIALSML.inputs, domain=None, range=Optional[Union[float, List[float]]])

slots.frequencies = Slot(uri=MATERIALSML.frequencies, name="frequencies", curie=MATERIALSML.curie('frequencies'),
                   model_uri=MATERIALSML.frequencies, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.reference_phases = Slot(uri=MATERIALSML.reference_phases, name="reference_phases", curie=MATERIALSML.curie('reference_phases'),
                   model_uri=MATERIALSML.reference_phases, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.cutoff_frequencies = Slot(uri=MATERIALSML.cutoff_frequencies, name="cutoff_frequencies", curie=MATERIALSML.curie('cutoff_frequencies'),
                   model_uri=MATERIALSML.cutoff_frequencies, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.harmonics = Slot(uri=MATERIALSML.harmonics, name="harmonics", curie=MATERIALSML.curie('harmonics'),
                   model_uri=MATERIALSML.harmonics, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.filters_orders = Slot(uri=MATERIALSML.filters_orders, name="filters_orders", curie=MATERIALSML.curie('filters_orders'),
                   model_uri=MATERIALSML.filters_orders, domain=None, range=Optional[Union[float, List[float]]])

slots.phase_settings = Slot(uri=MATERIALSML.phase_settings, name="phase_settings", curie=MATERIALSML.curie('phase_settings'),
                   model_uri=MATERIALSML.phase_settings, domain=None, range=Optional[Union[dict, PhaseSettings]])

slots.p_gain = Slot(uri=MATERIALSML.p_gain, name="p_gain", curie=MATERIALSML.curie('p_gain'),
                   model_uri=MATERIALSML.p_gain, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.i_gain = Slot(uri=MATERIALSML.i_gain, name="i_gain", curie=MATERIALSML.curie('i_gain'),
                   model_uri=MATERIALSML.i_gain, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.controller_on = Slot(uri=MATERIALSML.controller_on, name="controller_on", curie=MATERIALSML.curie('controller_on'),
                   model_uri=MATERIALSML.controller_on, domain=None, range=Optional[Union[bool, Bool]])

slots.frequency_shift = Slot(uri=MATERIALSML.frequency_shift, name="frequency_shift", curie=MATERIALSML.curie('frequency_shift'),
                   model_uri=MATERIALSML.frequency_shift, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.amplitude_settings = Slot(uri=MATERIALSML.amplitude_settings, name="amplitude_settings", curie=MATERIALSML.curie('amplitude_settings'),
                   model_uri=MATERIALSML.amplitude_settings, domain=None, range=Optional[Union[dict, AmplitudeSettings]])

slots.setpoint = Slot(uri=MATERIALSML.setpoint, name="setpoint", curie=MATERIALSML.curie('setpoint'),
                   model_uri=MATERIALSML.setpoint, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.excitation = Slot(uri=MATERIALSML.excitation, name="excitation", curie=MATERIALSML.curie('excitation'),
                   model_uri=MATERIALSML.excitation, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.output_settings = Slot(uri=MATERIALSML.output_settings, name="output_settings", curie=MATERIALSML.curie('output_settings'),
                   model_uri=MATERIALSML.output_settings, domain=None, range=Optional[Union[dict, OutputSettings]])

slots.output_active = Slot(uri=MATERIALSML.output_active, name="output_active", curie=MATERIALSML.curie('output_active'),
                   model_uri=MATERIALSML.output_active, domain=None, range=Optional[Union[bool, Bool]])

slots.output_add = Slot(uri=MATERIALSML.output_add, name="output_add", curie=MATERIALSML.curie('output_add'),
                   model_uri=MATERIALSML.output_add, domain=None, range=Optional[Union[bool, Bool]])

slots.pll_setup_settings = Slot(uri=MATERIALSML.pll_setup_settings, name="pll_setup_settings", curie=MATERIALSML.curie('pll_setup_settings'),
                   model_uri=MATERIALSML.pll_setup_settings, domain=None, range=Optional[Union[dict, PLLSetupSettings]])

slots.q_factor = Slot(uri=EMMO.EMMO_0658e7df_ffd9_4779_82fc_62efe0a7f3b1, name="q_factor", curie=EMMO.curie('EMMO_0658e7df_ffd9_4779_82fc_62efe0a7f3b1'),
                   model_uri=MATERIALSML.q_factor, domain=None, range=Optional[float])

slots.demodulation_bw_amp = Slot(uri=MATERIALSML.demodulation_bw_amp, name="demodulation_bw_amp", curie=MATERIALSML.curie('demodulation_bw_amp'),
                   model_uri=MATERIALSML.demodulation_bw_amp, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.demodulation_bw_pha = Slot(uri=MATERIALSML.demodulation_bw_pha, name="demodulation_bw_pha", curie=MATERIALSML.curie('demodulation_bw_pha'),
                   model_uri=MATERIALSML.demodulation_bw_pha, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.amplitude_excitation_ratio = Slot(uri=MATERIALSML.amplitude_excitation_ratio, name="amplitude_excitation_ratio", curie=MATERIALSML.curie('amplitude_excitation_ratio'),
                   model_uri=MATERIALSML.amplitude_excitation_ratio, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.amplitude_range = Slot(uri=MATERIALSML.amplitude_range, name="amplitude_range", curie=MATERIALSML.curie('amplitude_range'),
                   model_uri=MATERIALSML.amplitude_range, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.dwell_time = Slot(uri=MATERIALSML.dwell_time, name="dwell_time", curie=MATERIALSML.curie('dwell_time'),
                   model_uri=MATERIALSML.dwell_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.sample_temperature = Slot(uri=MATERIALSML.sample_temperature, name="sample_temperature", curie=MATERIALSML.curie('sample_temperature'),
                   model_uri=MATERIALSML.sample_temperature, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.recording_temperature = Slot(uri=MATERIALSML.recording_temperature, name="recording_temperature", curie=MATERIALSML.curie('recording_temperature'),
                   model_uri=MATERIALSML.recording_temperature, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.scan_dx = Slot(uri=MATERIALSML.scan_dx, name="scan_dx", curie=MATERIALSML.curie('scan_dx'),
                   model_uri=MATERIALSML.scan_dx, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.scan_z_min = Slot(uri=MATERIALSML.scan_z_min, name="scan_z_min", curie=MATERIALSML.curie('scan_z_min'),
                   model_uri=MATERIALSML.scan_z_min, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.scan_z_max = Slot(uri=MATERIALSML.scan_z_max, name="scan_z_max", curie=MATERIALSML.curie('scan_z_max'),
                   model_uri=MATERIALSML.scan_z_max, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.amplitude = Slot(uri=NCIT.NCIT_C70831, name="amplitude", curie=NCIT.curie('NCIT_C70831'),
                   model_uri=MATERIALSML.amplitude, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.resonance_frequency = Slot(uri=MATERIALSML.resonance_frequency, name="resonance_frequency", curie=MATERIALSML.curie('resonance_frequency'),
                   model_uri=MATERIALSML.resonance_frequency, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.p_model = Slot(uri=MATERIALSML.p_model, name="p_model", curie=MATERIALSML.curie('p_model'),
                   model_uri=MATERIALSML.p_model, domain=None, range=Optional[Union[str, "PModelEnum"]])

slots.wfms_uuid = Slot(uri=NCIT.NCIT_C54100, name="wfms_uuid", curie=NCIT.curie('NCIT_C54100'),
                   model_uri=MATERIALSML.wfms_uuid, domain=None, range=Optional[str])

slots.optimised = Slot(uri=MATERIALSML.optimised, name="optimised", curie=MATERIALSML.curie('optimised'),
                   model_uri=MATERIALSML.optimised, domain=None, range=Optional[Union[bool, Bool]])

slots.periodic_boundary_conditions = Slot(uri=MATERIALSML.periodic_boundary_conditions, name="periodic_boundary_conditions", curie=MATERIALSML.curie('periodic_boundary_conditions'),
                   model_uri=MATERIALSML.periodic_boundary_conditions, domain=None, range=Optional[Union[dict, PeriodicBoundaryConditions]])

slots.cell_vectors = Slot(uri=MATERIALSML.cell_vectors, name="cell_vectors", curie=MATERIALSML.curie('cell_vectors'),
                   model_uri=MATERIALSML.cell_vectors, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.atoms_positions = Slot(uri=MATERIALSML.atoms_positions, name="atoms_positions", curie=MATERIALSML.curie('atoms_positions'),
                   model_uri=MATERIALSML.atoms_positions, domain=None, range=Optional[Union[Union[dict, AtomsPositions], List[Union[dict, AtomsPositions]]]])

slots.atom_symbol = Slot(uri=EMMO.EMMO_eb77076b_a104_42ac_a065_798b2d2809ad, name="atom_symbol", curie=EMMO.curie('EMMO_eb77076b_a104_42ac_a065_798b2d2809ad'),
                   model_uri=MATERIALSML.atom_symbol, domain=None, range=Optional[str])

slots.k_points_conditions = Slot(uri=MATERIALSML.k_points_conditions, name="k_points_conditions", curie=MATERIALSML.curie('k_points_conditions'),
                   model_uri=MATERIALSML.k_points_conditions, domain=None, range=Optional[Union[dict, KPointsConditions]])

slots.bs_energies = Slot(uri=MATERIALSML.bs_energies, name="bs_energies", curie=MATERIALSML.curie('bs_energies'),
                   model_uri=MATERIALSML.bs_energies, domain=None, range=Optional[Union[dict, BSEnergy]])

slots.charge = Slot(uri=MATERIALSML.charge, name="charge", curie=MATERIALSML.curie('charge'),
                   model_uri=MATERIALSML.charge, domain=None, range=Optional[float])

slots.magnetic_momentum = Slot(uri=MATERIALSML.magnetic_momentum, name="magnetic_momentum", curie=MATERIALSML.curie('magnetic_momentum'),
                   model_uri=MATERIALSML.magnetic_momentum, domain=None, range=Optional[float])

slots.wfms_url = Slot(uri=MATERIALSML.wfms_url, name="wfms_url", curie=MATERIALSML.curie('wfms_url'),
                   model_uri=MATERIALSML.wfms_url, domain=None, range=Optional[str])

slots.state = Slot(uri=MATERIALSML.state, name="state", curie=MATERIALSML.curie('state'),
                   model_uri=MATERIALSML.state, domain=None, range=Optional[str])

slots.material = Slot(uri=SCHEMA.material, name="material", curie=SCHEMA.curie('material'),
                   model_uri=MATERIALSML.material, domain=None, range=Optional[str])

slots.face = Slot(uri=NCIT.NCIT_C126364, name="face", curie=NCIT.curie('NCIT_C126364'),
                   model_uri=MATERIALSML.face, domain=None, range=Optional[str])

slots.sample_plate = Slot(uri=EMMO_NANO['EMMO_69c892cc-b9c4-53cd-8e5e-216c57e2592f'], name="sample_plate", curie=EMMO_NANO.curie('EMMO_69c892cc-b9c4-53cd-8e5e-216c57e2592f'),
                   model_uri=MATERIALSML.sample_plate, domain=None, range=Optional[str])

slots.diameter = Slot(uri=NCIT.NCIT_C25285, name="diameter", curie=NCIT.curie('NCIT_C25285'),
                   model_uri=MATERIALSML.diameter, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.height = Slot(uri=NCIT.NCIT_C25347, name="height", curie=NCIT.curie('NCIT_C25347'),
                   model_uri=MATERIALSML.height, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.specifications = Slot(uri=SIO.SIO_000090, name="specifications", curie=SIO.curie('SIO_000090'),
                   model_uri=MATERIALSML.specifications, domain=None, range=Optional[str])

slots.reference_number = Slot(uri=NCIT.NCIT_99287, name="reference_number", curie=NCIT.curie('NCIT_99287'),
                   model_uri=MATERIALSML.reference_number, domain=None, range=Optional[str])

slots.stabilisation_time = Slot(uri=MATERIALSML.stabilisation_time, name="stabilisation_time", curie=MATERIALSML.curie('stabilisation_time'),
                   model_uri=MATERIALSML.stabilisation_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.deposition_time = Slot(uri=MATERIALSML.deposition_time, name="deposition_time", curie=MATERIALSML.curie('deposition_time'),
                   model_uri=MATERIALSML.deposition_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.substrate_temperature = Slot(uri=MATERIALSML.substrate_temperature, name="substrate_temperature", curie=MATERIALSML.curie('substrate_temperature'),
                   model_uri=MATERIALSML.substrate_temperature, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.molecule_temperature = Slot(uri=MATERIALSML.molecule_temperature, name="molecule_temperature", curie=MATERIALSML.curie('molecule_temperature'),
                   model_uri=MATERIALSML.molecule_temperature, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.evaporator_slot = Slot(uri=MATERIALSML.evaporator_slot, name="evaporator_slot", curie=MATERIALSML.curie('evaporator_slot'),
                   model_uri=MATERIALSML.evaporator_slot, domain=None, range=Optional[Union[dict, EvaporatorSlot]])

slots.evaporator_number = Slot(uri=MATERIALSML.evaporator_number, name="evaporator_number", curie=MATERIALSML.curie('evaporator_number'),
                   model_uri=MATERIALSML.evaporator_number, domain=None, range=Optional[int])

slots.details = Slot(uri=MATERIALSML.details, name="details", curie=MATERIALSML.curie('details'),
                   model_uri=MATERIALSML.details, domain=None, range=Optional[str])

slots.fermi_energy = Slot(uri=MATERIALSML.fermi_energy, name="fermi_energy", curie=MATERIALSML.curie('fermi_energy'),
                   model_uri=MATERIALSML.fermi_energy, domain=None, range=Optional[float])

slots.scf_convergence_threshold = Slot(uri=MATERIALSML.scf_convergence_threshold, name="scf_convergence_threshold", curie=MATERIALSML.curie('scf_convergence_threshold'),
                   model_uri=MATERIALSML.scf_convergence_threshold, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.xc_functional = Slot(uri=MATERIALSML.xc_functional, name="xc_functional", curie=MATERIALSML.curie('xc_functional'),
                   model_uri=MATERIALSML.xc_functional, domain=None, range=Optional[Union[str, "XCFunctionalEnum"]])

slots.spin_polarised_calculation = Slot(uri=MATERIALSML.spin_polarised_calculation, name="spin_polarised_calculation", curie=MATERIALSML.curie('spin_polarised_calculation'),
                   model_uri=MATERIALSML.spin_polarised_calculation, domain=None, range=Optional[Union[bool, Bool]])

slots.spin_multiplicity = Slot(uri=MATERIALSML.spin_multiplicity, name="spin_multiplicity", curie=MATERIALSML.curie('spin_multiplicity'),
                   model_uri=MATERIALSML.spin_multiplicity, domain=None, range=Optional[Union[bool, Bool]])

slots.initial_spin_guesses = Slot(uri=MATERIALSML.initial_spin_guesses, name="initial_spin_guesses", curie=MATERIALSML.curie('initial_spin_guesses'),
                   model_uri=MATERIALSML.initial_spin_guesses, domain=None, range=Optional[Union[float, List[float]]])

slots.absolute_magnetisation = Slot(uri=MATERIALSML.absolute_magnetisation, name="absolute_magnetisation", curie=MATERIALSML.curie('absolute_magnetisation'),
                   model_uri=MATERIALSML.absolute_magnetisation, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.total_magnetisation = Slot(uri=MATERIALSML.total_magnetisation, name="total_magnetisation", curie=MATERIALSML.curie('total_magnetisation'),
                   model_uri=MATERIALSML.total_magnetisation, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.total_spin_squared = Slot(uri=MATERIALSML.total_spin_squared, name="total_spin_squared", curie=MATERIALSML.curie('total_spin_squared'),
                   model_uri=MATERIALSML.total_spin_squared, domain=None, range=Optional[float])

slots.net_charge = Slot(uri=MATERIALSML.net_charge, name="net_charge", curie=MATERIALSML.curie('net_charge'),
                   model_uri=MATERIALSML.net_charge, domain=None, range=Optional[float])

slots.smearing = Slot(uri=MATERIALSML.smearing, name="smearing", curie=MATERIALSML.curie('smearing'),
                   model_uri=MATERIALSML.smearing, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.fermi_dirac_temperature = Slot(uri=MATERIALSML.fermi_dirac_temperature, name="fermi_dirac_temperature", curie=MATERIALSML.curie('fermi_dirac_temperature'),
                   model_uri=MATERIALSML.fermi_dirac_temperature, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.force_multiplicity = Slot(uri=MATERIALSML.force_multiplicity, name="force_multiplicity", curie=MATERIALSML.curie('force_multiplicity'),
                   model_uri=MATERIALSML.force_multiplicity, domain=None, range=Optional[Union[bool, Bool]])

slots.output_mulliken_population_analysis = Slot(uri=MATERIALSML.output_mulliken_population_analysis, name="output_mulliken_population_analysis", curie=MATERIALSML.curie('output_mulliken_population_analysis'),
                   model_uri=MATERIALSML.output_mulliken_population_analysis, domain=None, range=Optional[Union[dict, OutputPopulationAnalysis]])

slots.output_hirshfeld_population_analysis = Slot(uri=MATERIALSML.output_hirshfeld_population_analysis, name="output_hirshfeld_population_analysis", curie=MATERIALSML.curie('output_hirshfeld_population_analysis'),
                   model_uri=MATERIALSML.output_hirshfeld_population_analysis, domain=None, range=Optional[Union[dict, OutputPopulationAnalysis]])

slots.draft_type = Slot(uri=MATERIALSML.draft_type, name="draft_type", curie=MATERIALSML.curie('draft_type'),
                   model_uri=MATERIALSML.draft_type, domain=None, range=Optional[Union[str, "DraftTypeEnum"]])

slots.weight_before = Slot(uri=MATERIALSML.weight_before, name="weight_before", curie=MATERIALSML.curie('weight_before'),
                   model_uri=MATERIALSML.weight_before, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.weight_after = Slot(uri=MATERIALSML.weight_after, name="weight_after", curie=MATERIALSML.curie('weight_after'),
                   model_uri=MATERIALSML.weight_after, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.substance = Slot(uri=SCHEMA.Substance, name="substance", curie=SCHEMA.curie('Substance'),
                   model_uri=MATERIALSML.substance, domain=None, range=Optional[str])

slots.serial_number = Slot(uri=NCIT.NCIT_C73518, name="serial_number", curie=NCIT.curie('NCIT_C73518'),
                   model_uri=MATERIALSML.serial_number, domain=None, range=Optional[str])

slots.empa_id = Slot(uri=MATERIALSML.empa_id, name="empa_id", curie=MATERIALSML.curie('empa_id'),
                   model_uri=MATERIALSML.empa_id, domain=None, range=Optional[str])

slots.documentation_website = Slot(uri=NCIT.NCIT_C93821, name="documentation_website", curie=NCIT.curie('NCIT_C93821'),
                   model_uri=MATERIALSML.documentation_website, domain=None, range=Optional[str],
                   pattern=re.compile(r'(https?:\/\/(?:www\\.|(?!www))[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|www\\.[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|https?:\/\/(?:www\\.|(?!www))[a-zA-Z0-9]+\\.[^\\s]{2,}|www\\.[a-zA-Z0-9]+\\.[^\\s]{2,})'))

slots.dewar = Slot(uri=NCIT.NCIT_C43190, name="dewar", curie=NCIT.curie('NCIT_C43190'),
                   model_uri=MATERIALSML.dewar, domain=None, range=Optional[Union[str, "DewarEnum"]])

slots.force_convergence_threshold = Slot(uri=MATERIALSML.force_convergence_threshold, name="force_convergence_threshold", curie=MATERIALSML.curie('force_convergence_threshold'),
                   model_uri=MATERIALSML.force_convergence_threshold, domain=None, range=Optional[Union[dict, ForceConvergenceThreshold]])

slots.geometry_constraints = Slot(uri=MATERIALSML.geometry_constraints, name="geometry_constraints", curie=MATERIALSML.curie('geometry_constraints'),
                   model_uri=MATERIALSML.geometry_constraints, domain=None, range=Optional[Union[float, List[float]]])

slots.funder = Slot(uri=MATERIALSML.funder, name="funder", curie=MATERIALSML.curie('funder'),
                   model_uri=MATERIALSML.funder, domain=None, range=Optional[str])

slots.start_date = Slot(uri=SCHEMA.startDate, name="start_date", curie=SCHEMA.curie('startDate'),
                   model_uri=MATERIALSML.start_date, domain=None, range=Optional[Union[str, XSDDate]])

slots.end_date = Slot(uri=SCHEMA.endDate, name="end_date", curie=SCHEMA.curie('endDate'),
                   model_uri=MATERIALSML.end_date, domain=None, range=Optional[Union[str, XSDDate]])

slots.budget = Slot(uri=NCIT.NCIT_C60757, name="budget", curie=NCIT.curie('NCIT_C60757'),
                   model_uri=MATERIALSML.budget, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.project_id = Slot(uri=SCHEMA.identifier, name="project_id", curie=SCHEMA.curie('identifier'),
                   model_uri=MATERIALSML.project_id, domain=None, range=Optional[str])

slots.acknowledgement_sentence = Slot(uri=IAO.IAO_0000324, name="acknowledgement_sentence", curie=IAO.curie('IAO_0000324'),
                   model_uri=MATERIALSML.acknowledgement_sentence, domain=None, range=Optional[str])

slots.layers_2d = Slot(uri=MATERIALSML.layers_2d, name="layers_2d", curie=MATERIALSML.curie('layers_2d'),
                   model_uri=MATERIALSML.layers_2d, domain=None, range=Optional[Union[Union[dict, Layers2DDetails], List[Union[dict, Layers2DDetails]]]])

slots.number_layers = Slot(uri=MATERIALSML.number_layers, name="number_layers", curie=MATERIALSML.curie('number_layers'),
                   model_uri=MATERIALSML.number_layers, domain=None, range=Optional[int])

slots.dopants = Slot(uri=MATERIALSML.dopants, name="dopants", curie=MATERIALSML.curie('dopants'),
                   model_uri=MATERIALSML.dopants, domain=None, range=Optional[str])

slots.coverage = Slot(uri=MATERIALSML.coverage, name="coverage", curie=MATERIALSML.curie('coverage'),
                   model_uri=MATERIALSML.coverage, domain=None, range=Optional[str])

slots.substrates = Slot(uri=MATERIALSML.substrates, name="substrates", curie=MATERIALSML.curie('substrates'),
                   model_uri=MATERIALSML.substrates, domain=None, range=Optional[Union[Union[dict, SubstratesDetails], List[Union[dict, SubstratesDetails]]]])

slots.width = Slot(uri=SCHEMA.width, name="width", curie=SCHEMA.curie('width'),
                   model_uri=MATERIALSML.width, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.length = Slot(uri=NCIT.NCIT_C25334, name="length", curie=NCIT.curie('NCIT_C25334'),
                   model_uri=MATERIALSML.length, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.thickness = Slot(uri=NCIT.NCIT_C41145, name="thickness", curie=NCIT.curie('NCIT_C41145'),
                   model_uri=MATERIALSML.thickness, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.supplier_own_name = Slot(uri=SCHEMA.alternateName, name="supplier_own_name", curie=SCHEMA.curie('alternateName'),
                   model_uri=MATERIALSML.supplier_own_name, domain=None, range=Optional[str])

slots.geometry_constraints_increments = Slot(uri=MATERIALSML.geometry_constraints_increments, name="geometry_constraints_increments", curie=MATERIALSML.curie('geometry_constraints_increments'),
                   model_uri=MATERIALSML.geometry_constraints_increments, domain=None, range=Optional[Union[float, List[float]]])

slots.number_geometries = Slot(uri=MATERIALSML.number_geometries, name="number_geometries", curie=MATERIALSML.curie('number_geometries'),
                   model_uri=MATERIALSML.number_geometries, domain=None, range=Optional[int])

slots.method_type = Slot(uri=MATERIALSML.method_type, name="method_type", curie=MATERIALSML.curie('method_type'),
                   model_uri=MATERIALSML.method_type, domain=None, range=Optional[Union[str, "MethodTypeEnum"]])

slots.energy_barrier = Slot(uri=MATERIALSML.energy_barrier, name="energy_barrier", curie=MATERIALSML.curie('energy_barrier'),
                   model_uri=MATERIALSML.energy_barrier, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.energies = Slot(uri=MATERIALSML.energies, name="energies", curie=MATERIALSML.curie('energies'),
                   model_uri=MATERIALSML.energies, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.chemist = Slot(uri=OCCO.OCCO_19203100, name="chemist", curie=OCCO.curie('OCCO_19203100'),
                   model_uri=MATERIALSML.chemist, domain=None, range=Optional[Union[dict, Chemist]])

slots.atomic_selections = Slot(uri=CHEBI.CHEBI_33250, name="atomic_selections", curie=CHEBI.curie('CHEBI_33250'),
                   model_uri=MATERIALSML.atomic_selections, domain=None, range=Optional[Union[float, List[float]]])

slots.orbitals_selections = Slot(uri=MATERIALSML.orbitals_selections, name="orbitals_selections", curie=MATERIALSML.curie('orbitals_selections'),
                   model_uri=MATERIALSML.orbitals_selections, domain=None, range=Optional[Union[float, List[float]]])

slots.degauss_energy = Slot(uri=MATERIALSML.degauss_energy, name="degauss_energy", curie=MATERIALSML.curie('degauss_energy'),
                   model_uri=MATERIALSML.degauss_energy, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.plot_group_by = Slot(uri=MATERIALSML.plot_group_by, name="plot_group_by", curie=MATERIALSML.curie('plot_group_by'),
                   model_uri=MATERIALSML.plot_group_by, domain=None, range=Optional[Union[str, "PlotGroupByEnum"]])

slots.plot_contributions = Slot(uri=MATERIALSML.plot_contributions, name="plot_contributions", curie=MATERIALSML.curie('plot_contributions'),
                   model_uri=MATERIALSML.plot_contributions, domain=None, range=Optional[Union[str, "PlotContributionsEnum"]])

slots.energy = Slot(uri=EMMO.EMMO_31ec09ba_1713_42cb_83c7_b38bf6f9ced2, name="energy", curie=EMMO.curie('EMMO_31ec09ba_1713_42cb_83c7_b38bf6f9ced2'),
                   model_uri=MATERIALSML.energy, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.short_name = Slot(uri=FOAF.accountName, name="short_name", curie=FOAF.curie('accountName'),
                   model_uri=MATERIALSML.short_name, domain=None, range=Optional[str])

slots.work_status = Slot(uri=OPMI.OPMI_0000122, name="work_status", curie=OPMI.curie('OPMI_0000122'),
                   model_uri=MATERIALSML.work_status, domain=None, range=Optional[Union[bool, Bool]])

slots.mobile_phone = Slot(uri=NCIT.NCIT_C98193, name="mobile_phone", curie=NCIT.curie('NCIT_C98193'),
                   model_uri=MATERIALSML.mobile_phone, domain=None, range=Optional[str])

slots.abstract = Slot(uri=SCHEMA.abstract, name="abstract", curie=SCHEMA.curie('abstract'),
                   model_uri=MATERIALSML.abstract, domain=None, range=Optional[str])

slots.doi = Slot(uri=NCIT.NCIT_C71462, name="doi", curie=NCIT.curie('NCIT_C71462'),
                   model_uri=MATERIALSML.doi, domain=None, range=Optional[str])

slots.year = Slot(uri=NCIT.NCIT_C29848, name="year", curie=NCIT.curie('NCIT_C29848'),
                   model_uri=MATERIALSML.year, domain=None, range=Optional[int])

slots.url = Slot(uri=NCIT.NCIT_C42743, name="url", curie=NCIT.curie('NCIT_C42743'),
                   model_uri=MATERIALSML.url, domain=None, range=Optional[str],
                   pattern=re.compile(r'(https?:\/\/(?:www\\.|(?!www))[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|www\\.[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|https?:\/\/(?:www\\.|(?!www))[a-zA-Z0-9]+\\.[^\\s]{2,}|www\\.[a-zA-Z0-9]+\\.[^\\s]{2,})'))

slots.dataset_url = Slot(uri=DCAT.downloadURL, name="dataset_url", curie=DCAT.curie('downloadURL'),
                   model_uri=MATERIALSML.dataset_url, domain=None, range=Optional[str],
                   pattern=re.compile(r'(https?:\/\/(?:www\\.|(?!www))[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|www\\.[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|https?:\/\/(?:www\\.|(?!www))[a-zA-Z0-9]+\\.[^\\s]{2,}|www\\.[a-zA-Z0-9]+\\.[^\\s]{2,})'))

slots.repository_url = Slot(uri=APOLLO.APOLLO_SV_00000488, name="repository_url", curie=APOLLO.curie('APOLLO_SV_00000488'),
                   model_uri=MATERIALSML.repository_url, domain=None, range=Optional[str],
                   pattern=re.compile(r'(https?:\/\/(?:www\\.|(?!www))[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|www\\.[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|https?:\/\/(?:www\\.|(?!www))[a-zA-Z0-9]+\\.[^\\s]{2,}|www\\.[a-zA-Z0-9]+\\.[^\\s]{2,})'))

slots.discharge_voltage = Slot(uri=EMMO_ELECTROCHEMISTRY.electrochemistry_c7b26177_21bf_4787_b656_8e78edf27f88, name="discharge_voltage", curie=EMMO_ELECTROCHEMISTRY.curie('electrochemistry_c7b26177_21bf_4787_b656_8e78edf27f88'),
                   model_uri=MATERIALSML.discharge_voltage, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.e_min = Slot(uri=MATERIALSML.e_min, name="e_min", curie=MATERIALSML.curie('e_min'),
                   model_uri=MATERIALSML.e_min, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.e_max = Slot(uri=MATERIALSML.e_max, name="e_max", curie=MATERIALSML.curie('e_max'),
                   model_uri=MATERIALSML.e_max, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.de = Slot(uri=MATERIALSML.de, name="de", curie=MATERIALSML.curie('de'),
                   model_uri=MATERIALSML.de, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.fwhm = Slot(uri=NCIT.NCIT_C94903, name="fwhm", curie=NCIT.curie('NCIT_C94903'),
                   model_uri=MATERIALSML.fwhm, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.extrap_plane = Slot(uri=MATERIALSML.extrap_plane, name="extrap_plane", curie=MATERIALSML.curie('extrap_plane'),
                   model_uri=MATERIALSML.extrap_plane, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.constant_current = Slot(uri=MATERIALSML.constant_current, name="constant_current", curie=MATERIALSML.curie('constant_current'),
                   model_uri=MATERIALSML.constant_current, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.constant_height = Slot(uri=MATERIALSML.constant_height, name="constant_height", curie=MATERIALSML.curie('constant_height'),
                   model_uri=MATERIALSML.constant_height, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.experiment_type = Slot(uri=FAIRMAT_STS['nxsts-entry-experiment-type-field'], name="experiment_type", curie=FAIRMAT_STS.curie('nxsts-entry-experiment-type-field'),
                   model_uri=MATERIALSML.experiment_type, domain=None, range=Optional[Union[str, "ExperimentTypeEnum"]])

slots.acquisition_coordinate_x = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-position-x-field'], name="acquisition_coordinate_x", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-position-x-field'),
                   model_uri=MATERIALSML.acquisition_coordinate_x, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.acquisition_coordinate_y = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-position-y-field'], name="acquisition_coordinate_y", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-position-y-field'),
                   model_uri=MATERIALSML.acquisition_coordinate_y, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.acquisition_coordinate_z = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-position-z-field'], name="acquisition_coordinate_z", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-position-z-field'),
                   model_uri=MATERIALSML.acquisition_coordinate_z, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.final_z = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-final-z-field'], name="final_z", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-final-z-field'),
                   model_uri=MATERIALSML.final_z, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.filter_type = Slot(uri=OME_ONTO.FilterType, name="filter_type", curie=OME_ONTO.curie('FilterType'),
                   model_uri=MATERIALSML.filter_type, domain=None, range=Optional[Union[str, "FilterTypeEnum"]])

slots.filter_order = Slot(uri=MATERIALSML.filter_order, name="filter_order", curie=MATERIALSML.curie('filter_order'),
                   model_uri=MATERIALSML.filter_order, domain=None, range=Optional[int])

slots.filter_cutoff = Slot(uri=NCIT.NCIT_C94895, name="filter_cutoff", curie=NCIT.curie('NCIT_C94895'),
                   model_uri=MATERIALSML.filter_cutoff, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.num_pixel = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-num-pixel-field'], name="num_pixel", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-num-pixel-field'),
                   model_uri=MATERIALSML.num_pixel, domain=None, range=Optional[int])

slots.z_avg_time = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-z-avg-time-field'], name="z_avg_time", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-z-avg-time-field'),
                   model_uri=MATERIALSML.z_avg_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.first_settling_time = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-first-settling-time-field'], name="first_settling_time", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-first-settling-time-field'),
                   model_uri=MATERIALSML.first_settling_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.settling_time = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-settling-time-field'], name="settling_time", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-settling-time-field'),
                   model_uri=MATERIALSML.settling_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.integration_time = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-integration-time-field'], name="integration_time", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-integration-time-field'),
                   model_uri=MATERIALSML.integration_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.end_settling_time = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-end-settling-time-field'], name="end_settling_time", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-end-settling-time-field'),
                   model_uri=MATERIALSML.end_settling_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_control_time = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-control-time-field'], name="z_control_time", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-control-time-field'),
                   model_uri=MATERIALSML.z_control_time, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.max_slew_rate = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-max-slew-rate-field'], name="max_slew_rate", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-max-slew-rate-field'),
                   model_uri=MATERIALSML.max_slew_rate, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.backward_sweep = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-backward-sweep-field'], name="backward_sweep", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-backward-sweep-field'),
                   model_uri=MATERIALSML.backward_sweep, domain=None, range=Optional[Union[bool, Bool]])

slots.num_sweeps = Slot(uri=FAIRMAT_STS['nxsts-entry-resolution-indicators-number-of-sweeps-field'], name="num_sweeps", curie=FAIRMAT_STS.curie('nxsts-entry-resolution-indicators-number-of-sweeps-field'),
                   model_uri=MATERIALSML.num_sweeps, domain=None, range=Optional[int])

slots.channel_names = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-scan-control-circuit-channels-current-field'], name="channel_names", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-scan-control-circuit-channels-current-field'),
                   model_uri=MATERIALSML.channel_names, domain=None, range=Optional[Union[str, List[str]]])

slots.record_final_z = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-final-z-field'], name="record_final_z", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-final-z-field'),
                   model_uri=MATERIALSML.record_final_z, domain=None, range=Optional[Union[bool, Bool]])

slots.bias_spectroscopy_reset_bias = Slot(uri=MATERIALSML.bias_spectroscopy_reset_bias, name="bias_spectroscopy_reset_bias", curie=MATERIALSML.curie('bias_spectroscopy_reset_bias'),
                   model_uri=MATERIALSML.bias_spectroscopy_reset_bias, domain=None, range=Optional[Union[bool, Bool]])

slots.bias_spectroscopy_sweep_start = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-sweep-start-field'], name="bias_spectroscopy_sweep_start", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-sweep-start-field'),
                   model_uri=MATERIALSML.bias_spectroscopy_sweep_start, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.bias_spectroscopy_sweep_end = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-sweep-end-field'], name="bias_spectroscopy_sweep_end", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-environment-sweep-control-bias-spectroscopy-sweep-end-field'),
                   model_uri=MATERIALSML.bias_spectroscopy_sweep_end, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.bias_spectroscopy_z_offset = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-offset-field'], name="bias_spectroscopy_z_offset", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-offset-field'),
                   model_uri=MATERIALSML.bias_spectroscopy_z_offset, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.bias_spectroscopy_z_controller_hold = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-control-hold-field'], name="bias_spectroscopy_z_controller_hold", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-control-hold-field'),
                   model_uri=MATERIALSML.bias_spectroscopy_z_controller_hold, domain=None, range=Optional[Union[bool, Bool]])

slots.z_spectroscopy_reset_z = Slot(uri=MATERIALSML.z_spectroscopy_reset_z, name="z_spectroscopy_reset_z", curie=MATERIALSML.curie('z_spectroscopy_reset_z'),
                   model_uri=MATERIALSML.z_spectroscopy_reset_z, domain=None, range=Optional[Union[bool, Bool]])

slots.z_spectroscopy_initial_z_offset = Slot(uri=MATERIALSML.z_spectroscopy_initial_z_offset, name="z_spectroscopy_initial_z_offset", curie=MATERIALSML.curie('z_spectroscopy_initial_z_offset'),
                   model_uri=MATERIALSML.z_spectroscopy_initial_z_offset, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_spectroscopy_sweep_distance = Slot(uri=MATERIALSML.z_spectroscopy_sweep_distance, name="z_spectroscopy_sweep_distance", curie=MATERIALSML.curie('z_spectroscopy_sweep_distance'),
                   model_uri=MATERIALSML.z_spectroscopy_sweep_distance, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_spectroscopy_time_between_forward_backward = Slot(uri=MATERIALSML.z_spectroscopy_time_between_forward_backward, name="z_spectroscopy_time_between_forward_backward", curie=MATERIALSML.curie('z_spectroscopy_time_between_forward_backward'),
                   model_uri=MATERIALSML.z_spectroscopy_time_between_forward_backward, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_controller_setpoint = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-controller-setpoint-field'], name="z_controller_setpoint", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-controller-setpoint-field'),
                   model_uri=MATERIALSML.z_controller_setpoint, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_controller_p_gain = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-control-p-gain-field'], name="z_controller_p_gain", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-control-p-gain-field'),
                   model_uri=MATERIALSML.z_controller_p_gain, domain=None, range=Optional[float])

slots.z_controller_i_gain = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-control-i-gain-field'], name="z_controller_i_gain", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-control-i-gain-field'),
                   model_uri=MATERIALSML.z_controller_i_gain, domain=None, range=Optional[float])

slots.z_controller_time_constant = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-control-time-const-field'], name="z_controller_time_constant", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-control-time-const-field'),
                   model_uri=MATERIALSML.z_controller_time_constant, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_controller_tip_lift = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-control-tip-lift-field'], name="z_controller_tip_lift", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-control-tip-lift-field'),
                   model_uri=MATERIALSML.z_controller_tip_lift, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.z_controller_switch_off_delay = Slot(uri=FAIRMAT_STS['nxsts-entry-reproducibility-indicators-z-control-switchoff-delay-field'], name="z_controller_switch_off_delay", curie=FAIRMAT_STS.curie('nxsts-entry-reproducibility-indicators-z-control-switchoff-delay-field'),
                   model_uri=MATERIALSML.z_controller_switch_off_delay, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.lock_in_status = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-lock-in-status-field'], name="lock_in_status", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-lock-in-status-field'),
                   model_uri=MATERIALSML.lock_in_status, domain=None, range=Optional[Union[bool, Bool]])

slots.lock_in_modulated_signal = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-lock-in-modulation-signal-field'], name="lock_in_modulated_signal", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-lock-in-modulation-signal-field'),
                   model_uri=MATERIALSML.lock_in_modulated_signal, domain=None, range=Optional[str])

slots.lock_in_frequency = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-lock-in-modulation-frequency-field'], name="lock_in_frequency", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-lock-in-modulation-frequency-field'),
                   model_uri=MATERIALSML.lock_in_frequency, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.lock_in_amplitude = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-lock-in-modulation-amplitude-field'], name="lock_in_amplitude", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-lock-in-modulation-amplitude-field'),
                   model_uri=MATERIALSML.lock_in_amplitude, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.lock_in_hp_filter_cutoffs = Slot(uri=MATERIALSML.lock_in_hp_filter_cutoffs, name="lock_in_hp_filter_cutoffs", curie=MATERIALSML.curie('lock_in_hp_filter_cutoffs'),
                   model_uri=MATERIALSML.lock_in_hp_filter_cutoffs, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.lock_in_hp_filter_orders = Slot(uri=MATERIALSML.lock_in_hp_filter_orders, name="lock_in_hp_filter_orders", curie=MATERIALSML.curie('lock_in_hp_filter_orders'),
                   model_uri=MATERIALSML.lock_in_hp_filter_orders, domain=None, range=Optional[Union[int, List[int]]])

slots.lock_in_demodulated_signal = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-lock-in-demodulated-signal-field'], name="lock_in_demodulated_signal", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-lock-in-demodulated-signal-field'),
                   model_uri=MATERIALSML.lock_in_demodulated_signal, domain=None, range=Optional[Union[str, "LockInDemodulatedSignalEnum"]])

slots.lock_in_harmonics = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-lock-in-harmonic-order-n-field'], name="lock_in_harmonics", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-lock-in-harmonic-order-n-field'),
                   model_uri=MATERIALSML.lock_in_harmonics, domain=None, range=Optional[Union[int, List[int]]])

slots.lock_in_reference_phases = Slot(uri=FAIRMAT_STS['nxsts-entry-instrument-lock-in-ref-phase-n-field'], name="lock_in_reference_phases", curie=FAIRMAT_STS.curie('nxsts-entry-instrument-lock-in-ref-phase-n-field'),
                   model_uri=MATERIALSML.lock_in_reference_phases, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.lock_in_lp_filter_cutoffs = Slot(uri=MATERIALSML.lock_in_lp_filter_cutoffs, name="lock_in_lp_filter_cutoffs", curie=MATERIALSML.curie('lock_in_lp_filter_cutoffs'),
                   model_uri=MATERIALSML.lock_in_lp_filter_cutoffs, domain=None, range=Optional[Union[Union[dict, QuantityValue], List[Union[dict, QuantityValue]]]])

slots.lock_in_lp_filter_orders = Slot(uri=MATERIALSML.lock_in_lp_filter_orders, name="lock_in_lp_filter_orders", curie=MATERIALSML.curie('lock_in_lp_filter_orders'),
                   model_uri=MATERIALSML.lock_in_lp_filter_orders, domain=None, range=Optional[Union[int, List[int]]])

slots.lock_in_sync_filter = Slot(uri=MATERIALSML.lock_in_sync_filter, name="lock_in_sync_filter", curie=MATERIALSML.curie('lock_in_sync_filter'),
                   model_uri=MATERIALSML.lock_in_sync_filter, domain=None, range=Optional[Union[bool, Bool]])

slots.vacuum_system_status = Slot(uri=MATERIALSML.vacuum_system_status, name="vacuum_system_status", curie=MATERIALSML.curie('vacuum_system_status'),
                   model_uri=MATERIALSML.vacuum_system_status, domain=None, range=Optional[Union[dict, VacuumSystemStatusSettings]])

slots.principal_chamber = Slot(uri=MATERIALSML.principal_chamber, name="principal_chamber", curie=MATERIALSML.curie('principal_chamber'),
                   model_uri=MATERIALSML.principal_chamber, domain=None, range=Optional[Union[dict, PrincipalChamberSettings]])

slots.analysis_chamber = Slot(uri=MATERIALSML.analysis_chamber, name="analysis_chamber", curie=MATERIALSML.curie('analysis_chamber'),
                   model_uri=MATERIALSML.analysis_chamber, domain=None, range=Optional[Union[dict, AnalysisChamberSettings]])

slots.fel_chamber = Slot(uri=MATERIALSML.fel_chamber, name="fel_chamber", curie=MATERIALSML.curie('fel_chamber'),
                   model_uri=MATERIALSML.fel_chamber, domain=None, range=Optional[Union[dict, FELChamberSettings]])

slots.tsp_filament_used = Slot(uri=MATERIALSML.tsp_filament_used, name="tsp_filament_used", curie=MATERIALSML.curie('tsp_filament_used'),
                   model_uri=MATERIALSML.tsp_filament_used, domain=None, range=Optional[str])

slots.cryo_system_status = Slot(uri=MATERIALSML.cryo_system_status, name="cryo_system_status", curie=MATERIALSML.curie('cryo_system_status'),
                   model_uri=MATERIALSML.cryo_system_status, domain=None, range=Optional[Union[dict, CryoSystemStatusSettings]])

slots.liquid_helium = Slot(uri=MATERIALSML.liquid_helium, name="liquid_helium", curie=MATERIALSML.curie('liquid_helium'),
                   model_uri=MATERIALSML.liquid_helium, domain=None, range=Optional[Union[dict, QuantityValue]])

slots.liquid_nitrogen_refill_date = Slot(uri=MATERIALSML.liquid_nitrogen_refill_date, name="liquid_nitrogen_refill_date", curie=MATERIALSML.curie('liquid_nitrogen_refill_date'),
                   model_uri=MATERIALSML.liquid_nitrogen_refill_date, domain=None, range=Optional[str])

slots.tip_preparation_type = Slot(uri=MATERIALSML.tip_preparation_type, name="tip_preparation_type", curie=MATERIALSML.curie('tip_preparation_type'),
                   model_uri=MATERIALSML.tip_preparation_type, domain=None, range=Optional[Union[str, "TipPreparationTypeEnum"]])

slots.wire_material = Slot(uri=MATERIALSML.wire_material, name="wire_material", curie=MATERIALSML.curie('wire_material'),
                   model_uri=MATERIALSML.wire_material, domain=None, range=Optional[str])

slots.holder_id = Slot(uri=MATERIALSML.holder_id, name="holder_id", curie=MATERIALSML.curie('holder_id'),
                   model_uri=MATERIALSML.holder_id, domain=None, range=Optional[str])

slots.device_parameters = Slot(uri=NCIT.NCIT_80474, name="device_parameters", curie=NCIT.curie('NCIT_80474'),
                   model_uri=MATERIALSML.device_parameters, domain=None, range=Optional[str])

slots.model = Slot(uri=SCHEMA.model, name="model", curie=SCHEMA.curie('model'),
                   model_uri=MATERIALSML.model, domain=None, range=Optional[str])

slots.tip_sensor_type = Slot(uri=MATERIALSML.tip_sensor_type, name="tip_sensor_type", curie=MATERIALSML.curie('tip_sensor_type'),
                   model_uri=MATERIALSML.tip_sensor_type, domain=None, range=Optional[Union[str, "TipSensorTypeEnum"]])

slots.uhv_component_type = Slot(uri=MATERIALSML.uhv_component_type, name="uhv_component_type", curie=MATERIALSML.curie('uhv_component_type'),
                   model_uri=MATERIALSML.uhv_component_type, domain=None, range=Optional[Union[str, "UHVComponentTypeEnum"]])

slots.precursor_molecules = Slot(uri=SCHEMA.MolecularEntity, name="precursor_molecules", curie=SCHEMA.curie('MolecularEntity'),
                   model_uri=MATERIALSML.precursor_molecules, domain=None, range=Optional[Union[Union[dict, Molecule], List[Union[dict, Molecule]]]])

slots.supplier = Slot(uri=NCIT.NCIT_C43530, name="supplier", curie=NCIT.curie('NCIT_C43530'),
                   model_uri=MATERIALSML.supplier, domain=None, range=Optional[Union[dict, Supplier]])

slots.storage = Slot(uri=NCIT.NCIT_C16143, name="storage", curie=NCIT.curie('NCIT_C16143'),
                   model_uri=MATERIALSML.storage, domain=None, range=Optional[Union[dict, Storage]])

slots.room = Slot(uri=SCHEMA.Room, name="room", curie=SCHEMA.curie('Room'),
                   model_uri=MATERIALSML.room, domain=None, range=Optional[Union[dict, Room]])

slots.institution = Slot(uri=SCHEMA.Organization, name="institution", curie=SCHEMA.curie('Organization'),
                   model_uri=MATERIALSML.institution, domain=None, range=Optional[Union[dict, Institution]])

slots.annealing = Slot(uri=CHMO.CHMO_0001465, name="annealing", curie=CHMO.curie('CHMO_0001465'),
                   model_uri=MATERIALSML.annealing, domain=None, range=Optional[Union[dict, Annealing]])

slots.sputtering = Slot(uri=CHMO.CHMO_0001569, name="sputtering", curie=CHMO.curie('CHMO_0001569'),
                   model_uri=MATERIALSML.sputtering, domain=None, range=Optional[Union[dict, Sputtering]])

slots.deposition = Slot(uri=CHMO.CHMO_0001310, name="deposition", curie=CHMO.curie('CHMO_0001310'),
                   model_uri=MATERIALSML.deposition, domain=None, range=Optional[Union[dict, Deposition]])

slots.uhv_components = Slot(uri=MATERIALSML.uhv_components, name="uhv_components", curie=MATERIALSML.curie('uhv_components'),
                   model_uri=MATERIALSML.uhv_components, domain=None, range=Optional[Union[Union[dict, UHVComponent], List[Union[dict, UHVComponent]]]])

slots.instrument = Slot(uri=NCIT.NCIT_C16742, name="instrument", curie=NCIT.curie('NCIT_C16742'),
                   model_uri=MATERIALSML.instrument, domain=None, range=Optional[Union[dict, Instrument]])

slots.dft = Slot(uri=EMMO_DA['83a69cf2-00c1-58b8-915b-76cb3549890a'], name="dft", curie=EMMO_DA.curie('83a69cf2-00c1-58b8-915b-76cb3549890a'),
                   model_uri=MATERIALSML.dft, domain=None, range=Optional[Union[dict, DFT]])

slots.tb = Slot(uri=MATERIALSML.tb, name="tb", curie=MATERIALSML.curie('tb'),
                   model_uri=MATERIALSML.tb, domain=None, range=Optional[Union[dict, TightBinding]])

slots.crystal = Slot(uri=EMMO.EMMO_0bb3b434_73aa_428f_b4e8_2a2468648e19, name="crystal", curie=EMMO.curie('EMMO_0bb3b434_73aa_428f_b4e8_2a2468648e19'),
                   model_uri=MATERIALSML.crystal, domain=None, range=Optional[Union[dict, Crystal]])

slots.molecule = Slot(uri=SCHEMA.MolecularEntity, name="molecule", curie=SCHEMA.curie('MolecularEntity'),
                   model_uri=MATERIALSML.molecule, domain=None, range=Optional[Union[dict, Molecule]])

slots.person = Slot(uri=SCHEMA.Person, name="person", curie=SCHEMA.curie('Person'),
                   model_uri=MATERIALSML.person, domain=None, range=Optional[Union[dict, Person]])

slots.atomistic_models = Slot(uri=EMMO.EMMO_84cadc45_6758_46f2_ba2a_5ead65c70213, name="atomistic_models", curie=EMMO.curie('EMMO_84cadc45_6758_46f2_ba2a_5ead65c70213'),
                   model_uri=MATERIALSML.atomistic_models, domain=None, range=Optional[Union[Union[dict, AtomisticModel], List[Union[dict, AtomisticModel]]]])

slots.manufacturer = Slot(uri=NCIT.NCIT_C25392, name="manufacturer", curie=NCIT.curie('NCIT_C25392'),
                   model_uri=MATERIALSML.manufacturer, domain=None, range=Optional[Union[dict, Manufacturer]])

slots.tip_sensor = Slot(uri=EMMO_CHAMEO.Probe, name="tip_sensor", curie=EMMO_CHAMEO.curie('Probe'),
                   model_uri=MATERIALSML.tip_sensor, domain=None, range=Optional[Union[dict, TipSensor]])

slots.results = Slot(uri=EMMO.EMMO_0f6f0120_c079_4d95_bb11_4ddee05e530e, name="results", curie=EMMO.curie('EMMO_0f6f0120_c079_4d95_bb11_4ddee05e530e'),
                   model_uri=MATERIALSML.results, domain=None, range=Optional[Union[Union[dict, Result], List[Union[dict, Result]]]])

slots.empirical_modelling = Slot(uri=MATERIALSML.empirical_modelling, name="empirical_modelling", curie=MATERIALSML.curie('empirical_modelling'),
                   model_uri=MATERIALSML.empirical_modelling, domain=None, range=Optional[Union[dict, EmpiricalModelling]])

slots.grants = Slot(uri=SCHEMA.Grant, name="grants", curie=SCHEMA.curie('Grant'),
                   model_uri=MATERIALSML.grants, domain=None, range=Optional[Union[Union[dict, Grant], List[Union[dict, Grant]]]])

slots.drafts = Slot(uri=NCIT.NCIT_C85255, name="drafts", curie=NCIT.curie('NCIT_C85255'),
                   model_uri=MATERIALSML.drafts, domain=None, range=Optional[Union[Union[dict, Draft], List[Union[dict, Draft]]]])

slots.authors = Slot(uri=NCIT.NCIT_C42781, name="authors", curie=NCIT.curie('NCIT_C42781'),
                   model_uri=MATERIALSML.authors, domain=None, range=Optional[Union[Union[dict, Author], List[Union[dict, Author]]]])

slots.stms = Slot(uri=CHMO.CHMO_0000132, name="stms", curie=CHMO.curie('CHMO_0000132'),
                   model_uri=MATERIALSML.stms, domain=None, range=Optional[Union[Union[dict, STM], List[Union[dict, STM]]]])

slots.afms = Slot(uri=CHMO.CHMO_0000113, name="afms", curie=CHMO.curie('CHMO_0000113'),
                   model_uri=MATERIALSML.afms, domain=None, range=Optional[Union[Union[dict, AFM], List[Union[dict, AFM]]]])

slots.stss = Slot(uri=FAIRMAT_STS.nxsts, name="stss", curie=FAIRMAT_STS.curie('nxsts'),
                   model_uri=MATERIALSML.stss, domain=None, range=Optional[Union[Union[dict, STS], List[Union[dict, STS]]]])

slots.softwares = Slot(uri=SCHEMA.SoftwareSourceCode, name="softwares", curie=SCHEMA.curie('SoftwareSourceCode'),
                   model_uri=MATERIALSML.softwares, domain=None, range=Optional[Union[Union[dict, Software], List[Union[dict, Software]]]])

slots.geometry_optimisations = Slot(uri=MATERIALSML.geometry_optimisations, name="geometry_optimisations", curie=MATERIALSML.curie('geometry_optimisations'),
                   model_uri=MATERIALSML.geometry_optimisations, domain=None, range=Optional[Union[Union[dict, GeometryOptimisation], List[Union[dict, GeometryOptimisation]]]])

slots.layered_2d_material = Slot(uri=MATERIALSML.layered_2d_material, name="layered_2d_material", curie=MATERIALSML.curie('layered_2d_material'),
                   model_uri=MATERIALSML.layered_2d_material, domain=None, range=Optional[Union[dict, Layered2DMaterial]])

slots.wafer_substrate = Slot(uri=MATERIALSML.wafer_substrate, name="wafer_substrate", curie=MATERIALSML.curie('wafer_substrate'),
                   model_uri=MATERIALSML.wafer_substrate, domain=None, range=Optional[Union[dict, WaferSubstrate]])

slots.sample = Slot(uri=MATERIALSML.sample, name="sample", curie=MATERIALSML.curie('sample'),
                   model_uri=MATERIALSML.sample, domain=None, range=Optional[Union[dict, Sample]])

slots.cangasbottles = Slot(uri=MATERIALSML.cangasbottles, name="cangasbottles", curie=MATERIALSML.curie('cangasbottles'),
                   model_uri=MATERIALSML.cangasbottles, domain=None, range=Optional[Union[Union[dict, CanGasBottle], List[Union[dict, CanGasBottle]]]])

slots.maintenance = Slot(uri=NCIT.NCIT_C53297, name="maintenance", curie=NCIT.curie('NCIT_C53297'),
                   model_uri=MATERIALSML.maintenance, domain=None, range=Optional[Union[Union[dict, Maintenance], List[Union[dict, Maintenance]]]])

slots.errors = Slot(uri=NCIT.NCIT_C43369, name="errors", curie=NCIT.curie('NCIT_C43369'),
                   model_uri=MATERIALSML.errors, domain=None, range=Optional[Union[Union[dict, Errors], List[Union[dict, Errors]]]])

slots.location_before = Slot(uri=MATERIALSML.location_before, name="location_before", curie=MATERIALSML.curie('location_before'),
                   model_uri=MATERIALSML.location_before, domain=None, range=Optional[Union[dict, Room]])

slots.location_after = Slot(uri=MATERIALSML.location_after, name="location_after", curie=MATERIALSML.curie('location_after'),
                   model_uri=MATERIALSML.location_after, domain=None, range=Optional[Union[dict, Room]])

slots.container__afms = Slot(uri=MATERIALSML.afms, name="container__afms", curie=MATERIALSML.curie('afms'),
                   model_uri=MATERIALSML.container__afms, domain=None, range=Optional[Union[Union[dict, AFM], List[Union[dict, AFM]]]])

slots.container__annealings = Slot(uri=MATERIALSML.annealings, name="container__annealings", curie=MATERIALSML.curie('annealings'),
                   model_uri=MATERIALSML.container__annealings, domain=None, range=Optional[Union[Union[dict, Annealing], List[Union[dict, Annealing]]]])

slots.container__atomisticmodels = Slot(uri=MATERIALSML.atomisticmodels, name="container__atomisticmodels", curie=MATERIALSML.curie('atomisticmodels'),
                   model_uri=MATERIALSML.container__atomisticmodels, domain=None, range=Optional[Union[Union[dict, AtomisticModel], List[Union[dict, AtomisticModel]]]])

slots.container__authors = Slot(uri=MATERIALSML.authors, name="container__authors", curie=MATERIALSML.curie('authors'),
                   model_uri=MATERIALSML.container__authors, domain=None, range=Optional[Union[Union[dict, Author], List[Union[dict, Author]]]])

slots.container__band_structures = Slot(uri=MATERIALSML.band_structures, name="container__band_structures", curie=MATERIALSML.curie('band_structures'),
                   model_uri=MATERIALSML.container__band_structures, domain=None, range=Optional[Union[Union[dict, BandStructure], List[Union[dict, BandStructure]]]])

slots.container__can_gas_bottles = Slot(uri=MATERIALSML.can_gas_bottles, name="container__can_gas_bottles", curie=MATERIALSML.curie('can_gas_bottles'),
                   model_uri=MATERIALSML.container__can_gas_bottles, domain=None, range=Optional[Union[Union[dict, CanGasBottle], List[Union[dict, CanGasBottle]]]])

slots.container__chemists = Slot(uri=MATERIALSML.chemists, name="container__chemists", curie=MATERIALSML.curie('chemists'),
                   model_uri=MATERIALSML.container__chemists, domain=None, range=Optional[Union[Union[dict, Chemist], List[Union[dict, Chemist]]]])

slots.container__crystals = Slot(uri=MATERIALSML.crystals, name="container__crystals", curie=MATERIALSML.curie('crystals'),
                   model_uri=MATERIALSML.container__crystals, domain=None, range=Optional[Union[Union[dict, Crystal], List[Union[dict, Crystal]]]])

slots.container__depositions = Slot(uri=MATERIALSML.depositions, name="container__depositions", curie=MATERIALSML.curie('depositions'),
                   model_uri=MATERIALSML.container__depositions, domain=None, range=Optional[Union[Union[dict, Deposition], List[Union[dict, Deposition]]]])

slots.container__dfts = Slot(uri=MATERIALSML.dfts, name="container__dfts", curie=MATERIALSML.curie('dfts'),
                   model_uri=MATERIALSML.container__dfts, domain=None, range=Optional[Union[Union[dict, DFT], List[Union[dict, DFT]]]])

slots.container__dosings = Slot(uri=MATERIALSML.dosings, name="container__dosings", curie=MATERIALSML.curie('dosings'),
                   model_uri=MATERIALSML.container__dosings, domain=None, range=Optional[Union[Union[dict, Dosing], List[Union[dict, Dosing]]]])

slots.container__drafts = Slot(uri=MATERIALSML.drafts, name="container__drafts", curie=MATERIALSML.curie('drafts'),
                   model_uri=MATERIALSML.container__drafts, domain=None, range=Optional[Union[Union[dict, Draft], List[Union[dict, Draft]]]])

slots.container__empirical_modellings = Slot(uri=MATERIALSML.empirical_modellings, name="container__empirical_modellings", curie=MATERIALSML.curie('empirical_modellings'),
                   model_uri=MATERIALSML.container__empirical_modellings, domain=None, range=Optional[Union[Union[dict, EmpiricalModelling], List[Union[dict, EmpiricalModelling]]]])

slots.container__errors = Slot(uri=MATERIALSML.errors, name="container__errors", curie=MATERIALSML.curie('errors'),
                   model_uri=MATERIALSML.container__errors, domain=None, range=Optional[Union[Union[dict, Errors], List[Union[dict, Errors]]]])

slots.container__fillcryostats = Slot(uri=MATERIALSML.fillcryostats, name="container__fillcryostats", curie=MATERIALSML.curie('fillcryostats'),
                   model_uri=MATERIALSML.container__fillcryostats, domain=None, range=Optional[Union[Union[dict, FillCryostat], List[Union[dict, FillCryostat]]]])

slots.container__geoopts = Slot(uri=MATERIALSML.geoopts, name="container__geoopts", curie=MATERIALSML.curie('geoopts'),
                   model_uri=MATERIALSML.container__geoopts, domain=None, range=Optional[Union[Union[dict, GeometryOptimisation], List[Union[dict, GeometryOptimisation]]]])

slots.container__grants = Slot(uri=MATERIALSML.grants, name="container__grants", curie=MATERIALSML.curie('grants'),
                   model_uri=MATERIALSML.container__grants, domain=None, range=Optional[Union[Union[dict, Grant], List[Union[dict, Grant]]]])

slots.container__hydrogen_crackers = Slot(uri=MATERIALSML.hydrogen_crackers, name="container__hydrogen_crackers", curie=MATERIALSML.curie('hydrogen_crackers'),
                   model_uri=MATERIALSML.container__hydrogen_crackers, domain=None, range=Optional[Union[Union[dict, HydrogenCracker], List[Union[dict, HydrogenCracker]]]])

slots.container__institutions = Slot(uri=MATERIALSML.institutions, name="container__institutions", curie=MATERIALSML.curie('institutions'),
                   model_uri=MATERIALSML.container__institutions, domain=None, range=Optional[Union[Union[dict, Institution], List[Union[dict, Institution]]]])

slots.container__instruments = Slot(uri=MATERIALSML.instruments, name="container__instruments", curie=MATERIALSML.curie('instruments'),
                   model_uri=MATERIALSML.container__instruments, domain=None, range=Optional[Union[Union[dict, Instrument], List[Union[dict, Instrument]]]])

slots.container__layered2dmaterials = Slot(uri=MATERIALSML.layered2dmaterials, name="container__layered2dmaterials", curie=MATERIALSML.curie('layered2dmaterials'),
                   model_uri=MATERIALSML.container__layered2dmaterials, domain=None, range=Optional[Union[Union[dict, Layered2DMaterial], List[Union[dict, Layered2DMaterial]]]])

slots.container__maintenances = Slot(uri=MATERIALSML.maintenances, name="container__maintenances", curie=MATERIALSML.curie('maintenances'),
                   model_uri=MATERIALSML.container__maintenances, domain=None, range=Optional[Union[Union[dict, Maintenance], List[Union[dict, Maintenance]]]])

slots.container__manufacturers = Slot(uri=MATERIALSML.manufacturers, name="container__manufacturers", curie=MATERIALSML.curie('manufacturers'),
                   model_uri=MATERIALSML.container__manufacturers, domain=None, range=Optional[Union[Union[dict, Manufacturer], List[Union[dict, Manufacturer]]]])

slots.container__meps = Slot(uri=MATERIALSML.meps, name="container__meps", curie=MATERIALSML.curie('meps'),
                   model_uri=MATERIALSML.container__meps, domain=None, range=Optional[Union[Union[dict, MinimumEnergyPotential], List[Union[dict, MinimumEnergyPotential]]]])

slots.container__molecules = Slot(uri=MATERIALSML.molecules, name="container__molecules", curie=MATERIALSML.curie('molecules'),
                   model_uri=MATERIALSML.container__molecules, domain=None, range=Optional[Union[Union[dict, Molecule], List[Union[dict, Molecule]]]])

slots.container__notes = Slot(uri=MATERIALSML.notes, name="container__notes", curie=MATERIALSML.curie('notes'),
                   model_uri=MATERIALSML.container__notes, domain=None, range=Optional[Union[Union[dict, Notes], List[Union[dict, Notes]]]])

slots.container__pdoss = Slot(uri=MATERIALSML.pdoss, name="container__pdoss", curie=MATERIALSML.curie('pdoss'),
                   model_uri=MATERIALSML.container__pdoss, domain=None, range=Optional[Union[Union[dict, PDOS], List[Union[dict, PDOS]]]])

slots.container__persons = Slot(uri=MATERIALSML.persons, name="container__persons", curie=MATERIALSML.curie('persons'),
                   model_uri=MATERIALSML.container__persons, domain=None, range=Optional[Union[Union[dict, Person], List[Union[dict, Person]]]])

slots.container__protocols = Slot(uri=MATERIALSML.protocols, name="container__protocols", curie=MATERIALSML.curie('protocols'),
                   model_uri=MATERIALSML.container__protocols, domain=None, range=Optional[Union[Union[dict, Protocol], List[Union[dict, Protocol]]]])

slots.container__publications = Slot(uri=MATERIALSML.publications, name="container__publications", curie=MATERIALSML.curie('publications'),
                   model_uri=MATERIALSML.container__publications, domain=None, range=Optional[Union[Union[dict, Publication], List[Union[dict, Publication]]]])

slots.container__results = Slot(uri=MATERIALSML.results, name="container__results", curie=MATERIALSML.curie('results'),
                   model_uri=MATERIALSML.container__results, domain=None, range=Optional[Union[Union[dict, Result], List[Union[dict, Result]]]])

slots.container__rooms = Slot(uri=MATERIALSML.rooms, name="container__rooms", curie=MATERIALSML.curie('rooms'),
                   model_uri=MATERIALSML.container__rooms, domain=None, range=Optional[Union[Union[dict, Room], List[Union[dict, Room]]]])

slots.container__samples = Slot(uri=MATERIALSML.samples, name="container__samples", curie=MATERIALSML.curie('samples'),
                   model_uri=MATERIALSML.container__samples, domain=None, range=Optional[Union[Union[dict, Sample], List[Union[dict, Sample]]]])

slots.container__softwares = Slot(uri=MATERIALSML.softwares, name="container__softwares", curie=MATERIALSML.curie('softwares'),
                   model_uri=MATERIALSML.container__softwares, domain=None, range=Optional[Union[Union[dict, Software], List[Union[dict, Software]]]])

slots.container__sputterings = Slot(uri=MATERIALSML.sputterings, name="container__sputterings", curie=MATERIALSML.curie('sputterings'),
                   model_uri=MATERIALSML.container__sputterings, domain=None, range=Optional[Union[Union[dict, Sputtering], List[Union[dict, Sputtering]]]])

slots.container__srds = Slot(uri=MATERIALSML.srds, name="container__srds", curie=MATERIALSML.curie('srds'),
                   model_uri=MATERIALSML.container__srds, domain=None, range=Optional[Union[Union[dict, SRD], List[Union[dict, SRD]]]])

slots.container__stms = Slot(uri=MATERIALSML.stms, name="container__stms", curie=MATERIALSML.curie('stms'),
                   model_uri=MATERIALSML.container__stms, domain=None, range=Optional[Union[Union[dict, STM], List[Union[dict, STM]]]])

slots.container__storages = Slot(uri=MATERIALSML.storages, name="container__storages", curie=MATERIALSML.curie('storages'),
                   model_uri=MATERIALSML.container__storages, domain=None, range=Optional[Union[Union[dict, Storage], List[Union[dict, Storage]]]])

slots.container__stss = Slot(uri=MATERIALSML.stss, name="container__stss", curie=MATERIALSML.curie('stss'),
                   model_uri=MATERIALSML.container__stss, domain=None, range=Optional[Union[Union[dict, STS], List[Union[dict, STS]]]])

slots.container__suppliers = Slot(uri=MATERIALSML.suppliers, name="container__suppliers", curie=MATERIALSML.curie('suppliers'),
                   model_uri=MATERIALSML.container__suppliers, domain=None, range=Optional[Union[Union[dict, Supplier], List[Union[dict, Supplier]]]])

slots.container__tight_bindings = Slot(uri=MATERIALSML.tight_bindings, name="container__tight_bindings", curie=MATERIALSML.curie('tight_bindings'),
                   model_uri=MATERIALSML.container__tight_bindings, domain=None, range=Optional[Union[Union[dict, TightBinding], List[Union[dict, TightBinding]]]])

slots.container__tip_preparations = Slot(uri=MATERIALSML.tip_preparations, name="container__tip_preparations", curie=MATERIALSML.curie('tip_preparations'),
                   model_uri=MATERIALSML.container__tip_preparations, domain=None, range=Optional[Union[Union[dict, TipPreparation], List[Union[dict, TipPreparation]]]])

slots.container__tip_sensors = Slot(uri=MATERIALSML.tip_sensors, name="container__tip_sensors", curie=MATERIALSML.curie('tip_sensors'),
                   model_uri=MATERIALSML.container__tip_sensors, domain=None, range=Optional[Union[Union[dict, TipSensor], List[Union[dict, TipSensor]]]])

slots.container__transfers = Slot(uri=MATERIALSML.transfers, name="container__transfers", curie=MATERIALSML.curie('transfers'),
                   model_uri=MATERIALSML.container__transfers, domain=None, range=Optional[Union[Union[dict, Transfer], List[Union[dict, Transfer]]]])

slots.container__uhv_components = Slot(uri=MATERIALSML.uhv_components, name="container__uhv_components", curie=MATERIALSML.curie('uhv_components'),
                   model_uri=MATERIALSML.container__uhv_components, domain=None, range=Optional[Union[Union[dict, UHVComponent], List[Union[dict, UHVComponent]]]])

slots.container__wafer_substrates = Slot(uri=MATERIALSML.wafer_substrates, name="container__wafer_substrates", curie=MATERIALSML.curie('wafer_substrates'),
                   model_uri=MATERIALSML.container__wafer_substrates, domain=None, range=Optional[Union[Union[dict, WaferSubstrate], List[Union[dict, WaferSubstrate]]]])

slots.Molecule_name = Slot(uri=SCHEMA.name, name="Molecule_name", curie=SCHEMA.curie('name'),
                   model_uri=MATERIALSML.Molecule_name, domain=Molecule, range=str)
