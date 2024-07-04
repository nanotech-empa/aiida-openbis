# Auto generated from moleculeinfo.yaml by pythongen.py version: 0.0.1
# Generation date: 2024-07-02T14:01:06
# Schema: moleculeinfo
#
# id: https://w3id.org/linkml/examples/moleculeinfo
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
EMMO = CurieNamespace('emmo', 'https://w3id.org/emmo#')
LINKML = CurieNamespace('linkml', 'https://w3id.org/linkml/')
MOLECULEINFO = CurieNamespace('moleculeinfo', 'https://w3id.org/linkml/examples/moleculeinfo/')
NCIT = CurieNamespace('ncit', 'http://purl.obolibrary.org/obo/ncit.owl#')
OCCO = CurieNamespace('occo', 'http://purl.obolibrary.org/obo/occo.owl#')
SCHEMA = CurieNamespace('schema', 'http://schema.org/')
XSD = CurieNamespace('xsd', 'http://www.w3.org/2001/XMLSchema#')
DEFAULT_ = MOLECULEINFO


# Types
class CustomDatetime(str):
    type_class_uri = XSD["dateTime"]
    type_class_curie = "xsd:dateTime"
    type_name = "custom_datetime"
    type_model_uri = MOLECULEINFO.CustomDatetime


# Class references



@dataclass
class Chemist(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = OCCO["OCCO_19203100"]
    class_class_curie: ClassVar[str] = "occo:OCCO_19203100"
    class_name: ClassVar[str] = "Chemist"
    class_model_uri: ClassVar[URIRef] = MOLECULEINFO.Chemist

    name: str = None
    email: Optional[str] = None
    work_phone: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.name):
            self.MissingRequiredField("name")
        if not isinstance(self.name, str):
            self.name = str(self.name)

        if self.email is not None and not isinstance(self.email, str):
            self.email = str(self.email)

        if self.work_phone is not None and not isinstance(self.work_phone, str):
            self.work_phone = str(self.work_phone)

        super().__post_init__(**kwargs)


@dataclass
class Molecule(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["MolecularEntity"]
    class_class_curie: ClassVar[str] = "schema:MolecularEntity"
    class_name: ClassVar[str] = "Molecule"
    class_model_uri: ClassVar[URIRef] = MOLECULEINFO.Molecule

    name: str = None
    sum_formula: str = None
    smiles: str = None
    empa_number: int = None
    batch: str = None
    iupac_name: Optional[str] = None
    cas_number: Optional[str] = None
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
    amount: Optional[Union[dict, "QuantitativeValue"]] = None
    receive_date: Optional[Union[str, XSDDate]] = None
    comments: Optional[str] = None
    chemist: Optional[Union[dict, Chemist]] = None
    parent_molecules: Optional[Union[Union[dict, "Molecule"], List[Union[dict, "Molecule"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.name):
            self.MissingRequiredField("name")
        if not isinstance(self.name, str):
            self.name = str(self.name)

        if self._is_empty(self.sum_formula):
            self.MissingRequiredField("sum_formula")
        if not isinstance(self.sum_formula, str):
            self.sum_formula = str(self.sum_formula)

        if self._is_empty(self.smiles):
            self.MissingRequiredField("smiles")
        if not isinstance(self.smiles, str):
            self.smiles = str(self.smiles)

        if self._is_empty(self.empa_number):
            self.MissingRequiredField("empa_number")
        if not isinstance(self.empa_number, int):
            self.empa_number = int(self.empa_number)

        if self._is_empty(self.batch):
            self.MissingRequiredField("batch")
        if not isinstance(self.batch, str):
            self.batch = str(self.batch)

        if self.iupac_name is not None and not isinstance(self.iupac_name, str):
            self.iupac_name = str(self.iupac_name)

        if self.cas_number is not None and not isinstance(self.cas_number, str):
            self.cas_number = str(self.cas_number)

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

        if self.amount is not None and not isinstance(self.amount, QuantitativeValue):
            self.amount = QuantitativeValue(**as_dict(self.amount))

        if self.receive_date is not None and not isinstance(self.receive_date, XSDDate):
            self.receive_date = XSDDate(self.receive_date)

        if self.comments is not None and not isinstance(self.comments, str):
            self.comments = str(self.comments)

        if self.chemist is not None and not isinstance(self.chemist, Chemist):
            self.chemist = Chemist(**as_dict(self.chemist))

        if not isinstance(self.parent_molecules, list):
            self.parent_molecules = [self.parent_molecules] if self.parent_molecules is not None else []
        self.parent_molecules = [v if isinstance(v, Molecule) else Molecule(**as_dict(v)) for v in self.parent_molecules]

        super().__post_init__(**kwargs)


@dataclass
class QuantitativeValue(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SCHEMA["QuantitativeValue"]
    class_class_curie: ClassVar[str] = "schema:QuantitativeValue"
    class_name: ClassVar[str] = "QuantitativeValue"
    class_model_uri: ClassVar[URIRef] = MOLECULEINFO.QuantitativeValue

    value: Optional[float] = None
    unit: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.value is not None and not isinstance(self.value, float):
            self.value = float(self.value)

        if self.unit is not None and not isinstance(self.unit, str):
            self.unit = str(self.unit)

        super().__post_init__(**kwargs)


@dataclass
class EvaporationTemperature(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MOLECULEINFO["EvaporationTemperature"]
    class_class_curie: ClassVar[str] = "moleculeinfo:EvaporationTemperature"
    class_name: ClassVar[str] = "EvaporationTemperature"
    class_model_uri: ClassVar[URIRef] = MOLECULEINFO.EvaporationTemperature

    temperature: Optional[Union[dict, QuantitativeValue]] = None
    device: Optional[str] = None
    datetime: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.temperature is not None and not isinstance(self.temperature, QuantitativeValue):
            self.temperature = QuantitativeValue(**as_dict(self.temperature))

        if self.device is not None and not isinstance(self.device, str):
            self.device = str(self.device)

        if self.datetime is not None and not isinstance(self.datetime, str):
            self.datetime = str(self.datetime)

        super().__post_init__(**kwargs)


@dataclass
class Container(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = MOLECULEINFO["Container"]
    class_class_curie: ClassVar[str] = "moleculeinfo:Container"
    class_name: ClassVar[str] = "Container"
    class_model_uri: ClassVar[URIRef] = MOLECULEINFO.Container

    molecules: Optional[Union[Union[dict, Molecule], List[Union[dict, Molecule]]]] = empty_list()
    chemists: Optional[Union[Union[dict, Chemist], List[Union[dict, Chemist]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if not isinstance(self.molecules, list):
            self.molecules = [self.molecules] if self.molecules is not None else []
        self.molecules = [v if isinstance(v, Molecule) else Molecule(**as_dict(v)) for v in self.molecules]

        if not isinstance(self.chemists, list):
            self.chemists = [self.chemists] if self.chemists is not None else []
        self.chemists = [v if isinstance(v, Chemist) else Chemist(**as_dict(v)) for v in self.chemists]

        super().__post_init__(**kwargs)


# Enumerations


# Slots
class slots:
    pass

slots.name = Slot(uri=SCHEMA.name, name="name", curie=SCHEMA.curie('name'),
                   model_uri=MOLECULEINFO.name, domain=None, range=str)

slots.work_phone = Slot(uri=SCHEMA.telephone, name="work_phone", curie=SCHEMA.curie('telephone'),
                   model_uri=MOLECULEINFO.work_phone, domain=None, range=Optional[str])

slots.email = Slot(uri=SCHEMA.email, name="email", curie=SCHEMA.curie('email'),
                   model_uri=MOLECULEINFO.email, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$'))

slots.iupac_name = Slot(uri=SCHEMA.iupacName, name="iupac_name", curie=SCHEMA.curie('iupacName'),
                   model_uri=MOLECULEINFO.iupac_name, domain=None, range=Optional[str])

slots.sum_formula = Slot(uri=SCHEMA.molecularFormula, name="sum_formula", curie=SCHEMA.curie('molecularFormula'),
                   model_uri=MOLECULEINFO.sum_formula, domain=None, range=str)

slots.smiles = Slot(uri=SCHEMA.smiles, name="smiles", curie=SCHEMA.curie('smiles'),
                   model_uri=MOLECULEINFO.smiles, domain=None, range=str)

slots.cas_number = Slot(uri=EMMO.EMMO_d2a47cd8_662f_438f_855a_b4378eb992ff, name="cas_number", curie=EMMO.curie('EMMO_d2a47cd8_662f_438f_855a_b4378eb992ff'),
                   model_uri=MOLECULEINFO.cas_number, domain=None, range=Optional[str])

slots.empa_number = Slot(uri=MOLECULEINFO.empa_number, name="empa_number", curie=MOLECULEINFO.curie('empa_number'),
                   model_uri=MOLECULEINFO.empa_number, domain=None, range=int)

slots.batch = Slot(uri=NCIT.NCIT_C67073, name="batch", curie=NCIT.curie('NCIT_C67073'),
                   model_uri=MOLECULEINFO.batch, domain=None, range=str,
                   pattern=re.compile(r'^[a-z]$'))

slots.vial = Slot(uri=AFE.AFE_0000329, name="vial", curie=AFE.curie('AFE_0000329'),
                   model_uri=MOLECULEINFO.vial, domain=None, range=Optional[str],
                   pattern=re.compile(r'^(i|ii|iii|iv|v|vi|vii|viii|ix|x)$'))

slots.hazardous = Slot(uri=NCIT.NCIT_C73538, name="hazardous", curie=NCIT.curie('NCIT_C73538'),
                   model_uri=MOLECULEINFO.hazardous, domain=None, range=Optional[Union[bool, Bool]])

slots.hazardous_specification = Slot(uri=MOLECULEINFO.hazardous_specification, name="hazardous_specification", curie=MOLECULEINFO.curie('hazardous_specification'),
                   model_uri=MOLECULEINFO.hazardous_specification, domain=None, range=Optional[str])

slots.evaporation_temperatures = Slot(uri=MOLECULEINFO.evaporation_temperatures, name="evaporation_temperatures", curie=MOLECULEINFO.curie('evaporation_temperatures'),
                   model_uri=MOLECULEINFO.evaporation_temperatures, domain=None, range=Optional[Union[Union[dict, EvaporationTemperature], List[Union[dict, EvaporationTemperature]]]])

slots.fridge = Slot(uri=MOLECULEINFO.fridge, name="fridge", curie=MOLECULEINFO.curie('fridge'),
                   model_uri=MOLECULEINFO.fridge, domain=None, range=Optional[Union[bool, Bool]])

slots.no_light = Slot(uri=MOLECULEINFO.no_light, name="no_light", curie=MOLECULEINFO.curie('no_light'),
                   model_uri=MOLECULEINFO.no_light, domain=None, range=Optional[Union[bool, Bool]])

slots.dry = Slot(uri=MOLECULEINFO.dry, name="dry", curie=MOLECULEINFO.curie('dry'),
                   model_uri=MOLECULEINFO.dry, domain=None, range=Optional[Union[bool, Bool]])

slots.no_oxygen = Slot(uri=MOLECULEINFO.no_oxygen, name="no_oxygen", curie=MOLECULEINFO.curie('no_oxygen'),
                   model_uri=MOLECULEINFO.no_oxygen, domain=None, range=Optional[Union[bool, Bool]])

slots.other_storage_condition = Slot(uri=MOLECULEINFO.other_storage_condition, name="other_storage_condition", curie=MOLECULEINFO.curie('other_storage_condition'),
                   model_uri=MOLECULEINFO.other_storage_condition, domain=None, range=Optional[Union[bool, Bool]])

slots.other_storage_condition_specification = Slot(uri=MOLECULEINFO.other_storage_condition_specification, name="other_storage_condition_specification", curie=MOLECULEINFO.curie('other_storage_condition_specification'),
                   model_uri=MOLECULEINFO.other_storage_condition_specification, domain=None, range=Optional[str])

slots.chemist_molecule_name = Slot(uri=SCHEMA.alternateName, name="chemist_molecule_name", curie=SCHEMA.curie('alternateName'),
                   model_uri=MOLECULEINFO.chemist_molecule_name, domain=None, range=Optional[str])

slots.amount = Slot(uri=SCHEMA.QuantitativeValue, name="amount", curie=SCHEMA.curie('QuantitativeValue'),
                   model_uri=MOLECULEINFO.amount, domain=None, range=Optional[Union[dict, QuantitativeValue]])

slots.receive_date = Slot(uri=NCIT.NCIT_C164590, name="receive_date", curie=NCIT.curie('NCIT_C164590'),
                   model_uri=MOLECULEINFO.receive_date, domain=None, range=Optional[Union[str, XSDDate]])

slots.comments = Slot(uri=SCHEMA.comment, name="comments", curie=SCHEMA.curie('comment'),
                   model_uri=MOLECULEINFO.comments, domain=None, range=Optional[str])

slots.value = Slot(uri=MOLECULEINFO.value, name="value", curie=MOLECULEINFO.curie('value'),
                   model_uri=MOLECULEINFO.value, domain=None, range=Optional[float])

slots.unit = Slot(uri=MOLECULEINFO.unit, name="unit", curie=MOLECULEINFO.curie('unit'),
                   model_uri=MOLECULEINFO.unit, domain=None, range=Optional[str])

slots.device = Slot(uri=MOLECULEINFO.device, name="device", curie=MOLECULEINFO.curie('device'),
                   model_uri=MOLECULEINFO.device, domain=None, range=Optional[str])

slots.temperature = Slot(uri=MOLECULEINFO.temperature, name="temperature", curie=MOLECULEINFO.curie('temperature'),
                   model_uri=MOLECULEINFO.temperature, domain=None, range=Optional[Union[dict, QuantitativeValue]])

slots.datetime = Slot(uri=MOLECULEINFO.datetime, name="datetime", curie=MOLECULEINFO.curie('datetime'),
                   model_uri=MOLECULEINFO.datetime, domain=None, range=Optional[str])

slots.chemist = Slot(uri=OCCO.OCCO_19203100, name="chemist", curie=OCCO.curie('OCCO_19203100'),
                   model_uri=MOLECULEINFO.chemist, domain=None, range=Optional[Union[dict, Chemist]])

slots.parent_molecules = Slot(uri=SCHEMA.MolecularEntity, name="parent_molecules", curie=SCHEMA.curie('MolecularEntity'),
                   model_uri=MOLECULEINFO.parent_molecules, domain=None, range=Optional[Union[Union[dict, Molecule], List[Union[dict, Molecule]]]])

slots.container__molecules = Slot(uri=MOLECULEINFO.molecules, name="container__molecules", curie=MOLECULEINFO.curie('molecules'),
                   model_uri=MOLECULEINFO.container__molecules, domain=None, range=Optional[Union[Union[dict, Molecule], List[Union[dict, Molecule]]]])

slots.container__chemists = Slot(uri=MOLECULEINFO.chemists, name="container__chemists", curie=MOLECULEINFO.curie('chemists'),
                   model_uri=MOLECULEINFO.container__chemists, domain=None, range=Optional[Union[Union[dict, Chemist], List[Union[dict, Chemist]]]])
