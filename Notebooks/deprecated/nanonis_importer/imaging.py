#   Copyright ETH 2023 - 2024 Zürich, Scientific IT Services
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


from pybis import Openbis
import abc
import json
import base64
import requests
import os
import threading
from urllib.parse import urljoin


DEFAULT_SERVICE_NAME = "imaging"
IMAGING_CONFIG_PROP_NAME = "IMAGING_DATA_CONFIG".lower()
DEFAULT_DATASET_VIEW_PROP_NAME = "default_dataset_view"


def get_instance(url="http://localhost:8888/openbis"):
    openbis_instance = Openbis(
        url=url,
        verify_certificates=False,
        allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks=True
    )
    token = openbis_instance.login('admin', 'changeit')
    print(f'Connected to {url} -> token: {token}')
    return openbis_instance


class AtomicIncrementer:
    def __init__(self, value=0):
        self._value = int(value)
        self._lock = threading.Lock()

    def inc(self, d=1):
        with self._lock:
            self._value += int(d)
            return self._value


class AbstractImagingClass(metaclass=abc.ABCMeta):
    def to_json(self):

        c = AtomicIncrementer()

        def dictionary_creator(x):
            dictionary = x.__dict__
            val = c.inc()
            dictionary['@id'] = val
            return dictionary

        return json.dumps(self, default=dictionary_creator, sort_keys=True, indent=4)

    def __str__(self):
        return json.dumps(self.__dict__, default=lambda x: x.__dict__)

    def __repr__(self):
        return json.dumps(self.__dict__, default=lambda x: x.__dict__)


class AbstractImagingRequest(AbstractImagingClass, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def _validate_data(self):
        return


class ImagingDataSetPreview(AbstractImagingRequest):
    config: dict
    format: str
    bytes: str
    width: int
    height: int
    index: int
    show: bool
    metadata: dict
    comment: str
    tags: list

    def __init__(self, preview_format, config=None, metadata=None, index=0, comment="", tags=[]):
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetPreview"
        self.bytes = None
        self.format = preview_format
        self.config = config if config is not None else dict()
        self.metadata = metadata if metadata is not None else dict()
        self.index = index
        self.comment = comment
        self.tags = tags
        self._validate_data()

    def set_preview_image_bytes(self, width, height, bytes):
        self.width = width
        self.height = height
        self.bytes = bytes

    def _validate_data(self):
        assert self.format is not None, "Format can not be null"

    def save_to_file(self, file_path):
        assert self.bytes is not None, "There is no image information!"
        img_data = bytearray(self.bytes, encoding='utf-8')
        with open(file_path, "wb") as fh:
            fh.write(base64.decodebytes(img_data))

    @classmethod
    def from_dict(cls, data):
        if data is None:
            return None
        if "@id" in data:
            del data["@id"]
        preview = cls('', None, None)
        for prop in cls.__annotations__.keys():
            attribute = data.get(prop)
            preview.__dict__[prop] = attribute
        return preview

class ImagingDataSetExportConfig(AbstractImagingClass):
    archiveFormat: str
    imageFormat: str
    resolution: str
    include: list

    def __init__(self, archive_format, image_format, resolution, include=None):
        if include is None:
            include = ["IMAGE", "RAW_DATA"]
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetExportConfig"
        self.imageFormat = image_format
        self.archiveFormat = archive_format
        if resolution is None:
            resolution = "original"
        self.resolution = resolution
        self.include = include
        self._validate_data()

    def _validate_data(self):
        assert self.imageFormat is not None, "image format can not be null"
        assert self.archiveFormat is not None, "image format can not be null"

    @classmethod
    def from_dict(cls, data):
        if data is None:
            return None
        if "@id" in data:
            del data["@id"]
        preview = cls(None, None, None)
        for prop in cls.__annotations__.keys():
            attribute = data.get(prop)
            preview.__dict__[prop] = attribute
        return preview



class ImagingDataSetExport(AbstractImagingRequest):
    config: ImagingDataSetExportConfig
    metadata: dict

    def __init__(self, config, metadata=None):
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetExport"
        self.config = config
        self.metadata = metadata if metadata is not None else dict()
        self._validate_data()

    def _validate_data(self):
        assert self.config is not None, "Config can not be null"



class ImagingDataSetMultiExport(AbstractImagingRequest):
    permId: str
    imageIndex: int
    previewIndex: int
    config: ImagingDataSetExportConfig
    metadata: dict

    def __init__(self, permId, imageIndex, previewIndex, config, metadata=None):
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetMultiExport"
        self.permId = permId
        self.imageIndex = imageIndex
        self.previewIndex = previewIndex
        self.config = config
        self.metadata = metadata if metadata is not None else dict()
        self._validate_data()

    def _validate_data(self):
        assert self.permId is not None, "PermId can not be null"
        assert self.imageIndex is not None, "imageIndex can not be null"
        assert self.previewIndex is not None, "previewIndex can not be null"


class ImagingDataSetControlVisibility(AbstractImagingClass):
    label: str
    values: list[str]
    range: list[str]
    unit: str

    def __init__(self, label: str, values: list[str], values_range: list[str], unit: str = None):
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetControlVisibility"
        self.label = label
        self.values = values
        self.range = values_range
        self.unit = unit

    @classmethod
    def from_dict(cls, data):
        if data is None:
            return None
        if "@id" in data:
            del data["@id"]
        control = cls(None, None, None, None)
        for prop in cls.__annotations__.keys():
            attribute = data.get(prop)
            control.__dict__[prop] = attribute
        return control


class ImagingDataSetControl(AbstractImagingClass):
    label: str
    section: str
    type: str
    values: list[str]
    unit: str
    range: list[str]
    multiselect: bool
    playable: bool
    speeds: list[int]
    visibility: list[ImagingDataSetControlVisibility]
    metadata: dict

    def __init__(self, label: str, control_type: str, section: str = None, values: list[str] = None,
                 unit: str = None, values_range: list[str] = None, multiselect: bool = None,
                 playable: bool = False, speeds: list[int] = None,
                 visibility: list[ImagingDataSetControlVisibility] = None, metadata: dict = None):
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetControl"
        self.label = label
        self.type = control_type
        self.section = section
        self.unit = unit
        if control_type.lower() in ["slider", "range"]:
            self.range = values_range
        elif control_type.lower() in ["dropdown", "colormap"]:
            self.values = values
            if multiselect is None:
                self.multiselect = False
            else:
                self.multiselect = multiselect

        if playable is True:
            self.playable = True
            self.speeds = speeds
        self.visibility = visibility
        self.metadata = metadata

    @classmethod
    def from_dict(cls, data):
        if data is None:
            return None
        if "@id" in data:
            del data["@id"]
        control = cls(None, "", None, None)
        for prop in cls.__annotations__.keys():
            attribute = data.get(prop)
            if prop == 'visibility' and attribute is not None:
                attribute = [ImagingDataSetControlVisibility.from_dict(visibility) for visibility in attribute]
            control.__dict__[prop] = attribute
        return control


class ImagingDataSetConfig(AbstractImagingClass):
    adaptor: str
    version: float
    speeds: list[int]
    resolutions: list[str]
    playable: bool
    exports: list[ImagingDataSetControl]
    inputs: list[ImagingDataSetControl]
    metadata: dict

    def __init__(self, adaptor: str, version: float, resolutions: list[str], playable: bool,
                 speeds: list[int] = None, exports: list[ImagingDataSetControl] = None,
                 inputs: list[ImagingDataSetControl] = None, metadata: dict = None):
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetConfig"
        self.adaptor = adaptor
        self.version = version
        self.resolutions = resolutions
        self.playable = playable
        if playable:
            self.speeds = speeds
        self.exports = exports
        self.inputs = inputs
        self.metadata = metadata

    @classmethod
    def from_dict(cls, data):
        if data is None:
            return None
        if "@id" in data:
            del data["@id"]
        config = cls(None, None, None, None)
        for prop in cls.__annotations__.keys():
            attribute = data.get(prop)
            if prop in ['exports', 'inputs'] and attribute is not None:
                attribute = [ImagingDataSetControl.from_dict(control) for control in attribute]
            config.__dict__[prop] = attribute
        return config


class ImagingDataSetImage(AbstractImagingClass):
    config: ImagingDataSetConfig
    previews: list[ImagingDataSetPreview]
    image_config: dict
    index: int
    metadata: dict

    def __init__(self, config: ImagingDataSetConfig, image_config=None, previews=None, metadata=None, index=0):
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetImage"
        assert config is not None, "Config must not be None!"
        self.config = config
        self.image_config = image_config if image_config is not None else dict()
        self.previews = previews if previews is not None else [ImagingDataSetPreview("png")]
        self.metadata = metadata if metadata is not None else dict()
        self.index = index if index is not None else 0
        assert isinstance(self.previews, list), "Previews must be a list!"

    def add_preview(self, preview):
        self.previews += [preview]

    @classmethod
    def from_dict(cls, data):
        if data is None:
            return None
        if "@id" in data:
            del data["@id"]
        config = ImagingDataSetConfig.from_dict(data.get('config'))
        image = cls(config,None, None, None)
        for prop in cls.__annotations__.keys():
            attribute = data.get(prop)
            if prop == 'previews' and attribute is not None:
                attribute = [ImagingDataSetPreview.from_dict(preview) for preview in attribute]
            if prop not in ['config']:
                image.__dict__[prop] = attribute
        return image


class ImagingDataSetPropertyConfig(AbstractImagingClass):
    images: list[ImagingDataSetImage]
    metadata: dict

    def __init__(self, images: list[ImagingDataSetImage], metadata=None):
        self.__dict__["@type"] = "imaging.dto.ImagingDataSetPropertyConfig"
        self.images = images if images is not None else []
        self.metadata = metadata if metadata is not None else dict()

    @classmethod
    def from_dict(cls, data: dict):
        assert data is not None and any(data), "There is no property config found!"
        if "@id" in data:
            del data["@id"]
        attr = data.get('images')
        images = [ImagingDataSetImage.from_dict(image) for image in attr] if attr is not None else None
        metadata = data.get('metadata')
        return cls(images, metadata)

    def add_image(self, image: ImagingDataSetImage):
        if self.images is None:
            self.images = []
        self.images += [image]


class ImagingControl:

    def __init__(self, openbis_instance, service_name=DEFAULT_SERVICE_NAME):
        self._openbis = openbis_instance
        self._service_name = service_name

    def _execute_custom_dss_service(self, parameters):
        service_id = {
            "@type": "dss.dto.service.id.CustomDssServiceCode",
            "permId": self._service_name
        }
        options = {
            "@type": "dss.dto.service.CustomDSSServiceExecutionOptions",
            "parameters": parameters
        }
        request = {
            "method": "executeCustomDSSService",
            "params": [
                self._openbis.token,
                service_id,
                options
            ],
        }
        full_url = urljoin(self._openbis._get_dss_url(), self._openbis.dss_v3)
        return self._openbis._post_request_full_url(full_url, request)

    def make_preview(self, perm_id: str, index: int, preview: ImagingDataSetPreview) -> ImagingDataSetPreview:
        parameters = {
            "type": "preview",
            "permId": perm_id,
            "index": index,
            "error": None,
            "preview": preview.__dict__
        }
        service_response = self._execute_custom_dss_service(parameters)
        if service_response['error'] is None:
            if '@id' in service_response:
                del service_response['@id']
            if '@id' in service_response['preview']:
                del service_response['preview']['@id']
            preview.__dict__ = service_response["preview"]
            return preview
        else:
            raise ValueError(service_response['error'])

    def _get_export_url(self, perm_id: str, export: ImagingDataSetExport, image_index: int = 0) -> str:
        export_params = export.__dict__
        export_params["config"] = export_params["config"].__dict__
        parameters = {
            "type": "export",
            "permId": perm_id,
            "index": image_index,
            "error": None,
            "url": None,
            "export": export_params
        }
        service_response = self._openbis.execute_custom_dss_service(self._service_name, parameters)
        if service_response['error'] is None:
            return service_response['url']
        else:
            raise ValueError(service_response['error'])

    def _get_multi_export_url(self, exports: list[ImagingDataSetMultiExport]) -> str:
        export_params = [export.__dict__ for export in exports]
        for param in export_params:
            param["config"] = param["config"].__dict__
        parameters = {
            "type": "multi-export",
            "error": None,
            "url": None,
            "exports": [export.__dict__ for export in exports]
        }
        service_response = self._openbis.execute_custom_dss_service(self._service_name, parameters)
        if service_response['error'] is None:
            return service_response['url']
        else:
            raise ValueError(service_response['error'])

    def single_export_download(self, perm_id: str, export: ImagingDataSetExport, image_index: int = 0, directory_path=""):
        export_url = self._get_export_url(perm_id, export, image_index)
        self._download(export_url, directory_path)

    def multi_export_download(self, exports: list[ImagingDataSetMultiExport], directory_path=""):
        export_url = self._get_multi_export_url(exports)
        self._download(export_url, directory_path)

    def _download(self, url, directory_path=""):
        get_response = requests.get(url, stream=True, verify=self._openbis.verify_certificates)
        file_name = url.split("/")[-1]
        path = os.path.join(directory_path, file_name)
        with open(path, 'wb') as f:
            for chunk in get_response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)

    def get_property_config(self, perm_id: str) -> ImagingDataSetPropertyConfig:
        dataset = self._openbis.get_dataset(perm_id)
        imaging_property = json.loads(dataset.props[IMAGING_CONFIG_PROP_NAME])
        return ImagingDataSetPropertyConfig.from_dict(imaging_property)

    def update_property_config(self, perm_id: str, config: ImagingDataSetPropertyConfig):
        dataset = self._openbis.get_dataset(perm_id)
        dataset.props[IMAGING_CONFIG_PROP_NAME] = config.to_json()
        dataset.save()

    def create_imaging_dataset(self, dataset_type: str, config: ImagingDataSetPropertyConfig,
                               experiment: str, sample: str,
                               files: list[str], other_properties=None):
        if other_properties is None:
            other_properties = {}
        assert dataset_type is not None
        assert files is not None and len(files) > 0, "Files parameter must not be empty!"
        assert config is not None
        props = other_properties
        props[IMAGING_CONFIG_PROP_NAME] = config.to_json()
        props[DEFAULT_DATASET_VIEW_PROP_NAME] = 'IMAGING_DATASET_VIEWER'
        dataset = self._openbis.new_dataset(dataset_type, experiment=experiment, sample=sample, files=files, props=props)
        return dataset.save()





