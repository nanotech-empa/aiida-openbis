#   Copyright ETH 2023 Zürich, Scientific IT Services
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
#
import sys
import os

import imaging as imaging

import math
from pybis import Openbis
import numpy as np

from spmpy_terry import spm
import spmpy_terry as spmpy
from datetime import datetime

SXM_ADAPTOR = "ch.ethz.sis.openbis.generic.server.dss.plugins.imaging.adaptor.NanonisSxmAdaptor"
DAT_ADAPTOR = "ch.ethz.sis.openbis.generic.server.dss.plugins.imaging.adaptor.NanonisDatAdaptor"
VERBOSE = False
# DEFAULT_URL = "http://localhost:8080/openbis"
# DEFAULT_URL = "http://localhost:8888/openbis"
DEFAULT_URL = "https://local.openbis.ch/openbis"


def get_instance(url=None):
    if url is None:
        url = DEFAULT_URL
    openbis_instance = Openbis(
        url=url,
        verify_certificates=False,
        allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks=True
    )
    # token = openbis_instance.login('admin', 'changeit')
    token = openbis_instance.login('admin', '123456789')
    print(f'Connected to {url} -> token: {token}')
    return openbis_instance


def get_color_scale_range(img, channel):
    minimum = np.nanmin(img.get_channel(channel)[0])
    maximum = np.nanmax(img.get_channel(channel)[0])

    step = abs(round((maximum - minimum) / 100, 2))
    if step >= 1:
        step = 1
    elif step > 0:
        step = 0.01
    else:
        step = abs((maximum - minimum) / 100)
        step = math.log10(step)
        if math.isnan(step) or math.isinf(step):
            step = 0.01
        else:
            step = 10 ** math.floor(step)

    return [str(minimum), str(maximum), str(step)]


def create_sxm_dataset(openbis, experiment, file_path, sample=None):
    img = spm(file_path)
    channels = [x['ChannelNickname'] for x in img.SignalsList]

    imaging_control = imaging.ImagingControl(openbis)

    color_scale_visibility = [
        imaging.ImagingDataSetControlVisibility(
            "Channel",
            [channel],
            get_color_scale_range(img, channel),
            img.get_channel(channel)[1])
        for channel in channels]

    exports = [imaging.ImagingDataSetControl('include', "Dropdown", values=['image', 'raw data'], multiselect=True),
               imaging.ImagingDataSetControl('image-format', "Dropdown", values=['png', 'svg']),
               imaging.ImagingDataSetControl('archive-format', "Dropdown", values=['zip', 'tar']),
               imaging.ImagingDataSetControl('resolution', "Dropdown", values=['original', '150dpi', '300dpi'])]

    inputs = [
        imaging.ImagingDataSetControl('Channel', "Dropdown", values=channels),
        imaging.ImagingDataSetControl('X-axis', "Range", values_range=["0", str(img.get_param('width')[0]), "0.01"]),
        imaging.ImagingDataSetControl('Y-axis', "Range", values_range=["0", str(img.get_param('height')[0]), "0.01"]),
        imaging.ImagingDataSetControl('Color-scale', "Range", visibility=color_scale_visibility),
        imaging.ImagingDataSetControl('Colormap', "Colormap", values=['gray', 'YlOrBr', 'viridis', 'cividis', 'inferno', 'rainbow', 'Spectral', 'RdBu', 'RdGy']),
        imaging.ImagingDataSetControl('Scaling', "Dropdown", values=['linear', 'logarithmic']),
    ]

    imaging_config = imaging.ImagingDataSetConfig(
        SXM_ADAPTOR,
        1.0,
        ['original', '200x200', '2000x2000'],
        True,
        [1000, 2000, 5000],
        exports,
        inputs,
        {})

    images = [imaging.ImagingDataSetImage(previews=[imaging.ImagingDataSetPreview(preview_format="png")])]
    imaging_property_config = imaging.ImagingDataSetPropertyConfig(imaging_config, images)
    if VERBOSE:
        print(imaging_property_config.to_json())

    return imaging_control.create_imaging_dataset(
        dataset_type="IMAGING_DATA",
        config=imaging_property_config,
        experiment=experiment,
        sample=sample,
        files=[file_path])


def create_dat_dataset(openbis, folder_path, file_prefix='', sample=None, experiment=None):
    assert experiment is not None or sample is not None, "Either sample or experiment needs to be provided!"
    data = spmpy.importall(folder_path, file_prefix, 'spec')
    imaging_control = imaging.ImagingControl(openbis)

    for d in data:
        if d.type == 'scan':
            date = d.get_param('rec_date')
            time = d.get_param('rec_time')
            date_time = '%s %s' % (date, time)
            d.date_time = datetime.strptime(date_time, "%d.%m.%Y %H:%M:%S")

        if d.type == 'spec':
            date_time = d.get_param('Saved Date')
            d.date_time = datetime.strptime(date_time, "%d.%m.%Y %H:%M:%S") if date_time is not None else datetime.now()

    data.sort(key=lambda da: da.date_time)

    channels = list(set([(channel['ChannelNickname'], channel['ChannelUnit'], channel['ChannelScaling']) for spec in data for channel in spec.SignalsList]))
    print(channels)
    color_scale_visibility_x = []
    color_scale_visibility_y = []
    for (channel, unit, scaling) in channels:
        minimum = []
        maximum = []
        for spec in data:
            channel_in_signals_list = False
            for signal_settings in spec.SignalsList:
                if channel in signal_settings["ChannelNickname"]:
                    channel_in_signals_list = True
            if channel_in_signals_list:
                minimum += [np.nanmin(spec.get_channel(f'{channel}')[0])]
                maximum += [np.nanmax(spec.get_channel(f'{channel}')[0])]
        minimum = np.nanmin(minimum)
        maximum = np.nanmax(maximum)
        step = abs(round((maximum - minimum) / 100, 2))

        if step >= 1:
            step = 1
        elif step > 0:
            step = 0.01
        else:
            step = abs((maximum - minimum) / 100)
            step = math.log10(step)
            if math.isnan(step) or math.isinf(step):
                step = 0.01
            else:
                step = 10 ** math.floor(step)

        color_scale_visibility_x += [imaging.ImagingDataSetControlVisibility(
            "Channel X",
            [channel],
            [str(minimum), str(maximum), str(step)],
            unit
        )]

        color_scale_visibility_y += [imaging.ImagingDataSetControlVisibility(
            "Channel Y",
            [channel],
            [str(minimum), str(maximum), str(step)],
            unit
        )]

    exports = [imaging.ImagingDataSetControl('include', "Dropdown", values=['image', 'raw data'], multiselect=True),
               imaging.ImagingDataSetControl('image-format', "Dropdown", values=['png', 'svg']),
               imaging.ImagingDataSetControl('archive-format', "Dropdown", values=['zip', 'tar']),
               imaging.ImagingDataSetControl('resolution', "Dropdown", values=['original', '150dpi', '300dpi'])]

    inputs = [
        imaging.ImagingDataSetControl('Channel X', "Dropdown", values=[channel[0] for channel in channels]),
        imaging.ImagingDataSetControl('Channel Y', "Dropdown", values=[channel[0] for channel in channels]),
        imaging.ImagingDataSetControl('X-axis', "Range", visibility=color_scale_visibility_x),
        imaging.ImagingDataSetControl('Y-axis', "Range", visibility=color_scale_visibility_y),
        imaging.ImagingDataSetControl('Grouping', "Dropdown", values=[d.name for d in data], multiselect=True),
        imaging.ImagingDataSetControl('Colormap', "Colormap", values=['gray', 'YlOrBr', 'viridis', 'cividis', 'inferno', 'rainbow', 'Spectral', 'RdBu', 'RdGy']),
        imaging.ImagingDataSetControl('Scaling', "Dropdown", values=['lin-lin', 'lin-log', 'log-lin', 'log-log']),
    ]

    imaging_config = imaging.ImagingDataSetConfig(
        DAT_ADAPTOR,
        1.0,
        ['original', '200x200', '2000x2000'],
        True,
        [1000, 2000, 5000],
        exports,
        inputs,
        {})

    images = [imaging.ImagingDataSetImage()]
    imaging_property_config = imaging.ImagingDataSetPropertyConfig(imaging_config, images)
    if VERBOSE:
        print(imaging_property_config.to_json())

    return imaging_control.create_imaging_dataset(
        dataset_type="IMAGING_DATA",
        config=imaging_property_config,
        experiment=experiment,
        sample=sample,
        files=[d.path for d in data])


def create_preview(openbis, perm_id, config, preview_format="png", image_index=0):
    imaging_control = imaging.ImagingControl(openbis)
    preview = imaging.ImagingDataSetPreview(preview_format, config)
    preview = imaging_control.make_preview(perm_id, image_index, preview)
    return preview


def update_image_with_preview(openbis, perm_id, image_id, preview: imaging.ImagingDataSetPreview):
    imaging_control = imaging.ImagingControl(openbis)
    config = imaging_control.get_property_config(perm_id)
    image = config.images[image_id]
    if len(image.previews) > preview.index:
        image.previews[preview.index] = preview
    else:
        preview.index = len(image.previews)
        image.add_preview(preview)

    imaging_control.update_property_config(perm_id, config)


def export_image(openbis: Openbis, perm_id: str, image_id: int, path_to_download: str,
                 include=None, image_format='original', archive_format="zip", resolution='original'):
    if include is None:
        include = ['image', 'raw data']
    imaging_control = imaging.ImagingControl(openbis)
    export_config = {
        "include": include,
        "image-format": image_format,
        "archive-format": archive_format,
        "resolution": resolution
    }
    imaging_export = imaging.ImagingDataSetExport(export_config)
    imaging_control.single_export_download(perm_id, imaging_export, image_id, path_to_download)


def multi_export_images(openbis: Openbis, perm_ids: list[str], image_ids: list[int], preview_ids: list[int],
                        path_to_download: str, include=None, image_format='original',
                        archive_format="zip", resolution='original'):
    if include is None:
        include = ['image', 'raw data']
    imaging_control = imaging.ImagingControl(openbis)
    export_config = {
        "include": include,
        "image-format": image_format,
        "archive-format": archive_format,
        "resolution": resolution
    }
    imaging_multi_exports = []
    for i in range(len(perm_ids)):
        imaging_multi_exports += [imaging.ImagingDataSetMultiExport(perm_ids[i], image_ids[i],
                                                                   preview_ids[i], export_config)]
    imaging_control.multi_export_download(imaging_multi_exports, path_to_download)


def demo_sxm_flow(openbis, file_sxm, permId=None):

    perm_id = permId
    if perm_id is None:
        dataset_sxm = create_sxm_dataset(
            openbis=openbis,
            experiment='/CARBON_NANOMATERIALS/TRIANGULENE_SPIN_CHAINS/TRIANGULENE_SPIN_CHAINS_EXP_3',
            sample='/CARBON_NANOMATERIALS/TRIANGULENE_SPIN_CHAINS/TRIANGULENE_SPIN_CHAINS_EXP_3/TEMPLATE-SXM',
            file_path=file_sxm)
        perm_id = dataset_sxm.permId
        print(f'Created imaging .SXM dataset: {dataset_sxm.permId}')



def demo_dat_flow(openbis, folder_path, permId=None):

    print(f'Searching for .DAT files')
    perm_id = permId
    if permId is None:
        dataset_dat = create_dat_dataset(
            openbis=openbis,
            experiment='/CARBON_NANOMATERIALS/TRIANGULENE_SPIN_CHAINS/TRIANGULENE_SPIN_CHAINS_EXP_3',
            sample='/CARBON_NANOMATERIALS/TRIANGULENE_SPIN_CHAINS/TRIANGULENE_SPIN_CHAINS_EXP_3/TEMPLATE-DAT',
            folder_path=folder_path,
            file_prefix='')
        perm_id = dataset_dat.permId
        print(f'Created imaging .DAT dataset: {dataset_dat.permId}')



openbis_url = None
data_folder = 'data'

if len(sys.argv) > 2:
    openbis_url = sys.argv[1]
    data_folder = sys.argv[2]
else:
    print(f'Usage: python3 nanonis_importer.py <OPENBIS_URL> <PATH_TO_DATA_FOLDER>')
    print(f'Using default parameters')
    print(f'URL: {DEFAULT_URL}')
    print(f'Data folder: {data_folder}')

o = get_instance(openbis_url)

# sxm_files = [f for f in os.listdir(data_folder) if f.endswith('.sxm')]
# print(f'Found {len(sxm_files)} Nanonis .SXM files in {data_folder}')

# for sxm_file in sxm_files:
#     print(f"SXM file: {sxm_file}")
#     file_path = os.path.join(data_folder, sxm_file)
#     try:
#         demo_sxm_flow(o, file_path)
#     except:
#         print(f"Cannot upload {sxm_file}.")

demo_dat_flow(o, data_folder)

o.logout()
