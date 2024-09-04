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
import copy
import argparse

import imaging as imaging

import math
from pybis import Openbis
import numpy as np

from spmpy_terry import spm
import spmpy_terry as spmpy
from datetime import datetime
import json
import shutil
from collections import defaultdict

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
        step = np.log10(step)
        if np.isnan(step) or np.isinf(step):
            step = 0.01
        else:
            step = 10 ** np.floor(step)

    return [str(minimum), str(maximum), str(step)]

def reorder_sxm_channels(channels, header):
    """
    Lock-in>Lock-in status: ON > dIdV vs V
    Lock-in>Lock-in status: OFF:
        Z-Ctrl hold: TRUE > z vs V
        Z-Ctrl hold: FALSE:
            Oscillation Control>output off: TRUE > df vs V
            Oscillation Control>output off: FALSE > I vs V
    """
    channel_index = -1
    if header["lock-in>lock-in status"] == "ON":
        channel_index = channels.index("dIdV")
    else:
        if header["z-controller>controller status"] == "ON":
            channel_index = channels.index("z")
        else:
            if header["oscillation control>output off"] == "TRUE":
                channel_index = channels.index("df")
            else:
                channel_index = channels.index("I")
    
    # If the channel index is less than 0, it means the tree did not find the measurement type
    if channel_index >= 0:        
        channels[channel_index], channels[0] = channels[0], channels[channel_index]
        
    return channels

def reorder_dat_channels(channels, header):
    """
    Bias spectroscopy:
        Lock-in>Lock-in status: ON > dIdV vs V
        Lock-in>Lock-in status: OFF
            Z-Ctrl hold: TRUE > z vs V
            Z-Ctrl hold: FALSE
                Oscillation Control>output off: TRUE > df vs V
                Oscillation Control>output off: FALSE > I vs V

    Z spectroscopy:
        Lock-in>Lock-in status: ON > dIdV vs z
        Lock-in>Lock-in status: OFF
            Oscillation Control>output off: TRUE > df vs z
            Oscillation Control>output off: FALSE > I vs z
    """
    channels_x = copy.deepcopy(channels)
    channels_y = copy.deepcopy(channels)
    channel_x_index = -1
    channel_y_index = -1
    
    if header["Experiment"] == "bias spectroscopy":
        if header["Lock-in>Lock-in status"] == "ON":
            channel_x_index = channels_x.index(("V","V",1))
            channel_y_index = channels_y.index(("dIdV","pA",10**12))
        else:
            if header["Z-Ctrl hold"] == "FALSE":
                channel_x_index = channels_x.index(("V","V",1))
                channel_y_index = channels_y.index(("zspec","nm",10**9))
            else:
                if header["Oscillation Control>output off"] == "TRUE":
                    channel_x_index = channels_x.index(("V","V",1))
                    channel_y_index = channels_y.index(("df","Hz",1))
                else:
                    channel_x_index = channels_x.index(("V","V",1))
                    channel_y_index = channels_y.index(("I","pA",10**12))
    else:
        if header["Lock-in>Lock-in status"] == "ON":
            channel_x_index = channels_x.index(("zspec","nm",10**9))
            channel_y_index = channels_y.index(("dIdV","pA",10**12))
        else:
            if header["Oscillation Control>output off"] == "TRUE":
                channel_x_index = channels_x.index(("zspec","nm",10**9))
                channel_y_index = channels_y.index(("df","Hz",1))
            else:
                channel_x_index = channels_x.index(("zspec","nm",10**9))
                channel_y_index = channels_y.index(("I","pA",10**12))

    if channel_x_index >= 0:
        channels_x[channel_x_index], channels_x[0] = channels_x[0], channels_x[channel_x_index]
        channels_y[channel_y_index], channels_y[0] = channels_y[0], channels_y[channel_y_index]
    return channels_x, channels_y

def get_dat_type(header):
    """
    Bias spectroscopy:
        Lock-in>Lock-in status: ON > dIdV vs V
        Lock-in>Lock-in status: OFF
            Z-Ctrl hold: TRUE > z vs V
            Z-Ctrl hold: FALSE
                Oscillation Control>output off: TRUE > df vs V
                Oscillation Control>output off: FALSE > I vs V

    Z spectroscopy:
        Lock-in>Lock-in status: ON > dIdV vs z
        Lock-in>Lock-in status: OFF
            Oscillation Control>output off: TRUE > df vs z
            Oscillation Control>output off: FALSE > I vs z
    """
    measurement_type = ""
    
    if header["Experiment"] == "bias spectroscopy":
        if header["Lock-in>Lock-in status"] == "ON":
            measurement_type = "bias spectroscopy dIdV vs V"
        else:
            if header["Z-Ctrl hold"] == "TRUE":
                measurement_type = "bias spectroscopy z vs V"
            else:
                if header["Oscillation Control>output off"] == "TRUE":
                    measurement_type = "bias spectroscopy df vs V"
                else:
                    measurement_type = "bias spectroscopy I vs V"
    else:
        if header["Lock-in>Lock-in status"] == "ON":
            measurement_type = "z spectroscopy dIdV vs z"
        else:
            if header["Oscillation Control>output off"] == "TRUE":
                measurement_type = "z spectroscopy df vs z"
            else:
                measurement_type = "z spectroscopy I vs z"

    return measurement_type

def create_sxm_dataset(openbis, experiment, file_path, sample=None):
    img = spm(file_path)
    channels = [x['ChannelNickname'] for x in img.SignalsList]
    header = img.header
    
    # Select default channel according to the measurement type
    channels = reorder_sxm_channels(channels, header)

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
    channels_x, channels_y = reorder_dat_channels(channels, data[0].header) # All files inside data belong to the same measurement type. Thus, the header of the first file can be used for all of them.
    
    color_scale_visibility_x = []
    color_scale_visibility_y = []
    for idx, (channel_x, unit_x, scaling_x) in enumerate(channels_x):
        channel_y = channels_y[idx][0]
        unit_y = channels_y[idx][1]
        scaling_y = channels_y[idx][2]
        
        minimum_x, maximum_x = [], []
        minimum_y, maximum_y = [], []
        for spec in data:
            # # -------- Boolean flag was added to the code -------
            # channel_in_signals_list = False
            # for signal_settings in spec.SignalsList:
            #     if channel_x in signal_settings["ChannelNickname"]:
            #         channel_in_signals_list = True
            # if channel_in_signals_list:
            # # ---------------------------------------------------
            minimum_x += [np.nanmin(spec.get_channel(f'{channel_x}')[0])]
            maximum_x += [np.nanmax(spec.get_channel(f'{channel_x}')[0])]
            
            minimum_y += [np.nanmin(spec.get_channel(f'{channel_y}')[0])]
            maximum_y += [np.nanmax(spec.get_channel(f'{channel_y}')[0])]
        minimum_x = np.nanmin(minimum_x)
        maximum_x = np.nanmax(maximum_x)
        
        minimum_y = np.nanmin(minimum_y)
        maximum_y = np.nanmax(maximum_y)
        step_x = abs(round((maximum_x - minimum_x) / 100, 2))
        step_y = abs(round((maximum_y - minimum_y) / 100, 2))

        if step_x >= 1:
            step_x = 1
        elif step_x > 0:
            step_x = 0.01
        else:
            step_x = abs((maximum_x - minimum_x) / 100)
            step_x = np.log10(step_x)
            if np.isnan(step_x) or np.isinf(step_x):
                step_x = 0.01
            else:
                step_x = 10 ** np.floor(step_x)
                
        if step_y >= 1:
            step_y = 1
        elif step_y > 0:
            step_y = 0.01
        else:
            step_y = abs((maximum_y - minimum_y) / 100)
            step_y = np.log10(step_y)
            if np.isnan(step_y) or np.isinf(step_y):
                step_y = 0.01
            else:
                step_y = 10 ** np.floor(step_y)

        color_scale_visibility_x += [imaging.ImagingDataSetControlVisibility(
            "Channel X",
            [channel_x],
            [str(minimum_x), str(maximum_x), str(step_x)],
            unit_x
        )]

        color_scale_visibility_y += [imaging.ImagingDataSetControlVisibility(
            "Channel Y",
            [channel_y],
            [str(minimum_y), str(maximum_y), str(step_y)],
            unit_y
        )]

    exports = [imaging.ImagingDataSetControl('include', "Dropdown", values=['image', 'raw data'], multiselect=True),
               imaging.ImagingDataSetControl('image-format', "Dropdown", values=['png', 'svg']),
               imaging.ImagingDataSetControl('archive-format', "Dropdown", values=['zip', 'tar']),
               imaging.ImagingDataSetControl('resolution', "Dropdown", values=['original', '150dpi', '300dpi'])]
    
    inputs = [
        imaging.ImagingDataSetControl('Channel X', "Dropdown", values=[channel[0] for channel in channels_x]),
        imaging.ImagingDataSetControl('Channel Y', "Dropdown", values=[channel[0] for channel in channels_y]),
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

def demo_sxm_flow(openbis, file_sxm, experiment_permid, sample_permid, permId=None):

    perm_id = permId
    if perm_id is None:
        dataset_sxm = create_sxm_dataset(
            openbis=openbis,
            experiment=experiment_permid,
            sample=sample_permid,
            file_path=file_sxm)
        perm_id = dataset_sxm.permId
        print(f'Created imaging .SXM dataset: {dataset_sxm.permId}')

def demo_dat_flow(openbis, folder_path, experiment_permid, sample_permid, permId=None):

    print(f'Searching for .DAT files')
    perm_id = permId
    if permId is None:
        dataset_dat = create_dat_dataset(
            openbis=openbis,
            experiment=experiment_permid,
            sample=sample_permid,
            folder_path=folder_path,
            file_prefix='')
        perm_id = dataset_dat.permId
        print(f'Created imaging .DAT dataset: {dataset_dat.permId}')

def get_1D_measurement_props(full_dat_filepath):
    pass
    
def get_2D_measurement_props(full_sxm_filepath):
    img = spm(full_sxm_filepath)
    properties = {
        "$name": full_sxm_filepath.split("/")[-1],
        "start_time": datetime.strptime(f"{img.header['rec_date']} {img.header['rec_time']}", "%d.%m.%Y %H:%M:%S").strftime("%Y-%m-%d %H:%M:%S"),
        "duration": json.dumps({"has_value": float(img.header["acq_time"]), "has_unit": "http://qudt.org/vocab/unit/SEC"}),
        "bias_setpoint": json.dumps({"has_value": float(img.header["bias>bias (v)"]), "has_unit": "http://qudt.org/vocab/unit/V"}),
        "bias_calibration_factor": json.dumps({"has_value": float(img.header["bias>calibration (v/v)"]), "has_unit": "http://qudt.org/vocab/unit/V-PER-V"}),
        "bias_calibration_offset": json.dumps({"has_value": float(img.header["bias>offset (v)"]), "has_unit": "http://qudt.org/vocab/unit/V"}),
        "current_setpoint": json.dumps({"has_value": float(img.header["current>current (a)"]), "has_unit": "http://qudt.org/vocab/unit/A"}),
        "current_calibration_factor": json.dumps({"has_value": float(img.header["current>calibration (a/v)"]), "has_unit": "A/V"}),
        "current_calibration_offset": json.dumps({"has_value": float(img.header["current>offset (a)"]), "has_unit": "http://qudt.org/vocab/unit/A"}),
        "current_gain": json.dumps({"has_value": img.header["current>gain"], "has_unit": "None"})
    }
    return properties

parser = argparse.ArgumentParser(description = 'Upload Nanonis files into openBIS.')

# Define the arguments with flags
parser.add_argument('-o', '--openbis_url', type=str, help='OpenBIS URL', default = None)
parser.add_argument('-d', '--data_folder', type=str, help='Path to the data folder', default = 'data')
parser.add_argument('-c', '--collection_permid', type=str, help='Collection ID', default = '20240805121017676-1377')
parser.add_argument('-s', '--sample_permid', type=str, help='Sample ID', default = '20240809165355701-1394')

args = parser.parse_args()

openbis_url = args.openbis_url
data_folder = args.data_folder
collection_permid = args.collection_permid
sample_permid = args.sample_permid

if not openbis_url:
    print(f'Usage: python3 nanonis_importer.py -o <OPENBIS_URL> -d <PATH_TO_DATA_FOLDER> -c <COLLECTION_ID> -s <SAMPLE_ID>')
    print(f'Using default parameters')
    print(f'URL: {DEFAULT_URL}')
    print(f'Data folder: {data_folder}')
    print(f'Default collection ID: {collection_permid}')
    print(f'Default sample ID: {sample_permid}')

o = get_instance(openbis_url)

measurement_files = [f for f in os.listdir(data_folder)]
production = True

if production:
    measurement_datetimes = []

    for f in measurement_files:
        if f.endswith(".sxm"):
            img = spm(f"{data_folder}/{f}")
            img_datetime = datetime.strptime(f"{img.header['rec_date']} {img.header['rec_time']}", "%d.%m.%Y %H:%M:%S")
            measurement_datetimes.append(img_datetime)
            
        elif f.endswith(".dat"):
            img = spm(f"{data_folder}/{f}")
            img_datetime = datetime.strptime(img.header['Saved Date'], "%d.%m.%Y %H:%M:%S")
            measurement_datetimes.append(img_datetime)

    # Sort files by datetime
    sorted_measurement_files = [x for _, x in sorted(zip(measurement_datetimes, measurement_files))]

    # Dat files belonging to the same measurement session, i.e., that are consecutive, should be grouped into just one list of files.
    grouped_measurement_files = []
    group = []
    for i,f in enumerate(sorted_measurement_files):
        if f.endswith(".sxm"):
            if len(group) > 0:
                grouped_measurement_files.append(group)
                group = []
            grouped_measurement_files.append([f])
        elif f.endswith(".dat"):
            group.append(f)

    if len(group) > 0:
        grouped_measurement_files.append(group)
        group = []

    if collection_permid:
        for group in grouped_measurement_files:
            if group[0].endswith(".sxm"):
                print(f"SXM file: {group[0]}")
                twoD_measurement_sample = o.new_sample(
                    type = "2D_MEASUREMENT", 
                    experiment = collection_permid, 
                    props = get_2D_measurement_props(os.path.join(data_folder, group[0])),
                    parents = [sample_permid]
                )
                twoD_measurement_sample.save()
                file_path = os.path.join(data_folder, group[0])
                try:
                    demo_sxm_flow(o, file_path, collection_permid, twoD_measurement_sample.permId)
                except:
                    print(f"Cannot upload {group[0]}.")
            else:
                
                # Split the dat files by measurement type (e.g.: bias spec dI vs V in one list, bias spec z vs V in another list, etc.)
                dat_files_types = []
                for dat_file in group:
                    dat_data = spm(f"{data_folder}/{dat_file}")
                    dat_files_types.append(get_dat_type(dat_data.header))
                
                grouped = defaultdict(list)

                for item1, item2 in zip(group, dat_files_types):
                    grouped[item2].append(item1)

                dat_files_grouped_by_type = list(grouped.values())
                # ---------------
                
                for dat_files_group in dat_files_grouped_by_type:
                    oneD_measurement_sample = o.new_sample(
                        type = "1D_MEASUREMENT", 
                        experiment = collection_permid,
                        props = {"$name": "Experimental 1D Measurement"},
                        parents = [sample_permid]
                    )
                    oneD_measurement_sample.save()
                    
                    dat_files_directory = os.path.join(data_folder, "dat_files")
                    os.mkdir(dat_files_directory)
                    
                    for dat_file in dat_files_group:
                        shutil.copy(os.path.join(data_folder, dat_file), os.path.join(dat_files_directory, dat_file))
                    
                    demo_dat_flow(o, dat_files_directory, collection_permid, oneD_measurement_sample.permId)
                    
                    shutil.rmtree(dat_files_directory)
                    
    o.logout()

# If the openBIS microscopy app does not work as expected (you are not able to change plots inside openBIS),
# open openBIS app container (image: openbis/openbis-server:alpha) and install python together with some libraries.
# For the installation follow the following steps:
# 1. Open terminal of openBIS container
# 2. Run:
#    2.1. apt-get update
#    2.2. apt install software-properties-common -y
#    2.3. apt install python3-pip -y
#    2.4. pip3 install -r /etc/openbis/core-plugins/imaging-nanonis/1/scripts/python_requirements.txt
