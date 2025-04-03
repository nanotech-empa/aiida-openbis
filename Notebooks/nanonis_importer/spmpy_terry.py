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
#     # -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 15:13:13 2017

@author: Bruno Schuler
"""

##############################################
# spm_matrix class
##############################################

# This module provides functions to load and display Nanonis files

# requires packages: 
# pip install numpy
# pip install six 
# pip install nanonispy
# pip install spiepy

import os as os
import numpy as np
import nanonispy as nap
import spiepy
import warnings
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

class spm:

    #Dictionary Channels
    ChannelName = ['LI_Demod_1_X','LI_Demod_1_Y','Z','Current','Bias','Frequency_Shift','Amplitude','Excitation','Temperature_1',
                   'Bias (V)','Bias calc (V)', 'Bias [bwd] (V)', 'Current (A)','Current [bwd] (A)','Amplitude (m)',
                   'Amplitude [bwd] (m)', 'Excitation (V)', 'Excitation [bwd] (V)', 'Frequency Shift (Hz)', 'Frequency Shift [bwd] (Hz)',
                   'LI Demod 1 X (A)','LI Demod 1 X (A) [bwd] (A)','PMT (V)','Counter 1 (Hz)','Counter_1', 'Z rel (m)', 'Z (m)','Time (s)',
                   'Delay Sampling (s)','LI Demod 0 X (V)','LI Demod 0 Y (V)','LI Demod 3 X (A)','LI Demod 3 Y (A)',
                   'Position Phase1 (m)','Rotation1 (deg)','Rotation2 (deg)','Rotation (deg)','Index']
    ChannelNickname = ['dIdV','dIdV_Y','z','I','V','df','A','exc','T1',
                    'V', 'V', 'V_bw' ,'I','I_bw','A',
                    'A_bw', 'exc','exc_bw','df','df_bw',
                    'dIdV','dIdV_bw','PMT','counter','counter', 'zrel','zspec','t',
                    'Delay','EOS','EOS_Y','I_THz','I_THz_Y',
                    'Phase','Rot1','Rot2','Rot','Index']
    ChanneliScaling = [10**12,10**12,10**9,10**12,1,1,10**9,1,1,
                    1,1,1,10**12,10**12,10**9,
                    10**9,1,1,1,1,10**12,
                    1,1,1,1,10**12,10**9,1,
                    10**12,1,1,10**12,10**12,
                    10**3,1,1,1,1]
    ChannelUnit = ['pS','pS','nm','pA','V','Hz','nm','V','K',
                    'V','V','V','pA','pA','nm',
                    'nm','V','V','Hz','Hz','pA',
                    'a.u.','V','Hz','Hz','pm','nm','s',
                    'ps','V','V','pA','pA',
                    'mm','deg','deg','deg','N']

    global SignalsListReference
    SignalsListReference = []

    for (chName,chNickname,chScaling,chUnit) in zip(ChannelName,ChannelNickname,ChanneliScaling,ChannelUnit):
        SignalsListReference.append({'ChannelName': chName, 'ChannelNickname': chNickname , 'ChannelScaling': chScaling, 'ChannelUnit': chUnit})

    del ChannelName,ChanneliScaling,ChannelUnit    
    del chName,chNickname,chScaling,chUnit     

    #Dictionary Parameters
    ParamName = ['X (m)','Y (m)','Z (m)','bias','z-controller>setpoint','scan_angle','Comments','Comment01',
                'Lock-in>Amplitude','Lock-in>Reference phase D1 (deg)','Lock-in>Frequency (Hz)','Temperature 2>Temperature 2 (K)',
                'Current>Current (A)','Bias>Bias (V)','z-controller>tiplift (m)',
                'Parameter>Delay Sampling (s)','Parameter>Delay PP1 (s)','Parameter>Delay PP2 (s)','Parameter>Angle Rot1 (deg)','Parameter>Angle Rot2 (deg)','Parameter>Position Phase1 (m)','Parameter>Position_Phase1 (m)']
    ParamNickname = ['x','y','z','V','setpoint','angle','comment','comment_spec',
                'lockin_amplitude','lockin_phase','lockin_frequency','temperature',
                'setpoint_spec','V_spec','z_offset',
                'Sampling','PP1','PP2','Rot1','Rot2','Phase','Phase_depricated']
    ParamScaling = [10**9,10**9,10**9,1,10**12,1,'na','na',
                10**3,1,1,1,
                10**12,1,10**9,
                10**12,10**12,10**12,1,1,10**3,10**3] # na if not a numeric value
    ParamUnit = ['nm','nm','nm','V','pA','°','','',
                'mV','°','Hz','K',
                'pA','V','nm',
                'ps','ps','ps','deg','deg','mm','mm']


    global ParamListReference
    ParamListReference = []
    
    for (paName,paNickname,paScaling,paUnit) in zip(ParamName,ParamNickname,ParamScaling,ParamUnit):
        ParamListReference.append({'ParamName': paName, 'ParamNickname': paNickname, 'ParamScaling': paScaling, 'ParamUnit': paUnit})

    del ParamName,ParamScaling,ParamUnit    
    del paName,paNickname,paScaling,paUnit

    # constructor
    def __init__(self,path):
        #import os as os
        #import numpy as np
        #import nanonispy as nap
              
      
        # self.path = path.replace('//', '/')
        abspath = os.path.abspath(path)
        self.path = '/'.join(abspath.split('\\')[-4:])
        self.name = self.path.split('/')[-1]
        file_extension = os.path.splitext(path)[1]
        
        if file_extension == '.sxm':
            self.napImport = nap.read.Scan(path)
            self.type = 'scan'
        elif file_extension == '.dat':
            self.napImport = nap.read.Spec(path)
            self.type = 'spec'
        else:
            print('Datatype not supported.')
            return;
            
        ch = []    
        for key in self.napImport.signals:
            try:
                ch.append([d['ChannelName'] for d in SignalsListReference].index(key))
            except:
                pass;    
 
        self.SignalsList = [SignalsListReference[i] for i in ch] #List of all recorded channels
        self.channels = [c['ChannelNickname'] for c in self.SignalsList] 
        self.header = self.napImport.header
        
        
    def __repr__(self):
        return self.path
        
    #get channel
    def get_channel(self,channel,direction = 'forward', flatten = False, offset = False,zero = False):
        
        #import spiepy
        #import numpy as np
        
        if self.type == 'scan':
            chNum = [d['ChannelNickname'] for d in self.SignalsList].index(channel)
            im = self.napImport.signals[self.SignalsList[chNum]['ChannelName']][direction]
            im = im *self.SignalsList[chNum]['ChannelScaling']
            
            if flatten:
                if ~np.isnan(np.sum(im)):
                    im, _ = spiepy.flatten_xy(im)
                else:
                    m,n = np.shape(im)
                    i = np.argwhere(np.isnan(im))[0,0]
                    im_cut = im[:i-1,:]
                    im, _ = spiepy.flatten_xy(im_cut)
                    empty = np.full((m-i,n),np.nan)
                    im = np.vstack((im,empty))
                    
                    warnings.warn('NaN values in image')
                    #im-np.nanmean(im)
               
            if offset:
                im = im-np.mean(im)
            if zero:
                im = im+abs(np.min(im))

            unit = self.SignalsList[chNum]['ChannelUnit']
                
            return (im,unit)
        
        elif self.type == 'spec':
        
            if direction == 'backward':
                channel = channel + '_bw';
                #print(channel)
            
            chNum = [d['ChannelNickname'] for d in self.SignalsList].index(channel)
            data =  self.napImport.signals[self.SignalsList[chNum]['ChannelName']]
            data = data*self.SignalsList[chNum]['ChannelScaling']
            unit = self.SignalsList[chNum]['ChannelUnit']
            
            return (data,unit)
            
        return
    
    #get parameter            
    def get_param(self,param):
 
        if any(d['ParamNickname'] == param for d in ParamListReference):
            paNum = [d['ParamNickname'] for d in ParamListReference].index(param)

            if ParamListReference[paNum]['ParamScaling'] == 'na':
                return self.napImport.header[ParamListReference[paNum]['ParamName']]
            else:
                return (float(self.napImport.header[ParamListReference[paNum]['ParamName']])*ParamListReference[paNum]['ParamScaling'],ParamListReference[paNum]['ParamUnit'])
            
        
        elif param == 'width' or param == 'height':
            # height, width = self.get_param('scan_range')
            # height, width = [height*10**9, width*10**9]
            scanfield = self.get_param('scan>scanfield')
            width, height = float(scanfield.split(';')[2])*10**9, float(scanfield.split(';')[3])*10**9
             
            return (eval(param),'nm')
        
            
        else:
            if param in self.napImport.header.keys():
                return self.napImport.header[param]
            else:
                return        
        
              
        
    # print essential parameters for plotting  
    def print_params(self, show = True):
    
        # import numpy as np
        
        label = []
        
        if self.type == 'scan':
            fb_enable = self.get_param('z-controller>controller status')
            fb_ctrl = self.get_param('z-controller>controller name')
            bias = self.get_param('V')
            set_point = self.get_param('setpoint')
            height = self.get_param('height')
            width = self.get_param('width')
            angle = self.get_param('angle')
            z_offset = self.get_param('z_offset')
            comment = self.get_param('comments')
            
            
                         
            if fb_enable == 'OFF':
                label.append('constant height')
                label.append('z-offset: %.3f%s' % z_offset)
                
            if np.abs(bias[0])<0.1:
                bias = list(bias)
                bias[0] = bias[0]*1000
                bias[1] = 'mV'
                bias = tuple(bias)
                
            label.append('I: %.0f%s' % set_point)
            label.append('bias: %.2f%s' % bias)
            label.append('size: %.1f%s x %.1f%s (%.0f%s)' % (width+height+angle))
            label.append('comment: %s' % comment)
            
            
        elif self.type == 'spec':
            
            fb_enable = self.get_param('Z-Ctrl hold')
            set_point = self.get_param('setpoint_spec')
            bias = self.get_param('V_spec')
            #lockin_status = self.get_param('Lock-in>Lock-in status')
            lockin_amplitude = self.get_param('lockin_amplitude')
            lockin_phase= self.get_param('lockin_phase')
            lockin_frequency= self.get_param('lockin_frequency')
            comment = self.get_param('comment_spec')
            
                               
            #if lockin_status == 'ON':
            label.append('lockin: A = %.0f%s (θ = %.0f%s, f = %.0f%s)' % (lockin_amplitude+lockin_phase+lockin_frequency))
                 
            
            if fb_enable == 'FALSE':
                label.append('feedback: on')
                
            elif fb_enable == 'TRUE':
                label.append('feedback: off')
           
 
            label.append('setpoint: I = %.0f%s, V = %.1f%s' % (set_point+bias))    
            
            label.append('comment: %s' % comment)
    
        # label.append('path: %s' % self.path)
        label = '\n'.join(label)
        
        if show:
            print(label)
        
        return label
         
      
    
 # plot
    def plot(self, **params):
        import numpy as np
        #import matplotlib as ml
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        
        #cmaps = sorted(m for m in plt.cm.datad)
        # ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'jet', 'jet_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spectral', 'spectral_r', 'spring', 'spring_r', 'summer', 'summer_r', 'terrain', 'terrain_r', 'winter', 'winter_r']
        
        #ml.interactive(0)

        if self.type == 'scan':
            
            if 'channel' in params:
                channel = params['channel']
            else:
                channel = self.channels[0]
                
            if 'direction' in params:
                direction = params['direction']
            else:
                direction = 'forward'
                
            if 'flatten' in params:
                flatten = params['flatten']
            else:
                flatten = False
                
            if 'offset' in params:
                offset = params['offset']
            else:
                offset = False
                
            if 'cmap' in params:
                cmap = params['cmap']
            else:
                cmap = 'gray'
                
            if 'clim' in params:
                clim = params['clim']
            else:
                clim = False
                
            if 'log' in params:
                log = params['log']
            else:
                log = False    
                
            if 'show_params' in params:
                show_params = params['show_params']
            else:
                show_params = False

            if 'show' in params:
                show = params['show']
            else:
                show = True

            # PREMISE - specific params

            if 'hide' in params:
                hide = params['hide']
            else:
                hide = False

            if 'color_scale' in params:
                color_scale = params['color_scale']
            else:
                color_scale = False

            if 'x_axis' in params:
                x_axis = params['x_axis']
            else:
                x_axis = False

            if 'y_axis' in params:
                y_axis = params['y_axis']
            else:
                y_axis = False

            if 'colormap_scaling' in params:
                colormap_scaling = params['colormap_scaling']
            else:
                colormap_scaling = False
                
                
            
            (chData,chUnit) = self.get_channel(channel, direction = direction, flatten=flatten, offset=offset);
            
            if direction == 'backward':
                chData = np.fliplr(chData)
            
            width = self.get_param('width')
            height = self.get_param('height')
            pix_y,pix_x = np.shape(chData)

            # fig = plt.figure(figsize=(6, 4))
            fig = plt.figure()#figsize=(6,4))

            
            ImgOrigin = 'lower'
            if self.get_param('scan_dir') == 'down':
                ImgOrigin = 'upper'

            if color_scale:
                chData = np.clip(chData, color_scale[0], color_scale[1])

            if x_axis and y_axis:
                plt.ylim(y_axis[0], y_axis[1])
                plt.xlim(x_axis[0], x_axis[1])
                # plt.axis([x_axis[0], x_axis[1], y_axis[0], y_axis[1]])

            plot_params = dict(
                X = np.abs(chData) if log else chData,
                aspect = 'equal',
                extent = [0,width[0],0,height[0]],
                cmap = cmap,
                origin = ImgOrigin,
                interpolation = 'none'
            )
            if log:
                plot_params['norm'] = LogNorm()

            if colormap_scaling:
                min_x = x_axis[0]
                max_x = x_axis[1]
                min_y = y_axis[0]
                max_y = y_axis[1]

                shape = chData.shape
                minx = int(np.floor(shape[0] * min_x / width[0]))
                maxx = int(np.ceil(shape[0] * max_x / width[0]))
                miny = int(np.floor(shape[1] * min_y / height[0]))
                maxy = int(np.ceil(shape[1] * max_y / height[0]))

                mini = np.min(chData[minx:maxx, miny:maxy])
                maxi = np.max(chData[minx:maxx, miny:maxy])

                mini = int(mini*100)/100
                maxi = int(maxi*100)/100
                plot_params['vmin'] = mini
                plot_params['vmax'] = maxi

            im = plt.imshow(**plot_params)

            # if log:
            #     im = plt.imshow(np.abs(chData), aspect = 'equal', extent = [0,width[0],0,height[0]], cmap = cmap, norm=LogNorm(), origin = ImgOrigin)
            # else:
            #     im = plt.imshow(chData, aspect = 'equal',extent = [0,width[0],0,height[0]], cmap = cmap, origin = ImgOrigin)
            
            
            if clim:
                plt.clim(clim)
                
            
            # im.axes.set_xticks([0,width[0]])
            # im.axes.set_xticklabels([0,np.round(width[0],2)])
            # im.axes.set_yticks([0,height[0]])
            # im.axes.set_yticklabels([0,np.round(height[0],2)])
            
            if show_params:
                title = self.print_params(show = False);
            else:
                title = self.path  

            if not hide:
                # plt.title(title + '\n', loc='left')
                plt.xlabel('x (%s)' % width[1])
                plt.ylabel('y (%s)' % height[1])
                cbar = plt.colorbar(im)#,fraction=0.046, pad=0.02, format='%.2g',shrink = 0.5,aspect=10)
                cbar.set_label('channel:%s (%s)' % (channel,chUnit))
            
            if show:
                plt.show()
           # else:
           #     plt.close(fig)
            
            return fig 
            

           
        elif self.type == 'spec':
            
            
            if 'channelx' in params:
                channelx = params['channelx']
            else:
                channelx = self.channels[0]
            
            if 'channely' in params:
                channely = params['channely']
            else:
                channely = self.channels[1]
                
            if 'direction' in params:
                direction = params['direction']
            else:
                direction = direction = 'forward'
                
            if 'log' in params:
                log = params['log']
            else:
                log = False
                
            if 'loglog' in params:
                loglog = params['loglog']
            else:
                loglog = False
                
                
            if 'show_params' in params:
                show_params = params['show_params']
            else:
                show_params = False
                
            if 'show' in params:
                show = params['show']
            else:
                show = True
                
                
            (x_data,x_unit) = self.get_channel(channelx,direction=direction)
            (y_data,y_unit) = self.get_channel(channely,direction=direction)
            
            
            fig = plt.figure(figsize=(6,4))
            
            if log:
                plt.semilogy(x_data,np.abs(y_data))
                
            elif loglog:
                plt.loglog(np.abs(x_data),np.abs(y_data))
                
            else:
                plt.plot(x_data,y_data)          
                
           
            if show_params:
                title = self.print_params(show = False);
            else:
                title = self.path  
                
            plt.title(title + '\n', loc='left') 
                
            plt.xlabel('%s (%s)' % (channelx,x_unit))
            plt.ylabel('%s (%s)' % (channely,y_unit))
            
            if show:
                plt.show()
            else:
                plt.close(fig)
            
            return fig
                             
                             
 
            

##############################################
# functions for spm class
##############################################
        
        
# import all files in folder as list of spm objects
def importall(FilePath,FilePrefix = '',ImportOnly = '' ):
    
    from os import walk
    
    files = []
    NumSpec = 0
    NumScan = 0
    
    
    for root, dirs, filenames in walk(FilePath):
        for file in filenames:
            if file.endswith(".sxm") and file.startswith(FilePrefix):
                if not(ImportOnly == 'spec'):
                    files.append(spm(root + '/' + file))
                    NumScan = NumScan + 1
            elif file.endswith(".dat") and file.startswith(FilePrefix):
                if not(ImportOnly == 'scan'):
                    files.append(spm(root + '/' + file))
                    NumSpec = NumSpec + 1   
        
    print(str(len(files)) + ' files imported; ' + str(NumScan) + ' scan(s) and ' + str(NumSpec) + ' spectra')

    return files

# def importspecific(FilePath,FileList,FilePrefix = '',ImportOnly = '' ):
    
#     from os import walk
    
#     files = []
#     NumSpec = 0
#     NumScan = 0
    
    
#     for root, dirs, filenames in walk(FilePath):
#         for file in filenames:
#             if (file in FileList):
#                 if file.endswith(".sxm") and file.startswith(FilePrefix):
#                     if not(ImportOnly == 'spec'):
#                         files.append(spm(root + '/' + file))
#                         NumScan = NumScan + 1
#                 elif file.endswith(".dat") and file.startswith(FilePrefix):
#                     if not(ImportOnly == 'scan'):
#                         files.append(spm(root + '/' + file))
#                         NumSpec = NumSpec + 1   
        
#     print(str(len(files)) + ' files imported; ' + str(NumScan) + ' scan(s) and ' + str(NumSpec) + ' spectra')

#     return files




def specs_plot(specs,**params):
# plot spectra in list
    
    #import matplotlib.pyplot as plt
    
    if 'channelx' in params:
        channelx = params['channelx']
    else:
        channelx = specs[0].channels[0]
    
    if 'channely' in params:
        channely = params['channely']
    else:
        channely = specs[0].channels[1]
        
    if 'direction' in params:
        direction = params['direction']
    else:
        direction = list(specs[0].data.keys())[0]
    
    if 'color' in params:
        color = params['color']
    else:
        color = pl.cm.rainbow(np.linspace(0,1,len(specs)))
    if 'print_legend' in params:
        print_legend = params['print_legend']
    else:
        print_legend = False
    if 'offset' in params:
        offset = params['offset']
    else:
        offset = 0

    # PREMISE specific
    if 'show' in params:
        show = params['show']
    else:
        show = True

    if 'colormap' in params:
        colormap = params['colormap']
        cmap = plt.get_cmap(colormap)
        color = cmap(np.linspace(0,1,len(specs)))
    else:
        colormap = False

    if 'scaling' in params:
        scaling = params['scaling']
    else:
        scaling = False

    if 'x_axis' in params:
        x_axis = params['x_axis']
    else:
        x_axis = False

    if 'y_axis' in params:
        y_axis = params['y_axis']
    else:
        y_axis = False
    #######

    fig = plt.figure(figsize=(6,4))

    counter = 0
    for (s, c) in zip(specs, color):
        
        (x_data, x_unit) = s.get_channel(channelx, direction=direction)
        if x_axis:
            x_data = np.clip(x_data, x_axis[0], x_axis[1])

        (y_data, y_unit) = s.get_channel(channely, direction=direction)
        if y_axis:
            y_data = np.clip(y_data, y_axis[0], y_axis[1])

        if scaling:
            if scaling == 'lin-lin':
                plt.plot(x_data, y_data, color=c, label=s.name)
            elif scaling == "lin-log":
                plt.semilogy(x_data, np.abs(y_data), color=c, label=s.name)
            elif scaling == 'log-lin':
                plt.semilogx(np.abs(x_data), y_data, color=c, label=s.name)
            elif scaling == 'log-log':
                plt.loglog(np.abs(x_data), np.abs(y_data), color=c, label=s.name)
            else:
                plt.plot(x_data, y_data, color=c, label=s.name)
        else:
            plt.plot(x_data,y_data+counter*offset,color = c,label = s.name)
        counter = counter + 1
    
    plt.xlabel('%s (%s)' % (channelx,x_unit))
    plt.ylabel('%s (%s)' % (channely,y_unit))
    
    if print_legend:
        plt.legend()

    if show:
        plt.show()
    
    return fig


def relative_position(img,spec,**params):
    
    #width = ref.get_param('width')
    #height = ref.get_param('height')
    #[px_x,px_y] = ref.get_param('scan_pixels')

    [o_x,o_y] = img.get_param('scan_offset')
    width = img.get_param('width')[0]
    height = img.get_param('height')[0]
    [o_x,o_y] = [o_x*10**9,o_y*10**9]

    angle = float(img.get_param('scan_angle'))*-1* np.pi/180

    
    x_spec = spec.get_param('x')[0]
    y_spec = spec.get_param('y')[0]
    
    if angle != 0:
        #Transforming to relative coordinates with angle
        x_rel = (x_spec-o_x)*np.cos(angle) + (y_spec-o_y)*np.sin(angle)+width/2
        y_rel = -(x_spec-o_x)*np.sin(angle) + (y_spec-o_y)*np.cos(angle)+height/2
    else:
        x_rel = x_spec-o_x+width/2
        y_rel = y_spec-o_y+height/2
    return [x_rel,y_rel]

    

def ref_spec_plotting(ref_file,spec_files,fname_ref,fname_specs,**params):
    """
    Plots sxm reference image with spectra locations and multiple dat spectra

    Input:
    ref_file: filename of reference sxm image
    spec_files: list of filenames of dat spectra
    fname_ref: filename of saved reference image
    fname_specs: filename of saved spectra image

    Optional input:
    channelx_plot: channel to plot on x-axis of dat spectra (default: V)
    channely_plot: channel to plot on y-axis of dat spectra (default: dIdV)
    offset: offset between spectra (default: 0)
    annotate_spec: annotate spectra with numbers

    Output:
    specs_fig: figure of spectra
    """


    # Load modules
    import os

    # Optional input
    if 'channelx_plot' in params:
        channelx_plot = params['channelx_plot']
    else:
        channelx_plot = 'V'
    
    if 'channely_plot' in params:
        channely_plot = params['channely_plot']
    else:
        channely_plot = 'dIdV'

    if 'print_legend' not in params:
        params['print_legend'] = True

    if 'annotate_spec' in params:
        annotate_spec = params['annotate_spec']
    else:
        annotate_spec = False
      
    
    ## Import
    # Importing reference image
    ref = spm(ref_file)
    
    # Importing spec
    sp = []
    for s in spec_files:
        sp.append(spm(s))
    
    # Defining color scale
    col = pl.cm.rainbow(np.linspace(0,1,len(sp)))
    
    ## Plotting reference image with spec locations
    ref.plot(show_params = True, offset = False, show = False,channel = 'z',close_fig = False);

    # plot circle for each location
    for (s,c) in zip(sp,col):
        (x,y) = relative_position(ref,s)
        plt.plot(x,y,'o',color = c)
        if annotate_spec:
            plt.annotate(s.name, (x,y))

    plt.xlabel('x (nm)')
    plt.ylabel('y (nm)')
    
    figDir = os.path.abspath(os.path.join(fname_ref, os.pardir))
    if not os.path.exists(figDir):
        os.makedirs(figDir) 
    
    
    plt.savefig(fname_ref+'.png', dpi=500)
    plt.savefig(fname_ref+'.svg', dpi=500)

    # Plotting specs
    specs_fig = specs_plot(sp,channelx=channelx_plot,channely=channely_plot, direction = 'forward', color = col,**params);
    
    figDir = os.path.abspath(os.path.join(fname_specs, os.pardir))
    if not os.path.exists(figDir):
        os.makedirs(figDir) 
    
    specs_fig.savefig(fname_specs+'.png', dpi=500)
    specs_fig.savefig(fname_specs+'.svg', dpi=500)

    plt.show()
    
    return specs_fig
    
    
    
# returns dict list of fitted kpfm parabolas for list of spm objects
def kpfm(files,**params):
    
    if 'range' in params:
        range = params['range'];
    else:
        range = [0,len(files[0].get_channel('V')[0])-1]
    
    import numpy as np
    #sys.path.append('../')
    #import analyze as an
    
    data = {'V': [], 'df': [], 'V_max': [], 'df_max': [], 'p_fit': [], 'err_p': [], 'V_fit': [], 'df_fit': [], 'err_V_max': [], 'err_df_max': [], 'position': []}
    
    for f in files:
        if f.type == 'spec':
            if any(d['ChannelNickname'] == 'df' for d in f.SignalsList):
                
                data['V'].append(f.get_channel('V')[0])
                data['df'].append(f.get_channel('df')[0])
                
            if any(d['ChannelNickname'] == 'df_bw' for d in f.SignalsList):
                
                if not('df_bw' in data.keys()):
                    data['V_bw'] = [];
                    data['df_bw'] = [];
                    
                else:
                    data['V_bw'].append(f.get_channel('V_bw')[0])
                    data['df_bw'].append(f.get_channel('df_bw')[0])
                    
            x = f.get_param('x')[0]
            y = f.get_param('y')[0]
            
        data['position'].append([x,y])
         
                
    (p_fit,err_p,df_fit,V_max,err_V_max,df_max,err_df_max) = fit_parabola(np.array(data['V'])[:,range[0]:range[1]].tolist(),np.array(data['df'])[:,range[0]:range[1]].tolist(),single_spectrum=False,fitMin=False,fitMax=False)
    
    data['p_fit'] = p_fit;
    data['err_p'] = err_p;
    data['V_fit'] = np.array(data['V'])[:,range[0]:range[1]].tolist();
    data['df_fit'] = df_fit;
    data['V_max'] = V_max;
    data['err_V_max'] = err_V_max;
    data['df_max'] = df_max;
    data['err_df_max'] = err_df_max;
    
    if 'df_bw' in data.keys():
    
        (p_fit_bw,err_p_bw,df_fit_bw,V_max_bw,err_V_max_bw,df_max_bw,err_df_max_bw) = fit_parabola(np.array(data['V_bw'])[:,range[0]:range[1]].tolist(),np.array(data['df_bw'])[:,range[0]:range[1]].tolist(),single_spectrum=False,fitMin=False,fitMax=False)
        
        data['p_fit_bw'] = p_fit_bw;
        data['err_p_bw'] = err_p_bw;
        data['V_fit_bw'] = np.array(data['V_bw'])[:,range[0]:range[1]].tolist()
        data['df_fit_bw'] = df_fit_bw;
        data['V_max_bw'] = V_max_bw;
        data['err_V_max_bw'] = err_V_max_bw;
        data['df_max_bw'] = df_max_bw;
        data['err_df_max_bw'] = err_df_max_bw;
        
    return data
    

    
    
    
##############################################
# help functions
##############################################
    
def fit_parabola(
    bias,
    df,
    single_spectrum=False,
    fitMin=False,
    fitMax=False,
    ):
    
    #import numpy as np

    if single_spectrum:
        bias = np.array([bias])
        df = np.array([df])
    p_fit = np.zeros((len(bias), 3))
    df_fit = np.zeros((len(bias), len(bias[0])))
    err_p = np.zeros((len(bias), 3))
    bias_max = np.zeros(len(bias))
    err_bias_max = np.zeros(len(bias))
    df_max = np.zeros(len(bias))
    err_df_max = np.zeros(len(bias))

    for i in range(len(bias)):
        if not fitMin and not fitMax:
            (p, cov) = np.polyfit(bias[i], df[i], 2, cov=True)
        elif fitMin and not fitMax:

            fit_i = np.where(bias[i] > fitMin[i])[0]
            (p, cov) = np.polyfit(bias[i, fit_i], df[i, fit_i], 2,
                                  cov=True)
        p_fit[i] = p
        err1 = np.sqrt(abs(cov))
        err_p[i] = np.array([err1[0, 0], err1[1, 1], err1[2, 2]])

        df_fit[i] = np.polyval(p, bias[i])

        bias_max[i] = -p[1] / (2 * p[0])

        # if errors independent

        #err_bias_max[i] = abs(bias_max[i]) * np.sqrt((err1[1, 1]
        #        / p[1]) ** 2 + (err1[0, 0] / p[0]) ** 2)
        err_bias_max[i] = np.sqrt((err1[1, 1]/ (2*p[0]))** 2 + (err1[0, 0]*p[1] /(2*p[0]**2)) ** 2)

        # maximum error in any case (John Tayler Error analysis book)

        #err_bias_max[i] = abs(bias_max[i]) * (err1[1, 1] + err1[0, 0])

        df_max[i] = np.polyval(p, bias_max[i])
        err_df_max[i] = p[0] * bias_max[i] ** 2 * (err1[0, 0] + 2
                * err_bias_max[i] / abs(bias_max[i])) + p[1] \
            * bias_max[i] * (err1[1, 1] + err_bias_max[i]
                             / abs(bias_max[i])) + p[2] * err1[2, 2]
    return (
        p_fit,
        err_p,
        df_fit,
        bias_max,
        err_bias_max,
        df_max,
        err_df_max,
        )
    
# this function returns a rgb color code according to standard matlab colors    
def matlab_color(num):
    if num > 6:
        num = np.mod(num,7)

    colors = [[0,0.4470,0.7410],
            [0.8500,0.3250,0.0980],
            [0.9290,0.6940,0.1250],
            [0.4940,0.1840,0.5560],
            [0.4660,0.6740,0.1880],
            [0.3010,0.7450,0.9330],
            [0.6350,0.0780,0.1840]]
    
    return colors[np.round(num)]

