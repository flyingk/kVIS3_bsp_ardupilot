% kVIS3 Data Visualisation
%
% Copyright (C) 2012 - present  Kai Lehmkuehler, Matt Anderson and
% contributors
%
% Contact: kvis3@uav-flightresearch.com
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [BSP_Info] = BSP_ID()
% BSP ID function

% Name of BSP
BSP_Info.Name = 'ArduPilot';

BSP_Info.Version = '001';

% Data import function name
BSP_Info.importFcn = {'APM .bin File'     , 'import_APM_file'   ;
                      'APM .bin Folder'   , 'import_APM_folder' ;
                      'Export .csv Single', 'export_csv_single' ;
                      'Export .csv All'   , 'export_csv_all'   };

% a custom plot definition for a single axes plot to be shown on the main
% axes after import of data
BSP_Info.importPlot = '';

BSP_Info.armedChannel = '';

BSP_Info.flightModeChannel = '';

BSP_Info.mapChannels = {};

BSP_Info.aircraftVisualModelFile = '';

BSP_Info.customTabs = {};%{'Param','',''};

BSP_Info.addOns = {'Create SIDPAC file','','fill_fdata_Callback';
                   'Extract Event by Channel','','extract_event_by_channel_Callback'};

