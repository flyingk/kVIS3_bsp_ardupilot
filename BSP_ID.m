function [BSP_Info] = BSP_ID()
% BSP ID function

% Name of BSP
BSP_Info.Name = 'ArduPilot';

% Data import function name
BSP_Info.importFcn = {'ArduPilot mat file','import_APM'};

% a custom plot definition for a single axes plot to be shown on the main
% axes after import of data
BSP_Info.importPlot = '';

BSP_Info.armedChannel = '';

BSP_Info.flightModeChannel = '';

BSP_Info.mapChannels = {};

BSP_Info.aircraftVisualModelFile = '';

BSP_Info.customTabs = {};

BSP_Info.addOns = {'Create SIDPAC file','fill_fdata_Callback'};

