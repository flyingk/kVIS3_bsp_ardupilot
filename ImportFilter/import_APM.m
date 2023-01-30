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

function fds = import_APM(file)

debug = 1;

if (~nargin)
    file = './../Sample_Data/00000008.BIN';
    dbstop if error
end

% Constants
RTD = 180.0/pi;
t_start = inf;
t_end = -inf;

% Load file
if file==0
    disp('Error loading file.')
    return
end

tic

fprintf('Importing APM .bin file\n');
fprintf('\t%s\n\n',file);

% Import .bin file
log = Ardupilog(file);
log = log.getStruct();

% Get file name
[pathstr,name,ext] = fileparts(file);

%% get new fds structure
fds = kVIS_fdsInitNew();

fds.BoardSupportPackage = 'ArduPilot';

[fds, rootNode] = kVIS_fdsAddTreeBranch(fds, 0, 'APM_Data');

%% read data
data_stream_names = fieldnames(log);
storedGroupNames = {};

for ii = 1:numel(data_stream_names) % change to a 1 later
    
    field_name = data_stream_names{ii};

    % Exception list
    if ( strcmp(field_name,'fileName') || ...
            strcmp(field_name,'filePathName' ) || ...
            strcmp(field_name,'platform')  || ...
            strcmp(field_name,'version') || ...
            strcmp(field_name,'commit')  || ...
            strcmp(field_name,'bootTimeUTC') || ...
            strcmp(field_name,'msgsContained')  || ...
            strcmp(field_name,'totalLogMsgs') || ...
            strcmp(field_name,'msgFilter') || ...
            strcmp(field_name,'numMsgs') || ...
            ...
            strcmp(field_name,'FMT') || ...
            strcmp(field_name,'FMTU') || ...
            strcmp(field_name,'UNIT') || ...
            strcmp(field_name,'MULT') || ...
            strcmp(field_name,'MSG') || ...
            strcmp(field_name,'PARM') || ...
            ...
            strcmp(field_name,'ISBD') || ...
            ...
            isempty(field_name) )
        
        % skip these
        if (debug); fprintf('\t-Skip- %14s (Not a Data Set)\n',field_name); end
    else
        
        % Check to see if there is any time data.  If there isn't then the
        % file is empty and we don't need to store it
        if isempty(getfield(getfield(log,field_name),'TimeS'))
            if (debug); fprintf('\t-Skip- %14s (Empty Data Set)\n',field_name); end
            
        else
            % Looks like we have data, import it :)
            fprintf('\tImport %14s\n',field_name);
            sensorData = log.(field_name);

            % Work out how many instances
            n_instances = numel(sensorData);
            if n_instances > 1
                fprintf('\t\t%d Instances\n',n_instances);
                % Create a new branch for multiple-instance sensor data
                [fds, sensorNode] = kVIS_fdsAddTreeBranch(fds, rootNode, field_name);
            else
                % Use the root node
                sensorNode = rootNode;
            end

            % Loop through each instance
            for jj = 1:n_instances

                data = sensorData(jj);

                % Remove out what we don't need
                try data = rmfield(data,'typeNumID'); end
                try data = rmfield(data,'DatenumUTC'); end
                try data = rmfield(data,'TimeUS'); end
                
                % Get the number of points
                n_points = length(data.LineNo);
                data = rmfield(data,'LineNo');
    
                % Get the instance number (then remove the field)
                fields = fieldnames(data.fieldUnits);
                instance_number = data.MsgInstance;
                for kk = 1:numel(fields)
                    if strcmp(data.fieldUnits.(fields{kk}),'instance')
                        instance_number = data.(fields{kk})(1);
                    end
                end
                
                data = rmfield(data,'MsgInstance');
                
                % Get the units (then remove field)
                units_data = data.('fieldUnits');
                data = rmfield(data,'fieldUnits');
                
                % Get field multipliers
                multipliers_data = data.('fieldMultipliers');
                data = rmfield(data,'fieldMultipliers');
                
                % Create name for new struct (group name)
                groupName = data.('name'); ...
                data = rmfield(data,'name');
                if n_instances > 1

                    groupName = [groupName,'_',num2str(instance_number)];
                    nameUnique = ~any(strcmp(storedGroupNames,groupName));

                    if (nameUnique)
                        storedGroupNames{end+1,1} = groupName;
                    else
                        fprintf('\t\tGroup %s already exists, skipping\n',groupName);
                        continue
                    end 
                end
                
                % Get number of fields
                channel_names = fieldnames(data);
                n_channels = numel(channel_names);
                
                % List of channel names
                varNames = fieldnames(data); % reference name
                varNames{(strcmp(varNames,'TimeS'))} = 'Time'; % Change TimeS to Time for compatibility
                
                % List of channel units
                if (isempty(fieldnames(units_data)))
                    varUnits = repmat({'N/A'},n_channels,1);
                else
                    varUnits = struct2cell(units_data);
                    varUnits(cellfun('isempty',varUnits)) = {'N/A'};
                end
                
                % Reference frame of channel
                varFrames = repmat({''},n_channels,1);
                
                DAT=[];
                
                for kk = 1:n_channels
    
    %                 if (strcmp(channel_names{kk},'MsgInstance'))
    %                     % This is from the splitting up instances, we need to
    %                     % remove it otherwise things break...
    %                     a = 1;
    %                     varNames(kk) = [];
    %                     %varUnits(kk) = [];
    %                     varFrames(kk) = [];
    % 
    %                 else
                        channel_data = data.(channel_names{kk});
                        
                        % Add to fds.fdata if data is numeric, otherwise nan it
                        if (isnumeric(channel_data))
                            DAT(:,kk) = channel_data;
                        else
                            DAT(:,kk) = nan(size(channel_data,1),1);
                        end
    
    %                     a = 1;
    %                 end
                end
                
                % Check to see if the data is valid
                if (max(DAT(:,1)) > 1e6)
                    % Data is bad
                    fprintf('\t\t\t\t\t\t  Channel corruption detected, removing bad points\n');
                    locs = find(DAT(:,1) < 1e6);
                    DAT = DAT(locs,:);
                end
                
                % Add additional fields to certain groups
                if (strcmp(groupName,'NKF1') || ...
                        strcmp(groupName,'NKF6')   )
                    
                    varNames  = [ varNames', 'V', 'DistToHome']';
                    varUnits  = [ varUnits', 'm/s', 'm' ]';
                    varFrames = [ varFrames', {''}, {''}]';
                    
                    VN = DAT(:,strcmp(varNames,'VN'));
                    VE = DAT(:,strcmp(varNames,'VE'));
                    
                    PN = DAT(:,strcmp(varNames,'PN'));
                    PE = DAT(:,strcmp(varNames,'PE'));
                    
                    DAT = [ DAT, sqrt(VN.^2 + VE.^2), sqrt(PN.^2 + PE.^2)];
                    
                end
                
                if (strcmp(groupName,'BAT'))
                    
                    varNames  = [ varNames', 'Power']';
                    varUnits  = [ varUnits', 'W' ]';
                    varFrames = [ varFrames', {''}]';
                    
                    V = DAT(:,strcmp(varNames,'Volt'));
                    I = DAT(:,strcmp(varNames,'Curr'));
                    
                    DAT = [ DAT, V.*I];
                    
                end
                
                if (strcmp(groupName,'PSC'))
                    
                    varNames  = [ varNames', 'TV', 'V']';
                    varUnits  = [ varUnits', 'm/s', 'm/s']';
                    varFrames = [ varFrames', {''}, {''}]';
                    
                    TVX = DAT(:,strcmp(varNames,'TVX'));
                    TVY = DAT(:,strcmp(varNames,'TVY'));
                    VX  = DAT(:,strcmp(varNames,'VX'));
                    VY  = DAT(:,strcmp(varNames,'VY'));
                    
                    DAT = [ DAT, sqrt(TVX.^2 + TVY.^2), sqrt(VX.^2 + VY.^2)];
                    
                end
                
                if (strcmp(groupName,'BARO') || strcmp(groupName,'BAR2') || strcmp(groupName,'BAR3'))
                    
                    varNames  = [ varNames', 'Rho']';
                    varUnits  = [ varUnits', 'kg/m3']';
                    varFrames = [ varFrames', {''}]';
                    
                    P = DAT(:,strcmp(varNames,'Press'));
                    T = DAT(:,strcmp(varNames,'Temp'));
                    rho = calcDensity(P,T);
                    
                    DAT = [ DAT, rho ];
                    
                end
                
                % Add data to fds
                if ~isempty(DAT)
                    t_start = min(t_start,DAT(1,1));
                    t_end   = max(t_end  ,DAT(end,1));

                    % Replace nan data with 0
                    DAT(isnan(DAT))=0;

                    % Check len(varNames) == len(varUnits) ==
                    % len(varFrames) == len(DAT)
                    n_expected = size(DAT,2);
                    if (numel(varNames) ~= n_expected)
                        fprintf('\t\tName length mismatch, going to use default field names\n');
                        for kk = 1:n_expected
                            varNames{kk,1} = sprintf('field_%d',kk);
                        end
                        keyboard
                    end
                    if (numel(varUnits) ~= n_expected)
                        fprintf('\t\tUnit length mismatch, going to default to N/As\n');
                        varUnits = repmat({'N/A'},n_expected,1);
                    end
                    if (numel(varFrames) ~= n_expected)
                        fprintf('\t\tFrame length mismatch, going to default to N/As\n');
                        varFrames = repmat({'N/A'},n_expected,1);
                    end
                    
                    % Add a new leaf to the tree
                    fds = kVIS_fdsAddTreeLeaf(fds, groupName, varNames, varNames, varUnits, varFrames, DAT, sensorNode, false);
                end

            end
        end
    end

    if strcmp(field_name,'PARM') 
        % Create a special field for the params
        for jj = 1:length(log.(field_name).Name)
            param_name =  (log.(field_name).Name(jj,:));
            param_name = deblank(param_name);
            if (isvarname(param_name))
                fds.Param.(param_name) = (log.(field_name).Value(jj));   
            else
                fprintf('Couldn''t store param %s\n',param_name);
            end
        end
    end
    
    if strcmp(field_name,'MSG')
        % Create a special field for the messages
        fds.msg.Time = log.(field_name).TimeS/1e6;
        fds.msg.MSG = log.(field_name).Message;
    end
    
end

% Update the fds attributes
fds = kVIS_fdsUpdateAttributes(fds);

%% Break up sensor data that has an 'Id' feild
%fds = breakup_sensor_data(fds);

% Add vehicle data (if file found)
if exist([pathstr,'\',name,'.m'],'file')
    % Import data from a particular flight preferentially
    fds = add_vehicle_data([pathstr,'\',name,'.m'],fds);
    
elseif exist([pathstr,'\flight_info.m'],'file')
    % Import data from the folders
    fds = add_vehicle_data([pathstr,'\flight_info.m'],fds);
    
end

%% Sort the fdata fields alphabetically
% There needs to be a more intelligent way of doing this
% as it will break any references to parent nodes.  Currently
% files import alphabetically so doesn't really matter.
% [~,idx] = sort(fds.fdata(1,2:end));
% fds.fdata(:,2:end) = fds.fdata(:,idx+1);

%% Fix time stamping

fds.timeOffset = t_start; t_end = t_end - t_start;
for ii = 1:numel(fds.fdata(1,:))

    % Skip if this is a data header
    if isempty(fds.fdata{fds.fdataRows.data,ii})
        continue
    end

    % Fix up the time so starts at t = 0
    fds.fdata{fds.fdataRows.data,ii}(:,1) = fds.fdata{fds.fdataRows.data,ii}(:,1) - fds.timeOffset;

    % Fix anything time less than 0, was likely and nan before  
    idx = (fds.fdata{fds.fdataRows.data,ii}(:,1)) < 0.0;
    fds.fdata{fds.fdataRows.data,ii}(idx,1) = 0;

end

%% Add events based on flight mode changes
% We can use MODE.ModeNum to work this out
fprintf('Generating events based on MODE changes.\n');
modeTimes   = kVIS_fdsGetChannel(fds, 'MODE','Time');
modeNumbers = kVIS_fdsGetChannel(fds, 'MODE','ModeNum');
modeReasons = kVIS_fdsGetChannel(fds, 'MODE','Rsn');

modeTimes(1) = 0; % Force the first mode time to 0

% Combine short and same mode changes
if numel(modeTimes > 1)
    ii = 2;
    while ii < numel(modeNumbers)
        if modeNumbers(ii-1) == modeNumbers(ii)
            % Remove entry as mode didn't change
            modeNumbers(ii) = [];
            modeTimes(ii)   = [];
            modeReasons(ii) = [];
        else
            ii = ii+1;
        end
    end
end

% Add a marker for the mode at log's end
modeTimes(end+1) = t_end;
modeNumbers(end+1) = modeNumbers(end);
modeReasons(end+1) = modeReasons(end);

% Loop through modeTimes and store data
eventNumber = 0;
for ii = 1:numel(modeTimes)-1
    % Get info
    t_in  = modeTimes(ii);
    t_out = modeTimes(ii+1);
    
    % Check if change was valid
    if t_out-t_in > 1.0
        
        modeReason = modes_Reason(modeReasons(ii));
        
        % TODO:  Need to work out if using Copter/Plane/Rover etc
        %        Defaulting to COPTER for now though
        if contains(fds.msg.MSG(1,:),'ArduCopter')
            modeType = modes_ArduCopter(modeNumbers(ii));
        elseif contains(fds.msg.MSG(1,:),'ArduPlane')
            modeType = modes_ArduPlane(modeNumbers(ii));
        else
            modeType = sprintf('Mode %d',modeNumbers(ii));
        end
        
        % Fill out eList
        eventNumber = eventNumber+1;
        eList(eventNumber).type = modeType;
        eList(eventNumber).start= t_in;
        eList(eventNumber).end  = t_out;
        eList(eventNumber).description = modeReason;
        eList(eventNumber).plotDef='';
    end
        
end

% Add eList to eventList
if exist('eList','var')
   fds.eventList = eList;
end

%% All done!
fprintf('\nImport took %.2f s\n',toc);

return
end

function fds = add_vehicle_data(file,fds)
% Automatically adds information about the flight from a file
fprintf('\n\tAdding vechile data from %s\n',file);
run(file);

% Aircraft data
try; fds.aircraftData.acIdentifier = acIdentifier; end
try; fds.aircraftData.sRef = sRef; end
try; fds.aircraftData.cRef = cRef; end
try; fds.aircraftData.bRef = bRef; end
try; fds.aircraftData.mass = mass; end
try; fds.aircraftData.ixx = ixx; end
try; fds.aircraftData.iyy = iyy; end
try; fds.aircraftData.izz = izz; end
try; fds.aircraftData.ixz = ixz; end
try; fds.aircraftData.xCG = xCG; end
try; fds.aircraftData.yCG = yCG; end
try; fds.aircraftData.zCG = zCG; end

% Test Info
try; fds.testInfo.date = 0; end
try; fds.testInfo.startTime = 0; end
try; fds.testInfo.description = 0; end
try; fds.testInfo.pilot = 0; end
try; fds.testInfo.location = 0; end
try; fds.testInfo.weather = 0; end
try; fds.testInfo.windDir = 0; end
try; fds.testInfo.windSpeed = 0; end
try; fds.testInfo.ambientPressure = 0; end
try; fds.testInfo.ambientTemperature = 0; end

return;
end

function fds = breakup_sensor_data(fds)

% Certain data groups have an ID where multiple sensors are stored in the
% same file.  These need to be broken up
for ii = 1:size(fds.fdata,2)  
    % Check to see if data stream has an Id field
    if max((strcmp(fds.fdata{fds.fdataRows.varNames,ii},'Id')))
        % Find where the Id string is stored
        groupName_base = fds.fdata{1,ii};
        
        if strcmp(groupName_base,'D32') || ...
                strcmp(groupName_base,'DU32') || ...
                strcmp(groupName_base,'EV')
            
            % Skip these, they're meant to be ID
            
        else
            
            fprintf('Multiple sensors in %s\n',groupName_base)
            idx_Id = strcmp(fds.fdata{fds.fdataRows.varNames,ii},'Id');
            id_list = unique(fds.fdata{fds.fdataRows.data,ii}(:,idx_Id));
            
            for jj = 1:numel(id_list)
                % Work out which index we are using
                id = id_list(jj);
                idx = (fds.fdata{fds.fdataRows.data,ii}(:,idx_Id) == id);
                groupName = sprintf('%s[%d]',groupName_base,id);
                
                fprintf('\tAdding %s\n',groupName);
                
                % Assembly data
                varNames = fds.fdata{fds.fdataRows.varNames,ii};
                varUnits =  fds.fdata{fds.fdataRows.varUnits,ii};
                varFrames = fds.fdata{fds.fdataRows.varFrames,ii};
                parentNode = 1;
                
                DAT = fds.fdata{fds.fdataRows.data,ii}(idx,:);
                
                % Add leaf to tree
                fds = kVIS_fdsAddTreeLeaf(fds, groupName, varNames, varNames, varUnits, varFrames, DAT, parentNode, false);
                
            end
        end
    end
end

% Probably a good idea to remove the original file


% End of function
return
end

function rho = calcDensity(P,T)
% Calculates the air density (assuming dry air)
% P is pressure in Pascals
% T is temperature in Celcius

R = 287.05;      % Specific gas constant for dry air [ J/(kg*K) ]
T = T + 273.15;  % Convert temperature to Kelvin [ K ]

% Ideal gas law
rho = P ./ (R.*T);

return
end

function modeName = modes_ArduCopter(modeNumber)
% From https://github.com/ArduPilot/ardupilot/blob/master/ArduCopter/mode.h#L14

switch (modeNumber)
    case 0; modeName = 'STABILIZE';
    case 1; modeName = 'ACRO';
    case 2; modeName = 'ALT_HOLD';
    case 3; modeName = 'AUTO';
    case 4; modeName = 'GUIDED';
    case 5; modeName = 'LOITER';
    case 6; modeName = 'RTL';
    case 7; modeName = 'CIRCLE';
    case 9; modeName = 'LAND';
    case 11; modeName = 'DRIFT';
    case 13; modeName = 'SPORT';
    case 14; modeName = 'FLIP';
    case 15; modeName = 'AUTOTUNE';
    case 16; modeName = 'POSHOLD';
    case 17; modeName = 'BRAKE';
    case 18; modeName = 'THROW';
    case 19; modeName = 'AVOID_ADSB';
    case 20; modeName = 'GUIDED_NOGPS';
    case 21; modeName = 'SMART_RTL';
    case 22; modeName = 'FLOWHOLD';
    case 23; modeName = 'FOLLOW';
    case 24; modeName = 'ZIGZAG';
    case 25; modeName = 'SYSTEMID';
    case 26; modeName = 'AUTOROTATE';
    otherwise; modeName = ['UNKNOWN_(',num2str(modeNumber),')'];
end
return
end

function modeName = modes_ArduPlane(modeNumber)
% from https://github.com/ArduPilot/ardupilot/blob/master/ArduPlane/mode.h#L21

switch (modeNumber)
    case 0; modeName = 'MANUAL';
    case 1; modeName = 'CIRCLE';
    case 2; modeName = 'STABILIZE';
    case 3; modeName = 'TRAINING';
    case 4; modeName = 'ACRO';
    case 5; modeName = 'FLY_BY_WIRE_A';
    case 6; modeName = 'FLY_BY_WIRE_B';
    case 7; modeName = 'CRUISE';
    case 8; modeName = 'AUTOTUNE';
    case 10; modeName = 'AUTO';
    case 11; modeName = 'RTL';
    case 12; modeName = 'LOITER';
    case 13; modeName = 'TAKEOFF';
    case 14; modeName = 'AVOID_ADSB';
    case 15; modeName = 'GUIDED';
    case 16; modeName = 'INITIALISING';
    case 17; modeName = 'QSTABILIZE';
    case 18; modeName = 'QHOVER';
    case 19; modeName = 'QLOITER';
    case 20; modeName = 'QLAND';
    case 21; modeName = 'QRTL';
    case 22; modeName = 'QAUTOTUNE';
    case 23; modeName = 'QACRO';
    case 24; modeName = 'THERMAL';
    otherwise; modeName = ['UNKNOWN_(',num2str(modeNumber),')'];
end
return
end

function modeReason = modes_Reason(reasonNumber)
% from https://github.com/ArduPilot/ardupilot/blob/master/libraries/AP_Vehicle/ModeReason.h

switch (reasonNumber)
    case 0; modeReason = 'UNKNOWN';
    case 1; modeReason = 'RC_COMMAND';
    case 2; modeReason = 'GCS_COMMAND';
    case 3; modeReason = 'RADIO_FAILSAFE';
    case 4; modeReason = 'BATTERY_FAILSAFE';
    case 5; modeReason = 'GCS_FAILSAFE';
    case 6; modeReason = 'EKF_FAILSAFE';
    case 7; modeReason = 'GPS_GLITCH';
    case 8; modeReason = 'MISSION_END';
    case 9; modeReason = 'THROTTLE_LAND_ESCAPE';
    case 10; modeReason = 'FENCE_BREACHED';
    case 11; modeReason = 'TERRAIN_FAILSAFE';
    case 12; modeReason = 'BRAKE_TIMEOUT';
    case 13; modeReason = 'FLIP_COMPLETE';
    case 14; modeReason = 'AVOIDANCE';
    case 15; modeReason = 'AVOIDANCE_RECOVERY';
    case 16; modeReason = 'THROW_COMPLETE';
    case 17; modeReason = 'TERMINATE';
    case 18; modeReason = 'TOY_MODE';
    case 19; modeReason = 'CRASH_FAILSAFE';
    case 20; modeReason = 'SOARING_FBW_B_WITH_MOTOR_RUNNING';
    case 21; modeReason = 'SOARING_THERMAL_DETECTED';
    case 22; modeReason = 'SOARING_THERMAL_ESTIMATE_DETERIORATED';
    case 23; modeReason = 'VTOL_FAILED_TRANSITION';
    case 24; modeReason = 'VTOL_FAILED_TAKEOFF';
    case 25; modeReason = 'FAILSAFE';
    case 26; modeReason = 'INITIALISED';
    case 27; modeReason = 'SURFACE_COMPLETE';
    case 28; modeReason = 'BAD_DEPTH';
    case 29; modeReason = 'LEAK_FAILSAFE';
    case 30; modeReason = 'SERVOTEST';
    case 31; modeReason = 'STARTUP';
    case 32; modeReason = 'SCRIPTING';
    case 33; modeReason = 'UNAVAILABLE';
    case 34; modeReason = 'AUTOROTATION_START';
    case 35; modeReason = 'AUTOROTATION_BAILOUT';
    case 36; modeReason = 'SOARING_ALT_TOO_HIGH';
    case 37; modeReason = 'SOARING_ALT_TOO_LOW';
    case 38; modeReason = 'SOARING_DRIFT_EXCEEDED';
    case 39; modeReason = 'RTL_COMPLETE_SWITCHING_TO_VTOL_LAND_RTL';
    case 40; modeReason = 'RTL_COMPLETE_SWITCHING_TO_FIXEDWING_AUTOLAND';
    case 41; modeReason = 'MISSION_CMD';
    case 42; modeReason = 'FRSKY_COMMAND';
    otherwise; modeReason = 'UNKNOWN';
end

return
end

