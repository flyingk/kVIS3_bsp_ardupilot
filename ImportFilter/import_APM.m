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

% Constants
RTD = 180.0/pi;
t_start = inf;

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

[fds, parentNode] = kVIS_fdsAddTreeBranch(fds, 0, 'APM_Data');

%% read data
data_stream_names = fieldnames(log);

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
            data = log.(field_name);
            
            % Remove out what we don't need
            try data = rmfield(data,'typeNumID'); end
            try data = rmfield(data,'DatenumUTC'); end
            try data = rmfield(data,'TimeUS'); end
            
            % Get the number of points
            n_points = length(data.LineNo);
            data = rmfield(data,'LineNo');
            
            % Get the units (then remove field)
            units_data = data.('fieldUnits');
            data = rmfield(data,'fieldUnits');
            
            % Get field multipliers
            multipliers_data = data.('fieldMultipliers');
            data = rmfield(data,'fieldMultipliers');
            
            % Create name for new struct (group name)
            groupName = data.('name');
            data = rmfield(data,'name');
            
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
            
            for jj = 1:n_channels
                channel_data = data.(channel_names{jj});
                
                % Add to fds.fdata
                DAT(:,jj) = channel_data;
            end
            
            % Check to see if the data is valid
            if (max(DAT(:,1)) > 1e7)
                % Data is bad
                fprintf('\t\t\t\t\t\t  Channel corruption detected, removing bad points\n');
                locs = find(DAT(:,1) < 1e7);
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
            
            if (strcmp(groupName,'BARO'))
                
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
                
                % Add a new leaf to the tree
                fds = kVIS_fdsAddTreeLeaf(fds, groupName, varNames, varNames, varUnits, varFrames, DAT, parentNode, false);
            end
            
        end
    end
end

% Break up sensor data that has an 'Id' feild
fds = breakup_sensor_data(fds);

% Add vehicle data (if file found)
if exist([pathstr,'\',name,'.m'],'file')
    % Import data from a particular flight preferentially
    fds = add_vehicle_data([pathstr,'\',name,'.m'],fds);
    
elseif exist([pathstr,'\flight_info.m'],'file')
    % Import data from the folders
    fds = add_vehicle_data([pathstr,'\flight_info.m'],fds);
    
end


% Sort the fdata fields alphabetically
[~,idx] = sort(fds.fdata(1,2:end));
fds.fdata(:,2:end) = fds.fdata(:,idx+1);

% Fix up the time so starts at t = 0
fds.timeOffset = t_start;
for ii = 2:numel(fds.fdata(1,:))
    fds.fdata{7,ii}(:,1) = fds.fdata{7,ii}(:,1) - fds.timeOffset;
end

% All done!
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
    if max((strcmp(fds.fdata{2,ii},'Id')))
        % Find where the Id string is stored
        groupName_base = fds.fdata{1,ii};
        fprintf('Multiple sensors in %s\n',groupName_base)
        idx_Id = strcmp(fds.fdata{2,ii},'Id');
        id_list = unique(fds.fdata{7,ii}(:,idx_Id));
        
        for jj = 1:numel(id_list)
            % Work out which index we are using
            id = id_list(jj);
            idx = (fds.fdata{7,ii}(:,idx_Id) == id);
            groupName = sprintf('%s[%d]',groupName_base,id);
            
            % Assembly data
            varNames = fds.fdata{2,ii};
            varUnits =  fds.fdata{3,ii};
            varFrames = fds.fdata{4,ii};
            parentNode = 1;
            
            DAT = fds.fdata{7,ii}(idx,:);
            
            % Add leaf to tree
            fds = kVIS_fdsAddTreeLeaf(fds, groupName, varNames, varNames, varUnits, varFrames, DAT, parentNode, false);
            
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
