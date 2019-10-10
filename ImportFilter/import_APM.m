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

function [] = import_APM(hObject, ~)

% UAVmainframe v5 file
[file, pathname] = uigetfile('*.bin');

% Constants
RTD = 180.0/pi;
t_start = inf;

% Load file
if file==0
    disp('Error loading file.')
    return
else
    file = fullfile(pathname,file);
end


tic

fprintf('Importing APM .bin file\n');
fprintf('\t%s\n\n',file);

% Import .bin file
log = Ardupilog(file);
log = log.getStruct();

%% get new fds structure
fds = kVIS_fdsInitNew();

fds.BoardSupportPackage = 'ArduPilot';

[fds, parentNode] = kVIS_fdsAddTreeBranch(fds, 0, 'APM_Data');


%% read data
data_stream_names = fieldnames(log);

fds.fdata{ 2,1} = 'varNames';
fds.fdata{ 3,1} = 'varUnits';
fds.fdata{ 4,1} = 'varFrames';
fds.fdata{ 5,1} = 'varNamesDisp';
fds.fdata{ 6,1} = 'varLabelsTeX';
fds.fdata{ 7,1} = 'data';
fds.fdata{ 8,1} = 'treeParent';
fds.fdata{ 9,1} = 'treeGroupExpanded';
fds.fdata{10,1} = 'treeGroupSelected';

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
            strcmp(field_name,'MSG') || ...
            strcmp(field_name,'PARM') || ...
            ...
            isempty(field_name) )
        
        % skip these
        fprintf('\t-Skip- %14s\n',field_name);
    else
        
        % Check to see if there is any time data.  If there isn't then the
        % file is empty and we don't need to store it
        if isempty(getfield(getfield(log,field_name),'TimeS'))
            fprintf('\t-Skip- %14s (Empty Data Set)\n',field_name);
            
        else
            % Looks like we have data, import it :)
            fprintf('\tImport %14s\n',field_name);
            data = log.(field_name);
            
            % Remove out what we don't need
            data = rmfield(data,'typeNumID');
            data = rmfield(data,'DatenumUTC');
            data = rmfield(data,'TimeUS');
            
            % Get the number of points
            n_points = length(data.LineNo);
            data = rmfield(data,'LineNo');
            
            % Get the units (then remove field)
            units_data = data.('fieldUnits');
            data = rmfield(data,'fieldUnits');
            
            % Get field multipliers
            multipliers_data = data.('fieldMultipliers');
            data = rmfield(data,'fieldMultipliers');
            
            % Add a new column to fds.fdata for our new data
            fds.fdata = [fds.fdata, cell(size(fds.fdata,1),1)];
            idx_data = size(fds.fdata,2);
            
            % Create name for new struct (group name)
            fds.fdata( 1,idx_data) = {data.('name')};
            data = rmfield(data,'name');
            
            % Get number of fields
            channel_names = fieldnames(data);
            n_channels = numel(channel_names);
            
            % Group parent in tree
            fds.fdata( 8,idx_data) = {1};
            
            % Group node expanded
            fds.fdata( 9,idx_data) = {0};
            
            % Group node selected
            fds.fdata(10,idx_data) = {1};
            
            % List of channel names
            fds.fdata{ 2,idx_data} = fieldnames(data); % reference name
            fds.fdata{ 5,idx_data} = fieldnames(data); % display name
            fds.fdata{ 6,idx_data} = fieldnames(data); % tex name

            
            % List of channel units
            if (isempty(fieldnames(units_data)))
                fds.fdata{ 3,idx_data} = repmat({'N/A'},n_channels,1);
            else
                keyboard
            end
            
            % Reference frame of channel
            fds.fdata{ 4,idx_data} = repmat({'Unknown Frame'},n_channels,1);
            
            % Fill in the data
            fds.fdata{ 7,idx_data} = nan(n_points,n_channels);

            for jj = 1:n_channels
                channel_data = data.(channel_names{jj});
                
                % Add to fds.fdata
                fds.fdata{ 7,idx_data}(:,jj) = channel_data;   
            end
        end
    end
end

% Sort everything alphabetically
% Otherwise everything else is just random @_@


%
% % Fix up the time so starts at t = 0
% for ii = 3:size(file_data,2)
%     file_data{2,ii}(:,1) = file_data{2,ii}(:,1) - t_start;  % Start at t = 0
% end
%
% % Use a common time step (maybe)
%
% % Get user to fill in additional required fields
% prompt = { 'Version:', ...
%     'Activity:', ...
%     'Date:', ...
%     'Time:', ...
%     'Aircraft:', ...
%     'Pilot:', ...
%     'Sref:', ...
%     'cref:', ...
%     'bref:', ...
%     'mass:', ...
%     'Ixx:', ...
%     'Iyy:', ...
%     'Izz:', ...
%     'Comment:'  };
% title = 'APM Data';
% dims = [1 35];
% definput = {'APM:Plane', ...
%     'Testing', ...
%     date, ...
%     '12:00:00', ...
%     'Mini-Skywalker', ...
%     'Jag + Matt',...
%     '1', ...
%     '1', ...
%     '1', ...
%     '1', ...
%     '1', ...
%     '1', ...
%     '1', ...
%     ''   };
% answer = inputdlg(prompt,title,dims,definput);
%
% % Fill out fds data
% fds = evalin('base','fds');
%
% fds.aircraft  = answer{ 5};
% fds.ver       = answer{ 1};
% fds.maneuver  = answer{ 2};
% fds.test_date = [answer{ 3},' ',answer{ 4}];
% fds.pilot     = answer{ 6};
% fds.comment   = answer{14};
%
% % Fill out config data (h)
% h.sref     = answer{ 7};
% h.cref     = answer{ 8};
% h.bref     = answer{ 9};
% h.mass     = answer{10};
% h.ixx      = answer{11};
% h.iyy      = answer{12};
% h.izz      = answer{13};
%
fprintf('\nImport took %.2f s\n',toc);

%% Update KSID
fds = kVIS_fdsUpdateAttributes(fds);

fds = kVIS_fdsGenerateTexLabels(fds);

kVIS_addDataSet(hObject, fds, []);

return

