function [] = import_APM(hObject, ~)

% Constants
RTD = 180.0/pi;
t_start = inf;


% UAVmainframe v5 file
[file, pathname] = uigetfile('*.mat');

if file==0
    disp('Error loading file.')
    return
else
    file = fullfile(pathname,file);
end


tic

fprintf('Importing APM MATLAB file\n');
fprintf('\t%s\n\n',file);

% get valid file names inside mat file
aa = whos('-file', file);

for i = 1:size(aa, 1)
    
    if isstrprop(aa(i).name(1),'upper')  % Check if file is value (will all be capilitalised at first letter)
        APM_files{i} = aa(i).name; % Copy valid file over to APM_files struct
    end
    
end

%% get new fds structure
fds = kVIS_fdsInitNew();

fds.BoardSupportPackage = 'ArduPilot';

[fds, parentNode] = kVIS_fdsAddTreeBranch(fds, 0, 'APM_Data');


%% read data
for ii = 1:length(APM_files) % change to a 1 later
    
    % Exception list
    if strcmp(APM_files{ii},'MSG1') || ...
       strcmp(APM_files{ii},'FMT' ) || ...
       strcmp(APM_files{ii},'MSG')  || ...
       strcmp(APM_files{ii},'PARM') || ...
       strcmp(APM_files{ii},'Seen') || ...
       isempty(APM_files{ii})
    
    % skip these
    
    else
        
        % Loop through strut to find fields with data
        if isempty(strfind(APM_files{ii},'_label'))
            
            fprintf('Importing file %s\n',APM_files{ii});

            % data
            APM_data = load(file, '-mat', APM_files{ii});
            
            DAT = APM_data.(APM_files{ii});
            
            DAT = DAT(:,2:end);
            
            DAT(:,1) = DAT(:,1)/1e6;
            
            % labels
            
            APM_data = load(file, '-mat', [APM_files{ii} '_label']);
            
            varNames = APM_data.([APM_files{ii} '_label']);
            
            varNames = varNames(2:end);
            
            varNames{1} = 'Time';
            
            % units
            
            varUnits = cell(size(varNames));
            
            varUnits = cellfun(@(x) 'N/A', varUnits, 'UniformOutput', false);
            
            varFrames = cell(size(varUnits));
            varFrames = cellfun(@(x) '', varFrames, 'UniformOutput', false);
            
            fds = kVIS_fdsAddTreeLeaf(fds, APM_files{ii}, varNames, varNames, varUnits, varFrames, DAT, parentNode, false);
            
%             % Import data
%             file_data = [file_data ,cell(2,1)];  % Create new entry
%             file_data{1,end} = fields{ii}; 
%             temp = getfield(APM_data,fields{ii}); file_data{2,end} = temp(:,2:end);               
%             
%             % Create labels
%             % Unit labels are stored in UNITS
%             % Potentially will be able to import the unit automatically in the future
%             % See /libraries/DataFlash/LogStructure.h::const struct UnitStructure log_Units[] for more
%             varlabAA = [varlabAA, cell(1,1)];
%             temp = getfield(APM_data,[fields{ii},'_label']); varlabAA{end} = temp(2:end);
%             for jj = 1:length(varlabAA{end})
%                 varlabAA{end}{jj} = ['  ',varlabAA{end}{jj}];
%             end
%             
%             % Fake units
%             % Units are stored in FMTU and are multiplied by MULTI
%             % Potentially will be able to import the unit automatically in the future
%             % See /libraries/DataFlash/LogStructure.h::const struct UnitStructure log_Units[] for more
%             varunitsAA = [varunitsAA, cell(1,1)];
%             varunitsAA{end} = cell(size(file_data{2,end},2),1);
%             for jj = 1:length(varunitsAA{end})
%                 varunitsAA{end}{jj} = '  (??)';
%             end
%             
%             % Fix time vector
%             file_data{2,end}(:,1) = file_data{2,end}(:,1)/1e6;                    % Convert from us to s
%             varlabAA{end}{1} = '  Time';
%             varunitsAA{end}{1} = '  (s)';
%             
%             t_start = min(t_start,file_data{2,end}(1,1));
            
        end
    end
end
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

