function fill_fdata_Callback(hObject, ~)

% import package contents
cellfun(@import, KSID.import_package());



handles = guidata(hObject);

[fds, name] = get_current_fds(handles);

if isfield(fds, 'aircraftData')
    h = fds.aircraftData;
else
    errordlg('aircraftData not present :(')
    return
end

% add tree node
[fds, parentNode] = fdsAddTreeBranch(fds, 0, 'Post Flight Analysis');

% get names and units
[~, varNames, varUnits] = fdata_init_si_5_0;

% fill the array
[fdata0, fil] = fdataUAVmainframe(fds);

%45-52 geo
fdata0(:,45:52) = [fil*h.xCG, fil*h.yCG, fil*h.zCG, fil*h.mass, fil*h.ixx, fil*h.iyy, fil*h.izz, fil*h.ixz];

%77-79 geo
fdata0(:,77:79) = [fil*h.sRef, fil*h.bRef,fil*h.cRef];

% Computed quantities should always be assigned
% as fdata may have been filled with different data,
% e.g: raw->EKF, modifying EKF, etc.

fdata0 = compfmc_si(fdata0);

% alphadot, betadot
% fdata0(:,56) = deriv(fdata0(:,4),1/fds.SampleRate);
% fdata0(:,57) = deriv(fdata0(:,3),1/fds.SampleRate);

% Complete fds struct
fds = fdsAddTreeLeaf(fds, 'SIDPAC_fdata', varNames, varNames, varUnits, {}, fdata0, parentNode, false);

updateDataSet(hObject, fds, name);

end

function [fdata0, fil] = fdataUAVmainframe(fds)

% import package contents
cellfun(@import, KSID.import_package());

%% Compile SIDPAC data matrix
source = 1;%get(handles.source_popup,'Value');

R2D = 1*180/pi;
%A2G = 1/9.80665;

DAT = fds_get_group(fds,'Base');
VS = fds_get_group(fds,'Real Time Solution');

r=length(DAT(:,1));

dummy = zeros(r,1);
fil = ones(r,1);

% fill SIDPAC array
fdata0 = zeros(r,90);

% 1 time
fdata0(:,1) = DAT(:,1);

if source == 1 % fill with raw data (for wind tunnel or playing)
    
    %2-4 airdata
    if norm(DAT(:,21)) ~= 0
        fdata0(:,2:4) = DAT(:,[19,21,20]);
    else
        disp('filling fdata airdata from vehicle state...')
        fdata0(:,2:4) = VS(:,[24,26,25]);
    end
    
    %5-13 IMU
    fdata0(:,5:13) = [DAT(:,5)*R2D, DAT(:,6)*R2D, DAT(:,7)*R2D, DAT(:,11), DAT(:,12), DAT(:,13), DAT(:,2), DAT(:,3), DAT(:,4)];
    
    %27-30 airdata
    fdata0(:,27:30) = [DAT(:,23), dummy, DAT(:,22), DAT(:,18)];
    
    %80-81
    fdata0(:,80:81) = [DAT(:,21), DAT(:,20)];
    
elseif source == 2 % fill with post-processed EKF data
    
    EKF = fds_get_group(fds,'EKF_Post');
    
    %2-4 airdata
    fdata0(:,2:4) = EKF(:,[24,26,25]);
    
    %5-13 IMU
    fdata0(:,5:13) = [EKF(:,14:16), EKF(:,11:13), EKF(:,8:10)];
    
    %27-30 airdata
    fdata0(:,27:30) = [EKF(:,27), dummy, DAT(:,22), EKF(:,23)];
    
    %74-76 body axes speed components
    fdata0(:,74:76) = EKF(:,5:7);
    
end

%14-20 controls
fdata0(:,14:16) = [DAT(:,39), DAT(:,41), DAT(:,40)];
fdata0(:,38)    = DAT(:,42);

CTRL = fds_get_group(fds,'Flight Control 1');
%21-26 cmd
fdata0(:,21:24) = [CTRL(:,2), CTRL(:,3), CTRL(:,4), CTRL(:,5)];
fdata0(:,34)    = CTRL(:,2);

% 82 index
fdata0(:,82) = CTRL(:,8);

end

