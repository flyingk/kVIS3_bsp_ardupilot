function Create_flight_analysis_group(hObject, ~)
clc

[fds, name] = kVIS_getCurrentFds(hObject);

[root_name] = kVIS_fdsGetTreeRootLabel(fds);



[~, nodeIndex] = kVIS_fdsGetGroup(fds,'Post Flight Analysis');

if nodeIndex > 0
    % overwrite existing data
    parentNode = nodeIndex;
else
    % add tree node
    [fds, parentNode] = kVIS_fdsAddTreeBranch(fds, root_name,'Post Flight Analysis');
end

% Rotational
[varNames, varUnits, varData] = assembleGroupAPM(fds, 'IMU_0','Time', 'Time', {}, {}, []);
[varNames, varUnits, varData] = assembleGroupAPM(fds, 'IMU_0','GyrX', 'p_Sensor', varNames, varUnits, varData);
[varNames, varUnits, varData] = assembleGroupAPM(fds, 'IMU_0','GyrY', 'q_Sensor', varNames, varUnits, varData);
[varNames, varUnits, varData] = assembleGroupAPM(fds, 'IMU_0','GyrZ', 'r_Sensor', varNames, varUnits, varData);
% [varNames, varUnits, varData] = assembleGroupAPM(fds, 'AHR2','Roll', 'Roll', varNames, varUnits, varData);
% [varNames, varUnits, varData] = assembleGroupAPM(fds, 'AHR2','Pitch', 'Pitch', varNames, varUnits, varData);
% [varNames, varUnits, varData] = assembleGroupAPM(fds, 'AHR2','Yaw', 'Yaw', varNames, varUnits, varData);

% Translational
[varNames, varUnits, varData] = assembleGroupAPM(fds, 'IMU_0','AccX', 'Ax_Sensor', varNames, varUnits, varData);
[varNames, varUnits, varData] = assembleGroupAPM(fds, 'IMU_0','AccY', 'Ay_Sensor', varNames, varUnits, varData);
[varNames, varUnits, varData] = assembleGroupAPM(fds, 'IMU_0','AccZ', 'Az_Sensor', varNames, varUnits, varData);
% u
% v
% w
% Pn
% Pe
% Pd

% Airdata
% Vair
% AoA
% AoS
% AoG
% Density
% Temp

% units = cell(size(varNames));
% varUnits = cellfun(@(x) '', units, 'UniformOutput', false);

frames = cell(size(varNames));
varFrames = cellfun(@(x) '', frames, 'UniformOutput', false);

fds = kVIS_fdsAddDataGroup(fds, parentNode, 'State', varNames, varUnits, varFrames, varData);

fds = kVIS_fdsAddDataGroup(fds, parentNode, 'Controls',[],[],[],[]);

fds = kVIS_fdsAddDataGroup(fds, parentNode, 'Forces/Moments',[],[],[],[]);

fds = kVIS_fdsAddDataGroup(fds, parentNode, 'Autopilot',[],[],[],[]);

kVIS_updateDataSet(hObject, fds, name);

end


function [varNames, varUnits, varData] = assembleGroupAPM(fds, groupName, channel, newName, varNames, varUnits, varData)

[signal, signalMeta] = kVIS_fdsGetChannel(fds, groupName, channel);

varData = [varData, signal];
% varNames = [varNames; signalMeta.name];
varNames = [varNames; newName];
varUnits = [varUnits; signalMeta.unit];


end