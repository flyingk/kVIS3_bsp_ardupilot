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

function [] = export_csv_all(hObject, ~)

% Choose the output directory
root_directory = uigetdir(pwd,'Export .csv Files');

[~, ~, noOfEntries] = kVIS_dataSetListState(hObject);

% Loop through and save each open data set
for ii = 1:noOfEntries
    
    % Set the current fds
    h = guidata(hObject);
    h.uiTabData.dataSetList.Value = ii;
    guidata(hObject);

    % Get the name
    [~, name, ~] = kVIS_dataSetListState(hObject);
    
    % Load the fds file
    fds = kVIS_getCurrentFds(hObject);
    
    % Exprot to csv
    folder_name = [root_directory,'/',name];
    export_csv(fds,folder_name);

end

return

