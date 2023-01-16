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

function [] = export_csv(fds,folder_name)

% Exports the fds
mkdir(folder_name);

fprintf('Export csv files to %s ... ',folder_name);

% Loop through each of the sensor files in fdata and export them as a csv
for ii = 2:size(fds.fdata,2)
    % Initialise the table
    T = table();
        try
    % Read in the data
    sensor = fds.fdata{1,ii};
    varNames = fds.fdata{2,ii};
    units = fds.fdata{3,ii};
    data_out = fds.fdata{7,ii};
    
    % Fill out hte table
    variable_names = cell(numel(varNames),1);
    for jj = 1:numel(varNames)
        variable_names{jj} = [varNames{jj}, ' (', units{jj}, ')'];
        T.(varNames{jj}) = data_out(:,jj);
    end


    
    T.Properties.VariableUnits = units;
%     T.Properties.VariableNames = variable_names; % it doesn't like this
    
    % Write the table to a csv
    writetable(T,[folder_name,'/',sensor,'.csv'],'Delimiter',',') 
    catch
        fprintf('Skipping %s\n',sensor);
    end

end

% All done
fprintf('Done!\n');

return;
