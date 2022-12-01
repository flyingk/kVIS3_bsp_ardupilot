% Test Script to Make Sure Things are Working

close all
clear all
clc

addpath('./apm_converter');

data_dir = './../Sample_Data';

files = dir([data_dir,'\**\*.bin']);

for ii = 1:numel(files)
    filename = fullfile(files(ii).folder,files(ii).name);
    fprintf('Importing %s\n',filename);

    fds = Ardupilog(fullfile(files(ii).folder,files(ii).name)).getStruct();

end