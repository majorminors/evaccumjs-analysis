%% convert evaccumjs data from json to matlab format
% Dorian Minors
% Created: JAN21
%
% this takes bloody ages with loadjson, but because I am using jspsych
% components now, loadjson is the only one that works
% if i was not using components i would use jsondecode which is way faster
%
%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy
d = struct(); % set up a structure for the data we want to save

p.file_control = ones(1,27); % for each file, 1 = one participant, 0 = multiple, -1 = unknown, -2 = skip

% set up variables
rootdir = pwd; %% root directory - used to inform directory mappings

datadir = fullfile(rootdir,'data','behav_9'); % location of data

saveFileName = 'converted_data'; % what to save the converted data as

dataFilePattern = 'jatos_results_*'; % file pattern of input data

% map tools
addpath(genpath(fullfile(rootdir, 'lib'))); % add libraries to path

p.save_file = fullfile(datadir, saveFileName);

%% loop through subjects
fileInfo = dir(fullfile(datadir, dataFilePattern)); % find all the datafiles and get their info
allData = {};
skipCount = 0;
for thisFile = 1:length(fileInfo) % loop through the files
    fprintf(1, 'working with file %1.0f of %1.0f\n', thisFile, length(fileInfo)); % print that so you can check

    thisPath = fullfile(datadir, fileInfo(thisFile).name); % get the full path to the file
    fprintf(1, 'file path: %s\n', thisPath); % print that so you can check
    
    disp('loading file')
    loadedData = loadjson(thisPath); % load in the data
    %loadedData = jsondecode(fileread(thisPath));
    loadedData = loadedData';
    disp('file loaded')
    
    if p.file_control(thisFile) == -1
        disp(loadedData)
        disp('specified unknown for file control')
        thisPrompt = 'check entries loaded (above): if one participant (1), if not (0), if skip (-2): enter a value and press enter to continue';
        p.file_control(thisFile) = input(thisPrompt,'s');
    end
    
    if p.file_control(thisFile) == -2
        disp('file skipped')
        skipCount=skipCount+1;
        continue % skip this file, moving on to the next iteration of the for loop
    end
    
    if ~p.file_control
        datasetController = [1,length(loadedData(1,:))]; % variable to control the while loop - left side indicates what participant dataset in the file we're up to, right side indicates number of datasets
        while (datasetController(1) <= datasetController(2))
            if isempty(loadedData{1,datasetController(1)}) % if this dataset has nothing in it (i.e. a participant who started but didn't finish)
                loadedData(:,datasetController(1)) = []; % clear it and
                datasetController(2) = datasetController(2)-1; % reduce the number of datasets by one
            end
            datasetController(1) = datasetController(1)+1; % move on to the next dataset
        end
    end
    
    if p.file_control
        allData{thisFile-skipCount} = loadedData; % load that subject into one cell
    else
        allData = [allData,loadedData]; % concat those into one var, so each subject is a cell
    end
    
end

disp('all files imported')

d.alldata = allData; % save all the data

fprintf('saving as %s\n', p.save_file);
save(p.save_file,'p','d'); % save all data to a .mat file