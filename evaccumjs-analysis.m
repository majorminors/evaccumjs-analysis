%% deal with evaccumjs data
% Dorian Minors
% Created: JAN21
%
%
%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy
d = struct(); % set up a structure for the data info
t = struct(); % set up a structure for temp data

% set up variables
rootdir = pwd; %% root directory - used to inform directory mappings
datadir = fullfile(rootdir,'data/behav_1');
p.savefilename = 'processed_data';
p.datafilepattern = 'jatos_results_*';

% directory mapping
addpath(genpath(fullfile(rootdir, 'tools'))); % add tools folder to path (don't think we need this, but in case)

save_file = fullfile(datadir, p.savefilename);

%% loop through subjects
d.fileinfo = dir(fullfile(datadir, p.datafilepattern)); % find all the datafiles and get their info
t.alldata = {};
t.skip_this_dataset = 0;
for file = 1:length(d.fileinfo)
    t.path = fullfile(datadir, d.fileinfo(file).name); % get the full path to the file
    fprintf(1, 'working with %s\n', t.path); % print that so you can check
    
    t.load = loadjson(t.path); % load in the data
    
    if length(t.load) > 1
        dataset = [1,length(t.load(1,:))];
            disp(dataset);
        while (dataset(1) <= dataset(2))
            if isempty(t.load{1,dataset(1)})
                t.load(:,dataset(1)) = [];
                dataset(2) = dataset(2)-1;
                disp(dataset);
            end
            dataset(1) = dataset(1)+1;
        end
    end           

    t.alldata = [t.alldata,t.load]; % concat those into one var, so each subject is a cell
    
end
d.alldata = t.alldata; % save all the data
    
for subject = 1:length(t.alldata) % loop through each subject
    fprintf(1, 'working with subject %f\n', subject); % print that so you can check

    t.this_subj_data = t.alldata{subject};
    
    
end

fprintf('saving output from %s\n', mfilename);
save(save_file,'d'); % save all data to a .mat file

function accuracy = accthis(data)
    accuracy = sum(data)/length(data);
end
