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
addpath(genpath(fullfile(rootdir, 'lib'))); % add libraries to path

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
    fprintf(1, 'working with subject %f of %f\n', subject, length(t.alldata)); % print that so you can check
    
    t.this_subj_data = t.alldata{subject};
    
    t.id = t.this_subj_data{1}.unique_id; % since we generated unique ids we'll pull these in
    
    % data is all in a row, so we go through each col and pull the values we want
    t.coh_count = 0;
    t.rule_count = 0;
    t.exp_count = 0;
    for trial = 1:length(t.this_subj_data)
        clear t.current_trial;
        
        t.current_trial = t.this_subj_data{trial};
        
        % pull stimulus arrays from start screen trials
        if isfield(t.current_trial, 'coh_stim_array')
            t.stim_array = t.current_trial.coh_stim_array;
            t.coh.stim_array = [];
            for i = 1:length(t.stim_array)
                t.coh.stim_array = [t.coh.stim_array,t.stim_array{i}];
            end; clear i;
        elseif isfield(t.current_trial, 'rule_stim_array')
            t.stim_array = t.current_trial.rule_stim_array;
            t.rule.stim_array = [];
            for i = 1:length(t.stim_array)
                t.rule.stim_array = [t.rule.stim_array,t.stim_array{i}];
            end; clear i;
        elseif  isfield(t.current_trial, 'exp_stim_array')
            t.stim_array = t.current_trial.exp_stim_array;
            t.exp.stim_array = [];
            for i = 1:length(t.stim_array)
                t.exp.stim_array = [t.exp.stim_array,t.stim_array{i}];
            end; clear i;
        end % end stimulus array sorter
        
        % only deal with trials that are labelled with experiment part
        if isfield(t.current_trial, 'experiment_part')
            
            % filter rdk trials
            if strcmp(t.current_trial.experiment_part, 'cohtest_rdk')
                % just get the rdk trials for coherence test trials
                t.coh_count = t.coh_count+1;
                
                t.coh.rt(t.coh_count) = t.current_trial.rt;
                t.coh.correct(t.coh_count) = t.current_trial.correct;
                t.coh.direction(t.coh_count) = t.current_trial.coherent_direction;
                
            elseif strcmp(t.current_trial.experiment_part, 'ruletest_rdk')
                % just get the rdk trials for rule test trials
                t.rule_count = t.rule_count+1;
                
                t.rule.rt(t.rule_count) = t.current_trial.rt;
                t.rule.correct(t.rule_count) = t.current_trial.correct;
                t.rule.direction(t.rule_count) = t.current_trial.coherent_direction;
                
            elseif strcmp(t.current_trial.experiment_part, 'experiment_rdk')
                % just get the rdk trials for experimental trials
                t.exp_count = t.exp_count+1;
                
                t.exp.rt(t.exp_count) = t.current_trial.rt;
                t.exp.correct(t.exp_count) = t.current_trial.correct;
                t.exp.direction(t.exp_count) = t.current_trial.coherent_direction;
                
            end % end rdk trial filtering
            
        end % end experiment part filtering
        
    end % end trial loop
    
    % collate results
    d.subjects(subject).id = t.id;
    d.subjects(subject).coh = t.coh;
    d.subjects(subject).rule = t.rule;
    d.subjects(subject).exp = t.exp;
    
    disp('*coherence*')
    disp('accuracy')
    disp(accthis(d.subjects(subject).coh.correct))
    disp('num invalid')
    disp(length(find(d.subjects(subject).coh.rt == -1)))
    disp('num fast responses')
    disp(length(find(d.subjects(subject).coh.rt >=0 & d.subjects(subject).coh.rt < 500)))
    disp('*rule*')
    disp('accuracy')
    disp(accthis(d.subjects(subject).rule.correct))
    disp('num invalid')
    disp(length(find(d.subjects(subject).rule.rt == -1)))
    disp('num fast responses')
    disp(length(find(d.subjects(subject).rule.rt >=0 & d.subjects(subject).rule.rt < 500)))
    disp('*experiment*')
    disp('accuracy')
    disp(accthis(d.subjects(subject).exp.correct))
    disp('num invalid')
    disp(length(find(d.subjects(subject).exp.rt == -1)))
    disp('num fast responses')
    disp(length(find(d.subjects(subject).exp.rt >=0 & d.subjects(subject).exp.rt < 500)))
    
end % end subject loop

fprintf('saving output from %s\n', mfilename);
save(save_file,'d'); % save all data to a .mat file

function accuracy = accthis(data)
accuracy = sum(data)/length(data);
end
