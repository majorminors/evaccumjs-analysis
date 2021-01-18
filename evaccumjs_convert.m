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
datadir = fullfile(rootdir,'data/behav_1'); % location of data
p.savefilename = 'processed_data'; % savefile for all data
p.datafilepattern = 'jatos_results_*'; % file pattern of input data

% directory mapping
addpath(genpath(fullfile(rootdir, 'lib'))); % add libraries to path

p.save_file = fullfile(datadir, p.savefilename);

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
    fprintf(1, 'working with subject %1.0f of %1.0f\n', subject, length(t.alldata)); % print that so you can check
    
    t.this_subj_data = t.alldata{subject};
    
    t.id = t.this_subj_data{1}.unique_id; % since we generated unique ids we'll pull these in
    
    % lets get a code for button press
    t.button_condition = t.this_subj_data{1}.condition;
    % condition 1 : respkeys o,p
    % condition 2 : respkeys p,o
    % keypress is JS, so 79 is o and 80 is p
    if t.button_condition{2} == 1
        t.keycode = [79,80];
    elseif t.button_condition{2} == 2
        t.keycode = [80,79];
    end
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
                if isempty(find(t.current_trial.key_press == t.keycode))
                    t.coh.button(t.coh_count) = -1;
                else
                    t.coh.button(t.coh_count) = find(t.current_trial.key_press == t.keycode);
                end
                t.coh.correct(t.coh_count) = t.current_trial.correct;
                t.coh.direction(t.coh_count) = t.current_trial.coherent_direction;
                
            elseif strcmp(t.current_trial.experiment_part, 'ruletest_rdk')
                % just get the rdk trials for rule test trials
                t.rule_count = t.rule_count+1;
                
                t.rule.rt(t.rule_count) = t.current_trial.rt;
                if isempty(find(t.current_trial.key_press == t.keycode))
                    t.rule.button(t.rule_count) = -1;
                else
                    t.rule.button(t.rule_count) = find(t.current_trial.key_press == t.keycode);
                end
                t.rule.correct(t.rule_count) = t.current_trial.correct;
                t.rule.direction(t.rule_count) = t.current_trial.coherent_direction;
                
            elseif strcmp(t.current_trial.experiment_part, 'experiment_rdk')
                % just get the rdk trials for experimental trials
                t.exp_count = t.exp_count+1;
                
                t.exp.rt(t.exp_count) = t.current_trial.rt;
                if isempty(find(t.current_trial.key_press == t.keycode))
                    t.exp.button(t.exp_count) = -1;
                else
                    t.exp.button(t.exp_count) = find(t.current_trial.key_press == t.keycode);
                end
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
    
    t.prompt = 'Continue to process for LBA fit with this participant? y/n [y]: ';
    t.do_lba = input(t.prompt,'s');
    if isempty(t.do_lba); t.do_lba = 'y'; end
    
    d.subjects(subject).lba = t.do_lba;
    
end; clear subject; % end subject loop for initial checking

fprintf('saving output from %s\n', mfilename);
save(p.save_file,'d'); % save all data to a .mat file

fprintf('converting selected participants for lba\n');

ilba = 0; % initialise a counter
p.conditions = {'EcEr','EcHr','HcEr','HcHr'}; % 2x2 coherence and rule
p.conditioncodes = {1,2,3,4};
t.lbadata = {}; % init this
for subject = 1:length(d.subjects) % loop through subjects
    if d.subjects(subject).lba == 'y' % if subject has been approved for lba
        
        ilba = ilba+1;
        
        for trial = 1:length(d.subjects(subject).exp.rt)
            
            % create a row that gives you a number for each condition in your 2x2
            if d.subjects(subject).exp.stim_array{1,trial}.coh_difficulty == 1 && d.subjects(subject).exp.stim_array{1,trial}.match_difficulty == 1
                t.condition(trial,1) = string(p.conditions{1});
                t.conditioncode(trial,1) = p.conditioncodes{1};
            elseif d.subjects(subject).exp.stim_array{1,trial}.coh_difficulty == 1 && d.subjects(subject).exp.stim_array{1,trial}.match_difficulty == 2
                t.condition(trial,1) = string(p.conditions{2});
                t.conditioncode(trial,1) = p.conditioncodes{2};
            elseif d.subjects(subject).exp.stim_array{1,trial}.coh_difficulty == 2 && d.subjects(subject).exp.stim_array{1,trial}.match_difficulty == 1
                t.condition(trial,1) = string(p.conditions{3});
                t.conditioncode(trial,1) = p.conditioncodes{3};
            elseif d.subjects(subject).exp.stim_array{1,trial}.coh_difficulty == 2 && d.subjects(subject).exp.stim_array{1,trial}.match_difficulty == 2
                t.condition(trial,1) = string(p.conditions{4});
                t.conditioncode(trial,1) = p.conditioncodes{4};
            end
            
            % get trial code
            t.trialtype(trial,1) = d.subjects(subject).exp.stim_array{1,trial}.trial_cond_num;
            
        end; clear trial;
            
        
          % consolidate all that data
          t.consolidata(:,1) = num2cell(t.condition);
          t.consolidata(:,2) = num2cell(t.conditioncode);
          t.consolidata(:,3) = num2cell(d.subjects(subject).exp.button');
          t.consolidata(:,4) = num2cell(d.subjects(subject).exp.rt');
          t.consolidata(:,5) = num2cell(d.subjects(subject).exp.correct');
          t.consolidata(:,6) = num2cell(t.trialtype);

          % stack it up
         d.subjects(subject).lba = t.consolidata;
         d.lbadata{ilba} = t.consolidata;
        
    end % end lba approved if statement
    
end; clear subject ilba; % end subject loop for lba

fprintf('saving lba adjusted output from %s\n', mfilename);
save(p.save_file,'d'); % save all data to a .mat file
 

function accuracy = accthis(data)
accuracy = sum(data)/length(data);
end
