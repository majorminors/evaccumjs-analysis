%% move json into format for DMC LBA FIT
% Dorian Minors
% Created: NOC21
%
%
%
%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy
d = struct(); % set up a structure for the data we want to save
t = struct(); % set up a structure for temp data

% set up variables
rootdir = pwd; %% root directory - used to inform directory mappings
datadir = fullfile(rootdir,'data','behav_9'); % location of data
codeTableOnly = 1; % off to recode from converted json, on to just recode the table output

toolsdir = fullfile(rootdir, 'lib');
dataToProcess = 'converted_data'; % where is the converted data?
saveFileName = 'DMC_processed_data'; % what the name of the processed data is or should be

addpath(genpath(toolsdir)); % add libraries to path
p.save_file = fullfile(datadir, saveFileName);
p.table_save_file = fullfile(datadir,[saveFileName 'table.csv']);

if codeTableOnly; warning('using your saveFileName to load the data'); end

% grab simulated data for comparison
% t.pathToSimRTs = fullfile(toolsdir,'simulations','simulateRt','BoundaryRT.csv'); % 'AccrateRT.csv' or 'BoundaryRT.csv'
% t.simulatedRTData = readtable(t.pathToSimRTs);


if ~codeTableOnly
    
    %% cobble together what we need to play with the data
    theData = load(fullfile(datadir,dataToProcess)); % load the data
    t.alldata = theData.d.alldata; % let's just get what we really want and put it in the format the copied code below wants it in

    
    for subject = 1:length(t.alldata) % loop through each subject
        fprintf(1, 'working with subject %1.0f of %1.0f\n', subject, length(t.alldata)); % print that so you can check
        close all
        
        t.this_subj_data = t.alldata{subject};
        
        t.this_subj_trials = []; % init this so we can pull the trials out of each component
        % now we loop through the components
        for component = 1:numel(t.this_subj_data)
            %         1) JATOS set up: has app id and condition
            %         2) experiment set up with consent and demographics
            %         3) instructions
            %         4) coherence thresholdhing 1
            %         5) decision thresholding
            %         6) coherence thresholding 2
            %         7) experiment
            %         8) experiment finish (nothing really here)
            t.this_component = t.this_subj_data{component};
            
            % start off looking for prolific (or other) id
            if isfield(t.this_component,'app_identifier_string')
                if component == 1 % use the first component to make an id
                    t.id = t.this_component.app_identifier_string;
                end
                if isfield(t,'id')
                    if ~strcmp(t.id,t.this_component.app_identifier_string) % error if the ids aren't the same in every component
                        sprintf('t.id is %s and t.this_component.app_identifier_string is %s',t.id,t.this_component.app_identifier_string)
                        %error('you seem to have multiple app id strings (i.e. multiple participant results) confused')
                    end
                end
            end
            
            % lets get a code for button press
            if isfield(t.this_component,'condition')
                t.button_condition = t.this_component.condition;
                % then we need to get a keycode so make that too
                % condition 1 : respkeys o,p
                % condition 2 : respkeys p,o
                % keypress is JS, so 79 is o and 80 is p
                if t.button_condition{2} == 1
                    t.keycode = [79,80];
                elseif t.button_condition{2} == 2
                    t.keycode = [80,79];
                end
                % see if we can fix this?
                %             if isfield(t,'button_condition')
                %                 if t.button_condition{2} ~= t.this_component.condition{2} % check that the condition code numbers match for each component
                %                     sprintf('t.button_condition is %s and t.this_component.condition is %s',t.button_condition{2},t.this_component.condition{2})
                %                     error('you seem to have multiple button conditions (i.e. multiple participant results) confused') % else error
                %                 end
                %             end
            end
            
            % get values of stuff
            if isfield(t.this_component, 'rule_values')
                t.rule_values = t.this_component.rule_values;
            end
            %         if isfield(t.this_component, 'rule_easy_dots')
            %             t.easy_dots_rule_value = t.this_component.rule_easy_dots;
            %         end
            %         if isfield(t.this_component, 'rule_hard_dots')
            %             t.hard_dots_rule_value = t.this_component.rule_hard_dots;
            %         end
            if isfield(t.this_component, 'coherence_values')
                t.coherence_values = t.this_component.coherence_values;
            end
            if isfield(t.this_component, 'updated_coherence_values')
                t.updated_coherence_values = t.this_component.updated_coherence_values;
            end
            
            % now because some weird interaction between jspsych and jatos and
            % components makes the trials get saved in a completely
            % unintelligible field, we loop through and put all the trials in a
            % 'row'
            fn = fieldnames(t.this_component);
            for thisField = 1:numel(fn)
                if( contains(fn{thisField},'x0x') )
                    t.this_subj_trials = [t.this_subj_trials,{t.this_component.(fn{thisField})}];
                end
            end; clear thisField fn
        end
        
        
        % data is all in a row, so we go through each col and pull the values we want
        
        t.coh_count = 0;
        t.rule_count = 0;
        t.coh_ang_count = 0;
        t.exp_count = 0;
        for trial = 1:length(t.this_subj_trials)
            clear t.current_trial;
            
            t.current_trial = t.this_subj_trials{trial};
            
            % pull the psignifit arrays that were used
            if isfield(t.current_trial, 'coh_data_array')
                t.coh.data_array = t.current_trial.coh_data_array.data_array;
            end
            if isfield(t.current_trial, 'easy_dots_rule_data_array')
                t.rule.data_array_easy = t.current_trial.easy_dots_rule_data_array.data_array;
            end
            if isfield(t.current_trial, 'hard_dots_rule_data_array')
                t.rule.data_array_hard = t.current_trial.hard_dots_rule_data_array.data_array;
            end
            if isfield(t.current_trial, 'coh_data_array_updated')
                t.coh_ang.data_array = t.current_trial.coh_data_array_updated.data_array;
            end
            
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
            elseif isfield(t.current_trial, 'coherence_angle_array')
                t.stim_array = t.current_trial.coherence_angle_array;
                t.coh_ang.stim_array = [];
                for i = 1:length(t.stim_array)
                    t.coh_ang.stim_array = [t.coh_ang.stim_array,t.stim_array{i}];
                end; clear i;
            elseif  isfield(t.current_trial, 'exp_stim_array')
                t.stim_array = t.current_trial.exp_stim_array;
                t.exp.stim_array = [];
                for i = 1:length(t.stim_array)
                    t.exp.stim_array = [t.exp.stim_array,t.stim_array{i}];
                end; clear i;
                % so now we need to do something to sort out the erroneous match
                % difficulty sorting
                for i = 1:length(t.exp.stim_array)
                    t.match_distances(i) = t.exp.stim_array{1,i}.match_dist_cue_dir;
                end
                t.match_distances_rounded = unique(round(t.match_distances,1));
                if length(t.match_distances_rounded) ~= length(unique(t.match_distances))
                    warning('javascript has spawned multiple match distances - your matching difficulty coding will be wrong, because it was coded using min/max. recoding now')
                    disp('unique match distances:')
                    disp(unique(t.match_distances))
                    disp('rounded match distances:')
                    disp(t.match_distances_rounded)
                    for i = 1:length(t.exp.stim_array)
                        if round( t.exp.stim_array{1,i}.match_dist_cue_dir,1) == min(t.match_distances_rounded) || round( t.exp.stim_array{1,i}.match_dist_cue_dir,1) == max(t.match_distances_rounded)
                            t.exp.stim_array{1,i}.match_difficulty = 1;
                        else
                            t.exp.stim_array{1,i}.match_difficulty = 2;
                        end
                    end
                end
            end % end stimulus array sorter
            
            
            
            
            if isfield(t.current_trial, 'rule_values')
                t.rule_values = t.current_trial.rule_values;
            end
            %                 if isfield(t.current_trial, 'rule_easy_dots')
            %                     t.easy_dots_rule_value = t.current_trial.rule_easy_dots;
            %                 end
            %                 if isfield(t.current_trial, 'rule_hard_dots')
            %                     t.hard_dots_rule_value = t.current_trial.rule_hard_dots;
            %                 end
            if isfield(t.current_trial, 'coherence_values')
                t.coherence_values = t.current_trial.coherence_values;
            end
            if isfield(t.current_trial, 'updated_coherence_values')
                t.updated_coherence_values = t.current_trial.updated_coherence_values;
            end
            
            
            
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
                    
                elseif strcmp(t.current_trial.experiment_part, 'ruletest_rdk_easy') || strcmp(t.current_trial.experiment_part, 'ruletest_rdk_hard')
                    % just get the rdk trials for rule test trials
                    t.rule_count = t.rule_count+1;
                    
                    if strcmp(t.current_trial.experiment_part, 'ruletest_rdk_easy')
                        t.rule.coherence(t.rule_count) = 1;
                    elseif strcmp(t.current_trial.experiment_part, 'ruletest_rdk_hard')
                        t.rule.coherence(t.rule_count) = 2;
                    end
                    
                    t.rule.rt(t.rule_count) = t.current_trial.rt;
                    if isempty(find(t.current_trial.key_press == t.keycode))
                        t.rule.button(t.rule_count) = -1;
                    else
                        t.rule.button(t.rule_count) = find(t.current_trial.key_press == t.keycode);
                    end
                    t.rule.correct(t.rule_count) = t.current_trial.correct;
                    t.rule.direction(t.rule_count) = t.current_trial.coherent_direction;
                    
                elseif strcmp(t.current_trial.experiment_part, 'cohtest_angle_rdk')
                    % just get the rdk trials for coherence test trials
                    t.coh_ang_count = t.coh_ang_count+1;
                    
                    t.coh_ang.rt(t.coh_ang_count) = t.current_trial.rt;
                    if isempty(find(t.current_trial.key_press == t.keycode))
                        t.coh_ang.button(t.coh_ang_count) = -1;
                    else
                        t.coh_ang.button(t.coh_ang_count) = find(t.current_trial.key_press == t.keycode);
                    end
                    t.coh_ang.correct(t.coh_ang_count) = t.current_trial.correct;
                    t.coh_ang.direction(t.coh_ang_count) = t.current_trial.coherent_direction;
                    
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
        disp('finished coding trials')
        
        
        
        %% collate results
        
        disp('collating results')
        % collate results
        d.subjects(subject).id = t.id;
        d.subjects(subject).coh = t.coh;
        d.subjects(subject).rule = t.rule;
        d.subjects(subject).coh_ang = t.coh_ang;
        d.subjects(subject).exp = t.exp;
        d.subjects(subject).rule_values = t.rule_values;
        %     d.subjects(subject).hard_dots_rule_value = t.hard_dots_rule_value;
        %     d.subjects(subject).easy_dots_rule_value = t.easy_dots_rule_value;
        d.subjects(subject).coherence_values = t.coherence_values;
        d.subjects(subject).updated_coherence_values = t.updated_coherence_values;

        
    end; clear subject; % end subject loop
    
    fprintf('saving matlab vars as %s\n', p.save_file);
    save(p.save_file,'p','d'); % save all data to a .mat file
end

if codeTableOnly; load([p.save_file '.mat']); end
%% create output table
idx=0; subject=numel(d.subjects);
while subject
    thisSubjectData =  d.subjects(subject).exp;
    
    trial = numel(thisSubjectData.rt);
    
    while trial
        idx=idx+1;
        
        Var1(idx)=idx;
        s(idx)=subject;
        % we don't want this kind of granularity yet
        %         if thisSubjectData.stim_array{trial}.cue_dir == 1 || thisSubjectData.stim_array{trial}.cue_dir == 3
        %             S{idx}='s1';
        %         elseif thisSubjectData.stim_array{trial}.cue_dir == 2 || thisSubjectData.stim_array{trial}.cue_dir == 4
        %             S{idx}='s2';
        %         end
        if thisSubjectData.stim_array{trial}.coh_difficulty == 1
            Coh{idx}='cohE';
        elseif thisSubjectData.stim_array{trial}.coh_difficulty == 2
            Coh{idx}='cohH';
        end
        if thisSubjectData.stim_array{trial}.match_difficulty == 1
            Angle{idx}='angleE';
        elseif thisSubjectData.stim_array{trial}.match_difficulty == 2
            Angle{idx}='angleH';
        end
        if thisSubjectData.button(trial) == 1
            R{idx}='r1';
            % let's just code for correctness for our first modelling attempt
            if thisSubjectData.correct(trial) == 1
                S{idx} = 's1'; else; S{idx} = 's2';
            end
        elseif thisSubjectData.button(trial) == 2
            R{idx}='r2';
            % let's just code for correctness for our first modelling attempt
            if thisSubjectData.correct(trial) == 1
                S{idx} = 's2'; else; S{idx} = 's1';
            end
        end
        RT(idx)=thisSubjectData.rt(trial)/1000;
        
        trial = trial-1;
    end
    
    subject = subject-1;
end

Var1=Var1';s=s';S=S';Coh=Coh';Angle=Angle';R=R';RT=RT';
d.outputTable = table(Var1,s,S,Coh,Angle,R,RT);

fprintf('saving table as %s\n', p.table_save_file);
writetable(d.outputTable,p.table_save_file); % save to csv
