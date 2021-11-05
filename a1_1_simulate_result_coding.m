%% simulate results from converted data and check for suitability in LBA fit
% Dorian Minors
% Created: OCT21
%
% we will use some existing dataset
% since we don't care about RT, we won't worry about trying to match these
% we'll just simulate a curve using the default function of psychcurve
% and then select correct/incorrect based on the probability
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

toolsdir = fullfile(rootdir, 'lib');
datadir = fullfile(rootdir,'data','behav_9_simulation_accumulation'); % location of data
dataToProcess = 'converted_data'; % where is the converted data?
saveFileName = 'processed_data'; % what to save the processed data as

t.pathToSimRTs = fullfile(toolsdir,'simulations','simulateRt','AccrateRT.csv'); % 'AccrateRT.csv' or 'BoundaryRT.csv'


p.do_checks = 0;
p.plot_norms = 0;
p.skip_check_pp = 1;
p.plot_coh = 1;
p.check_coh = 1;
p.print_coh = 1; % will print values from the coherence thresholding rather than asking you to look online
p.plot_match = 1;
p.check_match = 1;
p.print_match = 1; % will print values from the coherence thresholding rather than asking you to look online
p.plot_coh_ang = 1;
p.check_coh_ang = 1;
p.print_coh_ang = 1; % will print values from the coherence thresholding rather than asking you to look online
p.skip_check_lba = 1;

% some quick checks
if p.check_coh == 1 && p.plot_coh == 0; error('you have to plot coherence to check coherence'); end
if p.check_coh_ang == 1 && p.plot_coh_ang == 0; error('you have to plot coherence thresholding 2 to check it'); end
if p.check_match == 1 && p.plot_match == 0; error('you have to plot matching to check matching'); end
if p.check_coh == 0 && p.print_coh == 1; error('you have to check coherence to print it'); end
if p.check_coh_ang == 0 && p.print_coh_ang == 1; error('you have to check coherence 2 to print it'); end
if p.check_match == 0 && p.print_match == 1; error('you have to check matching to print it'); end

%% grab simulated RTs
t.simulatedRTData = readtable(t.pathToSimRTs);
t.simulatedRTData = table2cell(t.simulatedRTData);
% since we don't know how to scale the right values in the simulation
% script, we'll scale them here
t.oldMax = max(cell2mat(t.simulatedRTData(:,7))); % grab the max value
t.oldMin = min(cell2mat(t.simulatedRTData(:,7))); % grab the min value
t.newMax = 0.4; % whats the lower bound of the scale you want?
t.newMin = 1.5; % what's the upper bound?
% we'll use these to scale value by value when we plug these in
simIdx = size(t.simulatedRTData,1); t.simEE=[];t.simEH=[];t.simHE=[];t.simHH=[];
while simIdx
    if strcmp(t.simulatedRTData{simIdx,4},'cohE') && strcmp(t.simulatedRTData{simIdx,5},'angleE')
        t.simEE = [t.simEE,t.simulatedRTData{simIdx,7}];
    elseif strcmp(t.simulatedRTData{simIdx,4},'cohE') && strcmp(t.simulatedRTData{simIdx,5},'angleH')
        t.simEH = [t.simEH,t.simulatedRTData{simIdx,7}];
    elseif strcmp(t.simulatedRTData{simIdx,4},'cohH') && strcmp(t.simulatedRTData{simIdx,5},'angleE')
        t.simHE = [t.simHE,t.simulatedRTData{simIdx,7}];
    elseif strcmp(t.simulatedRTData{simIdx,4},'cohH') && strcmp(t.simulatedRTData{simIdx,5},'angleH')
        t.simHH = [t.simHH,t.simulatedRTData{simIdx,7}];
    end
    simIdx = simIdx-1;
end

% cobble together what we need to play with the data and save it
theData = load(fullfile(datadir,dataToProcess)); % load the data
t.alldata = theData.d.alldata; % let's just get what we really want and put it in the format the copied code below wants it in
addpath(genpath(toolsdir)); % add libraries to path
figdir = fullfile(datadir,'figures'); % place to save figures
if ~exist(figdir,'dir')
    mkdir(figdir);
end
p.save_file = fullfile(datadir, saveFileName);

for subject = 1:5%length(t.alldata) % loop through each subject
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
    
    %% generate and insert simulated accuracies and RTs
    
    simulatedPoints = simulatecurve; % simulate a sigmoid use
    for i = 1:length(t.coh.correct)
        t.coh.correct(i) = idealObserver(t.coh.stim_array{1,i}.coherence_value,simulatedPoints);
    end
    for i = 1:length(t.coh_ang.correct)
        t.coh_ang.correct(i) = idealObserver(t.coh_ang.stim_array{1,i}.coherence_value,simulatedPoints);
    end
    for i = 1:length(t.rule.correct)
        % for match, we want to feed idealobserver `abs(angle./100-1)` to
        % invert value and convert to sub-integers
        t.rule.correct(i) = idealObserver(abs(t.rule.stim_array{1,i}.rule_value/100-1),simulatedPoints);
    end
    for i = 1:length(t.exp.correct)
        if t.exp.button(i) ~= -1 % if response was not invalid (invalid responses need to be stripped, or later scripts wont work (or we'd have to work out what was the actual button that should have been pressed for the lba))
            if t.exp.stim_array{1,i}.coh_difficulty == 1 && t.exp.stim_array{1,i}.match_difficulty == 1
                thisIntensity = 0.60;
                thisRT = datasample(t.simEE,1);
            elseif t.exp.stim_array{1,i}.coh_difficulty == 1 && t.exp.stim_array{1,i}.match_difficulty == 2
                thisIntensity = 0.50;
                thisRT = datasample(t.simEH,1);
            elseif t.exp.stim_array{1,i}.coh_difficulty == 2 && t.exp.stim_array{1,i}.match_difficulty == 1
                thisIntensity = 0.40;
                thisRT = datasample(t.simHE,1);
            elseif t.exp.stim_array{1,i}.coh_difficulty == 2 && t.exp.stim_array{1,i}.match_difficulty == 2
                thisIntensity = 0.30;
                thisRT = datasample(t.simHH,1);
            end
            t.exp.correct(i) = idealObserver(thisIntensity,simulatedPoints);
            
            % move thisRT into the right range of values
            t.oldRTPercentile = (thisRT-t.oldMin)/(t.oldMax-t.oldMin);
            thisRT = ((t.newMax-t.newMin)*t.oldRTPercentile)+t.newMin;
            t.exp.rt(i) = thisRT*1000; % then add it, but scale from s to ms
        end
    end
    
    
    %% move on to checks
    
    disp('collating results and checking')
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
    
    if p.do_checks
        %     disp('*coherence thresholding one*')
        %     disp('accuracy (low coh/hard threshold is .9, high coh/easy threshold is .7)')
        %     disp(accthis(d.subjects(subject).coh.correct))
        %     disp('% invalid')
        %     disp((length(find(d.subjects(subject).coh.rt == -1))/length(d.subjects(subject).coh.rt)))
        %     disp('% fast responses')
        %     disp((length(find(d.subjects(subject).coh.rt >=0 & d.subjects(subject).coh.rt < 400))/length(d.subjects(subject).coh.rt)))
        if p.plot_norms
            figure; normplot(d.subjects(subject).coh.rt)
            title('coherence thresholding')
            export_fig(fullfile(figdir,strcat(num2str(subject),'_coh_normplot.jpeg')),'-transparent')
        end
        %     disp('*rule*')
        %     disp('accuracy (hard threshold is .6)')
        %     disp(accthis(d.subjects(subject).rule.correct))
        %     disp('% invalid')
        %     disp((length(find(d.subjects(subject).rule.rt == -1))/length(d.subjects(subject).rule.rt)))
        %     disp('% fast responses')
        %     disp((length(find(d.subjects(subject).rule.rt >=0 & d.subjects(subject).rule.rt < 400))/length(d.subjects(subject).rule.rt)))
        if p.plot_norms
            figure; normplot(d.subjects(subject).rule.rt)
            title('match thresholding')
            export_fig(fullfile(figdir,strcat(num2str(subject),'_match_normplot.jpeg')),'-transparent')
        end
        
        %     disp('*coherence second thresholding*')
        %     disp('accuracy (low coh/hard threshold is .9, high coh/easy threshold is .7)')
        %     disp(accthis(d.subjects(subject).coh_ang.correct))
        %     disp('% invalid')
        %     disp((length(find(d.subjects(subject).coh_ang.rt == -1))/length(d.subjects(subject).coh_ang.rt)))
        %     disp('% fast responses')
        %     disp((length(find(d.subjects(subject).coh_ang.rt >=0 & d.subjects(subject).coh_ang.rt < 400))/length(d.subjects(subject).coh_ang.rt)))
        if p.plot_norms
            figure; normplot(d.subjects(subject).coh_ang.rt)
            title('coherence second thresholding')
            export_fig(fullfile(figdir,strcat(num2str(subject),'_coh_ang_normplot.jpeg')),'-transparent')
        end
        %     disp('*experiment*')
        %     disp('accuracy (.6 is a floor effect)')
        %     disp(accthis(d.subjects(subject).exp.correct))
        %     disp('$ invalid')
        %     disp((length(find(d.subjects(subject).exp.rt == -1))/length(d.subjects(subject).exp.rt)))
        %     disp('% fast responses')
        %     disp((length(find(d.subjects(subject).exp.rt >=0 & d.subjects(subject).exp.rt < 400))/length(d.subjects(subject).exp.rt)))
        if p.plot_norms
            figure; normplot(d.subjects(subject).exp.rt)
            title('experiment')
            export_fig(fullfile(figdir,strcat(num2str(subject),'_exp_normplot.jpeg')),'-transparent')
        end
        
        % put together a table for quick checking
        tmpTableTitles = {'exp_part' 'accuracy' 'percent_invalid' 'percent_fast'};
        tmpTable = {...
            'coh thresh 1 (70-90)' ...
            accthis(d.subjects(subject).coh.correct) ...
            (length(find(d.subjects(subject).coh.rt == -1))/length(d.subjects(subject).coh.rt))*100 ...
            (length(find(d.subjects(subject).coh.rt >=0 & d.subjects(subject).coh.rt < 400))/length(d.subjects(subject).coh.rt))*100;...
            'rule thresh (>60)' ...
            accthis(d.subjects(subject).rule.correct) ...
            (length(find(d.subjects(subject).rule.rt == -1))/length(d.subjects(subject).rule.rt))*100 ...
            (length(find(d.subjects(subject).rule.rt >=0 & d.subjects(subject).rule.rt < 400))/length(d.subjects(subject).rule.rt))*100;...
            'coh thresh 2 (70-90)' ...
            accthis(d.subjects(subject).coh_ang.correct) ...
            (length(find(d.subjects(subject).coh_ang.rt == -1))/length(d.subjects(subject).coh_ang.rt))*100 ...
            (length(find(d.subjects(subject).coh_ang.rt >=0 & d.subjects(subject).coh_ang.rt < 400))/length(d.subjects(subject).coh_ang.rt))*100;...
            'experiment (>60)' ...
            accthis(d.subjects(subject).exp.correct) ...
            (length(find(d.subjects(subject).exp.rt == -1))/length(d.subjects(subject).exp.rt))*100 ...
            (length(find(d.subjects(subject).exp.rt >=0 & d.subjects(subject).exp.rt < 400))/length(d.subjects(subject).exp.rt))*100;...
            };
        behavCheckTbl = cell2table(tmpTable);
        behavCheckTbl.Properties.VariableNames = tmpTableTitles;
        
        if p.skip_check_pp
            t.do_pp = 'y';
        else
            t.prompt = 'Continue to psychophys with this participant? y/n [y]: ';
            t.do_pp = input(t.prompt,'s');
            if isempty(t.do_pp); t.do_pp = 'y'; end
            close all
        end
        
        if t.do_pp == 'y'
            % check thresholds
            if p.plot_coh
                [t.coh_easy,t.coh_hard,~,t.psignifit_array] = coh_thresholding(d.subjects(subject).coh,figdir,'_1',subject,1);
                if p.check_coh
                    [t.jscoh_easy,t.jscoh_hard,~,t.jscoh_psignifit_array] = coh_thresholding(t.coh.data_array,figdir,'_thresholding_with_jsarray_1',subject,0,1);
                    tmp = psychcurve(t.coh.data_array,0.9);
                    t.correct_thresh_easy = tmp.Fit(1);
                    tmp = psychcurve(t.coh.data_array,0.7);
                    t.correct_thresh_hard = tmp.Fit(1); clear tmp;
                    % put together a table for quick checking
                    disp('check coherence threshold values')
                    tmpTableTitles = {'condition' 'coh_thresh_used' 'coh_thresh_matlab_jspsych_array_used' 'correctly_obtained_threshold' 'coh_thresh_matlab_array_generated' 'difference_with_same_array' 'difference_from_correct_threshold'};
                    tmpTable = {...
                        'easy' ...
                        t.coherence_values(1) ...
                        t.jscoh_easy ...
                        t.coh_easy ...
                        t.correct_thresh_easy ...
                        t.coherence_values(1)-t.jscoh_easy ...
                        t.correct_thresh_hard-t.jscoh_easy;...
                        'hard' ...
                        t.coherence_values(2) ...
                        t.jscoh_hard ...
                        t.coh_hard ...
                        t.correct_thresh_hard ...
                        t.coherence_values(2)-t.jscoh_hard ...
                        t.correct_thresh_hard-t.jscoh_hard;...
                        };
                    tmpCheckTbl = cell2table(tmpTable);
                    tmpCheckTbl.Properties.VariableNames = tmpTableTitles;
                    disp(tmpCheckTbl)
                    d.subjects(subject).coherence_vals_matlab_orig = [t.jscoh_easy,t.jscoh_hard];
                    d.subjects(subject).coherence_vals_matlab_correct = [t.correct_thresh_easy,t.correct_thresh_hard];
                    if ~p.print_coh
                        t.prompt = 'Press enter to continue';
                        input(t.prompt,'s');
                    else
                        figure;
                        vals = [d.subjects(subject).coherence_values; d.subjects(subject).coherence_vals_matlab_correct; d.subjects(subject).coherence_vals_matlab_orig];
                        bar(vals,'FaceColor',[0.0 0.502 0.502]);
                        set(gca,'XTickLabel',{'python' 'matlab correct' 'matlab orig'});
                        ylim([0, 1]);
                        hold on
                        for i = 1:size(vals,2)
                            plot([0 1],[1 1]*vals(i,1),'--g')
                            plot([0 1],[1 1]*vals(i,2),'--r')
                        end
                        hold off
                        export_fig(fullfile(figdir,strcat(num2str(subject),'_coh_python_v_matlab.jpeg')),'-transparent')
                    end
                    % check the arrays
                    disp('check arrays')
                    disp('array used')
                    t.coh.data_array
                    disp('array assembled by matlab')
                    t.psignifit_array
                    disp('array assembled by matlab from array used')
                    t.jscoh_psignifit_array
                    d.subjects(subject).coh_array_check = t.coh.data_array-t.jscoh_psignifit_array;
                    if ~p.print_coh
                        t.prompt = 'press enter to continue';
                        input(t.prompt,'s');
                    end
                end
            end
            if p.plot_match
                [t.match_easy,t.match_hard,t.match_summary,t.psignifit_array] = match_thresholding(d.subjects(subject).rule,figdir,p.save_file,subject,1);
                % threshold values (1) = easy coherence, (2) = hard coherence,
                % (3) = trials combined threshold, (4) = average threshold
                if p.check_match
                    tmp(:,:,1) = t.rule.data_array_easy;
                    tmp(:,:,2) = t.rule.data_array_hard;
                    [t.jsmatch_easy,t.jsmatch_hard,t.jsmatch_summary,t.jsmatch_psignifit_array] = match_thresholding(tmp,figdir,p.save_file,subject,0,1); clear tmp;
                    tmp2 = psychcurve(t.rule.data_array_easy,0.6);
                    t.correct_thresh_easy_dots = -tmp2.Fit(1);
                    tmp2 = psychcurve(t.rule.data_array_hard,0.6);
                    t.correct_thresh_hard_dots = -tmp2.Fit(1); clear tmp;
                    t.correct_thresh_ave(2) = (t.correct_thresh_easy_dots+t.correct_thresh_hard_dots)/2;
                    t.correct_thresh_ave(1) = 90-t.correct_thresh_ave(2);
                    % put together a table for quick checking
                    disp('check matching threshold values')
                    tmpTableTitles = {'condition' 'match_thresh_used' 'matlab_jspsych_array_easy_hard_combo_ave' 'matlab_array_generated_easy_hard_combo_ave' 'matlab_correct_thresh'};
                    tmpTable = {...
                        'easy' ...
                        t.rule_values(1) ...
                        [t.jsmatch_easy(1) t.jsmatch_easy(2) t.jsmatch_easy(3) t.jsmatch_easy(4)] ...
                        [t.match_easy(1) t.match_easy(2) t.match_easy(3) t.match_easy(4)] ...
                        t.correct_thresh_ave(1);...
                        'hard' ...
                        t.rule_values(2) ...
                        [t.jsmatch_hard(1) t.jsmatch_hard(2) t.jsmatch_hard(3) t.jsmatch_hard(4)] ...
                        [t.match_hard(1) t.match_hard(2) t.match_hard(3) t.match_hard(4)] ...
                        t.correct_thresh_ave(2);...
                        };
                    tmpCheckTbl = cell2table(tmpTable);
                    tmpCheckTbl.Properties.VariableNames = tmpTableTitles;
                    disp(tmpCheckTbl)
                    d.subjects(subject).rule_vals_matlab_orig = [t.jsmatch_easy(4),t.jsmatch_hard(4)];
                    d.subjects(subject).rule_vals_matlab_correct = t.correct_thresh_ave;
                    if ~p.print_match
                        t.prompt = 'Press enter to continue';
                        input(t.prompt,'s');
                    else
                        figure;
                        vals = [d.subjects(subject).rule_values; d.subjects(subject).rule_vals_matlab_correct; d.subjects(subject).rule_vals_matlab_orig];
                        bar(vals,'FaceColor',[0.0 0.502 0.502]);
                        set(gca,'XTickLabel',{'python' 'matlab correct' 'matlab orig'});
                        ylim([0, 90]);
                        hold on
                        for i = 1:size(vals,2)
                            plot([0 1],[1 1]*vals(i,1),'--g')
                            plot([0 1],[1 1]*vals(i,2),'--r')
                        end
                        hold off
                        export_fig(fullfile(figdir,strcat(num2str(subject),'_match_python_v_matlab.jpeg')),'-transparent')
                    end
                    % check the arrays
                    disp('check arrays')
                    disp('array used (easy; hard)')
                    t.rule.data_array_easy
                    t.rule.data_array_hard
                    disp('array assembled by matlab (easy; hard; combo)')
                    t.psignifit_array
                    disp('array assembled by matlab from array used (easy; hard; combo)')
                    t.jsmatch_psignifit_array
                    d.subjects(subject).match_array_check_easy = t.rule.data_array_easy-t.jsmatch_psignifit_array(:,:,1);
                    d.subjects(subject).match_array_check_hard = t.rule.data_array_hard-t.jsmatch_psignifit_array(:,:,2);
                    if ~p.print_match
                        t.prompt = 'press enter to continue';
                        input(t.prompt,'s');
                    end
                end
            end
            if p.plot_coh_ang
                [t.coh_easy,t.coh_hard,~,t.psignifit_array] = coh_thresholding(d.subjects(subject).coh_ang,figdir,'_2',subject,1);
                if p.check_coh_ang
                    [t.jscoh_easy,t.jscoh_hard,~,t.jscoh_psignifit_array] = coh_thresholding(t.coh_ang.data_array,figdir,'_thresholding_with_jsarray_2',subject,0,1);
                    tmp = psychcurve(t.coh_ang.data_array,0.9);
                    t.correct_thresh_easy = tmp.Fit(1);
                    tmp = psychcurve(t.coh_ang.data_array,0.7);
                    t.correct_thresh_hard = tmp.Fit(1); clear tmp;
                    % put together a table for quick checking
                    disp('check updated threshold values')
                    tmpTableTitles = {'condition' 'coh_thresh_used' 'coh_thresh_matlab_jspsych_array_used' 'correctly_obtained_threshold' 'coh_thresh_matlab_array_generated' 'difference_with_same_array' 'difference_from_correct_threshold'};
                    tmpTable = {...
                        'easy' ...
                        t.updated_coherence_values(1) ...
                        t.jscoh_easy ...
                        t.coh_easy ...
                        t.correct_thresh_easy ...
                        t.updated_coherence_values(1)-t.jscoh_easy ...
                        t.correct_thresh_hard-t.jscoh_easy;...
                        'hard' ...
                        t.updated_coherence_values(2) ...
                        t.jscoh_hard ...
                        t.coh_hard ...
                        t.correct_thresh_hard ...
                        t.updated_coherence_values(2)-t.jscoh_hard ...
                        t.correct_thresh_hard-t.jscoh_hard;...
                        };
                    tmpCheckTbl = cell2table(tmpTable);
                    tmpCheckTbl.Properties.VariableNames = tmpTableTitles;
                    disp(tmpCheckTbl)
                    d.subjects(subject).updated_coherence_vals_matlab_orig = [t.jscoh_easy,t.jscoh_hard];
                    d.subjects(subject).updated_coherence_vals_matlab_correct = [t.correct_thresh_easy,t.correct_thresh_hard];
                    if ~p.print_coh_ang
                        t.prompt = 'Press enter to continue';
                        input(t.prompt,'s');
                    else
                        figure;
                        vals = [d.subjects(subject).updated_coherence_values; d.subjects(subject).updated_coherence_vals_matlab_correct; d.subjects(subject).updated_coherence_vals_matlab_orig];
                        bar(vals,'FaceColor',[0.0 0.502 0.502]);
                        set(gca,'XTickLabel',{'python' 'matlab correct' 'matlab orig'});
                        ylim([0, 1]);
                        hold on
                        for i = 1:size(vals,2)
                            plot([0 1],[1 1]*vals(i,1),'--g')
                            plot([0 1],[1 1]*vals(i,2),'--r')
                        end
                        hold off
                        export_fig(fullfile(figdir,strcat(num2str(subject),'_coh_2_python_v_matlab.jpeg')),'-transparent')
                    end
                    % check the arrays
                    disp('check arrays')
                    disp('array used')
                    t.coh_ang.data_array
                    disp('array assembled by matlab')
                    t.psignifit_array
                    disp('array assembled by matlab from array used')
                    t.jscoh_psignifit_array
                    d.subjects(subject).updated_coh_array_check = t.coh_ang.data_array-t.jscoh_psignifit_array;
                    if ~p.print_coh_ang
                        t.prompt = 'press enter to continue';
                        input(t.prompt,'s');
                    end
                end
            end
            
            % put together a table for quick checking
            tmpTableTitles = {'exp_part' 'easy_threshold' 'hard_threshold'};
            tmpTable = {...
                'coh thresh 1 (percent coh)' ...
                t.coherence_values(1)*100 ...
                t.coherence_values(2)*100;...
                'match thresh (angle)' ...
                t.rule_values(1) ...
                t.rule_values(2);...
                'match thresh (mean rt)' ...
                mean(t.match_summary(4,:)) ...
                mean(t.match_summary(6,:));...
                'coh thresh 2 (percent coh)' ...
                t.updated_coherence_values(1)*100 ...
                t.updated_coherence_values(2)*100;...
                };
            psychphysCheckTbl = cell2table(tmpTable);
            psychphysCheckTbl.Properties.VariableNames = tmpTableTitles;
            
            fprintf(1,'participant id: %s\n',t.id);
            disp(behavCheckTbl)
            disp(psychphysCheckTbl)
            
            if p.skip_check_lba
                t.do_lba = 'y';
            else
                t.prompt = 'Continue to process for LBA fit with this participant? y/n [y]: ';
                t.do_lba = input(t.prompt,'s');
                if isempty(t.do_lba); t.do_lba = 'y'; end
                close all
            end
        else
            t.do_lba = 'n';
        end
    else
        t.do_lba = 'y';
    end
    d.subjects(subject).lba = t.do_lba;
    
end; clear subject; % end subject loop for initial checking

fprintf('saving as %s\n', p.save_file);
save(p.save_file,'p','d'); % save all data to a .mat file
