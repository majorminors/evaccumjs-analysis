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

p.file_control = [-2, -2, -2, 1]; % for each file, 1 = one participant, 0 = multiple, -1 = unknown, -2 = skip
p.plot_norms = 0;
p.skip_check_pp = 1;
p.plot_coh = 0;
p.check_coh = 0;
p.plot_match = 0;
p.check_match = 0;
p.plot_coh_ang = 0;
p.check_coh_ang = 0;
p.skip_check_lba = 0;
p.plot_rt_hist = 1;
p.plot_rts = 1;
p.plot_pc = 1;
p.skip_check_lbacont = 0;


% set up variables
rootdir = pwd; %% root directory - used to inform directory mappings
datadir = fullfile(rootdir,'data/behav_4'); % location of data
p.savefilename = 'processed_data'; % savefile for all data
figdir = fullfile(datadir,'figures'); % place to save figures
if ~exist(figdir,'dir')
    mkdir(figdir);
end
p.datafilepattern = 'jatos_results_*'; % file pattern of input data

% directory mapping
addpath(genpath(fullfile(rootdir, 'lib'))); % add libraries to path

if p.check_coh == 1 && p.plot_coh == 0; error('you have to plot coherence to check coherence'); end
if p.check_match == 1 && p.plot_match == 0; error('you have to plot matching to check matching'); end

p.save_file = fullfile(datadir, p.savefilename);

%% loop through subjects
d.fileinfo = dir(fullfile(datadir, p.datafilepattern)); % find all the datafiles and get their info
t.alldata = {};
t.skip_count = 0;
for file = 1:length(d.fileinfo) % loop through the files
    t.path = fullfile(datadir, d.fileinfo(file).name); % get the full path to the file
    fprintf(1, 'working with %s\n', t.path); % print that so you can check
    
    %t.load = loadjson(t.path); % load in the data
    t.load = jsondecode(fileread(t.path));
    t.load = t.load';
    if p.file_control(file) == -1
        disp(t.load)
        disp('specified unknown for file control')
        t.prompt = 'check entries loaded (above): if one participant (1), if not (0), if skip (-2): enter a value and press enter to continue';
        p.file_control(file) = input(t.prompt,'s');
    end
    
    if p.file_control(file) == -2
        disp('file skipped')
        t.skip_count=t.skip_count+1;
        continue % skip this file, moving on to the next iteration of the for loop
    end
        
    if ~p.file_control
        dataset = [1,length(t.load(1,:))]; % variable to control the while loop - left side indicates what participant dataset in the file we're up to, right side indicates number of datasets
        while (dataset(1) <= dataset(2))
            if isempty(t.load{1,dataset(1)}) % if this dataset has nothing in it (i.e. a participant who started but didn't finish)
                t.load(:,dataset(1)) = []; % clear it and
                dataset(2) = dataset(2)-1; % reduce the number of datasets by one
            end
            dataset(1) = dataset(1)+1; % move on to the next dataset
        end
    end
    
    if p.file_control
        t.alldata{file-t.skip_count} = t.load; % load that subject into one cell
    else
        t.alldata = [t.alldata,t.load]; % concat those into one var, so each subject is a cell
    end
    
end
d.alldata = t.alldata; % save all the data

for subject = 1:length(t.alldata) % loop through each subject
    fprintf(1, 'working with subject %1.0f of %1.0f\n', subject, length(t.alldata)); % print that so you can check
    
    t.this_subj_data = t.alldata{subject};
    
    warning('if you get an error on one of these next operations, you probably need to check whether you imported participants properly (e.g. turn on/off one participant in p.file_control)')
    t.id = t.this_subj_data{1}.unique_id; % since we generated unique ids we'll pull these in
    t.coh_psignifit_array = t.this_subj_data{1}.coh_data_array.data_array;
    t.rule_ecoh_psignifit_array = t.this_subj_data{1}.easy_dots_rule_data_array.data_array;
%     t.coh_psignifit_array_updated = t.this_subj_data{1}.coh_data_array_updated.data_array;
    t.rule_hcoh_psignifit_array = t.this_subj_data{1}.hard_dots_rule_data_array.data_array;
    
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
    t.coh_ang_count = 0;
    t.exp_count = 0;
    for trial = 1:length(t.this_subj_data)
        clear t.current_trial;
        
        t.current_trial = t.this_subj_data{trial};
        
        if isfield(t.current_trial, 'rule_values')
            t.rule_values = t.current_trial.rule_values;
        end
        if isfield(t.current_trial, 'rule_easy_dots')
            t.easy_dots_rule_value = t.current_trial.rule_easy_dots;
        end
        if isfield(t.current_trial, 'rule_hard_dots')
            t.hard_dots_rule_value = t.current_trial.rule_hard_dots;
        end
        if isfield(t.current_trial, 'coherence_values')
            t.coherence_values = t.current_trial.coherence_values;
        end
        if isfield(t.current_trial, 'updated_coherence_values')
            t.updated_coherence_values = t.current_trial.updated_coherence_values;
        end
        
        % pull stimulus arrays from start screen trials
        if isfield(t.current_trial, 'coh_stim_array')
            t.stim_array = t.current_trial.coh_stim_array;
            t.coh.stim_array = [];
            for i = 1:length(t.stim_array)
                t.coh.stim_array{i} = t.stim_array(i);
            end; clear i;
        elseif isfield(t.current_trial, 'rule_stim_array')
            t.stim_array = t.current_trial.rule_stim_array;
            t.rule.stim_array = [];
            for i = 1:length(t.stim_array(:,1))
                for ii = 1:length(t.stim_array)
                    tmp{ii} = t.stim_array(i,ii);
                end; clear ii;
                t.rule.stim_array = [t.rule.stim_array,tmp];
            end; clear i tmp;
        elseif isfield(t.current_trial, 'coherence_angle_array')
            t.stim_array = t.current_trial.coherence_angle_array;
            t.coh_ang.stim_array = [];
            for i = 1:length(t.stim_array)
                t.coh_ang.stim_array{i} = t.stim_array(i);
            end; clear i;
        elseif  isfield(t.current_trial, 'exp_stim_array')
            t.stim_array = t.current_trial.exp_stim_array;
            t.exp.stim_array = [];
            for i = 1:length(t.stim_array(:,1))
                for ii = 1:length(t.stim_array)
                    tmp{ii} = t.stim_array(i,ii);
                end; clear ii;
                t.exp.stim_array = [t.exp.stim_array,tmp];
            end; clear i tmp;
            
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
    
    % collate results
    d.subjects(subject).id = t.id;
    d.subjects(subject).coh = t.coh;
    d.subjects(subject).rule = t.rule;
    d.subjects(subject).coh_ang = t.coh_ang;
    d.subjects(subject).exp = t.exp;
%     d.subjects(subject).rule_values = t.rule_values;
%     d.subjects(subject).hard_dots_rule_value = t.hard_dots_rule_value;
%     d.subjects(subject).easy_dots_rule_value = t.easy_dots_rule_value;
%     d.subjects(subject).coherence_values = t.coherence_values;
%     d.subjects(subject).updated_coherence_values = t.updated_coherence_values;
%     d.subjects(subject).coh_psignifit_array = t.coh_psignifit_array;
%     d.subjects(subject).rule_easy_coh_psignifit_array = t.rule_ecoh_psignifit_array;
%     d.subjects(subject).rule_easy_coh_psignifit_array = t.rule_hcoh_psignifit_array;
    
    disp('*coherence*')
    disp('accuracy (low coh/hard threshold is .9, high coh/easy threshold is .7)')
    disp(accthis(d.subjects(subject).coh.correct))
    disp('num invalid')
    disp(length(find(d.subjects(subject).coh.rt == -1)))
    disp('num fast responses')
    disp(length(find(d.subjects(subject).coh.rt >=0 & d.subjects(subject).coh.rt < 400)))
    if p.plot_norms
        figure; normplot(d.subjects(subject).coh.rt)
        title('coherence thresholding')
        export_fig(fullfile(figdir,strcat(num2str(subject),'_coh_normplot.jpeg')),'-transparent')
    end
    disp('*rule*')
    disp('accuracy (hard threshold is .6)')
    disp(accthis(d.subjects(subject).rule.correct))
    disp('num invalid')
    disp(length(find(d.subjects(subject).rule.rt == -1)))
    disp('num fast responses')
    disp(length(find(d.subjects(subject).rule.rt >=0 & d.subjects(subject).rule.rt < 400)))
    if p.plot_norms
        figure; normplot(d.subjects(subject).rule.rt)
        title('match thresholding')
        export_fig(fullfile(figdir,strcat(num2str(subject),'_match_normplot.jpeg')),'-transparent')
    end
    disp('*coherence second thresholding*')
    disp('accuracy (low coh/hard threshold is .9, high coh/easy threshold is .7)')
    disp(accthis(d.subjects(subject).coh_ang.correct))
    disp('num invalid')
    disp(length(find(d.subjects(subject).coh_ang.rt == -1)))
    disp('num fast responses')
    disp(length(find(d.subjects(subject).coh_ang.rt >=0 & d.subjects(subject).coh_ang.rt < 400)))
    if p.plot_norms
        figure; normplot(d.subjects(subject).coh_ang.rt)
        title('coherence second thresholding')
        export_fig(fullfile(figdir,strcat(num2str(subject),'_coh_ang_normplot.jpeg')),'-transparent')
    end
    disp('*experiment*')
    disp('accuracy (.6 is a floor effect)')
    disp(accthis(d.subjects(subject).exp.correct))
    disp('num invalid')
    disp(length(find(d.subjects(subject).exp.rt == -1)))
    disp('num fast responses')
    disp(length(find(d.subjects(subject).exp.rt >=0 & d.subjects(subject).exp.rt < 400)))
    if p.plot_norms
        figure; normplot(d.subjects(subject).exp.rt)
        title('experiment')
        export_fig(fullfile(figdir,strcat(num2str(subject),'_exp_normplot.jpeg')),'-transparent')
    end
    
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
            [t.coh_easy,t.coh_hard,~,t.psignifit_array] = coh_thresholding(d.subjects(subject).coh,figdir,p.save_file,subject,1);
            if p.check_coh
                [t.jscoh_easy,t.jscoh_hard,~,t.jscoh_psignifit_array] = coh_thresholding(t.coh_psignifit_array,figdir,p.save_file,subject,0,1);
                disp('**checking js vs matlab psignifit for coherence**')
                disp('*easy:*')
                disp('matlab:')
                disp(t.coh_easy)
                disp('js:')
                disp(t.jscoh_easy)
                t.prompt = 'Press enter to continue';
                input(t.prompt,'s');
                disp('*hard:*')
                disp('matlab:')
                disp(t.coh_hard)
                disp('js:')
                disp(t.jscoh_hard)
                t.prompt = 'press enter to continue';
                input(t.prompt,'s');
                disp('*arrays:*')
                disp('matlab:')
                disp(t.psignifit_array)
                disp('js:')
                disp(t.jscoh_psignifit_array)
                t.prompt = 'press enter to continue';
                input(t.prompt,'s');
            else
                disp('*coherence thresholds*')
                disp('easy:')
                disp(t.coh_easy)
                disp('hard:')
                disp(t.coh_hard)
            end
        end
        if p.plot_match
            [t.match_easy,t.match_hard,t.match_summary,t.psignifit_array] = match_thresholding(d.subjects(subject).rule,figdir,p.save_file,subject,1);
            if p.check_match
                tmp(:,:,1) = t.rule_ecoh_psignifit_array;
                tmp(:,:,2) = t.rule_hcoh_psignifit_array;
                [t.jsmatch_easy,t.jsmatch_hard,t.jsmatch_summary,t.jsmatch_psignifit_array] = match_thresholding(tmp,figdir,p.save_file,subject,0,1); clear tmp;
                disp('**checking js vs matlab psignifit for matching**')
                disp('*easy (easy rule, hard rule):*')
                disp('matlab:')
                disp(t.match_easy)
                disp('js:')
                disp(t.jsmatch_easy)
                t.prompt = 'press enter to continue';
                input(t.prompt,'s');
                disp('*hard:*')
                disp('matlab:')
                disp(t.match_hard)
                disp('js:')
                disp(t.jsmatch_hard)
                t.prompt = 'press enter to continue';
                input(t.prompt,'s');
                disp('*arrays:*')
                disp('matlab:')
                disp(t.psignifit_array)
                disp('js:')
                disp(t.jsmatch_psignifit_array)
                t.prompt = 'press enter to continue';
                input(t.prompt,'s');
            else
                disp('*matching thresholds*')
                disp('easy:')
                disp(t.match_easy)
                disp('mean rt, correct trials, easy condition:')
                disp(mean(t.match_summary(4,:)))
                disp('hard:')
                disp(t.match_hard)
                disp('mean rt, correct trials, hard condition:')
                disp(mean(t.match_summary(6,:)))
            end
        end
        if p.plot_coh_ang
            [t.coh_easy,t.coh_hard,~,t.psignifit_array] = coh_thresholding(d.subjects(subject).coh_ang,figdir,p.save_file,subject,1);
            if p.check_coh_ang
                [t.jscoh_easy,t.jscoh_hard,~,t.jscoh_psignifit_array] = coh_thresholding(t.coh_psignifit_array_updated,figdir,p.save_file,subject,0,1);
                disp('**checking js vs matlab psignifit for coherence two**')
                disp('*easy:*')
                disp('matlab:')
                disp(t.coh_easy)
                disp('js:')
                disp(t.jscoh_easy)
                t.prompt = 'press enter to continue';
                input(t.prompt,'s');
                disp('*hard:*')
                disp('matlab:')
                disp(t.coh_hard)
                disp('js:')
                disp(t.jscoh_hard)
                t.prompt = 'press enter to continue';
                input(t.prompt,'s');
                disp('*arrays:*')
                disp('matlab:')
                disp(t.psignifit_array)
                disp('js:')
                disp(t.jscoh_psignifit_array)
                t.prompt = 'press enter to continue';
                input(t.prompt,'s');
            else
                disp('*coherence thresholds*')
                disp('easy:')
                disp(t.coh_easy)
                disp('hard:')
                disp(t.coh_hard)
            end
        end
        
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
        
        if p.plot_rt_hist
            % get some rt histograms
            disp('1 = EcEr, 2 = EcHr, 3 = HcEr, 4 = HcHr')
            titles = {'EcEr','EcHr','HcEr','HcHr'};
            all_conditions = cell2mat(t.consolidata(:,2));
            all_accuracies = cell2mat(t.consolidata(:,5));
            all_rts = cell2mat(t.consolidata(:,4));
            figure;
            for condition = 1:4
                condition_idx = all_conditions == condition;
                the_rts = condition_idx & all_accuracies; % condition index and accuracy is 1 (not 0)
                subplot(2,2,condition)
                h = histogram(all_rts(the_rts),'FaceColor',[0.0 0.502 0.502]);
                title(titles{condition});
                h.NumBins = 40;
                xlim([0,1500]);
                ylim([0, 40]);
            end; clear condition all_conditions all_accuracies all_rts condition_idx the_rts
            export_fig(fullfile(figdir,strcat('LBA_',num2str(subject),'_rt_hist.jpeg')),'-transparent')
        end
        
        if p.plot_rts
            % get some rt bars
            disp('1 = EcEr, 2 = EcHr, 3 = HcEr, 4 = HcHr')
            titles = {'EcEr','EcHr','HcEr','HcHr'};
            all_conditions = cell2mat(t.consolidata(:,2));
            all_accuracies = cell2mat(t.consolidata(:,5));
            all_rts = cell2mat(t.consolidata(:,4));
            figure;
            for condition = 1:4
                condition_idx = all_conditions == condition;
                the_rts = condition_idx & all_accuracies; % condition index and accuracy is 1 (not 0)
                mean_rts(condition) = mean(all_rts(the_rts),'omitnan');
                sem_rts(condition) = nansem(all_rts(the_rts));
            end; clear condition all_conditions all_accuracies all_rts accuracy_idx the_rts
            xvalues = categorical(titles);
            xvalues = reordercats(xvalues,titles);
            bar(xvalues,mean_rts,'FaceColor',[0.0 0.502 0.502]);
            hold on
            er = errorbar(xvalues,mean_rts,sem_rts);
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            hold off
            export_fig(fullfile(figdir,strcat('LBA_',num2str(subject),'_rts.jpeg')),'-transparent')
        end
        
        if p.plot_pc
            % get some accuracy bars
            disp('1 = EcEr, 2 = EcHr, 3 = HcEr, 4 = HcHr')
            titles = {'EcEr','EcHr','HcEr','HcHr'};
            all_conditions = cell2mat(t.consolidata(:,2));
            all_accuracies = cell2mat(t.consolidata(:,5));
            figure;
            for condition = 1:4
                condition_idx = all_conditions == condition;
                correct_accs = condition_idx & all_accuracies; % condition index and accuracy is 1 (not 0)
                condition_accs = all_accuracies(condition_idx);
                percent_correct(condition) = (sum(condition_accs)/length(condition_accs))*100;
            end; clear condition all_conditions all_accuracies all_rts accuracy_idx the_rts
            xvalues = categorical(titles);
            xvalues = reordercats(xvalues,titles);
            bar(xvalues,percent_correct,'FaceColor',[0.0 0.502 0.502]);
            export_fig(fullfile(figdir,strcat('LBA_',num2str(subject),'_percent_correct.jpeg')),'-transparent')
        end
        
        if p.skip_check_lbacont
            t.cont_lba = 'y';
        else
            t.prompt = 'Continue with this participant? y/n [y]: ';
            t.cont_lba = input(t.prompt,'s');
            if isempty(t.cont_lba); t.cont_lba = 'y'; end
            close all
        end
        
        if t.cont_lba
            % stack it up
            d.subjects(subject).lba = t.consolidata;
            d.lbadata{ilba} = t.consolidata;
        else
            d.subjects(subject).lba = 'n';
        end
        
    end % end lba approved if statement
    
end; clear subject ilba; % end subject loop for lba

fprintf('saving lba adjusted output from %s\n', mfilename);
save(p.save_file,'d'); % save all data to a .mat file


function accuracy = accthis(data)
accuracy = sum(data)/length(data);
end
