%% put together the processed data into a format for the lba fit
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

% set up variables
rootdir = pwd; %% root directory - used to inform directory mappings

datadir = fullfile(rootdir,'data','behav_9_test'); % location of data
dataToProcess = 'processed_data'; % where is the converted data?
saveFileName = 'lba_processed_data'; % what to save the processed data as

p.plot_rt_hist = 1;
p.plot_rts = 1;
p.plot_pc = 1;
p.skip_check_lbacont = 0;

p.conditions = {'EcEr','EcHr','HcEr','HcHr'}; % 2x2 coherence and rule
p.conditioncodes = {1,2,3,4};

% cobble together what we need to play with the data and save it
theData = load(fullfile(datadir,dataToProcess)); % load the data
d = theData.d; % let's just get what we really want and put it in the format the copied code below wants it in
addpath(genpath(fullfile(rootdir, 'lib'))); % add libraries to path
figdir = fullfile(datadir,'figures'); % place to save figures
if ~exist(figdir,'dir')
    mkdir(figdir);
end
p.save_file = fullfile(datadir, saveFileName);

fprintf('converting selected participants for lba\n');

ilba = 0; % initialise a counter
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
            disp('id:')
            disp(d.subjects(subject).id)
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
