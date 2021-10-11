%% do various sanity checks
%% Dorian Minors
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

datadir = fullfile(rootdir,'data','behav_9'); % location of data
dataToProcess = 'processed_data'; % where is the converted data?
saveFileName = 'are we gonna save this?'; % what to save the processed data as

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

for subject = 1:length(d.subjects) % loop through subjects
    close all
    
    thisSubject = d.subjects(subject);
    
    numTrials = numel(thisSubject.exp.rt);
    ecer=[];echr=[];hcer=[];hchr=[];
    for trial = 1:numTrials
        
        thisStimArray = thisSubject.exp.stim_array{trial};
        thisRT = thisSubject.exp.rt(trial);
        thisCorrect = thisSubject.exp.correct(trial);
        
        if thisStimArray.coh_difficulty == 1 && thisStimArray.match_difficulty == 1
            ecer = [ecer,[thisRT;thisCorrect]];
        elseif thisStimArray.coh_difficulty == 1 && thisStimArray.match_difficulty == 2
            echr = [echr,[thisRT;thisCorrect]];
        elseif thisStimArray.coh_difficulty == 2 && thisStimArray.match_difficulty == 1
            hcer = [hcer,[thisRT;thisCorrect]];
        elseif thisStimArray.coh_difficulty == 2 && thisStimArray.match_difficulty == 2
            hchr = [hchr,[thisRT;thisCorrect]];
        end
        
    end
    
    %% check some stoof
    summary(subject).ecer_num_trials = size(ecer,2);
    summary(subject).ecer_res = ecer;
    summary(subject).ecer_meanrt = mean(ecer(1,:),'omitnan');
    summary(subject).ecer_pc = accthis(ecer(2,:));
    summary(subject).echr_num_trials = size(echr,2);
    summary(subject).echr_res = echr;
    summary(subject).echr_meanrt = mean(echr(1,:),'omitnan');
    summary(subject).echr_pc = accthis(echr(2,:));
    summary(subject).hcer_num_trials = size(hcer,2);
    summary(subject).hcer_res = hcer;
    summary(subject).hcer_meanrt = mean(hcer(1,:),'omitnan');
    summary(subject).hcer_pc = accthis(hcer(2,:));
    summary(subject).hchr_num_trials = size(hchr,2);
    summary(subject).hchr_res = hchr;
    summary(subject).hchr_meanrt = mean(hchr(1,:),'omitnan');
    summary(subject).hchr_pc = accthis(hchr(2,:));
    
    
    %% let's look at first half vs second half
    titles = {'EcEr','EcHr','HcEr','HcHr'};
    % make this all indexable
    rts(:,:,1) = ecer;
    rts(:,:,2) = echr;
    rts(:,:,3) = hcer;
    rts(:,:,4) = hchr;
    figure;
    for condition = 1:size(rts,3)
        [mean_rts(1,condition),sem_rts(1,condition),mean_rts(2,condition),sem_rts(2,condition)] = meancomparison(rts(:,:,condition));
    end
    vals = [mean_rts(:,1)';mean_rts(:,2)';mean_rts(:,3)';mean_rts(:,4)'];
    err = [sem_rts(:,1)';sem_rts(:,2)';sem_rts(:,3)';sem_rts(:,4)'];
    bar(vals,'FaceColor',[0.0 0.502 0.502]);
    set(gca,'XTickLabel',titles);
    hold on
    [ngroups, nbars] = size(vals);
    % calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % set the position of each error bar in the centre of the main bar
    for i = 1:nbars
        % calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, vals(:,i), err(:,i), 'k', 'linestyle', 'none');
    end; clear x ngroups nbars groupwidth
    ylim([400, 1000]);
    hold off
    export_fig(fullfile(figdir,strcat(num2str(subject),'_mean_comparison.jpeg')),'-transparent')

    figure;
    for condition = 1:size(rts,3)
        [pc(1,condition),pc(2,condition)] = acccomparison(rts(:,:,condition));
    end
    vals = [pc(:,1)';pc(:,2)';pc(:,3)';pc(:,4)'];
    bar(vals,'FaceColor',[0.0 0.502 0.502]);
    set(gca,'XTickLabel',titles);
    ylim([50, 100]);
    hold off
    export_fig(fullfile(figdir,strcat(num2str(subject),'_pc_comparison.jpeg')),'-transparent')
    
    
    
end

all.ecerRT = mean(cell2mat({summary(:).ecer_meanrt}));
all.ecerPC = mean(cell2mat({summary(:).ecer_pc}));
all.echrRT = mean(cell2mat({summary(:).echr_meanrt}));
all.echrPC = mean(cell2mat({summary(:).echr_pc}));
all.hcerRT = mean(cell2mat({summary(:).hcer_meanrt}));
all.hcerPC = mean(cell2mat({summary(:).hcer_pc}));
all.hchrRT = mean(cell2mat({summary(:).hchr_meanrt}));
all.hchrPC = mean(cell2mat({summary(:).hchr_pc}));

% put together a table for quick checking
disp('check summary stats')
tmpTableTitles = {'condition' 'rt' 'pc'};
tmpTable = {...
    'ecer' ...
    all.ecerRT ...
    all.ecerPC; ...
    'echr' ...
    all.echrRT ...
    all.echrPC; ...
    'hcer' ...
    all.hcerRT ...
    all.hcerPC; ...
    'hchr' ...
    all.hchrRT ...
    all.hchrPC; ...
    };
tmpCheckTbl = cell2table(tmpTable);
tmpCheckTbl.Properties.VariableNames = tmpTableTitles;
disp(tmpCheckTbl)

all.ecerTrialNums = unique(cell2mat({summary(:).ecer_num_trials}));
all.echrTrialNums = unique(cell2mat({summary(:).echr_num_trials}));
all.hcerTrialNums = unique(cell2mat({summary(:).hcer_num_trials}));
all.hchrTrialNums = unique(cell2mat({summary(:).hchr_num_trials}));

% put together a table for quick checking
disp('check unique trial nums')
tmpTableTitles = {'condition' 'unique_trial_nums'};
tmpTable = {...
    'ecer' ...
    all.ecerTrialNums; ...
    'echr' ...
    all.echrTrialNums; ...
    'hcer' ...
    all.hcerTrialNums; ...
    'hchr' ...
    all.hchrTrialNums; ...
    };
tmpCheckTbl = cell2table(tmpTable);
tmpCheckTbl.Properties.VariableNames = tmpTableTitles;
disp(tmpCheckTbl)

function [firstHalf, firstHalfSem, secondHalf, secondHalfSem] = meancomparison(data)
% expects 2xn matrix. Row 1 is vals. Row 2 is correct/incorrect
% returns mean values for first and second half of vals

% split data in half
firstHalf = data(:,1:round(size(data,2)/2));
secondHalf = data(:,round(size(data,2)/2):end);

% reduce to correct values
firstHalf = firstHalf(1,firstHalf(1,:) & firstHalf(2,:));
secondHalf = secondHalf(1,secondHalf(1,:) & secondHalf(2,:));

% return sem
firstHalfSem = nansem(firstHalf);
secondHalfSem = nansem(secondHalf);

% return means
firstHalf = mean(firstHalf,'omitnan');
secondHalf = mean(secondHalf,'omitnan');

end

function [firstHalf, secondHalf] = acccomparison(data)

% split data in half
firstHalf = data(:,1:round(size(data,2)/2));
secondHalf = data(:,round(size(data,2)/2):end);

firstHalf = (sum(firstHalf(2,:))/length(firstHalf(2,:)))*100;
secondHalf = (sum(secondHalf(2,:))/length(secondHalf(2,:)))*100;

end

function semval = nansem(vector_data)
    % Recall that s.e.m. = std(x)/sqrt(length(x));
    nonan_std = nanstd(vector_data);
    nonan_len = length(vector_data(~isnan(vector_data)));
    % Plug in values
    semval = nonan_std / sqrt(nonan_len);
end