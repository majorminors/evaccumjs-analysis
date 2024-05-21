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
toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
datadir = fullfile(rootdir,'data','behav_9'); % location of data
figDir = fullfile(datadir,'behavFigs'); if ~exist(figDir,'dir'); mkdir(figDir); end

tmp = load(fullfile(datadir,'processed_data'));

teal = [0.2, 0.6, 0.7];
coral = [0.9, 0.4, 0.3];
lilac = [0.7, 0.5, 0.8];

n_colors = 3;
n_lightness = 10;

r = [linspace(teal(1), coral(1), n_lightness); linspace(coral(1), lilac(1), n_lightness); linspace(lilac(1), teal(1), n_lightness)];
g = [linspace(teal(2), coral(2), n_lightness); linspace(coral(2), lilac(2), n_lightness); linspace(lilac(2), teal(2), n_lightness)];
b = [linspace(teal(3), coral(3), n_lightness); linspace(coral(3), lilac(3), n_lightness); linspace(lilac(3), teal(3), n_lightness)];

colormap_matrix = reshape([r(:), g(:), b(:)], [], 3);

clear r g b n_colors n_lightness

% load and code

% since we coded this for behav_7, which has an 'lba' thing we must have
% used for initial result coding, we need to create that for behav_9
if ischar(tmp.d.subjects(1).lba)
    for s = 1:numel(tmp.d.subjects)
        tmp.d.subjects(s).lba = {};
        for t = 1:numel(tmp.d.subjects(s).exp.rt)
            
            
            if tmp.d.subjects(s).exp.stim_array{1,t}.coh_difficulty == 1 && tmp.d.subjects(s).exp.stim_array{1,t}.match_difficulty == 1
                x = 'EcEr';
            elseif tmp.d.subjects(s).exp.stim_array{1,t}.coh_difficulty == 1 && tmp.d.subjects(s).exp.stim_array{1,t}.match_difficulty == 2
                x = 'EcHr';
            elseif tmp.d.subjects(s).exp.stim_array{1,t}.coh_difficulty == 2 && tmp.d.subjects(s).exp.stim_array{1,t}.match_difficulty == 1
                x = 'HcEr';
            elseif tmp.d.subjects(s).exp.stim_array{1,t}.coh_difficulty == 2 && tmp.d.subjects(s).exp.stim_array{1,t}.match_difficulty == 2
                x = 'HcHr';
            end
            tmp.d.subjects(s).lba{t,1} = x;clear x%condlab
            tmp.d.subjects(s).lba{t,2} = tmp.d.subjects(s).exp.stim_array{1,t}.cue_dir;%cue
            tmp.d.subjects(s).lba{t,3} = tmp.d.subjects(s).exp.button(1,t);%button
            tmp.d.subjects(s).lba{t,4} = tmp.d.subjects(s).exp.rt(1,t);%rt
            tmp.d.subjects(s).lba{t,5} = tmp.d.subjects(s).exp.correct(1,t);%corr
            tmp.d.subjects(s).lba{t,6} = tmp.d.subjects(s).exp.stim_array{1,t}.trial_cond_num;%cond id
        end
    end
end

% make this callable, since we make one for all data, and one for ave data
varNames = {'subject', 'conditionCode', 'conditionLabel', 'rt', 'correct','cohDiff','catDiff'};
varTypes = {'double', 'double', 'string', 'double', 'double', 'string', 'string'};
makeTable = @() table('Size', [0 numel(varNames)],...
    'VariableTypes', varTypes,...
    'VariableNames', varNames);
clear varNames varTypes

allTotalTrials = []; allRejectedTrials = [];
allData = makeTable();

for s = 1:numel(tmp.d.subjects)
    
    subject = [];
    conditionCode = [];
    conditionLabel = [];
    rt = [];
    correct = [];
    cohDiff = [];
    catDiff = [];
    rejectedTrials = 0;
    
    thisSubjData = tmp.d.subjects(s).lba;
    numTrials = numel(tmp.d.subjects(s).exp.rt);
    
    for t = 1:numTrials
        
        
        if thisSubjData{t,3} == -1 || thisSubjData{t,4} < 300
            rejectedTrials = rejectedTrials+1;
            continue
        end
        
        if startsWith(thisSubjData(t,1),'Ec') % easy coh
            thisCoh = 'Easy Coherence';
            if endsWith(thisSubjData(t,1),'Er') % easy rule
                thisCat = 'Easy Categorisation';
                thisCond = 1;
                thisLab = 'EcEr';
            elseif endsWith(thisSubjData(t,1),'Hr') % hard rule
                thisCat = 'Hard Categorisation';
                thisCond = 2;
                thisLab = 'EcHr';
            end
        elseif startsWith(thisSubjData(t,1),'Hc') % easy coh
            thisCoh = 'Hard Coherence';
            if endsWith(thisSubjData(t,1),'Er') % easy rule
                thisCat = 'Easy Categorisation';
                thisCond = 3;
                thisLab = 'HcEr';
            elseif endsWith(thisSubjData(t,1),'Hr') % hard rule
                thisCat = 'Hard Categorisation';
                thisCond = 4;
                thisLab = 'HcHr';
            end
        end
        
        subject = [subject; s];
        conditionCode = [conditionCode; thisCond];
        conditionLabel = [conditionLabel; {thisLab}];
        rt = [rt; thisSubjData{t,4}];
        correct = [correct; thisSubjData{t,5}];
        cohDiff = [cohDiff; thisCoh];
        catDiff = [catDiff; thisCat];
        
        clear thisCond thisLab thisCoh thisCat
        
    end
    
    allTotalTrials = [allTotalTrials t]; clear t
    allRejectedTrials = [allRejectedTrials rejectedTrials]; clear rejectedTrials
    allData = [allData;table(subject,conditionCode,conditionLabel,rt,correct,cohDiff,catDiff)];
    
    clear conditionCode conditionLabel subject rt correct
end

numSubjects = s; clear s;
percentRejected = sum(allRejectedTrials)/sum(allTotalTrials) * 100; clear allRejectedTrials allTotalTrials

disp('averaging')

aveData = makeTable();

for subject = unique(allData.subject)'
    for conditionCode = unique(allData.conditionCode)'
        
        idx = allData.subject == subject & allData.conditionCode == conditionCode;
        
        conditionLabel = unique(allData.conditionLabel(idx));
        if numel(conditionLabel) > 1; error('something wrong with condition labels');end
                
        cohDiff = unique(allData.cohDiff(idx));
        if numel(cohDiff) > 1; error('something wrong with coh difficulty labels');end

        catDiff = unique(allData.catDiff(idx));
        if numel(catDiff) > 1; error('something wrong with cat difficulty labels');end
            
        rt = mean(allData.rt(idx));
        
        correct = sum(allData.correct(idx)>0)/numel(allData.correct(idx))*100;
        
        aveData = [aveData; table(subject,conditionCode,conditionLabel,rt,correct,cohDiff,catDiff)];
        
        clear idx correct conditionLabel rt
    end; clear conditionCode
end; clear subject

aveData.conditionLabel = cellstr(aveData.conditionLabel);
aveData.catDiff = cellstr(aveData.catDiff);
aveData.cohDiff = cellstr(aveData.cohDiff);

disp('done')

%% so now you have
% numSubjects
% percentRejected
% aveData



close all
clear g
g(1,1) = gramm('x',aveData.conditionLabel,'y',aveData.rt);
g(1,1).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,1).axe_property('YLim',[200 1000]);
g(1,1).set_names('x','Conditions','y','Mean RT (ms) with SEM');
g(1,1).set_title('a) Mean RTs for each condition');
g(1,1).set_color_options('map',teal);
g(1,2) = gramm('x',aveData.cohDiff,'y',aveData.rt);
g(1,2).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,2).axe_property('YLim',[200 1000]);
g(1,2).set_names('x','Averaged Across Catgorisation','y','Mean RT (ms) with SEM');
g(1,2).set_title('b) Mean RTs for Coherence Difficulty');
g(1,2).set_color_options('map',coral);
g(1,3) = gramm('x',aveData.catDiff,'y',aveData.rt);
g(1,3).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,3).axe_property('YLim',[200 1000]);
g(1,3).set_names('x','Averaged Across Coherence','y','Mean RT (ms) with SEM');
g(1,3).set_title('c) Mean RTs for Categorisation Difficulty');
g(1,3).set_color_options('map',lilac);

g(2,1) = gramm('x',aveData.conditionLabel,'y',aveData.correct);
g(2,1).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
% g(2,1).axe_property('YLim',[650 950]);
g(2,1).set_names('x','Conditions','y','Mean Accuracy (%) with SEM');
g(2,1).set_title('d) Mean Accuracy for each condition');
g(2,1).set_color_options('map',teal);
g(2,2) = gramm('x',aveData.cohDiff,'y',aveData.correct);
g(2,2).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
% g(2,2).axe_property('YLim',[650 950]);
g(2,2).set_names('x','Averaged Across Catgorisation','y','Mean Accuracy (%) with SEM');
g(2,2).set_title('e) Mean Accuracy for Coherence Difficulty');
g(2,2).set_color_options('map',coral);
g(2,3) = gramm('x',aveData.catDiff,'y',aveData.correct);
g(2,3).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
% g(2,3).axe_property('YLim',[650 950]);
g(2,3).set_names('x','Averaged Across Coherence','y','Mean Accuracy (%) with SEM');
g(2,3).set_title('f) Mean Accuracy for Categorisation Difficulty');
g(2,3).set_color_options('map',lilac);
g.draw();
f = gcf; f.Position = [10 10 1600 1600];
print([figDir filesep 'acc-rt.png'],'-dpng')

fileID = fopen([datadir filesep 'stats.txt'], 'w');
if fileID == -1
    error('File could not be opened.');
end
fprintf(fileID, 'N = %.0f\n', numSubjects);
fprintf(fileID, 'percent of trials rejected = %.2f\n', percentRejected);

disp('are easy and hard coherence different?')
[pval stats meanDiff] = dottest(aveData.rt(contains(aveData.cohDiff,'Easy')),aveData.rt(contains(aveData.cohDiff,'Hard')))
fprintf(fileID, 'rt coh diff stats\n');
fprintf(fileID, '  pval = %.0f\n', pval);
fprintf(fileID, '  tstat = %.0f\n', stats.tstat);
fprintf(fileID, '  df = %.0f\n', stats.df);
fprintf(fileID, '  sd = %.0f\n', stats.sd);
fprintf(fileID, '  mean diff = %.0f\n', meanDiff);

disp('are easy and hard categorisation different?')
[pval stats meanDiff] = dottest(aveData.rt(contains(aveData.catDiff,'Easy')),aveData.rt(contains(aveData.catDiff,'Hard')))
fprintf(fileID, 'rt cat diff stats\n');
fprintf(fileID, '  pval = %.0f\n', pval);
fprintf(fileID, '  tstat = %.0f\n', stats.tstat);
fprintf(fileID, '  df = %.0f\n', stats.df);
fprintf(fileID, '  sd = %.0f\n', stats.sd);
fprintf(fileID, '  mean diff = %.0f\n', meanDiff);

% disp('and the interaction?')
% dottest(...
%     aveData.rt(contains(aveData.cohDiff,'Easy'))-aveData.rt(contains(aveData.cohDiff,'Hard')),...
%     aveData.rt(contains(aveData.catDiff,'Easy'))-aveData.rt(contains(aveData.catDiff,'Hard'))...
%     )
disp('are easy and hard coherence different?')
[pval stats meanDiff] = dottest(aveData.correct(contains(aveData.cohDiff,'Easy')),aveData.correct(contains(aveData.cohDiff,'Hard')))
disp('are easy and hard categorisation different?')
fprintf(fileID, 'accuracy coh diff stats\n');
fprintf(fileID, '  pval = %.0f\n', pval);
fprintf(fileID, '  tstat = %.0f\n', stats.tstat);
fprintf(fileID, '  df = %.0f\n', stats.df);
fprintf(fileID, '  sd = %.0f\n', stats.sd);
fprintf(fileID, '  mean diff = %.0f\n', meanDiff);

[pval stats meanDiff] = dottest(aveData.correct(contains(aveData.catDiff,'Easy')),aveData.correct(contains(aveData.catDiff,'Hard')))
fprintf(fileID, 'accuracy cat diff stats\n');
fprintf(fileID, '  pval = %.0f\n', pval);
fprintf(fileID, '  tstat = %.0f\n', stats.tstat);
fprintf(fileID, '  df = %.0f\n', stats.df);
fprintf(fileID, '  sd = %.0f\n', stats.sd);
fprintf(fileID, '  mean diff = %.0f\n', meanDiff);

fclose(fileID);


writetable(makeTableWithNans(...
    {'EcEr' 'EcHr' 'HcEr' 'HcHr'},...
    aveData.rt(strcmp(aveData.conditionLabel,'EcEr')),...
    aveData.rt(strcmp(aveData.conditionLabel,'EcHr')),...
    aveData.rt(strcmp(aveData.conditionLabel,'HcEr')),...
    aveData.rt(strcmp(aveData.conditionLabel,'HcHr'))),...
    [datadir filesep 'rts-for-jasp.csv'])

writetable(makeTableWithNans(...
    {'EcEr' 'EcHr' 'HcEr' 'HcHr'},...
    aveData.correct(strcmp(aveData.conditionLabel,'EcEr')),...
    aveData.correct(strcmp(aveData.conditionLabel,'EcHr')),...
    aveData.correct(strcmp(aveData.conditionLabel,'HcEr')),...
    aveData.correct(strcmp(aveData.conditionLabel,'HcHr'))),...
    [datadir filesep 'accuracy-for-jasp.csv'])


%% subfunctions

function [b d m] = dottest(data1, data2)

[a b c d] = ttest(data1,data2);
fprintf('p: %.3f\n',b);
disp(d);
m = mean(data1-data2);
fprintf('meanDiff: %0.3f\n',m);

return
end
