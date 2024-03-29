%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy
d = struct(); % set up a structure for the data info
t = struct(); % set up a structure for temp data

% set up variables
rootdir = '/group/woolgar-lab/projects/Dorian/evaccum/evaccumjs-analysis'; % 'C:\Users\doria\Nextcloud\desiderata\desiderata\04 Research\05 Evidence Accumulation\01 EvAccum Code'; % % root directory - used to inform directory mappings

datadir = fullfile(rootdir,'data','behav_9_100optim_splitConds');
modelIdentifier = 'allConds';
modelNamePattern = ['Model_%s' '_' modelIdentifier '.mat'];

modeldir = fullfile(datadir,'lba_results'); % expects to find your modelling results here
toolsdir = fullfile(rootdir, 'lib'); % where are all your scripts/tools?

%p.datafilepattern =  'Model_*.mat';
%p.savefilename = '';

% directory mapping
addpath(genpath(toolsdir)); % add tools folder to path (don't think we need this, but in case)
%save_file = fullfile(modeldir, p.savefilename);

%% run the data2fit analysis again

lbadatadir = fullfile(datadir); % expects to find your data here and will save results in a sub-folder here
p.data_name = 'processed_data.mat'; % data file name

t.fileinfo = dir(fullfile(lbadatadir,p.data_name));
t.datapath = fullfile(lbadatadir,t.fileinfo.name);

% get the data
t.alldata = load(t.datapath);
t.data = t.alldata.d.lbadata;

num_subjs = length(t.data);%sum(~cellfun('isempty',{data.id})); % get the number of not empty arrays from the first field of the 'data' structure

%initialise cells
data2fit={};

for idxsubj = 1:num_subjs
    clear rt conds resp acc % clear vars from behavioural data to avoid reading mixed subject data
    fprintf('loading subject %1.0f of %1.0f \n',idxsubj,num_subjs);
    
    % pull the data
    subjdata = t.data{idxsubj}; % put the data in a readable variable
    dataValid = []; % now we'll strip invalid responses out
    validLabs = {};
    for i = 1:size(subjdata,1)
        if subjdata{i,4} >= 0 % if there's a valid response
            dataValid(end+1,:) = [subjdata{i,2} subjdata{i,3} subjdata{i,4} subjdata{i,5}]; % add the following rows to this new variable in order: condition code, response, rt, accuracy
            temp = subjdata{i,1}; % extract valid condition label (can't do in one step with multi-level structures and non-scalar indexing)
            validLabs{end+1} = temp; % add the valid condition label
        end
    end; clear i temp;
    
    % converting again to readable variables
    conds = dataValid(:,1);
    resp = dataValid(:,2);
    rt = dataValid(:,3);
    acc = dataValid(:,4);
    
    % do descriptive stats
    minRT(idxsubj) = min(rt(rt>=0.1)); % not used
    
    
    %%
    %Pool conditions according to design matrix & calculate stats
    unique_conds = unique(conds);
    for level = 1:length(unique_conds) % for all conditions
        
        trialn = conds == level;
        dataRaw{idxsubj,level}  = [resp(trialn) rt(trialn)];
        
        d.data2fit{idxsubj,level} = data_stats(dataRaw{idxsubj,level});
        
        % add some info to d.data2fit (the condition labels as a string, and
        % the minRT - in the case that the min rt is smaller than the
        % non-decision time - let's leave this for now, just in case)
        d.data2fit{idxsubj,level}.cond = validLabs{trialn};
        d.data2fit{idxsubj,level}.minRT = round(minRT(idxsubj),2);
        %parange(end,1) = round(minRT(idxsubj),2);%constrain the lower bound of T0 to the shortest RT
    end; clear level trialn
    
    clear rt conds resp acc dataRaw dataValid validLabs minRT % clear vars from behavioural data to avoid reading mixed subject data
end; clear idxsubj;

%% get model data
design_space={[1,3],[1,4],[1,3,4],[1,3,4,5],[1,2],[1,2,3],[1,2,4],[1,2,3,4],[1,2,3,4,5],[1,5],[1,3,5],[1,4,5],[1,2,5],[1,2,3,5],[1,2,4,5]};

for i = 1:length(design_space) % loop through each

    thisModel =  sprintf(modelNamePattern,num2str(i)); 
    thisModel = fullfile(modeldir,thisModel);    
    t.model_results = load(thisModel);
    
    %% plotting
    t.save_params = [];
    for idxsubj = 1:num_subjs
        [~,t.id] = min(t.model_results.bestval{idxsubj});
        t.bestparams = t.model_results.bestpar{idxsubj}(t.id,:);
        
        for level = 1:length(unique_conds)
            t.data2fit{level} = d.data2fit{idxsubj,level};
        end
        
        [~,parEE,parHE,parEH,parHH]=getModelParam_cell_RDK(t.model_results.settings.modfeat,2,t.bestparams);
        t.params(1) = parEE;
        t.params(2) = parEH;
        t.params(3) = parHE;
        t.params(4) = parHH;
        clear parEE parEH parHE parHH
        
        t.probmod = []; t.cumRT = [];
        for level = 1:length(unique_conds) % for all conditions
            [~,t.probmod(:,level),~,~,t.cumRT(:,level)]=mod_stats_sim('LBA_spec_FT3_c_even_template',t.params(level),1,t.data2fit{:,level}.allObs(:,1),2,1);
        end; clear level;
        
        % Plotting...
        if length(unique_conds) ~= length(t.data2fit); error('you dont appear to have as many conditions as sets of t.data2fit'); end % sanity check
        
        % hold onto params to plot later
        non_dec_time(idxsubj,:,i) = [mean(t.params(1).T0),mean(t.params(2).T0),mean(t.params(3).T0),mean(t.params(4).T0)];
        dec_bndry(idxsubj,:,i) = [mean(t.params(1).B),mean(t.params(2).B),mean(t.params(3).B),mean(t.params(4).B)];
        st_bias(idxsubj,:,i) = [mean(t.params(1).C0),mean(t.params(2).C0),mean(t.params(3).C0),mean(t.params(4).C0)];
        drift_rate(idxsubj,:,i) = [mean(t.params(1).Ame),mean(t.params(2).Ame),mean(t.params(3).Ame),mean(t.params(4).Ame)];
        
%         % plot parameters
%         %   t.params(1) = LL; t.params(2) = LH; t.params(3) = HL; t.params(4) = HH;
%         %   param names = B C0 Ame Astd T0
%         figure; hold on;
%         temp=[mean(t.params(1).C0),mean(t.params(2).C0),mean(t.params(3).C0),mean(t.params(4).C0)];
%         bar(temp');
%         ylim([min(temp)-.1 max(temp)+.1]);
%         %legend({'Decision Boundary'});
%         set(gca,'XTick',[1:4]);
%         set(gca,'XTickLabel',{'EE' 'EH' 'HE' 'HH'});
%         title('Start Bias');
%         
%         figure; hold on;
        temp=[mean(t.params(1).B),mean(t.params(2).B),mean(t.params(3).B),mean(t.params(4).B)];
%         bar(temp');
%         ylim([min(temp)-.1 max(temp)+.1]);
%         %legend({'Decision Boundary'});
%         set(gca,'XTick',[1:4]);
%         set(gca,'XTickLabel',{'EE' 'EH' 'HE' 'HH'});
%         title('Decision Boundary');
%         
%         figure; hold on;
%         temp=[mean(t.params(1).Ame),mean(t.params(2).Ame),mean(t.params(3).Ame),mean(t.params(4).Ame)];
%         bar(temp');
%         ylim([min(temp)-.1 max(temp)+.1]);
%         %legend({'Decision Boundary'});
%         set(gca,'XTick',[1:4]);
%         set(gca,'XTickLabel',{'EE' 'EH' 'HE' 'HH'});
%         title('Drift Rate');
% 
        t.save_params(idxsubj,:) = temp;
        
        for level = 1:length(unique_conds) % for all conditions

%             % plot quintiles
%             figure; hold on;
%             plot(t.data2fit{level}.allObs(1:end-1,1),t.data2fit{level}.allObs(1:end-1,2),'k','Marker','o');
%             plot(t.data2fit{level}.allObs(1:end-1,1),t.cumRT(1:end-1,level),'r','Marker','*');
%             legend({'data','model'},'Location','SouthEast');
%             % title({['Model params (B, A, Astd,To): [' num2str(t.model_results.bestpar(t.id,:),2) ']'] });
%             title('Chosen condition');
%             ylim([0 1]);
%             % xlim([0.2,0.6]);
%             ylabel('Cumulative probability');
%             xlabel('Response Time (s)');
%             
%             % plot choices
%             figure; hold on;
%             temp=[t.data2fit{level}.priorProb{1,1}(1),t.data2fit{level}.priorProb{1,2}(1)];
%             t.choiceProb=[temp;t.probmod(:,level)'];
%             bar(t.choiceProb');
%             legend({'data','model'});
%             set(gca,'XTick',[1:2]);
%             set(gca,'XTickLabel',{'Button 1' 'Button 2'});
%             title('Choice probability');
            
        end; clear level;
    end; clear idxsubj;
    


    
%     temp=[mean(t.save_params(:,1)),mean(t.save_params(:,2)),mean(t.save_params(:,3)),mean(t.save_params(:,4))];
%     figure; hold on;
%     bar(temp,'FaceColor',[0 .502 .502]);
%     ylim([min(temp)-0.05 max(temp)+0.05]);
%     set(gca,'XTick',[1:length(temp)]);
%     set(gca,'XTickLabel',{'EE' 'EH' 'HE' 'HH'});
%     title('Mean Decision Boundaries');
    
end; clear i;

for i = [1, 2, 5, 10]

    % plot T0
    tmpT0 = mean(non_dec_time(:,:,i),1);
    figure; hold on;
    bar(tmpT0);
    ylim([min(tmpT0)-0.5 max(tmpT0)+0.5]);
    %legend({'Decision Boundary'});
    set(gca,'XTick',[1:length(tmpT0)]);
    set(gca,'XTickLabel',{'EE' 'EH' 'HE' 'HH'});
    title(['Non-Decision Time, Model' num2str(i)]);
    
    % plot drift rate
    tmpAME = mean(drift_rate(:,:,i),1);
    figure; hold on;
    bar(tmpAME);
    ylim([min(min(tmpAME))-0.5 max(max(tmpAME))+0.5]);
    %legend({'Decision Boundary'});
    set(gca,'XTick',[1:length(tmpAME)]);
    set(gca,'XTickLabel',{'EE' 'EH' 'HE' 'HH'});
    title(['Drift Rate, Model' num2str(i)]);
    
    % plot decision boundary
    tmpB = mean(dec_bndry(:,:,i),1);
    figure; hold on;
    bar(tmpB);
    ylim([min(min(tmpB))-0.5 max(max(tmpB))+0.5]);
    %legend({'Decision Boundary'});
    set(gca,'XTick',[1:length(tmpB)]);
    set(gca,'XTickLabel',{'EE' 'EH' 'HE' 'HH'});
    title(['Decision Boundary, Model' num2str(i)]);
    
    % plot start bias
    tmpC0 = mean(st_bias(:,:,i),1);
    figure; hold on;
    bar(tmpC0);
    ylim([min(min(tmpC0))-0.5 max(max(tmpC0))+0.5]);
    %legend({'Decision Boundary'});
    set(gca,'XTick',[1:length(tmpC0)]);
    set(gca,'XTickLabel',{'EE' 'EH' 'HE' 'HH'});
    title(['Start Bias, Model' num2str(i)]);
end

