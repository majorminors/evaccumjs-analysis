function [easy_threshold,hard_threshold,summary] = match_thresholding(input_data,save_file)
% function [easy_value,hard_value] = coh_thresholding(p,d)
%
% matching threshold analysis
%
% finds matching rule value to achieve a specified percent correct for a
% participant
% 
% requires:
%    various test parameters (expects to find in structure 'p')
%    various saved outputs of test (expects to find in structure 'd')
%    the full path to the saved data (so it can use the same save directory
%       and savename conventions (save_file)
%
% specify in this function your desired percent correct for thresholding (low_threshold_pc,
%   high_threshold_pc)
%
% note: you're overwriting variables when you produce your plots. You
%       should address this at some point.
%
% produces:
%
% 'overview', both for 'low_coh' and 'hi_coh' which is three rows:
%   point condition per trial | correct/incorrect (1/0) per trial | matching angle (0-90 degrees) per trial
%       note: currently coded such that block one is low coherence and block
%           two is high coherence
%
% 'summary' which is six rows:
%   point condition | matching angle | percent correct for low coherence trials | mean rt for correct trials (low coh) | percent correct for high coherence trials | mean rt for correct trials (high coh)
%
%
%  will compute a sinusoid for the data, save parameters as 'sineparams'
%  and produce a figure
%    figure is x(top) - matching angle
%              y      - percent correct
%
% will save figs using the full path save_file

%% set up

% enter your thresholds
low_threshold_pc = 0.6; % percent correct for your hard threshold (uses inverse of this number for easy)
num_blocks = 2;
num_trials_per_block = length(input_data.stim_array)/num_blocks;

% makes a structure that looks like
% | match point | number correct | number of trials |
% and builds row-wise each iteration
all_accuracy = NaN(10,length(input_data.stim_array)); % make array as large as (num_points,num_trials)
all_rts = NaN(10,length(input_data.stim_array));
for i = 1:length(input_data.stim_array)
    % rows are rule point, cols are accuracy
    all_accuracy(input_data.stim_array{1,i}.rule_point_code,i) = input_data.correct(1,i);
    all_rule_point_values(i) = input_data.stim_array{1,i}.rule_value; % pull these for later
    if all_accuracy(input_data.stim_array{1,i}.rule_point_code,i) == 1
        all_rts(input_data.stim_array{1,i}.rule_point_code,i) = input_data.rt(1,i);
    end
end


for block = 1:num_blocks
    trials_start = num_trials_per_block*block-(num_trials_per_block-1);
    trials_end = num_trials_per_block*block;
    rule_point_values(:,:,block) = sort(unique(all_rule_point_values(:,trials_start:trials_end))); % cull to unique values only, and sort in order
    correct_rts(:,:,block) = mean(all_rts(:,trials_start:trials_end),2,'omitnan');
    accuracy(:,:,block) = all_accuracy(:,trials_start:trials_end);
    for i = 1:length(accuracy(:,1,block))
       psignifit_array(i,1,block) = rule_point_values(:,i,block);
       psignifit_array(i,2,block) = sum(accuracy(i,:,block),'omitnan');
    end
    psignifit_array(:,3,block) = sum(~isnan(accuracy(:,:,block)),2);
    
end

% Note: the psignifit tool only goes low to high, so if as in this case the
%       values of your stimulus level goes high to low, then you can flip
%       the intensity values (i.e. lie to matlab) or invert them to their
%       negative.
neg_psignifit_array(:,1,:) = -psignifit_array(:,1,:);
neg_psignifit_array(:,2:3,:) = psignifit_array(:,2:3,:);


    % summary is six rows:
%   1) point condition
%   2) matching angle
%   3) percent correct (low coherence)
%   4) mean rt for correct trials (low coherence)
%   5) percent correct (high coherence)
%   6) mean rt for correct trials (high coherence)
summary(1,:) = 1:length(rule_point_values(:,:,1));
summary(2,:) = rule_point_values(:,:,1)';
summary(3,:) = psignifit_array(:,2,1)./psignifit_array(:,3,1);
summary(4,:) = correct_rts(:,:,1);
summary(5,:) = psignifit_array(:,2,2)./psignifit_array(:,3,2);
summary(6,:) = correct_rts(:,:,2);
summary(1,:) = -summary(1,:);
summary(2,:) = -summary(2,:);

% prep data for psignifit



% make a sigmoid and put it on a figure
sigmoid_low = figure('visible','on');
[plotresult, plotline, plotdata] = psychcurve(neg_psignifit_array(:,:,2));
hold on
% find the x value for proportion correct:
[~, low_threshold_idx] = min(abs(plotline.YData-low_threshold_pc(1))); % first find the index of the value closest to the threshold pc
low_threshold = plotline.XData(low_threshold_idx); % then find the value using the index
high_threshold = -90-low_threshold; % this is the inverse of the low_threshold value - currently I've flipped everything into negative so it is oriented correctly in the psignifit tools
% add plot lines at the threshold value on y:
plot([-90 -0], [low_threshold_pc low_threshold_pc], '-', 'Color',[1 0 0])
[~, high_pc_idx] = min(abs(plotline.XData-high_threshold(1))); % then find the index of the value closest to the high_threshold
high_threshold_pc = plotline.YData(high_pc_idx); % then find the value using the index and print that in the command window
plot([-90 -0], [high_threshold_pc high_threshold_pc], '-', 'Color',[0 1 0])
% add plot lines at the threshold value on x:
plot([low_threshold low_threshold], [0.3 1], '-', 'Color',[1 0 0])
plot([high_threshold high_threshold], [0.3 1], '-', 'Color',[0 1 0])
savefig([save_file '_lowcohsigmoid']);
hold off
% diplay rts on a figure
rts_low = figure('visible','off');
plot(summary(2,:),summary(4,:),'ro:');
savefig([save_file '_low_coh_rts']);
% make a sigmoid and put it on a figure
sigmoid_hi = figure('visible','off');
[plotresult, plotline, plotdata] = psychcurve(neg_psignifit_array(:,:,1));
hold on
[~, low_pc_idx] = min(abs(plotline.XData-low_threshold(1))); % then find the index of the value closest to the high_threshold
low_threshold_pc = plotline.YData(low_pc_idx); % then find the value using the index and print that in the command window
[~, high_pc_idx] = min(abs(plotline.XData-high_threshold(1))); % then find the index of the value closest to the high_threshold
high_threshold_pc = plotline.YData(high_pc_idx); % then find the value using the index and print that in the command window
plot([-90 -0], [low_threshold_pc low_threshold_pc], '-', 'Color',[1 0 0])
plot([-90 -0], [low_threshold_pc low_threshold_pc], '-', 'Color',[1 0 0])
plot([low_threshold low_threshold], [0.3 1], '-', 'Color',[1 0 0])
plot([high_threshold high_threshold], [0.3 1], '-', 'Color',[0 1 0])
savefig([save_file '_hicohsigmoid']);
hold off
% diplay rts on a figure
rts_hi = figure('visible','off');
plot(summary(2,:),summary(6,:),'ro:');
savefig([save_file '_hi_coh_rts']);

low_threshold = abs(low_threshold);
high_threshold = abs(high_threshold);

% load those figures into variables
% sigmoid_low=hgload(fullfile(datadir, data_file, [data_file '_lowcohsigmoid.fig']));
% rts_low=hgload(fullfile(datadir, data_file, [data_file '_low_coh_rts.fig']));
% sigmoid_hi=hgload(fullfile(datadir, data_file, [data_file '_hicohsigmoid.fig']));
% rts_hi=hgload(fullfile(datadir, data_file, [data_file '_hi_coh_rts.fig']));

% prepare subplots
figure
visualise(1)=subplot(2,2,1);
visualise(2)=subplot(2,2,2);
visualise(3)=subplot(2,2,3);
visualise(4)=subplot(2,2,4);
% paste our figures on the subplots
copyobj(allchild(get(sigmoid_low,'CurrentAxes')),visualise(1));
copyobj(allchild(get(rts_low,'CurrentAxes')),visualise(2));
copyobj(allchild(get(sigmoid_hi,'CurrentAxes')),visualise(3));
copyobj(allchild(get(rts_hi,'CurrentAxes')),visualise(4));
% add a legend
t(1)=title(visualise(1),'percent correct low coh');
t(2)=title(visualise(2),'reaction time low coh');
t(3)=title(visualise(3),'percent correct hi coh');
t(4)=title(visualise(4),'reaction time hi coh');

% make the output a little less confusing to understand
easy_threshold = high_threshold;
hard_threshold = low_threshold;

return
