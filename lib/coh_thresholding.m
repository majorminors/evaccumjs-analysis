function [easy_threshold,hard_threshold,summary] = coh_thresholding(input_data,save_file)
% function [easy_value,hard_value] = coh_thresholding(p,d)
%
% coherence threshold analysis
%
% finds coherence value to achieve a specified percent correct for a
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
% produces:
%
% easy_threshold = value to achieve your higher percent correct
% hard_threshold = value to achieve your lower percent correct
%
%  will compute a sinusoid for the data, save parameters as 'sineparams'
%  and produce a figure
%    figure is x(top) - coherence
%              y      - percent correct
%
% will save figs using the full path save_file

%% set up

% enter your desired thresholds
low_threshold_pc = 0.7; % low coherence/lower percent correct = hard
high_threshold_pc = 0.9; % high coherence/higher percent correct = easy

% makes a structure that looks like
% | coherence point | number correct | number of trials |
% and builds row-wise each iteration
accuracy = NaN(10,length(input_data.stim_array)); % make array as large as (num_points,num_trials)
rts = NaN(10,length(input_data.stim_array));
for i = 1:length(input_data.stim_array)
    % rows are coh point, cols are accuracy
    accuracy(input_data.stim_array{1,i}.coh_point_code,i) = input_data.correct(1,i);
    coh_point_values(i) = input_data.stim_array{1,i}.coherence_value; % pull these for later
    if accuracy(input_data.stim_array{1,i}.coh_point_code,i) == 1
        rts(input_data.stim_array{1,i}.coh_point_code,i) = input_data.rt(1,i);
    end
end
coh_point_values = sort(unique(coh_point_values)); % cull to unique values only, and sort in order
correct_rts = mean(rts,2,'omitnan');
for i = 1:length(accuracy(:,1))
   psignifit_array(i,1) = coh_point_values(i);
   psignifit_array(i,2) = sum(accuracy(i,:),'omitnan');
end
psignifit_array(:,3) = sum(~isnan(accuracy),2);

% summary is four rows:
%   1) point condition
%   2) coherence value
%   3) percent correct
%   4) average rt for correct trials

summary(1,:) = psignifit_array(:,1);
summary(2,:) = coh_point_values';
summary(3,:) = psignifit_array(:,2)./psignifit_array(:,3);
summary(4,:) = correct_rts;

% make the sigmoid and put it on a figure
sigmoid = figure('visible','on');
[plotresult, plotline, plotdata] = psychcurve(psignifit_array);
hold on
% find the x value for proportion correct:
[~, low_threshold_idx] = min(abs(plotline.YData-low_threshold_pc(1))); % first find the index of the value closest to the threshold pc
[~, high_threshold_idx] = min(abs(plotline.YData-high_threshold_pc(1)));
low_threshold = plotline.XData(low_threshold_idx); % then find the values using the index
high_threshold = plotline.XData(high_threshold_idx);
% add plot lines at the threshold value on y:
plot([0 1], [low_threshold_pc low_threshold_pc], '-', 'Color',[1 0 0])
plot([0 1], [high_threshold_pc high_threshold_pc], '-', 'Color',[0 1 0])
% add plot lines at the threshold value on x:
plot([low_threshold low_threshold], [0.3 1], '-', 'Color',[1 0 0])
plot([high_threshold high_threshold], [0.3 1], '-', 'Color',[0 1 0])
savefig(sigmoid,[save_file '_sigmoid']);
hold off
% diplay rts on a figure
rts = figure('visible','off');
plot(summary(2,:),summary(4,:),'ro:')
savefig(rts,[save_file '_rts']);
% load those figures into variables
%sigmoid=hgload(fullfile(datadir, save_file, [data_file '_sigmoid.fig']));
%rts=hgload(fullfile(datadir, data_file, [data_file '_rts.fig']));

% prepare subplots
figure
visualise(1)=subplot(1,2,1);
visualise(2)=subplot(1,2,2);
% paste our figures on the subplots
copyobj(allchild(get(sigmoid,'CurrentAxes')),visualise(1));
copyobj(allchild(get(rts,'CurrentAxes')),visualise(2));
% add a legend
t(1)=title(visualise(1),'percent correct');
t(2)=title(visualise(2),'reaction time');

% make the output a little less confusing to understand
easy_threshold = high_threshold;
hard_threshold = low_threshold;

return

