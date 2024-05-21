close all;
clearvars;
clc;
% set up variables
rootdir = pwd; %% root directory - used to inform directory mappings

datadir = fullfile(rootdir,'data'); % location of data

x = readtable(fullfile(datadir,'behav_7','ages.csv'));

fprintf('for study with num subjs %.0f\n',numel(x.Var1))
fprintf('num females: ')
sum(strcmp(x.Var3,'f'))
fprintf('\nnum males: ')
sum(strcmp(x.Var3,'m'))
fprintf('\nmean age: %0.2f and std: %0.2f\n',mean(x.Var2),std(x.Var2))


y = readtable(fullfile(datadir,'behav_9','ages.csv'));

fprintf('for study with num subjs %.0f\n',numel(y.Var1))
fprintf('num females: ')
sum(strcmp(y.Var3,'f'))
fprintf('\nnum males: ')
sum(strcmp(y.Var3,'m'))
fprintf('\nmean age: %0.2f and std: %0.2f\n',mean(y.Var2),std(y.Var2))
