%% Simulate a matching threshold test for evidence accumulation task
% Dorian Minors
% Created: OCT19
% Last Edit: 18OCT19

% note that to simulate, you need to run simulatecurve.m outside the
% block/trial loops, and then run idealObserver in the trial loop

% trial settings all saved in 'p'

% data related information saved in 'd'
%   d.stim_mat_all contains trial condition matrices for each block
%       (d.stim_mat_all(:,:,[block number])
%   rows indicate trials
%   columns are explained when we define the matrix in 'define stimulus
%       parameters', below
%   in this script, the two blocks are two tests:
%       block 1 tests using easy coherence threshold.
%       block 2 tests using hard coherence threshold.


% other trial specific variables are in 't' in case something goes wrong
%   and we want to see them

% note: response key is swapped halfway through each block, and again at
%       the commencement of the subsequent block(s). this is not
%       counterbalanced within blocks, merely done as a training tool

% note: block 1 tests using easy coherence threshold.
%       block 2 tests using hard coherence threshold.
%       stored as d.stim_mat_all(:,:,[block number])

% note: for backwards compatibility (to 2016), subfunctions are disabled
%       in scripts (although functions can have subfunctions and scripts
%       instead use external functions.

% note: for compatibility with testing lab computers, KbQueueWait() only accepts
%       one argument - but from at least 2018, KbQueueWait can accept two,
%       which is commented out next to any call to KbQueueWait here.

%% subfunctions (at bottom):
%% pix = angle2pix(p,ang)
% calculates pixel size from visual angles, assuming isotropic (square) pixels

% requires:
% p.screen_distance           distance from screen in cm
% p.screen_width              width of screen in cm
% p.resolution                number of pixels of p in horizontal direction - this needs to be calculated after the window is opened for the dots

%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % est structure for parameter values
d = struct(); % est structure for trial data
t = struct(); % another structure for untidy trial specific floating variables that we might want to interrogate later if we mess up
points = simulatecurve; % simulate some points for idealObserver

% set up variables
rootdir = 'C:\Users\doria\Google Drive\04 Research\05 Evidence Accumulation\01 EvAccum Code'; % root directory - used to inform directory mappings
p.screen_num = 0; % screen to display experiment on (0 unless multiple screens)
p.fullscreen = 1; % 1 is full screen, 0 is whatever you've set p.window_size to
p.testing = 0; % change to 0 if not testing (1 skips PTB synctests and sets number of trials and blocks to test values) - see '% test variables' below
p.training = 0; % if 1, initiates training protocol (reduce dots presentation time from 'p.training_dots_duration' to 'p.dots_duration' by one 'p.training_reduction' every 'p.training_interval') - see '% training variables' below

% directory mapping
addpath(genpath(fullfile(rootdir, 'tools'))); % add tools folder to path (includes moving_dots function which is required for dot motion, as well as an external copy of subfunctions for backwards compatibility with MATLAB)
stimdir = fullfile(rootdir, 'stimuli');
datadir = fullfile(rootdir, 'data'); % will make a data directory if none exists
if ~exist(datadir,'dir')
    mkdir(datadir);
end

% training variables
p.training_dots_duration = 5; % duration (secs) dot stimulus presented for for training trials (will reduce trial by trial to p.dots_duration)
p.training_reduction = 1; % by how much do we reduce the duration during training (seconds)
p.training_interval = 2; % how many trials should we train on before reducing the dots presentation time
p.training_trials_per_level = 8; % how many trials to run each training level at (currently 3 levels, and first level will wait for 75%/p.training_trials_per_level while other levels will just do three and move on)
p.training_rule_lvl_easy = 10; % what matching rule level for easy (training level 1)
p.training_rule_lvl_mid = 45; % what matching rule level for mid (training level 2)
p.training_rule_lvl_hard = 75; % what matching rule level for hard (training level 3)

% test variables
p.num_test_trials = 4;
p.num_test_blocks = 2;
if p.testing == 1
    p.PTBsynctests = 1; % PTB will skip synctests if 1
    p.PTBverbosity = 1; % PTB will only display critical warnings with 1
elseif p.testing == 0
    p.PTBsynctests = 0;
    p.PTBverbosity = 3; % default verbosity for PTB
end
Screen('Preference', 'SkipSyncTests', p.PTBsynctests);
Screen('Preference', 'Verbosity', p.PTBverbosity);
if p.testing > 1; error('invalid value for p.testing'); end % check p.testing is a valid number, or error

% psychtoolbox setup
AssertOpenGL; % check Psychtoolbox (on OpenGL) and Screen() is working
KbName('UnifyKeyNames'); % makes key mappings compatible (mac/win)
rng('shuffle'); % seed rng using date and time

% set up participant info and save
t.prompt = {'enter participant number:','enter easy coherence threshold (fm 0-1, higher is easier)','enter hard coherence threshold (fm 0-1, lower is harder)'}; % prompt a dialog to enter subject info
t.prompt_defaultans = {num2str(99), num2str(0.75), num2str(0.25)}; % default answers corresponding to prompts
t.prompt_rsp = inputdlg(t.prompt, 'enter participant info', 1, t.prompt_defaultans); % save dialog responses
d.participant_id = str2double(t.prompt_rsp{1}); % add subject number to 'd'
d.easy_coherence = str2double(t.prompt_rsp{2}); % add participant coherence thresholds to 'd'
d.hard_coherence = str2double(t.prompt_rsp{3}); % add participant coherence thresholds to 'd'

% check participant info has been entered correctly for the script
if isnan(d.participant_id)
    error('no participant number entered');
elseif isnan(d.easy_coherence) || d.easy_coherence > 1
    error('invalid participant coherence threshold (easy)');
elseif isnan(d.hard_coherence) || d.hard_coherence > 1
    error('invalid participant coherence threshold (hard)')
end

% create a save file
save_file_name = [num2str(d.participant_id,'S%02d'),'_',mfilename];
save_file = fullfile(datadir, save_file_name);
if exist([save_file '.mat'],'file') % check if the file already exists and throw a warning if it does
    warning('the following save file already exists - overwrite? (y/n)\n %s.mat', save_file);
    while 1 % loop forever until y or n
        ListenChar(2);
        [secs,keyCode] = KbWait; % wait for response
        key_name = KbName(keyCode); % find out name of key that was pressed
        if strcmp(key_name, 'y')
            fprintf('instructed to overwrite:\n %s.mat\n overwriting and continuing with %s\n', save_file, mfilename)
            ListenChar(0);
            clear secs keyCode key_name
            break % break the loop and continue
        elseif strcmp(key_name, 'n')
            ListenChar(0);
            clear secs keyCode key_name
            error('instructed not to overwrite:\n %s.mat\n aborting %s\n', save_file, mfilename); % error out
        end
    end % end response loop
end % end check save file exist
save(save_file); % save all data to a .mat file

%% define experiment parameters

fprintf('defining exp params for %s\n', mfilename);

% define keys
p.quitkey = {'q'};
p.resp_keys = {'a','s'}; % only accepts two response options
% establish response keys for the trial based off whether participant id is odd or even
if mod(d.participant_id,2) % if pid not divisible by 2 (i.e. leaves a modulus after division)
    p.resp_keys = {p.resp_keys{1},p.resp_keys{2}}; % essentially do nothing - keep resp keys the same
elseif ~mod(d.participant_id,2) % if pid is divisible by 2 (i.e. does not leave a modulus after division)
    p.resp_keys = {p.resp_keys{2},p.resp_keys{1}}; % then swap the response keys
end

% define display info
p.bg_colour = [0 0 0]; % needs to be the same as the cue stimuli background colour (unless transparent)
p.text_colour= [255 255 255]; % colour of instructional text
p.cue_colour_blue = [121 181 240]; % colour of cue, for text formatting
p.cue_colour_orange = [240 181 121]; % colour of cue, for text formatting
p.text_size = 40; % size of text
p.window_size = [0 0 1200 800]; % size of window when ~p.fullscreen
p.screen_width = 35;   % Screen width in cm
p.screen_height = 50;    % Screen height in cm
p.screen_distance = 50; % Screen distance from participant in cm
p.visual_angle_cue = 15; % visual angle of the cue expressed as a decimal - determines size
p.visual_angle_dots = 0.15; % visual angle of the dots expressed as a decimal - determines size

% timing info
p.min_cue_time = 1; % minimum period to display cue (participants can't continue during this time)
p.dots_duration = 2; % seconds for the dot cloud to be displayed
p.min_resp_mapping_time = 1; % minimum period to display response mapping (participants can't continue during this time)
p.feedback_time = 0.5; % period to display feedback after response

% trial settings (*p.stim_mat* = parameter required to calculate stimulus condition matrix)
p.num_trials_per_block = 160; % *p.stim_mat* - must be divisible by p.num_cues && >= 15*p.num_points
p.num_blocks = 2; % if set to 1, keys will swap halfway through trials in each test
p.num_cues = 4; % *p.stim_mat*
p.start_coherence = 1; % starting coherence for the coherence test for training trials (0-1 - 1 is 100% coherence)
p.start_distance = 0; % start distance in degrees from the cue direction for the matching test for training trials
p.matching_cue_1 = 'BLUE'; % variable used to indicate response keys - this the upward arrow of the doublesided arrow cue in stimdir
p.matching_cue_2 = 'ORANGE'; % variable used to indicate response keys - this the downward arrow of the doublesided arrow cue in stimdir
p.cue_directions = 45:90:315; % *p.stim_mat* - refers to the direction of the upward arrow of the doublesided arrow cue in stimdir
p.num_points = 10; % *p.stim_mat* - number of points to test participants on for each test
p.rule_points = 0:10:90; %union([0:5:20],[70:5:90]); % *p.stim_mat* - length(p.rule_points) must == p.num_points
if length(p.rule_points) ~= p.num_points; error('number of coherence points is not equal to the number of testing points you specified'); end % check length(p.rule_points) == p.num_points, or error

%% define stimuli parameters

fprintf('defining stimuli params for %s\n', mfilename);

% read in stimulus file for the cue
p.cue = imread(fullfile(stimdir, 'arrows_cue_colours_with_line.png'));

% create matrix specifying stimulus conditions per trial:
%    1)  cue direction (1-4) - evenly allocates trials to cues
%    2)  blue cue direction in degrees for each trial - evenly adds cue
%        directions to trials in a similar manner to (1)
%    3)  coherence condition (1 = match blue, 2 = match orange) - evenly allocates
%        trials
%    4)  point condition (1-10) - each repeated
%        p.num_trials_per_block/p.num_points times
%    5)  matching points allocated to point conditions
%    6)  whether trial should add or subtract degrees from cues test (1 =
%        add, 2 = subtract) - currently 2 of each per cue (since four reps
%        of point conditions per cue)
%    7)  coherence direction from cued direction in degrees - calculated from cue direction and matching
%        point
%    8)  coherence direction in degrees incorporating blue match or orange match

% create matrix
p.stim_mat = zeros(p.num_trials_per_block,8);
% create columns
p.stim_mat(:,1) = sort(repmat(1:p.num_cues,[1,p.num_trials_per_block/p.num_cues]));
p.stim_mat(:,2) = p.cue_directions(p.stim_mat(:,1));
p.stim_mat(:,3) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/2])),[1,p.num_cues]);
p.stim_mat(:,4) = repmat(1:p.num_points,[1,p.num_trials_per_block/p.num_points]);
p.stim_mat(:,5) = repmat(p.rule_points,[1,p.num_trials_per_block/p.num_points]);
p.stim_mat(:,6) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/4])),[1,p.num_cues*2]);
add = (p.stim_mat(:,6)==1); % index all posns in (6) with 1s into 'add'
subtract = (p.stim_mat(:,6)==2); % index all posns in (6) with 2s into 'subtract'
p.stim_mat(add,7) = p.stim_mat(add,2)+p.stim_mat(add,5); % insert addition of matching point to cue direction into (7)
p.stim_mat(subtract,7) = p.stim_mat(subtract,2)-p.stim_mat(subtract,5); % insert subtraction of matching point to cue direction into (7)
temp = p.stim_mat(:,7)>360; % get all places where (7) > 360 degrees
p.stim_mat(temp,7) = p.stim_mat(temp,7)-360; % make any in (7) over 360 wrap around from 0 again
temp = p.stim_mat(:,7)<0; % get all places where (7) < 0 degrees
p.stim_mat(temp,7) = p.stim_mat(temp,7)+360; % make any in (7) under 0 wrap around from 360 again
temp = p.stim_mat(:,7)+180; % get 180 degrees on each cue
temp(temp>360) = temp(temp>360)-360; % make any over 360 wrap around from 0 again
p.stim_mat(:,8) = times(p.stim_mat(:,3)-1,temp); % make 1s and 2s from (3) into 0s and 1s, then times by temp (so becomes 0s and cue+180degs)
temp = (p.stim_mat(:,8)==0); % index all posns in (8) with 0s into 'temp'
p.stim_mat(temp,8) = p.stim_mat(temp,7); % insert coherence directions from (7) into (8) where there are 0s (using 'temp' index)
% clear floating variables
clear temp add subtract;

% shuffle the trial condition order for each block into 'd.stim_mat_all' then re-sort by cue
for block=1:p.num_blocks
    d.stim_mat_all(:,:,block) = p.stim_mat(Shuffle(1:p.num_trials_per_block),:);
    d.stim_mat_all(:,:,block) = sortrows(d.stim_mat_all(:,:,block),1);
end
% clear floating variables
clear block;

% dot cloud parameters for moving_dots function
%   note that for the following required params for moving_dots:
%       'dots.direction' is specified within the trial loop
%       'dots.coherence' is also specified withing the trial loop using outputs
%           from 'p.stim_mat' and 'd.easy_coherence'/'d.hard_coherence'
dots.aperture_size = [10,10]; % [width,height] in degrees of the rectangular aperture dots are displayed in
dots.centre = [0,0]; % [x,y] centre of the dot cloud
dots.num_dots = 100; % number of dots
dots.colour = [255,255,255]; % colour of the dots in [r,g,b]
dots.visual_angle = p.visual_angle_dots; % visual angle of the dots expressed as a decimal - determines size
dots.speed = 5; % speed of dots in degrees per second
dots.lifetime = 5; % number of frames dots live for

%% test start

fprintf('running matching threshold assessment (%s)\n', mfilename);

try
    %% block start
    
    % changes number of blocks to testing amount if testing
    if p.testing == 1
        p.act_block_num = p.num_test_blocks;
        fprintf('testing (p.testing set to 1) - will run %u blocks\n', p.num_test_blocks);
    else
        p.act_block_num = p.num_blocks;
    end
    
    % repeat trials for each block
    for block = 1:p.act_block_num
        fprintf('entering block %u of %u\n',block, p.act_block_num); %report block number to command window
        
        % pick up trial condition order for this block
        p.stim_mat = d.stim_mat_all(:,:,block);
        
        %% trials start
        
        % changes number of trials to testing amount for trial loop if testing
        if p.testing == 1
            p.act_trial_num = p.num_test_trials;
            fprintf('testing (p.testing set to 1) - will only run %u trials\n', p.num_test_trials);
        else
            p.act_trial_num = p.num_trials_per_block;
        end
        
        % open trial loop
        i = 0; % initialise trial index
        while i < p.act_trial_num
            i = i + 1;
            fprintf('trial %u of %u\n',i,p.act_trial_num); %report trial number to command window
            
            %set up a queue to collect response info
            t.queuekeys = [KbName(p.resp_keys{1}), KbName(p.resp_keys{2}), KbName(p.quitkey)]; % define the keys the queue cares about
            t.queuekeylist = zeros(1,256); % create a list of all possible keys (all 'turned off' i.e. zeroes)
            t.queuekeylist(t.queuekeys) = 1; % 'turn on' the keys we care about in the list (make them ones)
            KbQueueCreate([], t.queuekeylist); % initialises queue to collect response information from the list we made (not listening for response yet)
            KbQueueStart(); % starts delivering keypress info to the queue
            
            KbQueueFlush(); % flush the response queue from the continue screen presses
            
            % all presentation detail and running of movingdots removed
            
            %% deal with response and provide feedback (or abort if 'p.quitkey' pressed)
            
            % code response info
            t.firstPress = KbName(p.resp_keys{1}); % just put a value in t.firstPress
            d.resp_key_time(block,i) = 1; % just put a value in here % sum(t.firstPress); % get the timing info of the key used to respond
            d.rt(block,i) = 1; % just put a value in here %d.resp_key_time(block,i) - d.dots_onset(block,i); % rt is the timing of key info - time of dots onset (if you get minus values something's wrong with how we deal with nil/early responses)
            resp = idealObserver(abs(p.stim_mat(i,5)./100-1),points); % simulate response 1 0;
            
            % create variable for correct response
            if p.stim_mat(i,3) == 1 % if trial matches blue
                resp = (resp==0)*1+1;% 1 if correct 2 if wrong
                d.correct_resp(block,i) = p.resp_keys{1};
                d.incorrect_resp(block,i) = p.resp_keys{2};
            elseif p.stim_mat(i,3) == 2 % if trial matches orange
                resp = (resp==1)*1+1;% 2 if correct 1 if wrong
                d.correct_resp(block,i) = p.resp_keys{2};
                d.incorrect_resp(block,i) = p.resp_keys{1};
            end
            d.resp_key_name{block,i} = p.resp_keys{resp}; % get the name of the key used to respond
            
            % score response
            if strcmp(d.resp_key_name(block,i), d.correct_resp(block,i))
                d.correct(block,i) = 1; %correct trial
                t.feedback = 'correct';
            elseif strcmp(d.resp_key_name(block,i), d.incorrect_resp(block,i))
                d.correct(block,i) = 0; %incorrect trial
                t.feedback = 'incorrect';
            elseif strcmp(d.resp_key_name(block,i),p.quitkey)
                fclose('all');
                error('%s quit by user (p.quitkey pressed)\n', mfilename);
            else
                d.correct(block,i) = -1; % nil response
                d.rt(block,i) = 0;
                t.feedback = 'no valid input';
            end % end check correct
            
            %% post trial clean up
            
            % clear the response queue for the next trial and related floating variables
            KbQueueRelease();
            
            % save the trial data
            save(save_file); % save all data in '.mat' format
            
        end % trial loop
        
    end % end block loop
    
    %% wrap up
    
    fprintf('done running %s\n', mfilename);
    
catch err
    save(save_file);
    ShowCursor;
    KbQueueRelease(); %KbReleaseWait();
    sca; %Screen('Close',p.win);
    rethrow(err);
end

%% subfunctions

%% convert visual angles in degrees to pixels
% function pix = angle2pix(p,ang)
% pixSize = p.screen_width/p.resolution(1);   % cm/pix
% sz = 2*p.screen_distance*tan(pi*ang/(2*180));  %cm
% pix = round(sz/pixSize);   % pix
% return
% end