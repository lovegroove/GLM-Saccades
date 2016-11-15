function [glmData, glm, unit, stim] = loadSaccadeData(fileNameStim,fileNameSpikes,goodUnits,saveOutput,binSize)

% Script for dealing with Jake and Leor's saccade data (or function
% perhaps)

% NEEDS TONS OF CLEANING UP - THIS WAS A SUPER SLOPPY FIRST PASS UNDER DURESS

% load data:
filePath = '/Users/erichart/Documents/Data/LIP GLM saccade data (select cells from Jake and Leor)/just_the_good_stuff/'; %hard coded for now
% fileNameStim = 'pat120512saccademapping1518_stim.mat';
% fileNameSpikes = 'pat120512map1518_spikes.mat';
fullFileNameStim = fullfile(filePath,fileNameStim);
fullFileNameSpikes = fullfile(filePath,fileNameSpikes);


%regexp  % make a loop to go through all the files in the folder and
%appropriately read in the spikes and stim files, don't need to do this
%right now let's lone at a time

% POTENTIAL INPUTS
% Good units = [1 2 3];
% Win to count spikes over ???
%DEFAULTS
if nargin <5
    glm.binSize             = 1e-3;             % bin size (s)

else
    glm.binSize = binSize;
end

if nargin < 4
   saveOutput = 0; % default, don't save 
end
% if nargin <4
% win                 = [-.1 1.5];        % window (s) around event (e.g
% tTarg..) % npw defined below
% end
if nargin <3
    goodUnits = 1;
end



% load('pat120512saccademapping1518_stim.mat');
% load('pat120512map1518_spikes.mat');

stim            = load(fullFileNameStim,'-mat');
spikes          = load(fullFileNameSpikes,'-mat');
if isfield(spikes, 'spikes')
    spikes = spikes.spikes; % cause old versions are wonky
end

% Leor's code is in processData.m (ganking some of it to get the basics)

%% THE STIM STRUCT: 

% standardize across versions:
if min(size(stim.timing.targon))==2
    targon = stim.timing.targon(:,1);
else
    targon = stim.timing.targon;
end

% timing:
stim.tTarg   = (stim.timing.plxstart(:) + targon(:));    % get times of target on
stim.tSacc   = (stim.timing.plxstart(:) + stim.timing.choice(:));    % get times of saccades
stim.tFpon = (stim.timing.plxstart(:) + stim.timing.fpon(:,1));

% Vis or Mem guided saccades
%stim.task==1 %Vis
%stim.task==2 %Mem

% good trials indices:
stim.goodIdx     = stim.goodTrials & ~isnan(stim.timing.plxstart) & ~isnan(targon(:));
% % if isfield(ds.info, 'trialStart')
% %     stim.goodIdx(1:ds.info.trialStart) = false; % *** need to get trial start and end
% % end
% % if isfield(ds.info, 'trialEnd')
% %     stim.goodIdx((ds.info.trialEnd+1):end) = false;
% % end


%% THE UNIT STRUCT: (***have to switch out stuff Leor assinged in ds)

% prepare a new struct for the units:
unit = struct('ADfreq', [], 'channel', [], 'id', [], ...
            'snr', [], 'spikeTimes', [], ...
            'timeaxis', [], 'waveform', [], ...
            'rf', []);

% populate unit struct:
neuronId_list   = goodUnits; %*** I think this was just assigned from inspect, just passing stuff along to be conisistent 

nNeurons        = numel(neuronId_list);
for iN = 1:nNeurons
    neuronId = neuronId_list(iN);
    unit(iN).number     = goodUnits(iN);
    unit(iN).channel    = spikes.channel(neuronId);
    unit(iN).id         = neuronId;
    unit(iN).snr        = spikes.snr(neuronId);
    unit(iN).spikeTimes = spikes.time(spikes.id==neuronId);
    try unit(iN).ADfreq     = spikes.ADfreq; catch; end    % not all datafiles have ADfreq saved but it's ok cause spike times have already been computed.
    try unit(iN).timeaxis   = spikes.timeaxis; catch; end    % not all dataset have timeaxis but it's ok. It's just time for the wavefors.
    unit(iN).waveform   = spikes.waveform(spikes.id==neuronId,:);
end

%% estimate RF per unit:
% % 
% % % settting binsize to equal entire range, i.e. 1 bin within which to
% % % count spikes:
% % win         = [0 .300];         % window (s) around event (e.g tTarg..)
% % binSz       = range(win);       % bin size (s)
% % threshLevel = 0.4;              % threshold level for RF relative to gaussain fit peak.
% % xRange      = [-26 26];
% % yRange      = [-16 16];
% % 
% % % prealloc:
% % badUnit_rfTooBig    = false(nNeurons,1);
% % badUnit_rfTooSmall  = false(nNeurons,1);
% % 
% % % go:
% % for iN = 1:nNeurons
% %     neuronId = neuronId_list(iN);
% %     fprintf('Estimating RF for unit is %d', neuronId);
% %     sptimes = spikes.time(spikes.id==neuronId); % spike times
% %     % bin spikes:
% %     binnedSpikes    = pdsa.binSpTimes(sptimes, stim.tTarg, win, binSz);
% %     % get neuron RF matrix:
% %     [imageRf, thetaHat] = getNeuronRf(binnedSpikes, stim.targ1XY, xRange, yRange); % get rf matrix of neuron
% %     % threshold rf:
% %     imageThresh = imageRf > threshLevel;      % thresholded RF
% %     
% %     % CHECK IF RF IS DECENT. 
% %     % I'm doing something wierd here but it's good for now. 
% %     % any unit that has an imageThresh that takes over a higher percentage
% %     % of space than its (1-threshold) value is excluded. (or if too little)
% %     if mean(imageThresh) > (1-threshLevel)
% %         % i.e. rf is too big
% %         badUnit_rfTooBig(iN) = true;
% %         continue;
% %     elseif mean(imageThresh) < 0.1
% %         % i.e. rf is too small
% %         badUnit_rfTooSmall(iN) = true;
% %         continue;
% %     else
% %         badUnit_rfTooBig(iN)    = false;
% %         badUnit_rfTooSmall(iN)  = false;
% %     end
% %     
% %     % determine which trial was inRf and which was not:
% %     idxInRf    = isInRf(stim.targ1XY(:,1), stim.targ1XY(:,2), imageThresh);
% %     % store in unit:
% %     unit(iN).rf.binnedSpikes    = binnedSpikes;
% %     unit(iN).rf.thetaHatLabels  = {'amp', 'xMu', 'yMu', 'sigma'};
% %     unit(iN).rf.thetaHat        = thetaHat;
% %     unit(iN).rf.xRange          = xRange;
% %     unit(iN).rf.yRange          = yRange;
% %     unit(iN).rf.imageRf         = imageRf;
% %     unit(iN).rf.threshLevel     = threshLevel;
% %     unit(iN).rf.imageThresh     = imageThresh;
% %     unit(iN).rf.idxIn           = idxInRf;
% %     unit(iN).rf.idxOut          = ~idxInRf;
% %     fprintf('\tDone!\r')
% % end    
% % 
% % % If RF is too big/small, off with its head!
% % goodIdx             = ~badUnit_rfTooBig & ~badUnit_rfTooSmall;
% % goodUnits   = goodUnits(goodIdx); % ?????? THIS DANGEROUS?
% % unit(~goodIdx)      = [];
% % 
% % %% Enter Leor's insane pipeline (well already in it) (trying to extract the essence and only use what's necessary and make a simpler cleaner version)
% % % % %covName = {'history', 'targon','saccade'}; 
% % % % [glm] = run_glm_onUnit(unit, stim.tTarg, stim.tSacc, covName); % THIS IS THE FUNCTION THIS JUNK WAS IN
% % 
% % assert(numel(stim.tTarg)==numel(stim.tSacc), 'ERROR: number of elements doesnt match between tTarg and tSacc');
% % 
% % glm = struct('nTrials', [], 'bins', [], 'binnedSpikes', [], 'bTarg', [], 'bSacc', [], ...
% %     'glmData', [], 'n', [], 'nMetaCov', [], 'y', [], 'wml', [], 'wvar', []);
% % 
% % 
% % glm.nTrials         = size(stim.tTarg(:),1);
% % if glm.nTrials < 5
% %     return;
% % end
% % 
% % %win                 = [-.1 1.5];        % window (s) around event (e.g tTarg..)
% % glm.binSize             = 1e-3;             % bin size (s)
% % 
% % 
% % % go through units save each good unit to a different matl file***
% % for un = 1:length(goodUnits)
% % % % [spcnt, bc]         = pdsa.binSpTimes(unit(un).spikeTimes, stim.tTarg, win, binSize);  % ****** This win is only using taron specifications  % is this ok since we are jsut using this to get rid of bad trials and we plot things in teh otehr script with appropriate windows
% % % % %[spcnt, bc]         = pdsa.binSpTimes(unit.spikeTimes, tTarg, win, binSize);
% % % % 
% % % % glm.bins            = bc;
% % 
% % %glm.binnedSpikes = spcnt;
% % % get rid of NaNs (bad trials now - trials with no spikes or trials with no saccade)
% % %glm.binnedSpikes = spcnt(~any(isnan(spcnt),2),:); % GET RID OF NAN TRIALS
% % % % goodTrials = ~any(isnan(spcnt)==1,2); % ROWS WITHOUT NaNs could have used stim.goodIdx or stim.goodTrials !!!???)
% % % % glm.binnedSpikes = spcnt(goodTrials,:);
% % 
% % % WIP(********* % IMPORTANT WINDOW AND STIMULUS EVENTS: targon and saccade
% % win = [0 .5]; % WINDOW FOR targon
% % ev = stim.tTarg;
% % % stim event
% % [spcnt, bc]         = pdsa.binSpTimes(unit(un).spikeTimes,ev, win, glm.binSize);% what the hell to do about window for counting spikes here, just try to count for the whole trial
% % glm.binnedSpikes = spcnt(stim.goodIdx,:);
% % 
% % % don't have to put this stuff in a loop but whatever (for now)
% % glm.tTarg = stim.tTarg(stim.goodIdx);
% % glm.tSacc = stim.tSacc(stim.goodIdx);
% % glm.nTrials = length(glm.tSacc);
% % glm.task = stim.task(stim.goodIdx);
% % 
% % win = [-1.5 .5]; % WINDOW FOR saccade!
% % ev = stim.tSacc;
% % 
% % % % glm.inRF = unit(un).rf.idxIn(stim.goodIdx);
% % % % glm.outRF = unit(un).rf.idxOut(stim.goodIdx); % or just use indexing for In
% % % % 
% % 
% % % This has the new scrubbed trials now (THIS ALSO HAS A WIN THAT HAS TO DO
% % % WITH THE TARGET****
% % glm.bTarg           = repmat((0-win(1)) ./ glm.binSize, [glm.nTrials,1]);         % bin at which target came on
% % glm.bSacc           = round((glm.tSacc - glm.tTarg) ./ glm.binSize);   % bin at which saccade took place
% % 
% % 
% % % IMPORTANT: for each good unit, put binned spikes and events into trial struct for GLMing
% % % convert binned spikes from data into glmData struct:  
% % glmData = convert_binnedSpikes_to_glmData(glm.binnedSpikes, glm.bTarg, glm.bSacc);  
% % 
% % %**** need to save glmData.trial.sptrain for both sac and event locked
% % 
% % 
% % %save out the different units
% % % glmData.param.inRf = glm.inRF;
% % % glmData.param.outRf = glm.outRF;
% % glmData.param.monkey = 'pat'; % sloppy manual stuff for now
% % glmData.param.unit = un;
% % if saveOutput
% %     saveFile = strcat(fileNameStim,'_',num2str(glmData.param.unit),'.mat');
% %     save(saveFile,'glmData','glm','unit','stim')
% % end


%%  Get spike times into trial format for each good unit

assert(numel(stim.tTarg)==numel(stim.tSacc), 'ERROR: number of elements doesnt match between tTarg and tSacc');

%glm = struct('nTrials', [], 'bins', [], 'binnedSpikes', [], 'bTarg', [], 'bSacc', [], ...
%    'glmData', [], 'n', [], 'nMetaCov', [], 'y', [], 'wml', [], 'wvar', []);

glm.nTrials         = size(stim.tTarg(:),1);
if glm.nTrials < 5
    return;
end


glm.tTarg = stim.tTarg(stim.goodIdx);
glm.tSacc = stim.tSacc(stim.goodIdx);
glm.task = stim.task(stim.goodIdx);
glm.tFpon = stim.tFpon(stim.goodIdx);

% go through units save each good unit to a different matl file***
for un = 1:length(goodUnits)
%     win = [-.1 1.5]; % let's just make the window the entire trial
%     ev = stim.tTarg; % or glm.tTarg (NaNs already removed...) or stim.tTarg
%     win = [-0.1 2]; % tTarg
    % %%%%% IMPORTANT %%%%% make sure you know what event you are counting spikes around and what window and make sure you take thing into accoutn later
    glm.winModel = [-0.1 2.5]; % SUPER SUPER TRICKY - PICKING WINDOW AND EVENT TO ALIGN TO HERE - to count spikes over the period for which you stick into the GLM, this can totally change your kernels - I picked what seemed to capture the whole trial and be pretty stable - changing this now though, FYI trial lengths are about 2.5s
    ev = stim.tFpon; % THAT WINDOW CHOICE IS IMPORTANT - and remember later that you might want to carve out the interesting part of the trial for comparison to the conditioned PSTHs
    
    [spcnt, glm.bc]         = pdsa.binSpTimes(unit(un).spikeTimes,ev, glm.winModel, glm.binSize);
    glm.binnedSpikes = spcnt(stim.goodIdx,:);
    
    % NOTE: the win for bTarg below is dependent on the event used to count spikes is the targ on - these dependencies are terrible... ugh, I hate usign cobbled together code
    % THESE BIN CALCULATIONS DEPEND ON WHICH EVENT YOU CHOSE FOR REFERENCE ABOVE
    % FP ON IS NOW THE ZERO BIN
    
   
    
    glm.bFpon      = repmat((0-glm.winModel(1)) ./ glm.binSize, [glm.nTrials,1]); % bin at which the fixation is on
% %     %glm.bTarg           = repmat((0-glm.winModel(1)) ./ glm.binSize, [glm.nTrials,1]);         % bin at which target came on
% %     %glm.bSacc           = round((glm.tSacc - glm.tTarg) ./ glm.binSize);   % bin at which saccade took place
%     glm.bTarg           = round((glm.tTarg - glm.tFpon) ./ glm.binSize);          % bin at which target came on
%     glm.bSacc           = round((glm.tSacc - glm.tFpon) ./ glm.binSize);   % bin at which saccade took place
%     

    glm.binZero = (0-glm.winModel(1)) / glm.binSize; % THIS IS VITAL
    % IMPORTANT: adding in a offset incase you zeroBin isn't actually zero
    % (THIS IS A BIG DEAL AND FUCKED ME UP)
    glm.bTarg           = round((glm.tTarg - glm.tFpon) ./ glm.binSize) + glm.binZero;          % bin at which target came on
    glm.bSacc           = round((glm.tSacc - glm.tFpon) ./ glm.binSize) + glm.binZero;   % bin at which saccade took place
    
    
    glmData = convert_binnedSpikes_to_glmData(glm.binnedSpikes, glm.bTarg, glm.bSacc, glm.bFpon); % this is the important stuff that must be passed in
        
    
    
    
    % SAVING STUFF
    % also saving this out above glm.winModel so we have the overall window handy
    %glmData.param.monkey = 'pat'; % sloppy manual stuff for now
    glmData.param.monkey = 'pat';
    glmData.param.unit = un;
    glmData.param.fileNameStim = fileNameStim;
    if saveOutput
        saveFile = strcat(fileNameStim,'_unit',num2str(glmData.param.unit),'.mat');
        save(saveFile,'glmData','glm','unit','stim')
    end
    
end% end unit loop



end
% This should be ready for neuroGLM now (switch to GLMing Script)






