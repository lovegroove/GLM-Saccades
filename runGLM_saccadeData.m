%% Running the GLM with new shit (Round 2)

% load in .mat file, and go to town

% params
plotOn = 1;
saveOn = 1;
saveFigOn = 1;
% conds!!!!
conds = {'targ','sac'};


% Build Deisgn Matrix (Choose your adventure...)
%covName = {'history', 'targon','saccade'};
%covName = { 'targon','saccade','fpon'};
%covName = {'targon','saccade'};
covName = {'history', 'targon','saccade','fpon'};
 

% build design matrix:
[n, nMetaCov] = build_neuroGlmClass(glmData, covName); 

% These are run in build_neuroGlmClass 
%  n.removeConstantCols;
%  n.addBiasColumn;

% visualize_designMatrix(n, false, figDir)
%visualize_designMatrix(n)



%% Do the fitting

% Get the spike trains to regress against
y = n.getBinnedSpikeTrain(glmData.trial, 'sptrain', n.dm.trialIndices);

% fit glm to get maximumu likelihood weights (wml):
[wml, wvar] = fitGlm(n, y); % wrapper func

 % combine weights for particular covariates so we can plot the kernels
ws = n.combineWeights(wml);
wvar = n.combineWeights(wvar);

% % STORE EVERYTHING BACK IN glm:
% glm.glmData  = glmData;
% glm.n        = n;
% glm.nMetaCov = nMetaCov;
% glm.y        = y;
% glm.wml      = wml;
% glm.wvar     = wvar;

%% Plot kernels

fig = figure(2913); clf;
nCovar = numel(n.covar);
for kCov = 1:nCovar
    label = n.covar(kCov).label;
    subplot(nCovar, 1, kCov);
    plot(ws.(label).tr, ws.(label).data, ...
        ws.(label).tr, ws.(label).data+sqrt(wvar.(label).data), '--', ...
        ws.(label).tr, ws.(label).data-sqrt(wvar.(label).data), '--')
    title(label);
end

if saveFigOn
    saveFigName = strcat(glmData.param.fileNameStim,'_unit',num2str(glmData.param.unit),'_kernels','.fig');
    savefig(saveFigName)
end

%% prediction
X = n.dm.X; % get out of object
% size(wml)
% size(X)
yPred = exp(X * wml) / glm.binSize; % Divide by bin size to get Hz

% reconstruct trial format from collapsed form in design matrix
c = yPred;
nt = length(n.dm.trialIndices); %glmData.nTrials; % use number of trials from the design matrix just in case they drop one
nc = length(c(:))/nt;
yPredTrials = reshape(c,[nc,nt])'; % transpose... because reshape is SP00KY and you have to enter shit backwards
size(yPredTrials) % double check design matrix

%% Responses/predictions conditioned on different events (aka align stuff)

% choose your adventure... or loop through
%cond = 'targ'; % 'targ','sac'...


nConds = numel(conds);
for kCond = 1:nConds

    cond = conds{kCond};
    
%%%%% EVENTS %%%%% (glm.bTarg, glm.bSacc, glm.bFpon)
if strcmp(cond,'targ')
    alignedInd = glm.bTarg; % bins (ms if you chose it)
elseif strcmp(cond,'sac')
    alignedInd = glm.bSacc;
else
    error('No such condition')
end

%%%%% WINDOWS %%%%%
winTarg = [-100 1000]; 
winSac = [-1000 100]; % in ms (or whatever bin size you chose for the GLM data)
winFpon = [-100 2500];
%%%%%%%%%%%%%%%%%%%

if strcmp(cond,'targ')
    win = winTarg; % dumb but easy
elseif strcmp(cond,'sac')
    win = winSac; % dumb but easy
else
    error('No such condition')
end

% Get predictions actually aligned precisely to the events (not the most efficient way of course)
yPredAligned = cell(1,size(yPredTrials,1));
for iTrial = 1:size(yPredTrials,1) % maybe add in something like if alignedInd + 100 > 2500 for the outer boundery and just make it (end) if excedes that, but then will the dimensions be off?
    % or just throw away those few trials that went hella long - this seemed to fix it for most of them
    if alignedInd(iTrial) + winSac(2) > winFpon(2) % maybe need to add a lower bound as well...
        % do nothing, these trials are bogus
    else
        yPredAligned{iTrial} = yPredTrials(iTrial,((alignedInd(iTrial)+win(1)):alignedInd(iTrial)+win(2))); % IMPORTANT: indexing issue might be present if the saccade is so late in the trial taht when we add 100ms it's outside the the 2600ms im taking to stick in to model
    end
end
yPredAligned = cell2mat(yPredAligned'); % dumb quick fix,must flip because how I stuck it the cell array (OVERWRITING yPredAligned BE AWARE)

mPred = mean(yPredAligned);

xInds = win(1):win(2);


% Make data PSTHs

% data psth params
bs          = 1e-3;             % bin size (s) % might have to be careful with all this
skern       = ones(20,1);       %smoothing kernel

% THE UNIT!!!! %%%% BE CAREFUL - if you don't specify the unit, it will collapse all of the units into one and FUCK YOUR WORLD UP - this should be avoided now by loading in a mat file and every one is a different unit
%unit(un).spikeTimes
%%%%%%%%%%%%%%%%%%%%%%%%%
un = glmData.param.unit;% % DON'T FUCK UP - this should not change now because it just comes from the mat file you chose to load in via loadSaccadeData
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% WINDOWS %%%%%
winTarg = winTarg * glm.binSize;
winSac = winSac * glm.binSize; % converting ms(bin size above) to s, just for old school PSTH plotting 
winFpon = winFpon * glm.binSize;
%%%%%%%%%%%%%%%%%%%

% Note: These events come from the stim struct... and they are times, not the binned ones above in the glm struct
% stim.tTarg, stim.tSacc, stim.tFpon

ev = stim.tTarg; % targ
sptimes = unit(un).spikeTimes;
[m1,s,bc1] = pdsa.eventPsth(sptimes, ev, winTarg, bs, skern);

ev = stim.tSacc;  % sac
sptimes = unit(un).spikeTimes;
[m2,s,bc2] = pdsa.eventPsth(sptimes, ev, winSac, bs, skern);

ev = stim.tFpon; % fpon
sptimes = unit(un).spikeTimes;
[m3,s,bc3] = pdsa.eventPsth(sptimes, ev, winFpon, bs, skern);

% just for double checking stuff
% figure;
% subplot(311),plot(bc1,m1)
% title('Targ-on conditioned PSTH')
% subplot(312),plot(bc2,m2) 
% title('Saccade conditioned PSTH')
% subplot(313), plot(bc3,m3)
% title('FP-on conditioned PSTH')

% PLOTTING FOR REALZ

% smoothing (for the prediction)
span = 10; % 1 = no smoothing, 5 = default

% plot target aligned PSTHs
if plotOn && strcmp(cond,'targ')
figure
hold all
plot(bc1*1000,m1,'linewidth',2)
plot(xInds,smooth(mPred,span),'linewidth',2)
legend('Data','Model')
title('Target conditioned PSTHs')
ylabel('Firing Rate')
xlabel('Time (ms)')
hold off
saveFigName = strcat(glmData.param.fileNameStim,'_unit',num2str(glmData.param.unit),'_targ','.fig');
if saveFigOn
    savefig(saveFigName)
end
end

%plot saccade aligned PSTHs
if plotOn && strcmp(cond,'sac')
figure
hold all
plot(bc2*1000,m2,'linewidth',2)
plot(xInds,smooth(mPred,span),'linewidth',2)
legend('Data','Model')
title('Saccade conditioned PSTHs')
ylabel('Firing Rate')
xlabel('Time (ms)')
hold off
saveFigName = strcat(glmData.param.fileNameStim,'_unit',num2str(glmData.param.unit),'_sac','.fig');
if saveFigOn
    savefig(saveFigName)
end
end


% Saving

if saveOn && strcmp(cond,'targ')
saveFileName = strcat(glmData.param.fileNameStim,'_unit',num2str(glmData.param.unit),'_modeled_targ','.mat');

dataX = bc1*1000;
dataPSTH = m1;
modelX = xInds;
modelPSTH = mPred;

save(saveFileName,'dataPSTH','dataX','modelPSTH','modelX','ws')

end

if saveOn && strcmp(cond,'sac')
saveFileName = strcat(glmData.param.fileNameStim,'_unit',num2str(glmData.param.unit),'_modeled_sac','.mat');

dataX = bc2*1000;
dataPSTH = m2;
modelX = xInds;
modelPSTH = mPred;

save(saveFileName,'dataPSTH','dataX','modelPSTH','modelX','ws')

end




end % end Conds loop


%% TO DO
% add cross validation by resampling subsets of trials from the design matrix

%% Importing lots of units to look at means and plot everything together
% make sure you're in the right directory...

% targ files
targFiles = dir('*_modeled_targ.mat'); 
numFiles = length(targFiles);
targData = cell(1, numFiles);
targPSTH = cell(1, numFiles);
targModel = cell(1, numFiles);
targK = cell(1, numFiles);
sacK = cell(1, numFiles);
targDataX = cell(1, numFiles);
targModelX = cell(1, numFiles);

% change from cells or wait and to cell2mat
for k = 1:numFiles 
  targData{k} = importdata(targFiles(k).name); 
  
  targPSTH{k} = targData{k}.dataPSTH;
  targModel{k} = targData{k}.modelPSTH;
  targDataX{k} = targData{k}.dataX;
  targModelX{k} = targData{k}.modelX;
  
  targK{k} = targData{k}.ws.targon.data; % (ws struct doesn't matter if it comes from targ or sac data because it's not aligned to those, it's jsut the kernels) 
  sacK{k} = targData{k}.ws.saccade.data;
  histK{k} = targData{k}.ws.hist.data;
  
end
targPSTH = cell2mat(targPSTH'); % overwriting here... careful
targModel = cell2mat(targModel');
targK = cell2mat(targK)'; % be aware... different transposing in these
sacK = cell2mat(sacK)';
targDataX = cell2mat(targDataX'); % overwriting here... careful
targModelX = cell2mat(targModelX');
histK = cell2mat(histK)';

% sac files
sacFiles = dir('*_modeled_sac.mat'); 
numFiles = length(sacFiles);
sacData = cell(1, numFiles);
sacPSTH = cell(1, numFiles);
sacModel = cell(1, numFiles);
sacDataX = cell(1, numFiles);
sacModelX = cell(1, numFiles);

% change from cells or wait and to cell2mat
for k = 1:numFiles 
   sacData{k} = importdata(sacFiles(k).name); 

   sacPSTH{k} = sacData{k}.dataPSTH;
   sacModel{k} = sacData{k}.modelPSTH;
   sacDataX{k} = sacData{k}.dataX;
   sacModelX{k} = sacData{k}.modelX;

   % don't need to do the kernels again, thsoe are the same
   
end
sacPSTH = cell2mat(sacPSTH'); % overwriting here... careful
sacModel = cell2mat(sacModel');
sacDataX = cell2mat(sacDataX'); % overwriting here... careful
sacModelX = cell2mat(sacModelX');


%%
% all kernels
fig = figure;
hold all
plot(targData{k}.ws.targon.tr,targK) % x-axis: all tr should be the same if you were consistent so, picking the x range from any unit should be fine
plot(targData{k}.ws.targon.tr,mean(targK),'--k','linewidth',2)
title('Target On Kernels')
ylabel('Gain')
xlabel('Time (ms)')
hold off

fig = figure;
hold all
plot(targData{k}.ws.saccade.tr,sacK) % x-axis: all tr should be the same if you were consistent so, picking the x range from any unit should be fine
plot(targData{k}.ws.saccade.tr,mean(sacK),'--k','linewidth',2)
title('Saccade Kernels')
ylabel('Gain')
xlabel('Time (ms)')
hold off

fig = figure;
hold all
plot(targData{k}.ws.hist.tr,histK) % x-axis: all tr should be the same if you were consistent so, picking the x range from any unit should be fine
plot(targData{k}.ws.hist.tr,mean(histK),'--k','linewidth',2)
title('History Kernels')
ylabel('Gain')
xlabel('Time (ms)')
xlim([0 30]) % for better visualization
hold off

%%
%print(fig,'-dpdf','history_kernels','-bestfit')


%% mean data/prediction psths (across all units)

span = 10;

figure
hold all
plot(mean(targModelX),smooth(mean(targModel),span),'linewidth',2) % used mean x, but they should all be the same
plot(mean(targDataX),mean(targPSTH),'linewidth',2)
title('Target locked, population average')
ylabel('Firing rate (Hz)')
xlabel('Time (ms)')
xlim([-100 700])
legend('Model Prediction','Data','Location','Best')
hold off

figure
hold all
plot(mean(sacModelX),smooth(mean(sacModel),span),'linewidth',2) % used mean x, but they should all be the same
plot(mean(sacDataX),mean(sacPSTH),'linewidth',2)
title('Saccade locked, population average')
ylabel('Firing rate (Hz)')
xlabel('Time (ms)')
legend('Model Prediction','Data','Location','Best')
xlim([-700 100])
hold off


%% if you just want to inspect individual units
% to find what unit you want
% sacFiles(i).name or targFiles(i).name

% these are unsmoothed this far (model-wise)
span = 10;

% for plotting PSTHs
% x's depend on events
ylimMin = 0;
ylimMax = 45;

% for the kernels
ylimMinK = 0;
ylimMaxk = 2.5;

% Example 1 - (nancy 1424 - unit 4; i = 7), Example 2 - (nancy 1542 - unit 5; i = 9), Example 3 - (pat 1518 - unit 1; i = 19)
i = 9;
fig = figure;
subplot(221)
plot(targData{i}.ws.targon.tr,targK(i,:),'linewidth',2) % targ kernel
ylim([ylimMinK ylimMaxk])
xlabel('Time (ms)')
ylabel('Gain')
subplot(222)
plot(targData{i}.ws.saccade.tr,sacK(i,:),'linewidth',2) % saccade kernel
ylim([ylimMinK ylimMaxk])
xlabel('Time (ms)')
ylabel('Gain')
subplot(223)
hold all
plot(targModelX(i,:),smooth(targModel(i,:),span),'linewidth',2) % targ model
plot(targDataX(i,:),targPSTH(i,:),'linewidth',2) % targ psth
axis([-100 700 ylimMin ylimMax])
xlabel('Time (ms)')
ylabel('Firing Rate')
legend('Model Prediction','Data','Location','Best')
hold off
subplot(224)
hold all
plot(sacModelX(i,:),smooth(sacModel(i,:),span),'linewidth',2) % sac model
plot(sacDataX(i,:),sacPSTH(i,:),'linewidth',2) % sac psth
axis([-700 100 ylimMin ylimMax])
xlabel('Time (ms)')
ylabel('Firing Rate')
%legend('Model Prediction','Data','Location','Best')
hold off
%title(targFiles(i).name)

%%
print(fig,'-dpdf','example2_summary','-bestfit')

%% make and save summary figures for the whole set

% these are unsmoothed this far (model-wise)
span = 10;


for i = 1:numFiles


% for plotting PSTHs
% x's depend on events (just hard coding those)
ylimMin = 0;
maxY = max([targModel(i,:) sacModel(i,:) sacPSTH(i,:) targPSTH(i,:)]);
ylimMax = maxY + (maxY*.1);


% for the kernels
minyK = min([targK(i,:) sacK(i,:)]);
ylimMinK = minyK - (abs(minyK)*.1);
maxyK = max([targK(i,:) sacK(i,:)]);
ylimMaxk = maxyK + (abs(maxyK)*.1);


fig = figure;
subplot(221)
plot(targData{i}.ws.targon.tr,targK(i,:),'linewidth',2) % targ kernel
ylim([ylimMinK ylimMaxk])
xlabel('Time (ms)')
ylabel('Gain')
subplot(222)
plot(targData{i}.ws.saccade.tr,sacK(i,:),'linewidth',2) % saccade kernel
ylim([ylimMinK ylimMaxk])
xlabel('Time (ms)')
ylabel('Gain')
subplot(223)
hold all
plot(targModelX(i,:),smooth(targModel(i,:),span),'linewidth',2) % targ model
plot(targDataX(i,:),targPSTH(i,:),'linewidth',2) % targ psth
axis([-100 700 ylimMin ylimMax])
xlabel('Time (ms)')
ylabel('Firing Rate')
legend('Model Prediction','Data','Location','Best')
hold off
subplot(224)
hold all
plot(sacModelX(i,:),smooth(sacModel(i,:),span),'linewidth',2) % sac model
plot(sacDataX(i,:),sacPSTH(i,:),'linewidth',2) % sac psth
axis([-700 100 ylimMin ylimMax])
xlabel('Time (ms)')
ylabel('Firing Rate')
%legend('Model Prediction','Data','Location','Best')
hold off

%
fileName = strcat(targFiles(i).name(1:end-11),'_summary');
print(fig,'-dpdf',fileName,'-bestfit')


end




