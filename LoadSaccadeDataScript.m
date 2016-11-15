%% Load Saccade Data Script
% make some .mat files and then load them directly when you want to analyze
% them

%% Load data
%cd ../'GLM LIP'/'Data (processed)'/

% Pat

% fileNameStim = 'pat120512saccademapping1518_stim.mat';
% fileNameSpikes = 'pat120512map1518_spikes.mat';
% goodUnits = 1:3;% goodUnits come from info_lipglm
        
% fileNameStim = 'pat120412saccademapping1423_stim.mat';
% fileNameSpikes = 'pat120412map1423_spikes.mat';
% goodUnits = [1 2];% goodUnits come from info_lipglm

% *** Not sure why there is still an indexing issue with this one!!! (maybe the error is on the other side - the boundary for the target or soemthing
fileNameStim = 'pat111412saccademapping1351_stim.mat';
fileNameSpikes = 'pat111412map1351_spikes.mat';
goodUnits = [1 2];% goodUnits come from info_lipglm

% fileNameStim = 'pat110712saccademapping1402_stim.mat';
% fileNameSpikes = 'pat110712map1400_spikes.mat';
% goodUnits = 1:3;% goodUnits come from info_lipglm

% trial start = 277 (???? what)
% fileNameStim = 'pat110412saccademapping1247_stim.mat';
% fileNameSpikes = 'pat110412map1251_spikes.mat'; %'pat110412saccademapping1247_bp1000spikes.mat';
% goodUnits = 2:3;% goodUnits come from info_lipglm
% 

%% Nancy

% fileNameStim = 'nancy20150304saccademapping1353_stim';
% fileNameSpikes = 'n20150304sacMap_rcgFreeChoice1353_spikes';
% goodUnits = 3:5;% goodUnits come from info_lipglm

% *** NOT WORKING *** 
% fileNameStim = 'nancy20150305saccademapping1554_stim'; 
% fileNameSpikes = 'n20150305rcgfreechoice1353_lip8ch_mtch_d59_d84_spikes';
% goodUnits = [1 4 6];% goodUnits come from info_lipglm

% fileNameStim = 'nancy20150326saccademapping1424_stim';
% fileNameSpikes = '20150326_l1a4_d6300_l3a5_6200_t1257_spikes';
% goodUnits = 16:19;% goodUnits come from info_lipglm

% fileNameStim = 'nancy20150401saccademapping1542_stim';
% fileNameSpikes = '20150401_l1a4_d6300_l3a5_d6200_t1328_spikes';
% goodUnits = [6,9,13,16:20];% goodUnits come from info_lipglm


% %%%%%%% Notes & Parmeters 
% goodUnits are renamed 1...N 

saveOutput = 1; % MUST SAVE FOR NOW if you load without saving you will only have the the last unit's sptimes!!!!!! 
% just load the mat files directly to perform analysis

binSize = .001; % (s) this matters not only in what you load in the data at but also in neuroGLM

%[glmData, glm, unit, stim] = loadSaccadeData(fileNameStim,fileNameSpikes,goodUnits,saveOutput,binSize); %probably should quit spitting out output so we don't use that and get fucked up
% using only I saved mat file is legit
loadSaccadeData(fileNameStim,fileNameSpikes,goodUnits,saveOutput,binSize);


%glmData.param.monkey = 'pat';
%glmData.param.monkey = 'nancy'; % must change within loadSaccadeData for
%now, make this an input when you have time




% when I add in the cell with the RF computation shit fucks up plotting the
% kernels.... whataa ugh





