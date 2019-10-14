% script EEGfiltering

% script to create bandpass filtered EEG signals for event-field
% entrainment analysis

MyBase = 'MYFILEBASE';
nChs = 49; % the number of channels. This needs to be updated based on channel configuration
EEGch = 1; % the position of the EEG channel. This needs to be updated based on channel configuration
Fs = 1000; % sampling rate

%% params -- key parameters to define filtering bandwidthes
StartFreqs = [2 4 7 10 15 20 30 50 80 110];
FreqWidths = [2 3 3 5 5 10 20 30 30 40];
nFreqs = length(StartFreqs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load sleep scoring to know the length of the recording
% this section can be replaced depending on datasets.
load('sleepscoring_results.mat', 'Score');
FinishPoints = length(Score)*4*Fs; % 4 means 4-sec bin for sleep scoring.

%% load EEG signals
EEGpath = [MyBase, '.eeg'];
tmp = LoadEEGseg(EEGpath, nChs, 1, FinishPoints); % this load data as a 16bit vector
MyEEGorg = tmp(EEGch, :);
clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filtering - core
for f = 1:nFreqs
    % bandpass filtering of EEG signals
    MyF = StartFreqs(f);
    fprintf([num2str(MyF),'Hz\n']);
    
    if fopen([MyBase,'_EEG_',num2str(MyF),'Hz.mat']) == -1
        MyW = FreqWidths(f);
        %% create filter
        MyFilt = designfilt('bandpassfir','StopbandFrequency1',MyF-1,'PassbandFrequency1',MyF, ...
            'PassbandFrequency2',MyF+MyW,'StopbandFrequency2',MyF+MyW+1, ...
            'StopbandAttenuation1', 50,'PassbandRipple', 0.01,'StopbandAttenuation2', 50,...
            'SampleRate', Fs, 'DesignMethod', 'kaiserwin');
        
        %% filtering
        MyEEG = filtfilt(MyFilt, MyEEGorg);
        
        %% save
        FileName = [MyBase,'_EEG_',num2str(MyF),'Hz.mat'];
        save(FileName,'MyEEG','MyF','MyW');
    end
end

