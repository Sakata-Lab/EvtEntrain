% script EvtFieldEntrainment

% script for event (p-waves) - field entrainment analysis

clearvars;
close all

MyBase = 'MYFILEBASE';

%% params
StartFreqs = [2 4 7 10 15 20 30 50 80 110];
FreqWidths = [2 3 3 5 5 10 20 30 30 40];
nFreqs = length(StartFreqs);
nPhases = 4;
% phase bins ... consistent with Kayser et al JNS 2015
% 1, 135 ~ 180,-180 ~ -135 (EEG trough)
% 2. -135 ~ -45 (rising phase)
% 3. -45 ~ 45 (peak)
% 4. 45 ~ 135 (falling phase)

eWin = 4000; % window size for sleep scoring (in ms)

%% core
FileName = [MyBase,'_EEG_Phase.mat'];

%% loading data
%% sleep score
ScorePath = 'sleepscoring_results.mat'; % sleep scoring data
load(ScorePath, 'Score'); % 1, NREM; 2, REM; 3, AW

%% Pwav data
PwavPath = 'PwaveInfo.mat'; % P-wave information
load(PwavPath, 'DATA');
PwavTime = DATA.PwavTime; % in ms
PwavTime(PwavTime>length(Score)*eWin) = []; % remove P-waves during the period without sleep scores 
PwavState = Score(ceil(PwavTime/eWin)); % allocation of wake-sleep cycle for each P-wave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% phase analysis - preparation
BigMat = zeros(length(PwavTime), 2+nFreqs); % time, state, phases(x nFreqs)
BigMat(:,1) = PwavTime;
BigMat(:,2) = PwavState;

%% phase analysis - core
for f = 1:nFreqs
    MyF = StartFreqs(f);
    fprintf([num2str(MyF),'Hz\n']);

    %% load filtered signals
    FilterFile = [MyBase,'_EEG_',num2str(MyF),'Hz.mat'];
    load(FilterFile); % 'MyEEG','MyF','MyW'

    %% hilbert transform
    tmp = hilbert(MyEEG);
    Amp = real(tmp);
    Phase = circ_rad2ang(angle(tmp)); % in degree

    %% phase
    MyPhase = Phase(BigMat(:,1));
    BigMat(:,2+f) = MyPhase;
end

%% phase analysis by states
Pmat = ones(3, nFreqs); % p value of Rayleigh test
PhaseMat = zeros(nFreqs, nPhases, 3);

for s = 1:3
    idx = find(BigMat(:,2) == s);
    if ~isempty(idx)
        for f = 1:nFreqs
            MyPhase = BigMat(idx, 2+f);

            %% Rayleigh test
            MyPhaseInRadian = circ_ang2rad(MyPhase);
            Pmat(s, f) = circ_rtest(MyPhaseInRadian);

            %% phase-mod summary
            PhaseMat(f,1,s) = length(find(MyPhase>135 | MyPhase<=-135));
            PhaseMat(f,2,s) = length(find(MyPhase>-135 & MyPhase<=-45));
            PhaseMat(f,3,s) = length(find(MyPhase>-45 & MyPhase<=45));
            PhaseMat(f,4,s) = length(find(MyPhase>45 & MyPhase<=135));
        end
    end
end

%% save
save(FileName, 'BigMat', 'Pmat', 'PhaseMat', 'StartFreqs', 'FreqWidths');



