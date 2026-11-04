%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      calculating CD and ERPs                                  %
% Author:        Anna-Lena Schubert, Mareike Hülsemann, Henrike Jungeblut %
% Last Modified: 09/2026                                                  %
% Note:          this script is not reproducible as preprocessed EEG data %
%                is not provided in the code package                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------- Set Up ------------------------------------------------------

% define experiment and task specific settings
marker_start  = 'S  6'; marker_start_num = str2double(marker_start(3:end));
event_indices = {'S 51' 'S 52' 'S 53' 'S 54'}; % Cue: {'S 31' 'S 32' 'S 33' 'S 34'};
epoch_limits    = [-2.0 2.0]; % sec
trial_limits    = [-750 1000]; % ms
baseline_limits = [-750 -250]; % ms

% define PATHs
% notebook
% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});

PATH_EEGLAB    = fullfile(PATH_MAIN, 'Code_package', 'Tools', 'eeglab2025.0.0');

PATH_IN        = fullfile(PATH_MAIN, 'preprocessed');
PATH_OUT       = fullfile(PATH_MAIN, 'tf_power');
PATH_OUT_STATS = fullfile(PATH_MAIN, 'stats');
PATH_SUB_CONNY_MATS = fullfile(PATH_MAIN, 'subject_conny_matrices/');
PATH_CONNYDEGREE_TRACES = fullfile(PATH_MAIN, 'Code_package', 'subject_connydegree_traces/');
PATH_ERP = fullfile(PATH_MAIN, 'Code_package', 'erp/');

% define FILE names
FILE_MAIN = ['_epoched_target_' num2str(epoch_limits(1)) '_' num2str(epoch_limits(2))];


% load list of preprocessed EEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not reproducible as preprocessed EEG data is not provided               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allSets = dir(fullfile(PATH_IN, '*.set'));

% initialize preprocessing stats
artefact_varnames = {...
    'id', ...
    'task', ...
    'srate_for_TF', ...
    'nb_trials', ...
    'pnts_trial', ...
    'nb_trials_clean', ...
    'nb_trials_rejected', ...
    'prct_trials_kept', ...
    'prct_trials_rejected', ...
    'nb_trials_rejected_total', ...
    'nb_trials_rejected_thre', ...
    'nb_trials_rejected_prob', ...
    'nb_trials_rejected_kurt', ...
    'nb_trials_rejected_unique_thre', ...
    'nb_trials_rejected_unique_prob', ...
    'nb_trials_rejected_unique_kurt', ...
    'nb_trials_rejected_inTP_notK', ...
    'nb_trials_rejected_inTK_notP', ...
    'nb_trials_rejected_inPK_notT', ...
    'nb_trials_rejected_inTPK'};
artefact_epochs_varnames = {'trials_rejected_unique_thre', ...
    'trials_rejected_unique_prob', ...
    'trials_rejected_unique_kurt', ...
    'trials_rejected_inTP_notK', ...
    'trials_rejected_inTK_notP', ...
    'trials_rejected_inPK_notT', ...
    'trials_rejected_inTPK', ...
    'trials_rejected_all'};



%% ---------- start EEGLAB ------------------------------------------------
addpath(PATH_EEGLAB);
eeglab; %close(gcf);
pop_editoptions('option_single', 0, ...
    'option_savetwofiles', 0, ... % new default: FALSE
    'option_computeica',   1, ... % one only has EEG.icaact if TRUE
    'option_checkversion', 0);


% ---------- Loop Over Subjects ------------------------------------------
n = randi(length(allSets));
for n=2:length(allSets)
    disp(['**************** Started Set ' num2str(n) ' of ' num2str(length(allSets)) '; ' allSets(n).name '; ' datestr(datetime()) ' ****************'])

    % initialise
    artefact_overview  = nan(1,  numel(artefact_varnames));
    artefact_epochs    = cell(1, numel(artefact_epochs_varnames));


    % ---------- Loading data  --------------------------------------------
    % load data set
    EEG = pop_loadset('filename',allSets(n).name,'filepath',PATH_IN);


    % ---------- Preprocessing for Time-Frequency-Analysis ----------------

    % save stats of rejection
    col=    1; artefact_overview(1,col) = str2double(allSets(n).name(1:7));
    col=col+1; artefact_overview(1,col) = marker_start_num;
    col=col+1; artefact_overview(1,col) = EEG.srate;

    % extract trials (no baseline removal)
    EEG = pop_epoch(EEG, event_indices, epoch_limits);

    % detect artefactual epochss
    EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan ,-500,500,0,EEG.xmax,0,0);
    EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 5, 2, 0, 0, 0, [], 0); % joint probability, according to Makotos recommendations and Delorme et al. 2007 (5SD; more strict for EEG)
    EEG = pop_rejkurt(  EEG, 1, 1:EEG.nbchan, 5, 2, 0, 0, 0, [], 0);          
    REJ_THRE = logical(EEG.reject.rejthresh);
    REJ_JP   = EEG.reject.rejjp;
    REJ_KURT = EEG.reject.rejkurt;
    REJ_ALL  = logical(REJ_KURT + REJ_JP + REJ_THRE);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1); %marks trials for rejection

    % save stats of rejection
    col=col+1; artefact_overview(1,col) = EEG.trials;
    col=col+1; artefact_overview(1,col) = EEG.pnts;
    col=col+1; artefact_overview(1,col) = EEG.trials - sum(EEG.reject.rejglobal);
    col=col+1; artefact_overview(1,col) = sum(EEG.reject.rejglobal);
    col=col+1; artefact_overview(1,col) = 1 - sum(EEG.reject.rejglobal)/EEG.trials;
    col=col+1; artefact_overview(1,col) = sum(EEG.reject.rejglobal)/EEG.trials;

    col=col+1; artefact_overview(1,col) = numel(unique([find(REJ_THRE) find(REJ_JP) find(REJ_KURT)]));
    col=col+1; artefact_overview(1,col) = sum(REJ_THRE);
    col=col+1; artefact_overview(1,col) = sum(REJ_JP);
    col=col+1; artefact_overview(1,col) = sum(REJ_KURT);
    col=col+1; artefact_overview(1,col) = numel(setdiff(find(REJ_THRE), unique([find(REJ_JP)   find(REJ_KURT)])));
    col=col+1; artefact_overview(1,col) = numel(setdiff(find(REJ_JP),   unique([find(REJ_THRE) find(REJ_KURT)])));
    col=col+1; artefact_overview(1,col) = numel(setdiff(find(REJ_KURT), unique([find(REJ_THRE) find(REJ_JP)  ])));
    col=col+1; artefact_overview(1,col) = numel(setdiff(intersect(find(REJ_THRE),find(REJ_JP)),   find(REJ_KURT)));
    col=col+1; artefact_overview(1,col) = numel(setdiff(intersect(find(REJ_THRE),find(REJ_KURT)), find(REJ_JP)));
    col=col+1; artefact_overview(1,col) = numel(setdiff(intersect(find(REJ_JP),  find(REJ_KURT)), find(REJ_THRE)));
    col=col+1; artefact_overview(1,col) = numel(intersect(intersect(find(REJ_THRE),find(REJ_JP)),find(REJ_KURT)));

    % document data set statistics
    col=    1; artefact_epochs(1,col) = {(setdiff(find(REJ_THRE), unique([find(REJ_JP)   find(REJ_KURT)])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(find(REJ_JP),   unique([find(REJ_THRE) find(REJ_KURT)])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(find(REJ_KURT), unique([find(REJ_THRE) find(REJ_JP)  ])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_THRE),find(REJ_JP)),   find(REJ_KURT)))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_THRE),find(REJ_KURT)), find(REJ_JP)))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_JP),  find(REJ_KURT)), find(REJ_THRE)))};
    col=col+1; artefact_epochs(1,col) = {(intersect(intersect(find(REJ_THRE),find(REJ_JP)),find(REJ_KURT)))};
    col=col+1; artefact_epochs(1,col) = {unique([find(REJ_THRE) find(REJ_JP) find(REJ_KURT)])};    

    % save preprocessing stats
    if not(isfolder(PATH_OUT_STATS))
            mkdir(PATH_OUT_STATS)
    end
    csvwrite(fullfile(PATH_OUT_STATS, [allSets(n).name(1:7) FILE_MAIN '_artefact_overview.csv']), artefact_overview)
    out1 = matfile(fullfile(PATH_OUT_STATS, [allSets(n).name(1:7) FILE_MAIN '_artefact_epochs.mat']), 'writable', true);
    out1.artefact_epochs          = artefact_epochs;
    out1.artefact_epochs_varnames = artefact_epochs_varnames;

    % reject trials with EEG artifacts
    EEG = pop_rejepoch(EEG, find(EEG.reject.rejglobal), 0); 

    % ---------- Remove incorrect trials from EEG.event ------------------

    %find epochs that contain S250 or S255 or S251
    badEpochs = [EEG.event(ismember({EEG.event.type}, {'S250','S255','S251'})).epoch];

    % remove all events whose epoch is in badEpochs
    EEG.event = EEG.event(~ismember([EEG.event.epoch], badEpochs));
    EEG = pop_select(EEG, 'notrial', badEpochs);
    % ---------- Define conditions and RT outliers in EEG.event ----------

    for m = 1: length(EEG.event)
        if strcmp(EEG.event(m).type,'S 31')
           EEG.event(m+2).Condition = 'LM-R';
        elseif strcmp(EEG.event(m).type,'S 32')
            EEG.event(m+2).Condition = 'LM-S';
        elseif strcmp(EEG.event(m).type,'S 33')
            EEG.event(m+2).Condition = 'OE-R';
        elseif strcmp(EEG.event(m).type,'S 34')
            EEG.event(m+2).Condition = 'OE-S';
        end

        if strcmp(EEG.event(m).type,'S160')
           EEG.event(m-1).RT = EEG.event(m).latency - EEG.event(m-1).latency;
        elseif strcmp(EEG.event(m).type,'S150')
           EEG.event(m-1).RT = EEG.event(m).latency - EEG.event(m-1).latency;
        end

    end

    % ---------- Reject trials with no answer ---------------------
    for m=1:length(EEG.event)
       if strcmp(EEG.event(m).Condition, 'LM-R') | strcmp(EEG.event(m).Condition, 'LM-S') | strcmp(EEG.event(m).Condition, 'OE-R') | strcmp(EEG.event(m).Condition, 'OE-S') 
           if isempty(EEG.event(m).RT)
               EEG.event(m).response = 'noResponse';
           else
               EEG.event(m).response = 'yesResponse';
           end
       end
    end

    noAnswerEpochs = [EEG.event(strcmp({EEG.event.response}, 'noResponse')).epoch];

    % remove all events whose epoch is in badEpochs
    EEG.event = EEG.event(~ismember([EEG.event.epoch], noAnswerEpochs));
    EEG = pop_select(EEG, 'notrial', noAnswerEpochs);

    % ---------- Reject trials with RT < 150 ms ---------------------
    %find epochs that contain RT < 150
    toofastEpochs = [EEG.event([EEG.event.RT] < 150).epoch];

    % remove all events whose epoch is in badEpochs
    EEG.event = EEG.event(~ismember([EEG.event.epoch], toofastEpochs));
    EEG = pop_select(EEG, 'notrial', toofastEpochs);

    % ---------- Reject trials with logarithmized RTs +/- 3 SD of
    % calculate logarithmized RTs ---------------------
    for m=1:length(EEG.event)
        if strcmp(EEG.event(m).Condition, 'LM-R') | strcmp(EEG.event(m).Condition, 'LM-S') | strcmp(EEG.event(m).Condition, 'OE-R') | strcmp(EEG.event(m).Condition, 'OE-S') 
            EEG.event(m).RTlog = log(EEG.event(m).RT);
        end
    end

    % compute condition-wise means and SDs
    LMR = [];
    LMS = [];
    OER = [];
    OES = [];
    for m=1:length(EEG.event)
        if strcmp(EEG.event(m).Condition, 'LM-R') 
            LMR = [LMR, EEG.event(m).RTlog];
        elseif strcmp(EEG.event(m).Condition, 'LM-S') 
            LMS = [LMS, EEG.event(m).RTlog];
        elseif strcmp(EEG.event(m).Condition, 'OE-R') 
            OER = [OER, EEG.event(m).RTlog];
        elseif strcmp(EEG.event(m).Condition, 'OE-S') 
            OES = [OES, EEG.event(m).RTlog];
        end
    end
    mean_LMR = mean(LMR);
    std_LMR = std(LMR);
    mean_LMS = mean(LMS);
    std_LMS = std(LMS);
    mean_OER = mean(OER);
    std_OER = std(OER);
    mean_OES = mean(OES);
    std_OES = std(OES);

    %find epochs that contain RTlog exceeding +/- 3 SDs
    for m=1:length(EEG.event)
        if strcmp(EEG.event(m).Condition, 'LM-R') 
            if EEG.event(m).RTlog < mean_LMR-3*std_LMR || EEG.event(m).RTlog > mean_LMR+3*std_LMR
                EEG.event(m).RToutlier = 'outlier';
            else
                EEG.event(m).RToutlier = 'okay';
            end
        elseif strcmp(EEG.event(m).Condition, 'LM-S') 
            if EEG.event(m).RTlog < mean_LMS-3*std_LMS || EEG.event(m).RTlog > mean_LMS+3*std_LMS
                EEG.event(m).RToutlier = 'outlier';
            else
                EEG.event(m).RToutlier = 'okay';
            end
        elseif strcmp(EEG.event(m).Condition, 'OE-R') 
            if EEG.event(m).RTlog < mean_OER-3*std_OER || EEG.event(m).RTlog > mean_OER+3*std_OER
                EEG.event(m).RToutlier = 'outlier';
            else
                EEG.event(m).RToutlier ='okay';
            end
        elseif strcmp(EEG.event(m).Condition, 'OE-S') 
            if EEG.event(m).RTlog < mean_OES-3*std_OES || EEG.event(m).RTlog > mean_OES+3*std_OES
                EEG.event(m).RToutlier = 'outlier';
            else
                EEG.event(m).RToutlier = 'okay';
            end
        end
    end

      %find epochs that contain S250 or S260
    outlierEpochs = [EEG.event(strcmp({EEG.event.RToutlier}, 'outlier')).epoch];

    % remove all events whose epoch is in badEpochs
    EEG.event = EEG.event(~ismember([EEG.event.epoch], outlierEpochs));
    EEG = pop_select(EEG, 'notrial', outlierEpochs);


    %% TF analysis

    for cond = {'LM-R','LM-S','OE-R','OE-S'}%1:4

        DAT = EEG;

        % Get all channel combinations
        channel_combs = nchoosek(1:DAT.nbchan, 2);
        channel_name_combs = cell(size(channel_combs));

        for i = 1:size(channel_combs, 1)
            channel_name_combs{i, 1} = DAT.chanlocs(channel_combs(i, 1)).labels;
            channel_name_combs{i, 2} = DAT.chanlocs(channel_combs(i, 2)).labels;
        end


        % Set convolution parameters
        wavelet_time = -2 : 1 / DAT.srate : 2;
        frequencies = 6;
        minn = EEG.trials;
        n_data = DAT.pnts *minn;
        n_kernel = length(wavelet_time);
        n_convolution = n_data + n_kernel - 1;
        half_wavelet = (length(wavelet_time) - 1) / 2;
        n_cycles = 4;

        % Save frequencies and time
        if not(isfolder(PATH_SUB_CONNY_MATS))
            mkdir(PATH_SUB_CONNY_MATS)
        end
        dlmwrite([PATH_SUB_CONNY_MATS char(cond) '_tf_time.csv'], DAT.times);
        dlmwrite([PATH_SUB_CONNY_MATS char(cond) '_tf_frqs.csv'], frequencies);

        % Create wavelets
        cmwX = zeros(length(frequencies), n_convolution);
        for f = 1 : length(frequencies)

            % Create time-domain wavelet
            s2 = n_cycles(f) / (2 * pi * frequencies(f));
            cmw = exp(2 * 1i * pi * frequencies(f) .* wavelet_time) .* exp((-wavelet_time .^ 2) / (2 * s2 .^2));

            % Compute fourier coefficients of wavelet and normalize
            cmwX(f, :) = fft(cmw, n_convolution);
            cmwX(f, :) = cmwX(f, :) ./ max(cmwX(f, :));
        end


        % Select condition
        idx = [EEG.event(strcmp({EEG.event.Condition}, cond)).epoch];

        % Iterate chans
        for ch = 1:size(channel_combs, 1)
            try
            fprintf('\nsub (%i/%i) chans (%i/%i)', n, length(allSets), ch, size(channel_combs, 1));

            chanIdx1 = find(strcmpi({DAT.chanlocs.labels}, channel_name_combs{ch,1}));
            chanIdx2 = find(strcmpi({DAT.chanlocs.labels}, channel_name_combs{ch,2}));

            if isempty(chanIdx1) || isempty(chanIdx2)
                warning('Skipping ch %d: channel not found (%s, %s)', ch, channel_name_combs{ch,1}, channel_name_combs{ch,2});
                continue;
            end

            chandata_X1 = fft(reshape(double(squeeze(DAT.data(chanIdx1, :, idx))), 1, []), n_convolution);
            chandata_X2 = fft(reshape(double(squeeze(DAT.data(chanIdx2, :, idx))), 1, []), n_convolution);

            mat_pli = zeros(length(frequencies), DAT.pnts);

            for f = 1:length(frequencies)
                as1 = ifft(cmwX(f, :) .* chandata_X1);
                as1 = as1(half_wavelet + 1:end - half_wavelet);
                as1 = reshape(as1, DAT.pnts, minn);

                as2 = ifft(cmwX(f, :) .* chandata_X2);
                as2 = as2(half_wavelet + 1:end - half_wavelet);
                as2 = reshape(as2, DAT.pnts, minn);

                phasediffs = exp(1i * (angle(as1) - angle(as2)));
                mat_pli(f, :) = abs(mean(sign(imag(phasediffs)), 2));
            end

            dlmwrite([PATH_SUB_CONNY_MATS allSets(n).name(1:7) '_' char(cond) '_pli_' ...
                      char(channel_name_combs{ch,1}) '_' char(channel_name_combs{ch,2}) '.csv'], mat_pli);

        catch ME
            warning('Error at ch %d (%s, %s): %s', ch, channel_name_combs{ch,1}, channel_name_combs{ch,2}, ME.message);
        end
    end
    end % End condition loop



end

%% ======================= CONNECTIVITY DEGREE =================================================================================================

% define base measure. 
base_measure = 'pli';

% Get frequencies and time
tf_time = DAT.times;
tf_frqs = frequencies;

% Define frequency of interest
connyfreq = 6;
[~, frq_idx] = min(abs(tf_frqs - connyfreq));
actual_frq = tf_frqs(frq_idx);

% Result matrices
connytrace = [];

% Iterating conditions
for cond = {'LM-R','LM-S','OE-R','OE-S'}
    % Iterating subject list
    for s = 1: length(allSets)

        % Init subject connycubes
        connymat = zeros(EEG.nbchan, EEG.nbchan, length(tf_time));

        % Iterate chancombs
        for ch = 1 : size(channel_name_combs, 1)

            % Talk
            fprintf('\nreading sub (%i/%i) chans (%i/%i)', s, length(allSets), ch, size(channel_name_combs, 1));

            % Load conny matrices
            stim = dlmread([PATH_SUB_CONNY_MATS allSets(s).name(1:7) '_' char(cond) '_' base_measure '_' char(channel_name_combs(ch, 1)) '_' char(channel_name_combs(ch, 2)) '.csv']);

            % Fill connycubes
            connymat(channel_combs(ch, 1), channel_combs(ch, 2), :) = stim(frq_idx, :);
            connymat(channel_combs(ch, 2), channel_combs(ch, 1), :) = stim(frq_idx, :);


        end % End chancomb iteration

        % Init thresholded subject connycubes
        connythresh = zeros(EEG.nbchan, EEG.nbchan, length(tf_time));

        % Thresholding connymats at each time point
        for t = 1 : size(connymat, 3)

            % Determine thresh at t
            tmp = nonzeros(tril(squeeze(connymat(:, :, t))));
            thresh = median(tmp) + std(tmp);

            % Binarize mats
            connyslice = zeros(EEG.nbchan, EEG.nbchan);
            connyslice(squeeze(connymat(:, :, t)) > thresh) = 1;
            connythresh(:, :, t) = connyslice;

        end

        % Iterate chans and sum connections
        for ch = 1 : EEG.nbchan

            % Get chan-slice and average to connytrace
            connytrace(s, ch, :) = sum(squeeze(connythresh(ch, :, :)), 1);

        end % End chanloop

    end % End subject loop 

   % Save sub specific connydegree traces
    if not(isfolder(PATH_CONNYDEGREE_TRACES))
         mkdir(PATH_CONNYDEGREE_TRACES)
    end

    for ch = 1 : EEG.nbchan
        dlmwrite([PATH_CONNYDEGREE_TRACES 'connydegree_traces_chan_' EEG.chanlocs(ch).labels '_' base_measure '_' char(cond) '.csv'], squeeze(connytrace(:, ch, :)));
    end % End chanloop

end


%% ======================= PART4: ERPs =================================================================================================

for s = 1:length(allSets)

    % Load data
    subject = allSets(s).name(1:7);

    EEG = pop_loadset('filename',allSets(s).name,'filepath',PATH_IN);

    EEG = pop_eegfiltnew(EEG, 'locutoff', 0.1, 'hicutoff', 16);


    % Recode event types for ERPlab
    for j = 1:length(EEG.event)
        if strcmp(EEG.event(j).type, 'S 51')
            EEG.event(j).type = 1;
        elseif strcmp(EEG.event(j).type, 'S 52') 
            EEG.event(j).type = 2;
        elseif strcmp(EEG.event(j).type, 'S 53') 
            EEG.event(j).type = 3;
        elseif strcmp(EEG.event(j).type, 'S 54') 
            EEG.event(j).type = 4;
        end

    end



        % Remove response markers
        EEG = eeg_checkset(EEG, 'eventconsistency');

        % Create EventList for ERPLab
        EEG  = pop_creabasiceventlist(EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' } );

        % Create Bins for different conditions
        EEG  = pop_binlister(EEG , 'BDF', [PATH_MAIN '\binlister_switching.txt'], 'IndexEL',  1, 'SendEL2', 'Workspace&EEG', 'UpdateEEG', 'on', 'Voutput', 'EEG' );
        EEG.setname = ['switching_' subject '_bin'];

        % Epoch data
        EEG = pop_epochbin( EEG , [-200.0  1000.0],  [ -200 0]);
        EEG.setname = ['switching_' subject '_bin_epoch'];

        % Compute ERPs for the different bins
        ERP = pop_averager( EEG , 'Criterion', 'all', 'ExcludeBoundary', 'off', 'SEM', 'on' );
        ERP.erpname = ['switching_' subject '_bin_average'];

        % Save ERP set to
        if not(isfolder(PATH_ERP))
         mkdir(PATH_ERP)
        end

        ERP = pop_savemyerp(ERP, 'erpname', ERP.erpname, 'filename', ['switching_target_locked_' subject '-2.erp'], 'filepath', PATH_ERP);

end