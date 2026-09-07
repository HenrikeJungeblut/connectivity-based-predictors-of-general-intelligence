%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      general EEG preprocessing                                %
% Author:        Mareike Hülsemann, Henrike Jungeblut                     %
% Last Modified: 09/2026                                                  %
% Note:          this script is not reproducible as raw EEG data is not   %
%                provided in the code package                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ---------- Necessary Functions -----------------------------------------
% %
% % EEGLAB version eeglab2025.0.0 with the following plugins is used:
% % - "Cleanline" v2.1
% % - "ICLabel" v1.6
% % - "clean_rawdata" v2.10
% % - "firfilt" v2.8
% % The EEGLAB toolbox with all plugins installed is provided with this
% % code within the folder "Tools".
% %


%% ---------- Specifiying Paths -------------------------------------------
% %
% % If the provided folder structure is used, no paths need to be adjusted.
% % This script is saved in the subfolder "code" and automatically
% % references to the subfolders "data", where raw data is provided, and
% % "derivatives", where single-subject output files are exported to.
% % The EEGLAB toolbox with all plugins installed is provided with this
% % code within the folder "software".
% %


%% ---------- Set Up ------------------------------------------------------

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});

% load list of raw EEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not reproducible as raw EEG data is not provided in the code package    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allFiles = dir(fullfile(PATH_MAIN, 'segmented_task_data', 'sub*', 'sub*'));

% define experiment and task specific settings: preprocessing
DATA_SET = '_task-switch_eeg';
elecs2remove = {'G32'};
online_ref    = struct('labels',{'Cz'}, ... % G31 ist online ref
    'sph_radius', {1}, 'sph_theta', {0}, 'sph_phi', {90}, ...
    'theta', {'0'}, 'radius', {'0'}, 'X',{'0'}, 'Y',{'0'}, 'Z',{'1'}, ...
    'type',{'EEG'},'ref',{''},'urchan',{[]});
line_noise = [50 100];


%% ---------- start EEGLAB ------------------------------------------------
addpath(fullfile(PATH_MAIN, 'Code_package', 'Tools', 'eeglab2025.0.0'));
eeglab; close(gcf);
pop_editoptions('option_single', 0, 'option_computeica', 1);


%% ---------- Loop over Subjects ------------------------------------------

for n= 1:numel(allFiles)

    % ---------- set up ---------------------------------------------------

    % process only files that have not yet been processed
    if isfile(fullfile(PATH_MAIN, 'preprocessed', ...
            [allFiles(n).name(1:7) DATA_SET '_preprocessed.set']))
        disp([allFiles(n).name ' is already preprocessed.'])
    else
        disp(['*** Started Set ' num2str(n) ' of ' num2str(length(allFiles)) '; ' allFiles(n).name ' ***'])

        % initialise for reporting preprocessing stats and results
        stats = table;
       
        
        % ---------- preparing data segment -------------------------------

        % load data
        EEG = pop_loadset('filename', [allFiles(n).name(1:7) DATA_SET '.set'], ...
            'filepath', fullfile(PATH_MAIN, 'segmented_task_data', allFiles(n).name(1:7))); 
        % convert EEG.data to double precision
        EEG.data = double(EEG.data);

        % remove unused channels
        EEG = pop_select(EEG, 'nochannel', elecs2remove);

        % save urchanlocs structure without online reference
        orig_chanlocs = EEG.chanlocs;

        % save subject name
        EEG.subject = allFiles(n).name(1:7);

        % document data set statistics
        stats.id            = EEG.subject; % ist leer, sollte das gefüllt sein?
        stats.srate_raw     = EEG.srate;
        stats.pnts_raw      = EEG.pnts;
        stats.nb_chan_raw   = EEG.nbchan;
        stats.nb_events_raw = numel(EEG.event);

        % ---------- cleaning data segment --------------------------------

        % filter data temporarily for cleanline and noisy channel detection
        EEGfilt = pop_eegfiltnew(EEG, 'locutoff', 1);

        % detect line noise
        EEGclean = pop_cleanline(EEGfilt, ...
            'bandwidth',2, ...              % 'Bandwidth',2,
            'chanlist',1:EEGfilt.nbchan, ...% 'ChanCompIndices',[1:EEG.nbchan],
            'computepower',0, ...           % 'ComputeSpectralPower',true,
            'linefreqs',line_noise, ...     % 'LineFrequencies',[60 120],
            'newversion',1, ...
            'normSpectrum',0, ...           % 'NormalizeSpectrum',false,
            'p',0.01, ...                   % 'LineAlpha',0.01,
            'pad',2, ...                    % 'PaddingFactor',2,
            'plotfigures',0, ...            % 'PlotFigures',false,
            'scanforlines',1, ...           % 'ScanForLines',true,
            'sigtype','Channels', ...       % 'SignalType','Channels',
            'taperbandwidth',2, ...         % setting for newversion
            'tau',100, ...                  % 'SmoothingFactor',100,
            'verb',0, ...                   % 'VerboseOutput',1,
            'winsize',4, ...                % 'SlidingWinLength',EEGfilt.pnts/EEGfilt.srate,
            'winstep',1);                   % 'SlidingWinStep',EEG.pnts/EEG.srate

        % remove line noise without comitting to filter strategy
        EEG.data = EEG.data - EEGfilt.data + EEGclean.data;

        % detect noisy channels
        EEGclean = pop_clean_rawdata(EEGclean, ...
            'FlatlineCriterion', 'off', ... default: 5 ('off')
            'ChannelCriterion', 0.8, ... default: 0.8 ('off')
            'LineNoiseCriterion', 'off', ... default: 4 ('off')
            'Highpass', 'off', ...
            'BurstCriterion', 'off', ...
            'WindowCriterion', 'off', ...
            'BurstRejection', 'off', ...
            'Distance','Euclidian');
        % save noisy channels
        removed_chans = ~ismember({orig_chanlocs.labels}, {EEGclean.chanlocs.labels});
        removed_chans_table = array2table(removed_chans, ...
            'VariableNames',{orig_chanlocs.labels});
        % remove noisy channels from relevant data set
        EEG = pop_select( EEG, 'rmchannel',{orig_chanlocs(removed_chans).labels});

        % document data set statistics
        stats.nb_chan_removed = sum(removed_chans);

        
        % save noisy channels
        removed_corr_chans = ~ismember({orig_chanlocs.labels}, {EEGclean.chanlocs.labels});
        removed_corr_chans_table = array2table(removed_corr_chans, ...
            'VariableNames',{orig_chanlocs.labels});

        % insert zero-channel for online-reference (according to fullRankAveRef() from Makoto Miyakoshi)
        EEG.nbchan = EEG.nbchan+1;
        EEG.data(end+1,:) = zeros(1, EEG.pnts);
        EEG.chanlocs(1,EEG.nbchan).labels = 'onlineReference';
        % re-reference to average reference
        EEG = pop_reref(EEG, []);
        % discard one channel to adjust rank
        EEG = pop_select(EEG,'nochannel',{'onlineReference'});

        % tidy up workspace
        clear EEGfilt EEGclean


        % ---------- removing artifacts with ICA --------------------------

        % create data set for ICA
        ICA = EEG;

        % filter data for ICA (1Hz high-pass)
        ICA = pop_eegfiltnew(ICA, 'locutoff', 1);

        % adjust for deficient rank due to re-ref
        dataRank = sum(eig(cov(double(ICA.data'))) > 1E-7);

        % document data set statistics
        stats.dataRank = dataRank;

        % epoch data for automatic artifact rejection
        ICA = eeg_regepochs(ICA, 'recurrence', 2, 'limits', [0 2], 'rmbase', NaN);

        % document data set statistics
        stats.pnts_epochs = ICA.pnts;
        stats.nb_epochs   = ICA.trials;

        % artifact rejection
        ICA = pop_eegthresh(ICA,1,1:ICA.nbchan,-1000,1000,0,ICA.xmax,0,0);
        ICA = pop_jointprob(ICA,1,1:ICA.nbchan, 5, 2, 0, 0, 0, [], 0); % joint probability, according to Makotos recommendations and Delorme et al. 2007 (5SD; more strict for ICA)
        ICA = pop_rejkurt(  ICA,1,1:ICA.nbchan, 5, 2, 0, 0, 0, [], 0);
        REJ_THRE = logical(ICA.reject.rejthresh);
        REJ_JP   = ICA.reject.rejjp;
        REJ_KURT = ICA.reject.rejkurt;
        REJ_ALL  = logical(REJ_KURT + REJ_JP + REJ_THRE);
        ICA = eeg_rejsuperpose(ICA, 1, 1, 1, 1, 1, 1, 1, 1);
        ICA = pop_rejepoch(ICA, find(ICA.reject.rejglobal), 0);

        % document data set statistics
        stats.nb_epochs_ICA_clean         = ICA.trials;
        stats.nb_epochs_ICA_rejected      = sum(REJ_ALL);
        stats.nb_epochs_ICA_rejected_thre = sum(REJ_THRE);
        stats.nb_epochs_ICA_rejected_prob = sum(REJ_JP);
        stats.nb_epochs_ICA_rejected_kurt = sum(REJ_KURT);
        stats.nb_epochs_ICA_rejected_thre_unique = sum(REJ_THRE & ~ (REJ_JP   | REJ_KURT));
        stats.nb_epochs_ICA_rejected_prob_unique = sum(REJ_JP   & ~ (REJ_THRE | REJ_KURT));
        stats.nb_epochs_ICA_rejected_kurt_unique = sum(REJ_KURT & ~ (REJ_THRE | REJ_JP  ));

        % run ICA
        ICA = pop_runica(ICA, ...
            'icatype', 'runica', 'extended', 1, ...
            'lrate', 1e-5, 'maxsteps', 2000, ...
            'verbose', 'on');

        % run ICLabel
        ICA = pop_iclabel(ICA, 'default');


        % ---------- flag ICs manually for rejection ----------------------

        % find max probabilities and their indices
        [max_prob, max_idx] = ...
            max(ICA.etc.ic_classification.ICLabel.classifications, [], 2);

        % determine which indices meet the rejection criteria
        IC_muscle = (max_idx == 2 & max_prob >= 0.9);
        IC_eye    = (max_idx == 3 & max_prob >= 0.25);
        IC_heart  = (max_idx == 4 & max_prob >= 0.6);

        % set rejection flags
        ICA.reject.gcompreject = IC_muscle | IC_eye | IC_heart;

        % document data set statistics
        stats.nb_rejIC_eye    = sum(IC_eye);
        stats.nb_rejIC_muscle = sum(IC_muscle);
        stats.nb_rejIC_heart  = sum(IC_heart);
        stats.nb_rejIC_total  = sum(ICA.reject.gcompreject);


        % ---------- removing bad ICs -------------------------------------

        % copy IC weights to EEG data set
        EEG.icachansind              = ICA.icachansind;
        EEG.icaweights               = ICA.icaweights;
        EEG.icasphere                = ICA.icasphere;
        EEG.icawinv                  = ICA.icawinv;
        EEG.etc.ic_classification    = ICA.etc.ic_classification;
        EEG.etc.icaweights_beforerms = ICA.etc.icaweights_beforerms;
        EEG.etc.icasphere_beforerms  = ICA.etc.icasphere_beforerms;
        EEG.reject.gcompreject       = ICA.reject.gcompreject;

        % remove components previously marked for rejection
        EEG = pop_subcomp( EEG, [], 0); % remove IC saved in EEG.reject.gcompreject

        % interpolate removed channels and online reference
        orig_chanlocs(size(orig_chanlocs,2)+1) = online_ref;
        EEG = pop_interp(EEG, orig_chanlocs, 'spherical');


        % ---------- saving files -----------------------------------------
        % save continuous preprocessed EEG dataset
        EEG.setname = [EEG.subject DATA_SET '_preprocessed'];
        if not(isfolder(fullfile(PATH_MAIN, 'preprocessed')))
            mkdir(fullfile(PATH_MAIN, 'preprocessed'))
        end
        EEG = pop_saveset(EEG, 'filename', EEG.setname, ...
            'filepath', fullfile(PATH_MAIN, 'preprocessed'));

        % save preprocessing stats
        writetable(stats, ...
            fullfile(PATH_MAIN, 'preprocessed', ...
            [EEG.setname '_stats.csv']))

        % save removed channels
        writetable(removed_chans_table, ...
            fullfile(PATH_MAIN, 'preprocessed', ...
            [EEG.setname '_removed_chans.csv']))


    end
end

