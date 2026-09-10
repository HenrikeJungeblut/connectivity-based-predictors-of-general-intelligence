%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      calculate whole-brain SC-FC coupling                     %
% Author:        Henrike Jungeblut, University of Mainz                   %
% Last Modified: 09/2026                                                  %
% Note:                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------- Necessary Functions -----------------------------------------
% %
% % Brain Connectivity Toolbox
% % The Brain Connectivity Toolbox is provided with this
% % code within the folder "Tools".
% %

%% load brain connectivity toolbox and set directories
filepath   = cd;
parts = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});
PATH_IN_dti = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'dwi');
addpath(fullfile(PATH_MAIN, 'Code_package', 'Tools', '2019_03_03_BCT'));

%% set whether res/rawfunc data should be used
% the script has to be run twice: Once for res (background FC) and once for
% rawfunc (task-evoked FC). Only the variable data_type needs to be
% adjusted.
% for data_type = ["rest", "rawfunc", "rest"] %res: background connectivity/rawfunc: task-evoked
    if strcmp(data_type, 'res')
        PATH_IN = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'func', 'timeseries_per_region', 'res_');
        PATH_OUT = fullfile(PATH_MAIN, 'Code_package', 'Analysis_data', 'coupling_comm_measures.csv');
    elseif strcmp(data_type, 'rawfunc')
        PATH_IN = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'func', 'rawfunc_timeseries_per_region', 'rawfunc_');
        PATH_OUT = fullfile(PATH_MAIN, 'Code_package', 'Analysis_data', 'rawfunc_coupling_comm_measures.csv');
    elseif strcmp(data_type, 'rest')
        PATH_IN = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'func', 'rest_timeseries_per_region', 'rest_');
        PATH_OUT = fullfile(PATH_MAIN, 'Code_package', 'Analysis_data', 'rest_coupling_comm_measures.csv');
    end

    %% loading subject list:
    fileID = fopen(fullfile(PATH_MAIN, 'Code_package', 'Data', 'subject_list.txt'), 'r');

    subjectNames = textscan(fileID, '%s');
    fclose(fileID);
    subjectNames = subjectNames{1};

    %remove outliers who perform only around chance level
    to_remove = {'sub-032', 'sub-062', 'sub-076', 'sub-095', 'sub-096', 'sub-129', ...
        'sub-147', 'sub-151', 'sub-160', 'sub-165', 'sub-172', 'sub-174', 'sub-180', ...
        'sub-184', 'sub-188', 'sub-190', 'sub-194', 'sub-195', 'sub-205', 'sub-211', ...
        'sub-212', 'sub-213', 'sub-214', 'sub-215'};

    % Filter out the outliers
    subjectNames_filtered = subjectNames(~ismember(subjectNames, to_remove));

    %% initialize vectors for SC-FC coupling for all measures
    Ediff_all = [];
    SI_all = [];
    T_all = [];
    Eres_all = [];
    F_all = [];
    corrs_all = [];
    sub_name = [];
    size_func = [];
    size_dwi = [];

    %% whole-brain parcels
    % for whole-brain analysis,we use all 360 parcels
    wb_areas_orig = 1:360;

    %% loop through subject list to calculate SC-FC coupling

    for sub = 1:length(subjectNames_filtered)
        subject = subjectNames_filtered{sub};

        %% Subject-specific parcel removals
        % For all subjects, we remove parcels 120 and 300 because this
        % bilateral parcel disappeared for most participants as a consequence
        % of resampling the MRI data.

        % Four additional subjects need additional treatment because they are
        % mising one or two additional parcels.

        if strcmp(subject, 'sub-019') % functional data coverage
            removed = [118 120 300];
        elseif strcmp(subject, 'sub-021') % functional data coverage
            removed = [120 172 300 352];
        elseif strcmp(subject, 'sub-025') || strcmp(subject, 'sub-083') % zero nodes in SC matrix
            removed = [120 300 332];
        elseif strcmp(subject, 'sub-197') % functional data coverage
            removed = [120 183 300];
        else
            removed = [120 300]; % resampling consequence
        end

        % get parcels in wb_areas_orig that are not in removed ->
        % exclude all parcels that shall be removed
        wb_areas = setdiff(wb_areas_orig, removed);

        %% load func data
        func_lh = fullfile([PATH_IN, subject,'_lh_filtered_timeseries_per_glasser_region.txt']);
        func_rh = fullfile([PATH_IN, subject,'_rh_filtered_timeseries_per_glasser_region.txt']);

        % relevant if data_type == 'rest' because not all subjects from the
        % main analaysis also have rest data available
        if ~isfile(func_lh) || ~isfile(func_rh)
            fprintf('Skipping %s: missing file(s)\n', subject);
            continue
        end

        sub_name = [sub_name; subject];

        % In most cases, lh_mat and rh_mat will have 180 or 179 columns. There
        % are 180 columns if no parcel got lost during preprocessing.
        % lh_mat and rh_mat have 540 rows that corresponds to the number of
        % timepoints in  the fMRI timeseries.
        lh_mat = readmatrix(func_lh, 'Delimiter', ' ', 'DecimalSeparator', '.');
        rh_mat = readmatrix(func_rh, 'Delimiter', ' ', 'DecimalSeparator', '.');

        %% Adjust func data
        % We fill in missing columns (missing parcels with zero columns to
        % ensure that all matrices start out with 360 x 360 as dimensions

        % Subject-specific handling of subjects who are missing an additional
        % parcel
        if strcmp(subject, 'sub-019') %parcel 118
            lh_mat = [lh_mat(:, 1:117), zeros(size(lh_mat, 1), 1), lh_mat(:, 118:end)];
        elseif strcmp(subject, 'sub-021') %parcels 172 in both hemispheres
            lh_mat = [lh_mat(:, 1:171), zeros(size(lh_mat, 1), 1), lh_mat(:, 172:end)];
            rh_mat = [rh_mat(:, 1:171), zeros(size(rh_mat, 1), 1), rh_mat(:, 172:end)];
        elseif strcmp(subject, 'sub-197') %parcel 183
            rh_mat = [rh_mat(:, 1:2), zeros(size(rh_mat, 1), 1), rh_mat(:, 3:end)];
        end

        % in case, parcel 120 is missing as resampling consequence, fill matrix
        % with zero column to account for size mismatch
        if size(lh_mat,2) == 179
            lh_mat = [lh_mat(:, 1:119), zeros(size(lh_mat, 1), 1), lh_mat(:, 120:end)];
        end
        if size(rh_mat,2) == 179
            rh_mat = [rh_mat(:, 1:119), zeros(size(rh_mat, 1), 1), rh_mat(:, 120:end)];
        end

        % sanity check to see whether there are 180 clumns
        fprintf('%s: Dimension of lh: %i', subject, size(lh_mat,2));
        fprintf('%s: Dimension of rh: %i', subject, size(rh_mat,2));

        %% get FC matrix and apply Fisher to Z transformation
        func_mat = [lh_mat, rh_mat];

        func_mat = func_mat(:,wb_areas);

        % Fisher-Z correlation matrix
        Z_func_corr_mat = atanh(corrcoef(func_mat));
        Z_func_corr_mat(isinf(Z_func_corr_mat)) = NaN;
        %imagesc(Z_func_corr_mat);
        size_func = [size_func; length(Z_func_corr_mat)];

        %% load dti data
        % the SC matrix W always has 360 x 360 as dimensions.
        dti = fullfile(PATH_IN_dti, [subject,'_mean_FA_connectome.csv']);
        W = readmatrix(dti, 'Delimiter', ',', 'DecimalSeparator', '.');

        W = W(wb_areas, wb_areas);
        size_dwi = [size_dwi; length(W)];

        % sanity check: there should be no SC matrices with a degree of zero
        % because subjects who had zero nodes were already handled above
        deg = sum(W>0, 2);        % node degree (number of nonzero connections per node)
        zero_nodes = find(deg==0);
        if ~isempty(zero_nodes)
            fprintf('Subject %s: %d zero nodes. Nodes: %s\n', subject, numel(zero_nodes), mat2str(zero_nodes));
        else
            fprintf('Subject %s: no zero nodes\n', subject);
        end

        %% Calculate SC-FC
        % general pipeline
        % 1) Communication matrix (diffusion efficiency, search information,
        % path transitivity, resource efficiency, communicability) is
        % calculated on the SC matrix.
        % 2) Make communication matrix symmetric if needed and invert matrix if
        % needed (the latter is only done for search information).
        % 3) Correlate FC matrix and communication matrix row for row
        % 4) average across all row-wise SC-FC scores to get a whole-brain
        % SC-FC value per participant

        %% diffusion efficiency (Ediff)

        % --- Quality checks for structural matrix W ---
        if rank(W) == 0 || any(sum(W)==0) || any(sum(W,2)==0)
            disp(['Skipping diffusion efficiency for subject ', subject, ...
                ' (invalid or disconnected W matrix)']);
            Ediff_all = [Ediff_all; NaN]; % keep table lengths aligned
        else
            % Compute diffusion efficiency normally
            [GEdiff,Ediff] = diffusion_efficiency(W);
            EdiffSymm = (Ediff + transpose(Ediff))/2;

            % initialize empty correlation vector
            Ediff_corrs = [];
            for n=1:length(Z_func_corr_mat)
                correlation = corr(Z_func_corr_mat(:,n), EdiffSymm(:,n), 'Rows', 'complete');
                Ediff_corrs = [Ediff_corrs, correlation];
            end
            Ediff_all = [Ediff_all; mean(Ediff_corrs)];
        end

        %% search information (SI)
        L = weight_conversion(W, 'lengths');
        SI = search_information(W, L);
        SISymm = -1*((SI + transpose(SI))/2); % make symmetric and switch sign to reflect ease of communication

        %imagesc(SISymm);

        % initialize empty correlation vector
        SI_corrs = [];
        for n=1:length(Z_func_corr_mat)
            correlation = corr(Z_func_corr_mat(:,n), SISymm(:,n), 'Rows', 'complete');
            SI_corrs = [SI_corrs, correlation];
        end
        SI_all = [SI_all; mean(SI_corrs)];

        %% path transitivity (T)
        T = path_transitivity(L); % no transform needs to be supplied if input is a length matrix
        %imagesc(T);

        % initialize empty correlation vector
        T_corrs = [];
        for n=1:length(Z_func_corr_mat)
            correlation = corr(Z_func_corr_mat(:,n), T(:,n), 'Rows', 'complete');
            T_corrs = [T_corrs, correlation];
        end
        T_all = [T_all; mean(T_corrs)];

        %% resource efficiency (Eres)
        % thresholding to get unweighted matrix (UW)
        p = 0.5;
        W_thr = threshold_proportional(W, p);
        UW = W_thr > 0;

        % calculating resource efficiency
        lambda = .5; % in Goni et al., 2013 lambda = 0.5
        [Eres,prob_SPL] = resource_efficiency_bin(UW,lambda);

        EresSymm = (Eres + transpose(Eres))/2; % make symmetric

        %imagesc(EresSymm);

        % initialize empty correlation vector
        Eres_corrs = [];
        for n=1:length(Z_func_corr_mat)
            correlation = corr(Z_func_corr_mat(:,n), EresSymm(:,n), 'Rows', 'complete');
            Eres_corrs = [Eres_corrs, correlation];
        end
        Eres_all = [Eres_all; mean(Eres_corrs)];

        %% communicability (F)
        % Code von https://github.com/brain-networks/local_scfc/blob/main/fcn/communicability_wei.m
        N = size(W,1);
        B = sum(W')';
        C = diag(B);
        D = C^(-(1/2));
        E = D * W * D;
        F = expm(E);
        F = F.*~eye(N);
        %imagesc(F);

        % initialize empty correlation vector
        F_corrs = [];
        for n=1:length(Z_func_corr_mat)
            correlation = corr(Z_func_corr_mat(:,n), F(:,n), 'Rows', 'complete');
            F_corrs = [F_corrs, correlation];
        end
        F_all = [F_all; mean(F_corrs)];

        %% simple correlation
        % initialize empty correlation vector
        corrs = [];
        for n=1:length(Z_func_corr_mat)
            correlation = corr(Z_func_corr_mat(:,n), W(:,n), 'Rows', 'complete');
            corrs = [corrs, correlation];
        end
        corrs_all = [corrs_all; mean(corrs)];


    end

    %% combine data to table and save
    dataTable = table(sub_name, Ediff_all, SI_all, T_all, Eres_all, F_all);
    writetable(dataTable, PATH_OUT)
end