%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      calculate reliability of whole-brain SC-FC coupling      %
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

%% set paths
filepath   = cd;
parts = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});
PATH_IN_dti = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'dwi');
addpath(fullfile(PATH_MAIN, 'Code_package', 'Tools', '2019_03_03_BCT'));

%% set whether res/rawfunc/rest data should be used
data_type = 'res'; %res: background connectivity/rawfunc: task-evoked
if strcmp(data_type, 'res')
    PATH_IN = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'func', 'timeseries_per_region', 'res_'); 
elseif strcmp(data_type, 'rawfunc')
    PATH_IN = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'func', 'rawfunc_timeseries_per_region', 'rawfunc_'); 
end

fileID = fopen(fullfile(PATH_MAIN, 'Code_package', 'Data', 'subject_list.txt'), 'r');

subjectNames = textscan(fileID, '%s');
fclose(fileID);
subjectNames = subjectNames{1};

%% remove outliers who perform only around chance level
to_remove = {'sub-032', 'sub-062', 'sub-076', 'sub-095', 'sub-096', 'sub-129', ...
    'sub-147', 'sub-151', 'sub-160', 'sub-165', 'sub-172', 'sub-174', 'sub-180', ...
    'sub-184', 'sub-188', 'sub-190', 'sub-194', 'sub-195', 'sub-205', 'sub-211', ...
    'sub-212', 'sub-213', 'sub-214', 'sub-215'};

% Filter out the outliers
subjectNames_filtered = subjectNames(~ismember(subjectNames, to_remove));

%initialize vectors for SC-FC coupling for all measures
Ediff_all_L = [];
SI_all_L = [];
T_all_L = [];
Eres_all_L = [];
F_all_L = [];
Ediff_all_R = [];
SI_all_R = [];
T_all_R = [];
Eres_all_R = [];
F_all_R = [];
sub_name = [];
size_func = [];
size_dwi = [];

left_areas = 1:180;
right_areas = 181:360;

%% loop through subject list
for sub = 1:length(subjectNames_filtered)
    subject = subjectNames_filtered{sub};

    % Subject-specific removals
    if strcmp(subject, 'sub-019')
        removed_l = [118 120];
        removed_r = [300];
    elseif strcmp(subject, 'sub-021')
        removed_l = [120 172];
        removed_r = [300 352];
    elseif strcmp(subject, 'sub-025') || strcmp(subject, 'sub-083') % zero nodes
        removed_l = [120 300];
        removed_r = [332];
    elseif strcmp(subject, 'sub-197')
        removed_l = [120];
        removed_r = [183 300];
    else
        removed_l = [120]; %resampling consequence
        removed_r = [300]; %resampling consequence
    end
    
    % get parcels in wb_areas_orig that are not in removed ->
    % exclude all parcels that shall be removed
    wb_areas_l = setdiff(left_areas, removed_l);
    wb_areas_r = setdiff(right_areas, removed_r);
    
    %% load func data
    func_lh = fullfile([PATH_IN, subject,'_lh_filtered_timeseries_per_glasser_region.txt']);
    func_rh = fullfile([PATH_IN, subject,'_rh_filtered_timeseries_per_glasser_region.txt']);
    
    sub_name = [sub_name; subject];

    lh_mat = readmatrix(func_lh, 'Delimiter', ' ', 'DecimalSeparator', '.');
    rh_mat = readmatrix(func_rh, 'Delimiter', ' ', 'DecimalSeparator', '.');
    
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

    % in case, parcel 120 is ssing as resampling consequence, fill matrix
    % with zero column to account for size mismatch
    if size(lh_mat,2) == 179 
        lh_mat = [lh_mat(:, 1:119), zeros(size(lh_mat, 1), 1), lh_mat(:, 120:end)];
    end
    if size(rh_mat,2) == 179
        rh_mat = [rh_mat(:, 1:119), zeros(size(rh_mat, 1), 1), rh_mat(:, 120:end)];
    end
    %fprintf('%s: Dimension of lh: %i', subject, size(lh_mat,2));
    %fprintf('%s: Dimension of rh: %i', subject, size(rh_mat,2));

    func_mat = [lh_mat, rh_mat];

    func_mat_L = func_mat(:,wb_areas_l);
    func_mat_R = func_mat(:,wb_areas_r);

    % Fisher-Z correlation matrix
    Z_func_corr_mat_L = atanh(corrcoef(func_mat_L));
    Z_func_corr_mat_L(isinf(Z_func_corr_mat_L)) = NaN;

    Z_func_corr_mat_R = atanh(corrcoef(func_mat_R));
    Z_func_corr_mat_R(isinf(Z_func_corr_mat_R)) = NaN;

    %% load dti data
    dti = fullfile(PATH_IN_dti, [subject,'_mean_FA_connectome.csv']);
    W = readmatrix(dti, 'Delimiter', ',', 'DecimalSeparator', '.');
    
    W_L = W(wb_areas_l, wb_areas_l);
    W_R = W(wb_areas_r, wb_areas_r);

    

    %% diffusion efficiency (Ediff)

    % Compute diffusion efficiency normally
    [GEdiff_L,Ediff_L] = diffusion_efficiency(W_L);
    EdiffSymm_L = (Ediff_L + transpose(Ediff_L))/2;

    % initialize empty correlation vector
    Ediff_corrs_L = [];
    for n=1:length(Z_func_corr_mat_L)
        correlation = corr(Z_func_corr_mat_L(:,n), EdiffSymm_L(:,n), 'Rows', 'complete');
        Ediff_corrs_L = [Ediff_corrs_L, correlation];
    end
    Ediff_all_L = [Ediff_all_L; mean(Ediff_corrs_L)];

    % Compute diffusion efficiency normally
    [GEdiff_R,Ediff_R] = diffusion_efficiency(W_R);
    EdiffSymm_R = (Ediff_R + transpose(Ediff_R))/2;

    % initialize empty correlation vector
    Ediff_corrs_R = [];
    for n=1:length(Z_func_corr_mat_R)
        correlation = corr(Z_func_corr_mat_R(:,n), EdiffSymm_R(:,n), 'Rows', 'complete');
        Ediff_corrs_R = [Ediff_corrs_R, correlation];
    end
    Ediff_all_R = [Ediff_all_R; mean(Ediff_corrs_R)];


    %% search information (SI)
    L_L = weight_conversion(W_L, 'lengths');
    SI_L = search_information(W_L, L_L);
    SISymm_L = -1*((SI_L + transpose(SI_L))/2); % make symmetric and switch sign to reflect ease of communication

    % initialize empty correlation vector
    SI_corrs_L = [];
    for n=1:length(Z_func_corr_mat_L)
        correlation = corr(Z_func_corr_mat_L(:,n), SISymm_L(:,n), 'Rows', 'complete');
        SI_corrs_L = [SI_corrs_L, correlation];
    end
    SI_all_L = [SI_all_L; mean(SI_corrs_L)];


    L_R = weight_conversion(W_R, 'lengths');
    SI_R = search_information(W_R, L_R);
    SISymm_R = -1*((SI_R + transpose(SI_R))/2); % make symmetric and switch sign to reflect ease of communication

    % initialize empty correlation vector
    SI_corrs_R = [];
    for n=1:length(Z_func_corr_mat_R)
        correlation = corr(Z_func_corr_mat_R(:,n), SISymm_R(:,n), 'Rows', 'complete');
        SI_corrs_R = [SI_corrs_R, correlation];
    end
    SI_all_R = [SI_all_R; mean(SI_corrs_R)];

    %% path transitivity (T)
    T_L = path_transitivity(L_L); % no transform needs to be supplied if input is a length matrix

    % initialize empty correlation vector
    T_corrs_L = [];
    for n=1:length(Z_func_corr_mat_L)
        correlation = corr(Z_func_corr_mat_L(:,n), T_L(:,n), 'Rows', 'complete');
        T_corrs_L = [T_corrs_L, correlation];
    end
    T_all_L = [T_all_L; mean(T_corrs_L)];

    T_R = path_transitivity(L_R); % no transform needs to be supplied if input is a length matrix

    % initialize empty correlation vector
    T_corrs_R = [];
    for n=1:length(Z_func_corr_mat_R)
        correlation = corr(Z_func_corr_mat_R(:,n), T_R(:,n), 'Rows', 'complete');
        T_corrs_R = [T_corrs_R, correlation];
    end
    T_all_R = [T_all_R; mean(T_corrs_R)];

    %% resource efficiency (Eres)
    % thresholding to get unweighted matrix (UW)
    p = 0.5;
    W_thr_L = threshold_proportional(W_L, p);
    UW_L = W_thr_L > 0;

    % calculating resource efficiency 
    lambda = .5; % in Goni et al., 2013 lambda = 0.5
    [Eres_L,prob_SPL_L] = resource_efficiency_bin(UW_L,lambda);

    EresSymm_L = (Eres_L + transpose(Eres_L))/2; % make symmetric

    % initialize empty correlation vector
    Eres_corrs_L = [];
    for n=1:length(Z_func_corr_mat_L)
        correlation = corr(Z_func_corr_mat_L(:,n), EresSymm_L(:,n), 'Rows', 'complete');
        Eres_corrs_L = [Eres_corrs_L, correlation];
    end
    Eres_all_L = [Eres_all_L; mean(Eres_corrs_L)];

    
    W_thr_R = threshold_proportional(W_R, p);
    UW_R = W_thr_R > 0;

    % calculating resource efficiency 
    [Eres_R,prob_SPL_R] = resource_efficiency_bin(UW_R,lambda);

    EresSymm_R = (Eres_R + transpose(Eres_R))/2; % make symmetric

    % initialize empty correlation vector
    Eres_corrs_R = [];
    for n=1:length(Z_func_corr_mat_R)
        correlation = corr(Z_func_corr_mat_R(:,n), EresSymm_R(:,n), 'Rows', 'complete');
        Eres_corrs_R = [Eres_corrs_R, correlation];
    end
    Eres_all_R = [Eres_all_R; mean(Eres_corrs_R)];

    %% communicability (F)
    % Code von https://github.com/brain-networks/local_scfc/blob/main/fcn/communicability_wei.m
    N_L = size(W_L,1);
    B_L = sum(W_L')';
    C_L = diag(B_L);
    D_L = C_L^(-(1/2));
    E_L = D_L * W_L * D_L;
    F_L = expm(E_L);
    F_L = F_L.*~eye(N_L);

    % initialize empty correlation vector
    F_corrs_L = [];
    for n=1:length(Z_func_corr_mat_L)
        correlation = corr(Z_func_corr_mat_L(:,n), F_L(:,n), 'Rows', 'complete');
        F_corrs_L = [F_corrs_L, correlation];
    end
    F_all_L = [F_all_L; mean(F_corrs_L)];


    N_R = size(W_R,1);
    B_R = sum(W_R')';
    C_R = diag(B_R);
    D_R = C_R^(-(1/2));
    E_R = D_R * W_R * D_R;
    F_R = expm(E_R);
    F_R = F_R.*~eye(N_R);

    % initialize empty correlation vector
    F_corrs_R = [];
    for n=1:length(Z_func_corr_mat_R)
        correlation = corr(Z_func_corr_mat_R(:,n), F_R(:,n), 'Rows', 'complete');
        F_corrs_R = [F_corrs_R, correlation];
    end
    F_all_R = [F_all_R; mean(F_corrs_R)];


end

%% combine data to table
dataTable = table(sub_name, Ediff_all_L,Ediff_all_R, SI_all_L, SI_all_R, T_all_L, T_all_R, Eres_all_L, Eres_all_R, F_all_L, F_all_R);

%% compute split-half reliabilities
rel_Ediff = corr(dataTable{:,2},dataTable{:,3});
rel_Ediff_sh = (2*rel_Ediff)/(1+rel_Ediff);

rel_SI = corr(dataTable{:,4},dataTable{:,5});
rel_SI_sh = (2*rel_SI)/(1+rel_SI);

rel_T = corr(dataTable{:,6},dataTable{:,7});
rel_T_sh = (2*rel_T)/(1+rel_T);

rel_Eres = corr(dataTable{:,8},dataTable{:,9});
rel_Eres_sh = (2*rel_Eres)/(1+rel_Eres);

rel_F = corr(dataTable{:,10},dataTable{:,11});
rel_F_sh = (2*rel_F)/(1+rel_F);

min([rel_Ediff_sh,rel_SI_sh,rel_T_sh,rel_Eres_sh,rel_F_sh])
max([rel_Ediff_sh,rel_SI_sh,rel_T_sh,rel_Eres_sh,rel_F_sh])

