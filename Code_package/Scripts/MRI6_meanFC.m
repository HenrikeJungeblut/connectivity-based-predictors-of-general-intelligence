%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      calculate average and absolute FC                        %
% Author:        Henrike Jungeblut, University of Mainz                   %
% Last Modified: 09/2026                                                  %
% Note:          This analysis is EXPLORATORY                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set paths
filepath   = cd;
parts = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});
PATH_IN_dti = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'dwi');

%% loop through data_type
for data_type = ["rest", "res", "rawfunc"]    %res: background connectivity/rawfunc: task-evoked


    if strcmp(data_type, 'rest')
        PATH_IN = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'func', 'rest_timeseries_per_region', 'rest_');
        PATH_OUT = fullfile(PATH_MAIN, 'Code_package', 'Analysis_data', 'averageFC_rest.csv');
    elseif strcmp(data_type, 'res')
        PATH_IN = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'func', 'timeseries_per_region', 'res_');
        PATH_OUT = fullfile(PATH_MAIN, 'Code_package', 'Analysis_data', 'averageFC_res.csv');
    elseif strcmp(data_type, 'rawfunc')
        PATH_IN = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'func', 'rawfunc_timeseries_per_region', 'rawfunc_');
        PATH_OUT = fullfile(PATH_MAIN, 'Code_package', 'Analysis_data', 'averageFC_rawfunc.csv');
    end

    % loading subject list:
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

    %% frontoparietal parcels
    frontoparietal_areas_orig = [14 15 29 63 73 77 80 82 83 85 89 91 ...
        92 97 98 111 133 144 145 149 170 171 ...
        194 195 209 238 242 243 253 254 257 260 261 263 265 ...
        269 271 272 273 277 278 291 313 ...
        324 325 329 342 350 351 357];

    overall_FC =  [];
    overall_FC_abs =  [];
    sub_name = [];



    %% loop through subject list

    for sub = 1:length(subjectNames_filtered)
        if isempty(subjectNames_filtered{sub})
            continue
        end

        subject = subjectNames_filtered{sub};

        % Subject-specific removals
        if strcmp(subject, 'sub-019')
            removed = [118 120 300];
        elseif strcmp(subject, 'sub-021')
            removed = [120 172 300 352];
        elseif strcmp(subject, 'sub-025') || strcmp(subject, 'sub-083') % zero nodes
            removed = [120 300 332];
        elseif strcmp(subject, 'sub-197')
            removed = [120 183 300];
        else
            removed = [120 300]; %resampling consequence
        end

        % get parcels in frontoparietal_areas_orig that are not in removed ->
        % exclude all parcels that shall be removed
        frontoparietal_areas = setdiff(frontoparietal_areas_orig, removed);

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

        func_mat = [lh_mat, rh_mat];
        func_mat = func_mat(:,frontoparietal_areas);

        % Fisher-Z correlation matrix
        Z_func_corr_mat = atanh(corrcoef(func_mat));
        Z_func_corr_mat(isinf(Z_func_corr_mat)) = NaN;


        %% mean FC (average of unique connections)
        % FC = NxN functional connectivity matrix
        N = size(Z_func_corr_mat,1);
        % Extract upper triangle without diagonal
        upper = Z_func_corr_mat(triu(true(N),1));
        % Mean overall FC
        overall_FC = [overall_FC; mean(upper)];

        %% mean absolute FC (overall connection strength
        N = size(Z_func_corr_mat,1);
        upper = Z_func_corr_mat(triu(true(N),1));
        overall_FC_abs = [overall_FC_abs; mean(abs(upper))];



    end

    %% combine data to table
    dataTable = table(sub_name, overall_FC, overall_FC_abs);
    writetable(dataTable, PATH_OUT);
end

%rest (overall_FC, overall_FC_abs):
%0.2023 0.2696

%res (overall_FC, overall_FC_abs):
%0.1999 0.2432

%rawfunc (overall_FC, overall_FC_abs):
%0.2378 0.2962