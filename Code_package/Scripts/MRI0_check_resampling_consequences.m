%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      check resampling consequences                            %
% Author:        Henrike Jungeblut, University of Mainz                   % 
% Last Modified: 09/2026                                                  %
% Note:                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define paths
filepath   = cd;
parts = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});
PATH_IN_labels = fullfile(PATH_MAIN, 'Code_package', 'atlas_labels'); 

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

% set counter to count how many times areas disappear due to resampling
rest_left_counter = 0;
rest_right_counter = 0;
task_left_counter = 0;
task_right_counter = 0;

%% loop through subject list
for sub = 1:length(subjectNames_filtered)
    if isempty(subjectNames_filtered{sub})
        continue
    end
    subject = subjectNames_filtered{sub};

    %% load label data
    label_rest_lh = fullfile(PATH_IN_labels, [subject,'_restdata_labels_lh.txt']);
    label_rest_rh = fullfile(PATH_IN_labels, [subject,'_restdata_labels_rh.txt']);
    label_task_lh = fullfile(PATH_IN_labels, [subject,'_taskdata_labels_lh.txt']);
    label_task_rh = fullfile(PATH_IN_labels, [subject,'_taskdata_labels_rh.txt']);
    rest_lh = readmatrix(label_rest_lh, 'Delimiter', ',', 'DecimalSeparator', '.');
    rest_rh = readmatrix(label_rest_rh, 'Delimiter', ',', 'DecimalSeparator', '.');
    task_lh = readmatrix(label_task_lh, 'Delimiter', ',', 'DecimalSeparator', '.');
    task_rh = readmatrix(label_task_rh, 'Delimiter', ',', 'DecimalSeparator', '.');
   
    % the label files contain all intensity values contained in the
    % respective images. It always starts with zero (background intensity)
    % and then with 1001 (first atlas parcel) to 1180 (last atlas parcel).
    % Consequently, if no area was compromized by resampling procedure,
    % each file should have 181 lines. Here, we check whether this is true
    % and if not, we investigate which line is missing. This is done by
    % checking if each difference between lines is the increment of 1. 
    if length(rest_lh) < 181
        num_areas = length(rest_lh)-1;
        missing = find(diff(rest_lh) ~= 1); %this corresponds to the missing atlas parcel, as the file starts with zero and only then lists the parcel numbers
        missing_str = sprintf('%d ', missing(2:end));
        fprintf('Rest left hemisphere %s: area(s) %sare missing (%d areas in total)\n', ...
        subject, missing_str, num_areas);
        rest_left_counter = rest_left_counter + 1;
    end

    if length(rest_rh) < 181
        num_areas = length(rest_rh)-1;
        missing = find(diff(rest_rh) ~= 1); 
        missing_str = sprintf('%d ', missing(2:end));
        fprintf('Rest right hemisphere %s: area(s) %sare missing (%d areas in total)\n', ...
        subject, missing_str, num_areas);
        rest_right_counter = rest_right_counter + 1;
    end

    if length(task_lh) < 181
        num_areas = length(task_lh)-1;
        missing = find(diff(task_lh) ~= 1); 
        missing_str = sprintf('%d ', missing(2:end));
        fprintf('Task left hemisphere %s: area(s) %sare missing (%d areas in total)\n', ...
        subject, missing_str, num_areas);
        task_left_counter = task_left_counter + 1;
    end

    if length(task_rh) < 181
        num_areas = length(task_rh)-1;
        missing = find(diff(task_rh) ~= 1); 
        missing_str = sprintf('%d ', missing(2:end));
        fprintf('Task right hemisphere %s: area(s) %sare missing (%d areas in total)\n', ...
        subject, missing_str, num_areas);
        task_right_counter = task_right_counter + 1;
    end   
  
% the script finds that on the right hemisphere, area 120 is disappearing
% due to resampling in 23 cases (in both task and rest data). On the left
% hemisphere, area 120 is disappearing in 84 cases (in both task and rest
% data). Therefore, area 120 is excluded from analysis and all durther
% calculations will rely on 358x358 matrices for func as well as DTI data.
end

% percentages of missing parcels
task_left_counter/length(subjectNames_filtered)
task_right_counter/length(subjectNames_filtered)
rest_left_counter/length(subjectNames_filtered)
rest_right_counter/length(subjectNames_filtered)