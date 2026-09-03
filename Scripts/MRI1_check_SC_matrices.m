%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      check SC matrices for zero nodes                         %
% Author:        Henrike Jungeblut, University of Mainz                   % 
% Last Modified: 09/2026                                                  %
% Note:                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set directories
filepath   = cd;
parts = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});
PATH_IN_dti = fullfile(PATH_MAIN, 'Code_package', 'MRI', 'dwi');

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

sub_name = [];

%% loop through subject list
for sub = 1:length(subjectNames_filtered)
    subject = subjectNames_filtered{sub};
    sub_name = [sub_name; subject];

   
    %% load dti data
    %filename = 'V:/PsySchubert/projects/SC_FC_coupling_connectivity_HJ_MH_ALS/4_derivative/structural_connectomes/sub-001_mean_FA_connectome.csv';
    dti = fullfile(PATH_IN_dti, [subject,'_mean_FA_connectome.csv']);
    W = readmatrix(dti, 'Delimiter', ',', 'DecimalSeparator', '.');
    
    % remove area 120 on right and left side (120 and 180+120 = 300) due to
    % findings from script M0.
    W([120,300],:) = [];
    W(:,[120,300]) = [];

    % --- check for zero-degree nodes in this subject ---
    deg = sum(W > 0, 2);            % number of connections per node
    zero_nodes = find(deg == 0);   % indices of nodes with no connections
    
    % store result for later across-subjects analysis
    zero_node_list{sub} = zero_nodes;
    
    if ~isempty(zero_nodes)
        fprintf('Subject %s: %d zero nodes -> %s\n', subject, ...
            numel(zero_nodes), mat2str(zero_nodes));
    else
        fprintf('Subject %s: no zero nodes\n', subject);
    end


    Nnodes = size(W, 1);  % or set equal to your parcellation size
    all_zero_counts = zeros(1, Nnodes);
    
    for s = 1:length(zero_node_list)
        zn = zero_node_list{s};
        all_zero_counts(zn) = all_zero_counts(zn) + 1;
    end
    
    % which parcels are often zero?
    often_zero = find(all_zero_counts > 1);  % threshold: appears in >1 subject
    fprintf('\nNodes often disconnected across subjects:\n');
    disp(often_zero);

end