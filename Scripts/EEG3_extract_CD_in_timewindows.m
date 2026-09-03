%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      generating analysis data for CD modeling and create      %
%                Figure CD time course                                    %
% Author:        Henrike Jungeblut, University of Mainz                   % 
% Last Modified: 09/2026                                                  %
% Note:                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ---------- Necessary Functions -----------------------------------------
% %
% % EEGLAB version eeglab2025.0.0 with the following plugins is used:
% % The EEGLAB toolbox with all plugins installed is provided with this
% % code within the folder "Tools".


%% define paths
filepath   = cd;
parts = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});
PATH_IN_ERP = fullfile(PATH_MAIN, 'Code_package', 'erp');
PATH_IN_CD = fullfile(PATH_MAIN, 'Code_package', 'subject_connydegree_traces');
addpath(fullfile(PATH_MAIN, 'Code_package', 'Tools', 'eeglab2025.0.0'));

%% load erp data
erp = readtable(fullfile(PATH_IN_ERP, 'erp_ALS.csv'),'PreserveVariableNames',true);
erp([41 162 186 187],:) = []; % remove subjects sub-041, sub-164, sub-189, sub-190 who did not adhere to task instructions or display low overall accuracies
%% define parameters
% average P3 latency in the for conditions
P3_OE_R = table2array(mean(erp(:,"OE_R"), "omitnan"));
P3_OE_S = table2array(mean(erp(:,"OE_S"), "omitnan")); 
P3_LM_R = table2array(mean(erp(:,"LM_R"), "omitnan"));
P3_LM_S = table2array(mean(erp(:,"LM_S"), "omitnan")); 

% SD of the mean P3 in the four conditions
SD_P3_OE_R = std(table2array(erp(:,"OE_R")), "omitnan");
SD_P3_OE_S = std(table2array(erp(:,"OE_S")), "omitnan");
SD_P3_LM_R = std(table2array(erp(:,"LM_R")), "omitnan");
SD_P3_LM_S = std(table2array(erp(:,"LM_S")), "omitnan");

% time windows
E_OE_R = round([1, (P3_OE_R - 1.5*SD_P3_OE_R) - 1]);
E_OE_S = round([1, (P3_OE_S - 1.5*SD_P3_OE_S) - 1]);
E_LM_R = round([1, (P3_LM_R - 1.5*SD_P3_LM_R) - 1]);
E_LM_S = round([1, (P3_LM_S - 1.5*SD_P3_LM_S) - 1]);

L_OE_R = round([P3_OE_R - 1.5*SD_P3_OE_R, P3_OE_R + 1.5*SD_P3_OE_R]);
L_OE_S = round([P3_OE_S - 1.5*SD_P3_OE_S, P3_OE_S + 1.5*SD_P3_OE_S]);
L_LM_R = round([P3_LM_R - 1.5*SD_P3_LM_R, P3_LM_R + 1.5*SD_P3_LM_R]);
L_LM_S = round([P3_LM_S - 1.5*SD_P3_LM_S, P3_LM_S + 1.5*SD_P3_LM_S]);


%% load CD data

G_channels = "G" + compose("%02d", 1:30);
Y_channels = "Y" + compose("%02d", 1:30);
channels = [G_channels, Y_channels, "Cz"];

for ch = channels

    filenameOER = "connydegree_traces_chan_" + ch + "_pli_OE-R.csv";
    filenameOES = "connydegree_traces_chan_" + ch + "_pli_OE-S.csv";
    filenameLMR = "connydegree_traces_chan_" + ch + "_pli_LM-R.csv";
    filenameLMS = "connydegree_traces_chan_" + ch + "_pli_LM-S.csv";

    CD_OE_R = table2array(readtable(fullfile(PATH_IN_CD, filenameOER)));
    CD_OE_R(:,1:2000) = []; % wo genau ist der Target Onset?
    CD_OE_R([41 160 183 184],:) = []; % remove subjects sub-041, sub-162, sub-186, sub-187 who scored below chance level in the cognitive control task
    CD_OE_S = table2array(readtable(fullfile(PATH_IN_CD, filenameOES)));
    CD_OE_S(:,1:2000) = []; 
    CD_OE_S([41 160 183 184],:) = []; % remove subjects sub-041, sub-162, sub-186, sub-187 who scored below chance level in the cognitive control task
    CD_LM_R = table2array(readtable(fullfile(PATH_IN_CD, filenameLMR)));
    CD_LM_R(:,1:2000) = []; 
    CD_LM_R([41 160 183 184],:) = []; % remove subjects sub-041, sub-162, sub-186, sub-187 who scored below chance level in the cognitive control task
    CD_LM_S = table2array(readtable(fullfile(PATH_IN_CD, filenameLMS)));
    CD_LM_S(:,1:2000) = []; 
    CD_LM_S([41 160 183 184],:) = []; % remove subjects sub-041, sub-162, sub-186, sub-187 who scored below chance level in the cognitive control task


    %% loop through subjects and CD doubles

    % initialize vectors
    E_CD_OE_R = [];
    E_CD_OE_S = [];
    E_CD_LM_R = [];
    E_CD_LM_S = [];
    L_CD_OE_R = [];
    L_CD_OE_S = [];
    L_CD_LM_R = [];
    L_CD_LM_S = [];

    for n = 1:size(CD_OE_R,1)
        E_CD_OE_R = [E_CD_OE_R; mean(CD_OE_R(n, E_OE_R(1):E_OE_R(2)))];
        E_CD_OE_S = [E_CD_OE_S; mean(CD_OE_S(n, E_OE_S(1):E_OE_S(2)))];
        E_CD_LM_R = [E_CD_LM_R; mean(CD_LM_R(n, E_LM_R(1):E_LM_R(2)))];
        E_CD_LM_S = [E_CD_LM_S; mean(CD_LM_S(n, E_LM_S(1):E_LM_S(2)))];
        L_CD_OE_R = [L_CD_OE_R; mean(CD_OE_R(n, L_OE_R(1):L_OE_R(2)))];
        L_CD_OE_S = [L_CD_OE_S; mean(CD_OE_S(n, L_OE_S(1):L_OE_S(2)))];
        L_CD_LM_R = [L_CD_LM_R; mean(CD_LM_R(n, L_LM_R(1):L_LM_R(2)))];
        L_CD_LM_S = [L_CD_LM_S; mean(CD_LM_S(n, L_LM_S(1):L_LM_S(2)))];
    end

    % append to erp table
    erp.E_CD_OE_R  = E_CD_OE_R;
    erp.E_CD_OE_S  = E_CD_OE_S;
    erp.E_CD_LM_R  = E_CD_LM_R;
    erp.E_CD_LM_S  = E_CD_LM_S;

    erp.L_CD_OE_R  = L_CD_OE_R;
    erp.L_CD_OE_S  = L_CD_OE_S;
    erp.L_CD_LM_R  = L_CD_LM_R;
    erp.L_CD_LM_S  = L_CD_LM_S;

    outfile = ch + "_CD_in_timewindows_4subs_removed.csv";
    writetable(erp,fullfile(PATH_MAIN, 'Code_package', 'Analysis_data', outfile));
end


%% plotting
OES = table2array(readtable(fullfile(PATH_IN_CD, 'connydegree_traces_chan_G04_pli_OE-S.csv'), 'PreserveVariableNames',true));
OES_smoothed = smoothdata(OES,2,'Gaussian',200);

OER = table2array(readtable(fullfile(PATH_IN_CD, 'connydegree_traces_chan_G04_pli_OE-R.csv'), 'PreserveVariableNames',true));
OER_smoothed = smoothdata(OER,2,'Gaussian',200);

LMS = table2array(readtable(fullfile(PATH_IN_CD, 'connydegree_traces_chan_G04_pli_LM-S.csv'), 'PreserveVariableNames',true));
LMS_smoothed = smoothdata(LMS,2,'Gaussian',200);

LMR = table2array(readtable(fullfile(PATH_IN_CD, 'connydegree_traces_chan_G04_pli_LM-R.csv'), 'PreserveVariableNames',true));
LMR_smoothed = smoothdata(LMR,2,'Gaussian',200);

figure;
h1 = plot(-200:1000, mean(OES_smoothed(:,1800:3000),1), 'LineWidth', 2, 'Color', '#bae4bc');
ylim([12.2 14.4]);
xlim([-200 1000]);
hold on
h2 = plot(-200:1000, mean(OER_smoothed(:,1800:3000),1), 'LineWidth', 2, 'Color', '#7bccc4');
h3 = plot(-200:1000, mean(LMS_smoothed(:,1800:3000),1), 'LineWidth', 2, 'Color', '#43a2ca');
h4 = plot(-200:1000, mean(LMR_smoothed(:,1800:3000),1), 'LineWidth', 2, 'Color', '#0868ac');
xline(0, '--k', 'LineWidth', 1.5);
xlabel('Time in ms', 'FontSize', 14);
ylabel('Average CD', 'FontSize', 14);

patch([E_OE_S(1) E_OE_S(2) E_OE_S(2) E_OE_S(1)], ...
      [14.35 14.35 14.4 14.4], [186 228 188]/255, 'FaceAlpha', 1, 'EdgeColor', 'none');
patch([L_OE_S(1) L_OE_S(2) L_OE_S(2) L_OE_S(1)], ... 
      [14.35 14.35 14.4 14.4], [186 228 188]/255, 'FaceAlpha', 1, 'EdgeColor', 'none');

patch([E_OE_R(1) E_OE_R(2) E_OE_R(2) E_OE_R(1)], ...
      [14.3 14.3 14.35 14.35], [123 204 196]/255, 'FaceAlpha', 1, 'EdgeColor', 'none');
patch([L_OE_R(1) L_OE_R(2) L_OE_R(2) L_OE_R(1)], ... 
      [14.3 14.3 14.35 14.35], [123 204 196]/255, 'FaceAlpha', 1, 'EdgeColor', 'none');

patch([E_LM_S(1) E_LM_S(2) E_LM_S(2) E_LM_S(1)], ...
      [14.25 14.25 14.3 14.3], [67 162 202]/255, 'FaceAlpha', 1, 'EdgeColor', 'none');
patch([L_LM_S(1) L_LM_S(2) L_LM_S(2) L_LM_S(1)], ... 
      [14.25 14.25 14.3 14.3], [67 162 202]/255, 'FaceAlpha', 1, 'EdgeColor', 'none');

patch([E_LM_R(1) E_LM_R(2) E_LM_R(2) E_LM_R(1)], ...
      [14.2 14.2 14.25 14.25], [8 104 172]/255, 'FaceAlpha', 1, 'EdgeColor', 'none');
patch([L_LM_R(1) L_LM_R(2) L_LM_R(2) L_LM_R(1)], ... 
      [14.2 14.2 14.25 14.25], [8 104 172]/255, 'FaceAlpha', 1, 'EdgeColor', 'none');

legend([h1 h2 h3 h4], {'OE-S', 'OE-R', 'LM-S', 'LM-R'}, 'Location', 'best');

saveas(gcf,fullfile(PATH_MAIN, 'Code_package', 'Figures', 'EEG3_CD_timecourse.png'));