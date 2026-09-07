%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         Connectivity-based predictors of intelligence            %
% Subtitle:      plot latent correlations CD and intelligence             %
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
filepath    = cd;
parts       = regexp(filepath,filesep,'split');
PATH_MAIN   = fullfile(parts{1:end-2});
EEG         = load(fullfile(PATH_MAIN, 'Code_package', 'Data', 'chanlocs.mat'));
PATH_IN_cors = fullfile(PATH_MAIN, 'Code_package','Analysis_data');
addpath(fullfile(PATH_MAIN, 'Code_package', 'Tools'));
addpath(fullfile(PATH_MAIN, 'Code_package', 'Tools', 'eeglab2025.0.0'));


% plotting
latent_cors_late =  table2array(readtable(fullfile(PATH_IN_cors, 'latent_cors_late.csv')));

figure;
topoplot(latent_cors_late, EEG.x, ...
    'style', 'straight', ...
    'electrodes','labels');
title({'CD in the Late', 'Time Window', 'and Intelligence'},FontSize=16);
colormap(flipud(brewermap(256,'RdBu')))   % or any ColorBrewer map
clim([-0.3 0.3]);
h = colorbar;
h.Position = [0.92 0.11 0.02 0.815];   % moves the colorbar to the right side
h.Label.String = 'Latent Correlation';
h.Label.FontSize = 12;

saveas(gcf,fullfile(PATH_MAIN, 'Code_package', 'Figures', 'EEG4_CD_latent_corrs.png'));