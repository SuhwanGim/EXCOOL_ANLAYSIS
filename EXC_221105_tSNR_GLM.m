%% ADDpath
addpath(genpath('/Volumes/homeo/github/canlab/CanlabCore'))
addpath(genpath('/Volumes/homeo/github/cocoanlab/cocoanCORE'));
addpath(genpath('/Volumes/homeo/github/cocoanlab/humanfmri_preproc_bids'));
addpath(genpath('/Volumes/homeo/dropbox/resources/spm12'));
addpath(genpath('/Volumes/homeo/github/canlab/preprocess'))
%% Building HRF Build model: Onsets (with durations) convolve with HRF (need SPM on path)
basedir = '/Users/WIWIFH/Dropbox/Projects/EXCOOL/data'; % suhwan macbook
% basedir = '/Volumes/homeo/EXCOOL'; % homeo
%% Load data
sub_i = 2; % subject number
loadnames = filenames(fullfile(basedir,[sprintf('EXC_task_Sub-EX%03d_run*',sub_i) '.csv']));
T = []; 
for run_i = 1:length(loadnames) % run number 
    T{run_i} = readtable(loadnames{run_i});
end
%% FILES 
niifiles = []; 
TRlen = []; 
% basedir = '/Volumes/homeo/EXCOOL'; % homeo
subjdir = fullfile(basedir,'preprocessed','func',sprintf('sub-EXCOOL%03d',sub_i));
rp = []; %realignment parameters 
for run_i = 1:length(loadnames)
    
    niifiles{run_i}{1,:} = filenames(fullfile(subjdir, sprintf('func_task-task_run-%02d_bold',run_i), sprintf('swdcrfunc_task-task_run-%02d*_e2.nii',run_i)));             % single-echo
    niifiles{run_i}{2,:} = filenames(fullfile(subjdir, sprintf('func_task-task_run-%02d_bold',run_i), sprintf('swdcrfunc_task-task_run-%02d*_optcom.nii',run_i)));         % optimal-combined
    niifiles{run_i}{3,:} = filenames(fullfile(subjdir, sprintf('func_task-task_run-%02d_bold',run_i), sprintf('swdcrfunc_task-task_run-%02d_*optcomDenoised.nii',run_i))); % ME_ICA
    
    
    rp{run_i}=importdata(filenames(fullfile(subjdir,sprintf('func_task-task_run-%02d_bold',run_i) ,'*e2.par')));    
    TRlen(run_i,:)=length(spm_vol(niifiles{run_i}{1})); % Get TR length
end

%% Make event regressors single-trial cue model
% This model treats each cue event as a distinct event regressor.
%   - 18 event regressors for each cue event
%   - 1 event regressor for all rating events
X = []; 
% set parameters
TR = 1;                            % repeation time in pulse sequence
ratings = []; 
for run_i = 1:length(loadnames)
    numberOfTR = TRlen(run_i,:);   % it could be differed accross participant and/or run
    len = TR * numberOfTR;
    % set
    events = [];
    events_names = [];
    for i = 1:length(T{run_i}.onsets_cue)
        events{i} = [T{run_i}.onsets_cue(i) T{run_i}.durations_cue(i)];          % 18 cue event regressors
        
        % regressor info
        %   
        targets = [];
        if T{run_i}.Tcondition(i) == 1
            targets = 'self';
        elseif T{run_i}.Tcondition(i) == 2
            targets = 'friend';
        elseif T{run_i}.Tcondition(i) == 3
            targets = 'Celeb';
        end
            
        vals =[];
        if T{run_i}.Wvalence(i) == 1
            vals = 'Pos'; % Positive word
        elseif T{run_i}.Wvalence(i) == -1
            vals = 'neg'; % Negative word
        end
        
        ratings{run_i}(i,:) = T{run_i}.ratings(i); % ratings from -0.5 to 0.5        
        events_names{i} = {sprintf('Run%02d_Trial%02d_cue_%s_%s',run_i, i,targets,vals)};
    end
    events{end+1} = [T{run_i}.onsets_ratings T{run_i}.durations_ratings];        % ALL ratings events
    events_names{end+1} = {'Ratings_ALL'};
    X{run_i} = onsets2fmridesign(events, 1, len, spm_hrf(1));      % Making HRF convolved event regressors
end
%% PLOT
create_figure('X 2nd plot');
h = plot_matrix_cols(zscore(X{run_i}(:, 1:size(X,2)-1)), 'horiz');
%% Calculate tSNR
% check to be sure:
tSNR = []; 
for run_i = 1:length(niifiles)    
    for cond_i = 1:3 % condition 
        % cond_i == 1: single-echo
        %        == 2: optimally combined
        %        == 3: ME-ICA 
        tdat = fmri_dat(niifiles{run_i}{cond_i});        
        tdat.dat(isnan(tdat.dat)) = 0;
        m = mean(tdat.dat',1)';     % Mean values of each voxel
        s = std(tdat.dat',1)';      % Std of each voxel
        d = m./s;                   % Temporal Signal-to-Noise Ratio 
                                    % -> It is calculated by dividing the mean
                                    % of a time series by its standard deviation
        d(m == 0 | s == 0) = 0;
        
        % tSNR map 
        tSNR{run_i}{cond_i} = tdat;
        tSNR{run_i}{cond_i}.dat = d;        
    end    
end
clear tdat
run_i = 1;
% see tSNR examples
orthviews_multiple_objs({tSNR{run_i}{1} tSNR{run_i}{2} tSNR{run_i}{3}});
%% Extracting CSF AND WM (from single-echo image and optimally combined)
%
% :We will assume that the influence of CSF and WM is denoised by ME-ICA
% and extract the signals of CSF and WM from other two conditions.
%
for run_i = 1:length(niifiles)
    for cond_i = 1:2  
        % cond_i == 1: single-echo
        %        == 2: optimally combined
        tdat = fmri_data(niifiles{run_i}{cond_i});
        [~, components] = extract_gray_white_csf(tdat, 'masks', ...
        {'gray_matter_mask.nii', 'canonical_white_matter_thrp5_ero1.nii', ...
        'canonical_ventricles_thrp5_ero1.nii'});
        wm_csf{run_i}{cond_i} = [scale(double(components{2})) scale(double(components{3}))];
    end    
end
clear tdat
%% REGRESS 
% * Construct regressor matrix. 
%   : A regressor matrix consists of 
%       1) event regressors 
%       2) motion alignment 
%       3) CSF and WM
%       4) linear drift 
%   : Additionally, high-bass filtering will be applied for both matrix and
%   fMRI data. 
%
% * Regressor combinations
%   cond_i == 1: single-echo         with [Event regressors + Motion parameters + CSF and WM + linear drift]
%          == 2: optimally combined  with [Event regressors + Motion parameters + CSF and WM + linear drift]
%          == 3: ME-ICA              with [Event regressors + linear drift]
%

% set parameter for filtering
filt(1) = 1/100; % 100 secs (see ref: Chavez and Wager, 2020, JPSP)
filt(2) = Inf; % for highpassfilter

for run_i = 1:length(niifiles)    
    for cond_i = 1:3
        % 0. load fMRI data
        tdat = fmri_data(niifiles{run_i}{cond_i});
        
        % 1. reconstruct regressor matrix 
        R = []; comb_R = []; 
        if cond_i == 3
            R = [X{run_i} zscores((1:TRlen(run_i,:))')]; 
        else % for cond_i = 1, 2            
            % motion (realign parameters)
            temp_mv = [[rp{run_i} rp{run_i}.^2 ...
                    [zeros(1,6); diff(rp{run_i})] [zeros(1,6); diff(rp{run_i})].^2]];            
                
            R = [X{run_i} temp_mv wm_csf{cond_i} zscores((1:TRlen(run_i,:))')]; 
            
        end
        comb_R = R;
        nEvents = size(X{run_i},2); % the number of events 
        nBetas = size(comb_R,2);    % the number of all regressors 
        
        
        % 2. temporal filtering 
        fprintf('::: temporal filtering brain data... \n');
        % conn_filter does mean-center by default.
        % We need mean back in
        
        % Pass-filter brain data
        % -------------------------------------------
        % TR = 1 secs [repetion time in fMRI pulse sequence]
        tdat.dat = conn_filter(TR, filt, tdat.dat', 'full')' ...
            + repmat(mean(tdat.dat,2), 1, size(tdat.dat,2));
        
        % Pass-filter regressors         
        % -------------------------------------------
        outlier_idx  =[];
        for rr_i = 1:size(comb_R,2)
            if unique(comb_R(:,rr_i)) == 1
                outlier_idx = [outlier_idx rr_i];
            end
        end
        outliers = comb_R(:,outlier_idx);
        comb_R(:,outlier_idx) = [];
        
        comb_R = conn_filter(TR, filt, comb_R, 'full') ...
            + repmat(mean(comb_R), size(comb_R,1), 1);
        
        % add outlier indices back in
        comb_R = [comb_R outliers];
        
        
        % 3. REGRESSION
        % regression
        fprintf('::: Regression brain data... \n');
                
        beta_dat = dat;
        beta_dat.dat = [];
        beta_dat.dat = dat.dat * pinv(comb_R)';
        beta_dat.covariates = comb_R; 
        %     beta_dat.dat = beta_dat.dat(:,1:nEvents);
        
        
        
        % 4. save
        if cond_i == 3
            foldernames = 'MEICA';
        elseif cond_i == 2
            foldernames = 'OptComb';
        elseif cond_i == 1
            foldernames = 'SingleEcho';
        end
        % basedir = '/Volumes/homeo/EXCOOL'; % homeo
        subjsavdir = fullfile(basedir,'first_level','MODEL1_single_trial_cues',foldernames,sprintf('sub-EXCOOL%03d',sub_i));
        
        if ~exist(subjsavdir); mkdir(subjsavdir); end
        fprintf('::: save brain data... \n');
        save(fullfile(subjsavdir,'betas_all_obj.mat'),'beta_dat');        
        fprintf('::: donez... \n');
        
        
        beta_dat.dat = beta_dat.dat(:,1:nEvents);
        %events_names;
        for tri_i = 1:length(nEvents)
            templist{tri_i} = fullfile(subjsavdir, [events_names{tri_i} '.nii']);
            tempdat = beta_dat.get_wh_image(tri_i);
            tempdat.fullpath = char(templist{tri_i});
            write(tempdat);
        end
    end
end
%% 

% DONE
%% Visualization 
% basedir = '/Volumes/homeo/EXCOOL'; % homeo
if cond_i == 3
    foldernames = 'MEICA';
elseif cond_i == 2
    foldernames = 'OptComb';
elseif cond_i == 1
    foldernames = 'SingleEcho';
end
sub_i = 1;
run_i = 1;
subjsavdir = filenames(fullfile(basedir,'first_level','MODEL1_single_trial_cues','*','Run*trial*cue*.nii'));

% see tSNR examples
orthviews_multiple_objs({tSNR{run_i}{1} tSNR{run_i}{2} tSNR{run_i}{3}});
