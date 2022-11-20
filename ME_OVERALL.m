%% Multi-echo preprocessing scripts
% by Suhwan 
clear; close;
%% ADD PATH
%% SET PATH
username = char(java.lang.System.getProperty('user.name'));
project_name = 'EXCOOL';
basedir = '/Volumes/homeo/EXCOOL';
%% Subject ID:
subj_number = 2; % from 2 to 7
sub_i = subj_number;
%% 1. SET DIRECTORY
rawdir = fullfile(basedir,'raw');
prodir = fullfile(basedir,'preprocessed',sprintf('sub-EXCOOL%03d',sub_i));
preproc_func_dir = fullfile(prodir, 'func');
if ~exist(preproc_func_dir, 'dir'), mkdir(preproc_func_dir); end
preproc_mean_func_dir = fullfile(prodir, 'mean_func');
if ~exist(preproc_mean_func_dir, 'dir'), mkdir(preproc_mean_func_dir); end
preproc_anat_dir = fullfile(prodir, 'anat');
if ~exist(preproc_anat_dir, 'dir'), mkdir(preproc_anat_dir); end
preproc_fmap_dir = fullfile(prodir, 'fmap');
if ~exist(preproc_fmap_dir, 'dir'), mkdir(preproc_fmap_dir); end
qcdir = fullfile(prodir, 'qc_images');
if ~exist(qcdir, 'dir'), mkdir(qcdir); end
%% 2. COPY FILES FROM RAW TO PREPROCESSED
subjanatdir = [];
for i = 1:length(anadir{sub_i})
    [~,b] = fileparts(anadir{sub_i}{i});
    copyfile([fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'anat',b) '/*'], fullfile(preproc_anat_dir,b));
    subjanatdir{sub_i}{i} = fullfile(preproc_anat_dir,b);
end
subjfmapdir=[];
for i = 1:length(fmapdir{sub_i})
    [~,b] = fileparts(fmapdir{sub_i}{i});
    copyfile([fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'fmap',b) '/*'], fullfile(preproc_fmap_dir,b));
    subjfmapdir{sub_i}{i} = fullfile(preproc_fmap_dir,b);
end
subjfmapdir{sub_i} = subjfmapdir{sub_i}';

subjdir = [];
for i = 1:length(fundir{sub_i})
    [~,b] = fileparts(fundir{sub_i}{i});
    copyfile([fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'func',b) '/*'], fullfile(preproc_func_dir,b));
    subjdir{sub_i}{i} = fullfile(preproc_func_dir,b);
end
subjdir{sub_i} = subjdir{sub_i}';
%% 3. spike detection, slice time correction, and realignment

%% 4. TEDANA

%% 5. distortion correction, Coregistrations, normalization, and smoothing (5 mm)

%% OTHERS
% Check the affine matrix (?)
temp_t = [];
for i =1:4
    %temp_t{i} = spm_vol([func_bold_files{i} ',2']); % raw images
    %temp_t{i} = spm_vol([afunc_bold_files{i} ',2']); % after slice-timing correction
    %temp_t{i} = spm_vol([rafunc_bold_files{i} ',2']); % after realinment
end