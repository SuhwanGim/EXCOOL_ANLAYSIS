function outputs = STEP2_ME_preprocesse(subj_number, basedir)
% subject directory
rawdir = fullfile(basedir,'raw');

fundir = [];
anadir = [];
fmapdir = [];
for sub_i = 1:7
    fundir{sub_i} = filenames(fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'dicom','func','func*'));
    anadir{sub_i} = filenames(fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'dicom','anat','T1*'));
    fmapdir{sub_i} = filenames(fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'dicom','fmap','*ISO*'));
end
%% DICOM 2 NIIFT for anat
sub_i = 1;
for i = 1:length(anadir{sub_i})
    %addpath
    [~,b] = fileparts(anadir{sub_i}{i});
    outputdir = fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'anat',b);
    if ~isfolder(outputdir); mkdir(outputdir); end
    system(sprintf('/usr/local/bin/dcm2niix -o %s -z n %s ',outputdir,anadir{sub_i}{i})); % on sein
    %system(sprintf('/Users/admin/dcm2niix/build/bin/dcm2niix -o %s -z n %s ',outputdir,dcmSource{i})); % on homeo
end
%% DICOM 2 NIIFT for fmap
sub_i = 1;
for i = 1:length(fmapdir{sub_i})
    %addpath
    [~,b] = fileparts(fmapdir{sub_i}{i});
    outputdir = fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'fmap',b);
    if ~isfolder(outputdir); mkdir(outputdir); end
    system(sprintf('/usr/local/bin/dcm2niix -o %s -z n %s ',outputdir,fmapdir{sub_i}{i})); % on sein
    %system(sprintf('/Users/admin/dcm2niix/build/bin/dcm2niix -o %s -z n %s ',outputdir,dcmSource{i})); % on homeo
end
%% DICOM 2 NIIFT for func
% dcm2niix /Users/suhwan/Dropbox/Projects/sync_SEP/data/coil_test/20210803_multi_echo/COCOAN_SUHWAN_GIM_20210803_111531_408000/1_2_3ME_HALF_CMRR_2_7ISO_PA_64CH+FLEX_MB4_IA2_0
sub_i = 1;
for i = 1:length(fundir{sub_i})
    %addpath
    [~,b] = fileparts(fundir{sub_i}{i});
    outputdir = fullfile(rawdir,sprintf('sub-EXCOOL%03d',sub_i),'func',b);
    if ~isfolder(outputdir); mkdir(outputdir); end
    system(sprintf('/usr/local/bin/dcm2niix -o %s -z n %s ',outputdir,fundir{sub_i}{i})); % on sein
    %system(sprintf('/Users/admin/dcm2niix/build/bin/dcm2niix -o %s -z n %s ',outputdir,dcmSource{i})); % on homeo
end
end