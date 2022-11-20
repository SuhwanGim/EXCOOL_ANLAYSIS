

function outputs = ME_SETP1(subj_number, subjdir)
sub_i = subj_number;
%% SIMPLE PREPROCESSTING
mask = []; dat = [];
spike_covariates =[];

subjdat=subjdir;
subjdat(contains(subjdat,'sbref')) = [];

for i = 1:length(subjdat) % each task
    
    temp_list = filenames(fullfile(subjdat{i},'func_*.nii'));
    mask = [];
    dat = [];
    for e_i = 1:length(temp_list) % each echo
        clf;
        [a,b]=fileparts(temp_list{e_i});
        % implicit_mask
        [~, ~, ~, ~, outputname] = fmri_mask_thresh_canlab(char(temp_list{e_i}),...
            fullfile(a, sprintf('implicit_mask_e%02d.nii',e_i)));
        implicit_mask_file = outputname;
        mask{e_i} = fmri_data(implicit_mask_file,implicit_mask_file);
        dat{e_i} = fmri_data(temp_list{e_i}, implicit_mask_file);
        dat{e_i}.images_per_session = size(dat{e_i}.dat,2);
        % spike id
        diary(fullfile(a, ['qc_diary_' b '.txt']));
        dat{e_i} = preprocess(dat{e_i}, 'outliers', 'plot');  % Spike detect and globals by slice
        subplot(5, 1, 5);
        dat{e_i} = preprocess(dat{e_i}, 'outliers_rmssd', 'plot');  % RMSSD Spike detect
        diary off;
        
        sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
        set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);
        drawnow;
        
        qcspikefilename = fullfile(a, ['qc_spike_plot_' b '.png']); % Scott added some lines to actually save the spike images
        saveas(gcf,qcspikefilename);
        spike_covariates{e_i} = dat{e_i}.covariates;
    end
    save(fullfile(a,'spike_covariates_all.mat','spike_covariates'));
    
end

%% Slice time correction
subjdat=subjdir{sub_i};
subjdat(contains(subjdat,'sbref')) = [];

for i = 1:length(subjdat) % each task
    
    temp_list = filenames(fullfile(subjdat{i},'func_*.nii'));
    a = [];
    b = [];
    for e_i = 1:length(temp_list) % each echo
        
        slice_timing_job = [];
        json_read = []; json_file = [];
        % Read Json file
        [a,b] = fileparts(temp_list{e_i});
        json_file = fullfile(a,[b '.json']);
        fid = fopen(json_file);
        raw = fread(fid, inf);
        str = char(raw');
        fclose(fid);
        json_read = jsondecode(str);
        
        % set parameter
        slice_time = json_read.SliceTiming;
        tr = json_read.RepetitionTime;
        mbf = json_read.MultibandAccelerationFactor;
        %% DATA
        slice_timing_job{1}.spm.temporal.st.scans{1} = spm_select('expand', temp_list(e_i)); % individual 4d images in cell str
        %% 1. nslices
        Vfirst_vol = spm_vol([temp_list{e_i} ',1']);
        num_slices = Vfirst_vol(1).dim(3);
        slice_timing_job{1}.spm.temporal.st.nslices = num_slices; % number of slices
        %% 2. tr
        slice_timing_job{1}.spm.temporal.st.tr = tr;
        %% 3. ta: acquisition time
        slice_timing_job{1}.spm.temporal.st.ta = tr - tr * mbf / num_slices; % if not multi-band, mbf = 1;
        %% 4. so: Slice order
        
        slice_timing_job{1}.spm.temporal.st.so = slice_time;
        %     if ~exist('custom_slice_timing', 'var')
        %         slice_timing_job{1}.spm.temporal.st.refslice = find(slice_time==0, 1, 'first');
        %     else
        %         slice_timing_job{1}.spm.temporal.st.refslice = find(slice_time==min(slice_time), 1, 'first');
        %     end
        if min(slice_time) >= 0 && max(slice_time) <= tr % Time-based
            slice_timing_job{1}.spm.temporal.st.refslice = min(slice_time);
        end
        slice_timing_job{1}.spm.temporal.st.prefix = 'a';
        
        %% Saving slice time correction job
        %josb_slice_timing_job = slice_timing_job{1};
        %% RUN
        spm('defaults','fmri');
        spm_jobman('initcfg');
        spm_jobman('run', slice_timing_job);
        %spm_jobman('interactive', slice_timing_job);
    end
end
%% REALIGNMENT: Motion correction
% The consensus is to do 1) estimate realigment parameter using
% before-slice timing correction images and 2) alignment using after-slice
% timing correction images
%
% By J.J

subjdat=subjdir{sub_i};
subjdat(contains(subjdat,'sbref')) = [];


for i = 1:length(subjdat) % each task
    afunc_bold_files = filenames(fullfile(subjdat{i},'afunc_*e*.nii'));
    a = [];
    b = [];
    for e_i = 1:length(afunc_bold_files) % each echo
        [a,b]=fileparts(afunc_bold_files{e_i});
        def = spm_get_defaults('realign');
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions = def.estimate;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions = def.write;
        
        % change a couple things
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % do not register to mean (twice as long)
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0; % do not mask (will set data to zero at edges!)
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 0]; % do not output mean image
        
        data = [];
        data_all = afunc_bold_files(e_i);
        matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = data_all;
        
        %% RUN
        spm('defaults','fmri');
        spm_jobman('initcfg');
        spm_jobman('run', {matlabbatch});
        
        
        %% Save realignment parameter
        
        [d, f] = fileparts(data_all{1});
        
        tempcpfile = fullfile(d, ['rp_' f '.txt']);
        temp_mvmt = textread(tempcpfile);
        
        %temp_mvmt(1:(images_per_session+kk-1),:) = [];
        
        % save plot
        create_figure('mvmt', 2, 1)
        subplot(2,1,1);
        plot(temp_mvmt(:,1:3));
        legend('x', 'y', 'z');
        
        subplot(2,1,2);
        plot(temp_mvmt(:,4:6));
        legend('pitch', 'roll', 'yaw');
        
        sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
        set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);
        drawnow;
        
        %[~,a] = fileparts(afunc_bold_files{e_i});
        [a,b]=fileparts(afunc_bold_files{e_i});
        mvmt_qcfile = fullfile(a, ['qc_mvmt_' b '.png']); % Scott added some lines to actually save the spike images
        saveas(gcf,mvmt_qcfile);
        close all;
    end
end
end