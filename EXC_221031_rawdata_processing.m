
basedir = '/Users/WIWIFH/Dropbox/Projects/EXCOOL/resources';
datdir = fullfile(basedir,'data','task');

%% Find matlist
matlist = []; subjnames = [];
k=1; % EXC_task_Sub-SJY_run01.mat: SUBJECT001
subjnames{k} = 'SJY';
matlist{k} = filenames(fullfile(datdir, 'EXC_task_Sub-SJY*.mat'));
k=k+1;
%EXC_task_Sub-EX002_run01.mat
for sub_i = 1:10
    temp_files = filenames(fullfile(datdir, sprintf('EXC_task_Sub-EX%03d_run*.mat',sub_i)));    
    if ~contains(temp_files,'no matches')
        matlist{k}=temp_files;
        subjnames{k,1} = sprintf('EX%03d',sub_i);
        k=k+1;
        
    end
end
% [~,b] = fileparts(temp_files{end})

%% Load
datall = []; 
tbl = []; 
for i = 1:length(matlist)
    temp_dat = []; 
    for ii=1:length(matlist{i})
        temp_dat = load(matlist{i}{ii});
        %temp_dat.dat.ts.table
        
        [ratings, conds,vals,onsets_cue,onsets_ratings,durations_cue,durations_ratings]=deal([]);
        words = {};
        for trial_i=1:length(temp_dat.dat.dat)
            if length(temp_dat.dat.dat{trial_i}.norms_ratings) < 60 
                ratings(trial_i,:) = NaN;
            else
                ratings(trial_i,:) = mean(temp_dat.dat.dat{trial_i}.norms_ratings(end-50:end))-0.5; % ratings 
            end
            conds(trial_i,:) = find(contains(temp_dat.dat.ts.target_names,temp_dat.dat.dat{trial_i}.target)); %condition 
            % => 1 = self; 2=friends; 3=cerebrities
            
            vals(trial_i,:) = temp_dat.dat.dat{trial_i}.valence;  % word valence 
            words{trial_i,1} = temp_dat.dat.dat{trial_i}.words;
            
            
            % onsets and durations 
            onsets_cue(trial_i,:) = temp_dat.dat.dat{trial_i}.ITI_EndTime - temp_dat.dat.runscan_starttime;% cue (3secs)
            onsets_ratings(trial_i,:) = temp_dat.dat.dat{trial_i}.ISI_EndTime - temp_dat.dat.runscan_starttime;% trait and ratings (7 secs)
            durations_cue(trial_i,:) = 3;
            durations_ratings(trial_i,:) = 7;                        
        end
        
        %array2table        - Convert homogeneous array to table.
        %cell2table         - Convert cell array to table.
        %struct2table       
        
        %temp_run_tbl=array2table(1:length(temp_dat.dat.dat),'TrialN',conds,'Tcondition', vals,'Wvalence', ratings,'ratings');
        temp_run_tbl=array2table([ones(length(temp_dat.dat.dat),1).*ii, (1:length(temp_dat.dat.dat))',conds, vals, ratings,onsets_cue, durations_cue, onsets_ratings, durations_ratings],'VariableNames',{'RunN','TrialN','Tcondition','Wvalence','ratings','onsets_cue','durations_cue','onsets_ratings','durations_ratings'});        
        datall{i}{ii} = temp_run_tbl;
        [~,b] = fileparts(matlist{i}{ii});
        savedir = fullfile('/Users/WIWIFH/Dropbox/Projects/EXCOOL/data',[b '.csv']);
        writetable(temp_run_tbl,savedir);
    end
end
