%% Trace Combine
%script to combine the traces extracted with Suite2P with the aux and stim
%info so they can be analyzed
%% Clean up
clearvars
close all
%% Select the target folders

%define the general path
ca_path = 'J:\Drago Guggiana Nilo\Data\';

%bring up the file selection to pick target folders
tar_path = uipickfiles('FilterSpec',ca_path);
%% Loop through the folders

%get the number of folders picked
fold_num = length(tar_path);
%define the frame averaging to use (need to see if it's contained in any
%file)
frame_averaging = 2;

%define the type of experiment (gratings 0, stage 1)
exp_type = 0;

%for all the folders
for folds = 1:fold_num
    
    %identify the number of experiments in the current folder
    %get the string corresponding to the experiments
    exp_string = strsplit(tar_path{folds},'\');
    exp_string = strsplit(exp_string{end},'_');
    exp_num = length(exp_string);
    
    %% Load the ca traces
    
    %get the files in the folder
    ca_file = dir(strcat(tar_path{folds},'\suite2P\*_proc_clust.mat'));
    %load the file contents
    ca_data = load(fullfile(ca_file.folder,ca_file.name));
    ca_data = ca_data.dat;
    %for all the experiments in the analysis block
    for exps = 1:exp_num
        
        %assemble this file's root path
        single_path = strsplit(tar_path{folds},'\');
        single_path = fullfile(single_path{1:end-1},exp_string{exps});
        
        %allocate memory for both channels and ROI or neuropil (chan,ROI_neuro)
        ca_cell = cell(2,2);
        %get the matrix of ROIs ID'd as traces
        ca_cell{1,1} = ca_data.Fcell{exps}(cat(1,ca_data.stat(:).iscell)==1,:);
        ca_cell{1,2} = ca_data.FcellNeu{exps}(cat(1,ca_data.stat(:).iscell)==1,:);
        ca_cell{2,1} = ca_data.Fcell_red{exps}(cat(1,ca_data.stat(:).iscell)==1,:);
        ca_cell{2,2} = ca_data.FcellNeu_red{exps}(cat(1,ca_data.stat(:).iscell)==1,:);
        neurop_coeff = cat(1,ca_data.stat(cat(1,ca_data.stat(:).iscell)==1).neuropilCoefficient)';
        %get the number of imaging frames
        im_frames = size(ca_cell{1,1},2);
        %get the number of galvo scans, based on the frame averaging and
        %imaging frames
        num_scans_im = frame_averaging * im_frames;
        %% Load the lvd file
        %lvd_data rows: 1)visual stim, 2)ball sensor, 3)galvo, 4)shutter, 5) laser
        %power, 6) photodiode on screen
        
        %get the path to the lvd file
        lvd_file = dir(strcat(single_path,'\*.lvd'));
        %load the file contents
        [lvd_data, scanrateA, numchannels, timestamp, inputrange, status] = load_lvd(fullfile(lvd_file.folder,lvd_file.name));
        %% Load the eye camera files
        
        %get the path to the eye1 and eye2 files file
        eye1_file = dir(strcat(single_path,'\*.eye1'));
        eye2_file = dir(strcat(single_path,'\*.eye2'));

        %load the file contents
        [~,eye1_data] = load_eye_monitor_data_eye1UD(fullfile(eye1_file.folder,eye1_file.name),0);
        switch exp_type
            case 0
                [~,eye2_data] = load_eye_monitor_data_eye1UD(fullfile(eye2_file.folder,eye2_file.name),0);
            case 1
                [eye2_data] = load_eye_monitor_data_UD(fullfile(eye2_file.folder,eye2_file.name));
        end

        %get the frame times
        [iframe1_times,~] = get_iframe_times(eye1_data,lvd_data(4,:));
        [iframe2_times,~] = get_iframe_times(eye2_data,lvd_data(4,:)); 
        %% if a stage experiment, load the protocol names
        
        if exp_type == 1
            %get the file info and name
            prot_file = dir(strcat(single_path,'\*.txt'));
            %load the file and store
            fid = fopen(fullfile(prot_file.folder,prot_file.name));
            stimprot_ids = textscan(fid,'%s');
            stimprot_ids = stimprot_ids{1};
            fclose(fid);
            
            %get the stim_id
            stim_id = stimprot_ids{1};
            %add the spontaneous activity at the beginning
            stimprot_ids = cat(2,cat(1,{'Spont_act'},stimprot_ids(2:end)),num2cell(1:length(stimprot_ids))');
        end
        %% Trim the trace to include only the full frames of the mirror, starting from the end 
        %(since weird stuff happens at the beginning)

        close all
        
        %scale the galvo trace to +/- 1
        lvd_data(3,:) = lvd_data(3,:)./max(lvd_data(3,:));
        %take the derivative of the galvo trace
        diff_vector = diff(lvd_data(3,:));
        %grab the vector of sign transitions in the derivative, and pad a 0 at
        %the beginning and equate to 1 to make an index vector for the original
        %trace
        peak_vector = [0 diff_vector(1:end-1)>0&diff_vector(2:end)<0]==1;
        trough_vector = [0 diff_vector(1:end-1)<0&diff_vector(2:end)>0]==1;
        %find the positions in time of the peaks
        peak_pos = find(peak_vector);
        trough_pos = find(trough_vector);
        %get the values of the peaks from the raw vector
        peak_values = lvd_data(3,peak_vector);
        trough_values = lvd_data(3,trough_vector);
        %filter the peak values to only the ones right at the actual peak and
        %above a defined amount of separation from the neighbors
        valid_peaks = peak_values>0.75;
        %get the valid peak positions
        valid_peakpos = peak_pos(valid_peaks);
        %for the valid troughs, use the trough immediately preceding every
        %valid peak
        %allocate memory for the valid troughs
        valid_troughs = zeros(size(valid_peakpos));
        %for all the valid peak positions
        for peaks = 1:length(valid_peakpos)
            %find the closest trough position immediately behind it
            valid_troughs(peaks) = find(valid_peakpos(peaks)>trough_pos,1,'last');
        end

        %get the positions of the valid peaks in time. The size of the vector
        %should now match the number of frames in the imaging traces * the frame
        %averaging.
        peak_vector = peak_pos(valid_peaks);
        trough_vector = trough_pos(valid_troughs);

        %assemble a vector to look at the separation between frames
        sep_vector = [peak_vector(1:end-1);trough_vector(2:end)];
        sep_vector = [diff(sep_vector)>1 0]==1;
        %extend the frames with more than 1 separation to make it gapless
        %(since the relevant component is the mode of when the stimulus
        %happened, not the exact position of the galvo)
        peak_vector(sep_vector) = peak_vector(sep_vector) + 1;

%         %plot the peak selection
%         x_vector = 1:size(lvd_data,2);
%         plot(x_vector,lvd_data(3,:))
%         hold('on')
%         plot(x_vector(peak_vector),lvd_data(3,peak_vector),'r*')
%         plot(x_vector(trough_vector),lvd_data(3,trough_vector),'k*')

        %trim this vector from the end to match the number of imaging frames
        frame_vector = cat(1,trough_vector,peak_vector);
        frame_vector = frame_vector(:,end-num_scans_im+1:end);
        
        %based on the position of the stimuli, determine whether the file
        %might've been cut at the beginning or the end
        %find the delta between the trigger and the last stim start. if
        %it's less than 100 frames, flag as possible early termination
        trigger_off = find(diff(lvd_data(4,:))<-0.5,1,'last');
        last_stim = find(diff(lvd_data(1,:))<-2,1,'last');
        %also find out if the first stimulus starts before the first frame
        first_stim = find(diff(lvd_data(1,:))>2,1,'first');
        if (trigger_off-last_stim) > 100
            cut_flag = 2;  
        elseif first_stim<frame_vector(1,1)
            cut_flag = 1;
        else
            cut_flag = 0;
        end

        %now trim the lvd data from the beginning of the first frame to the end
        %of the last frame
        lvd_data = lvd_data(:,frame_vector(1,1):frame_vector(2,end));

        %get the frame intervals (gotta calculate from peak to peak)
    %     frame_intervals = diff(frame_vector,1,1);
        frame_intervals = [diff(frame_vector(1,:),1,2),frame_vector(2,end)+1-frame_vector(1,end)];
        %and the median frame time in s
        med_frame_time = median(frame_intervals)/scanrateA*frame_averaging;
        %finally, get the median sample rate
        sample_rate = 1./med_frame_time;
        %% Process stim file if it's a visual stimulation experiment 
        
        if exp_type == 0
  
            %get the path to the stim file
            stim_file = dir(strcat(single_path,'\InfoStim*.mat'));
            %load the file contents
            stim_data = load(fullfile(stim_file.folder,stim_file.name));
            %get a vector that indicates the particular stimulus protocol, stimulus
            %rep, and post stim interval for each trial (in the order in which they were
            %presented)
            
            
            %allocate memory for this vector
            trial_info = zeros(size(stim_data.Param.stimSeq,2),4);
            %load the vector of stimulus protocols used
            trial_info(:,1) = stim_data.Param.stimSeq';
            %identify the number of stimulus protocols used
            stimprot_num = length(fieldnames(stim_data.Param.StimProtocols));
            %allocate memory to store the ids of the different protocols
            stimprot_ids = cell(stimprot_num,2);
            %get the field names inside stim_data
            field_names = fieldnames(stim_data);
            %exclude Param
            stimprot_ids(:,1) = field_names(~strcmp(field_names,'Param'));
            %recalculate the number of protocols
            stimprot_num = length(stimprot_ids);
            %for all of the protocols
            for protocols = 1:stimprot_num
                %also load the id
                stimprot_ids{protocols,2} = stim_data.(stimprot_ids{protocols,1}).stim_id;
            end
            %get a vector with the sorted protocol numbers
            sorted_protocols = sort(cat(1,stimprot_ids{:,2}))';
            
            %for all of the protocols
            for protocols = sorted_protocols
                %get the field name of the corresponding protocol
                tar_name = stimprot_ids{cat(1,stimprot_ids{:,2})==protocols,1};
                
                switch tar_name
                    case 'DG'
%                         %fill in the particular stimulus (might have to use a switch, since
%                         %the names of the fields can vary depending on the stimulus. For
%                         %now I'll leave it as just seq_directions)
                        
                        %get the stimulus numbers in a unique sequence,
                        %combining all the conditions
                        tmp = stim_data.(tar_name).paramorder;
                        unique_seq = sub2ind([length(unique(tmp(1,:,1))),length(unique(tmp(1,:,2))),...
                            length(unique(tmp(1,:,3)))],tmp(:,:,1),tmp(:,:,2),tmp(:,:,3));
                        trial_info(trial_info(:,1)==protocols,2) = reshape(unique_seq',[],1);
                        
%                         trial_info(trial_info(:,1)==protocols,2) = reshape(stim_data.(tar_name).seqdirections',[],1);
                        %load the rep for each trial
                        rep_num = stim_data.(tar_name).n_reps;
%                         stim_num = stim_data.(tar_name).directions;
                        stim_num = length(unique(tmp(1,:,1)))*length(unique(tmp(1,:,2)))*...
                            length(unique(tmp(1,:,3)));
                        rep_idx = repmat(1:rep_num,stim_num,1);
                        trial_info(trial_info(:,1)==protocols,3) = rep_idx(:);
                        %load the post-stim interval
                        trial_info(trial_info(:,1)==protocols,4) = stim_data.(tar_name).poststim_time;
                    case 'RFM'
                        %get the stimulus numbers in a unique sequence,
                        %combining all the conditions
                        unique_seq = stim_data.(tar_name).stimseq;
                        trial_info(trial_info(:,1)==protocols,2) = reshape(unique_seq',[],1);
%                         unique_seq = sub2ind([length(unique(tmp(1,:,1))),length(unique(tmp(1,:,2)))],...
%                         tmp(:,:,1),tmp(:,:,2));
                        %load the rep for each trial
                        rep_num = stim_data.(tar_name).n_reps;
                        stim_num = size(unique_seq,2);
                        rep_idx = repmat(1:rep_num,stim_num,1);
                        trial_info(trial_info(:,1)==protocols,3) = rep_idx(:);
                        %load the post-stim interval
                        trial_info(trial_info(:,1)==protocols,4) = stim_data.(tar_name).interpatch_time;
                    case 'LO'
                        %get the stimulus numbers in a unique sequence,
                        %combining all the conditions
                        tmp = stim_data.(tar_name).paramorder;
                        unique_seq = sub2ind([length(unique(tmp(1,:,1))),length(unique(tmp(1,:,2)))],...
                        tmp(:,:,1),tmp(:,:,2));
                        trial_info(trial_info(:,1)==protocols,2) = reshape(unique_seq',[],1);
                        %load the rep for each trial
                        rep_num = stim_data.(tar_name).n_reps;
                        stim_num = length(unique(tmp(1,:,1)))*length(unique(tmp(1,:,2)));
                        rep_idx = repmat(1:rep_num,stim_num,1);
                        trial_info(trial_info(:,1)==protocols,3) = rep_idx(:);
                        %load the post-stim interval
                        trial_info(trial_info(:,1)==protocols,4) = stim_data.(tar_name).poststim_time;
                        
                end
            end
            %get the total number of trials
            trial_num = size(trial_info,1);
        end
        %% Process the ca traces (down to de-trended R(t))
        %get the sampling rate
        %design the filter to low pass the data (based on the sampling rate)
        lpFilt = designfilt('lowpassfir','PassbandFrequency',0.8, ...
            'StopbandFrequency',1,'StopbandAttenuation',65,'SampleRate',sample_rate);

        %allocate memory for the filtered, neuropil subtracted data
        f_cell = cell(2,1);
        %for both channels
        for chan = 1:2
            %filter the data
            filt_data = filtfilt(lpFilt,double(ca_cell{chan,1}'));
            filt_npdata = filtfilt(lpFilt,double(ca_cell{chan,2}'));

            %Subtract neuropil
            f_cell{chan} = filt_data-neurop_coeff.*filt_npdata+neurop_coeff.*median(filt_npdata,1);
        end

        %calculate R(t)
        R_t = f_cell{1}./f_cell{2};

        %detrend
        %define the length in time of the window (in s)
        window_time = 14;
        %define the anonymous function to use in the moving window
        function_handle = @(x) prctile(x,8,1);

        %perform the detrending
        R_detrend = R_t - moving_window(R_t,round(sample_rate*window_time),function_handle,0);
        %% Align each frame with the lvd info
        %get a trace of lvd_data binned (via mode) to match the imaging frames
        %generate the label vector to bin the frame times according to the
        %frame averaging
        label_vec = round((1:length(frame_intervals))/frame_averaging);
        %calculate the edges for the trace
        edge_vector = [0,cumsum(splitapply(@sum,frame_intervals,label_vec))];
        %generate the index vector for each frame
        idx_vec = discretize(1:length(lvd_data),edge_vector);
        
        %check for gaps in the idx vector. if any, interpolate
        gap_list = find(diff(idx_vec)>1)+1;
        
        if ~isempty(gap_list)
            while ~isempty(gap_list)                
                nanc = gap_list(1);
                %add an extra value in the index vector and the data vector
                %(via interpolation)
                idx_vec = cat(2,idx_vec(1:nanc-1),...
                    idx_vec(nanc-1)+1,idx_vec(nanc:end));
                %generate the interpolated data point
                interp_point = interp1(1:2,lvd_data(:,nanc-1:nanc)',1.5)';
                lvd_data = cat(2,lvd_data(:,1:nanc-1),...
                    interp_point,lvd_data(:,nanc:end));
                
%                 delta_cam1frame = diff(idx_vec)>1;
                
                gap_list = find(diff(idx_vec)>1)+1;
%                 last_nan1 = last_nan1 + 1;

            end
        end

        %if visual stimulation, get the periods with stimulation from the
        %lvd file
        if exp_type == 0
            %use the index vector to take the mode of the samples for each frame
            stim_perframe = splitapply(@mode,lvd_data(1,:),idx_vec);
        end
        %% Align each frame with the eye cam info

        %label the eye data frames with the corresponding microscope frame
        idx_eye1vec = discretize(iframe1_times,edge_vector);
        idx_eye2vec = discretize(iframe2_times,edge_vector);
        
        %depending on the type of eye file, select the relevant columns
        if size(eye1_data,1)<10
            %select the columns with the X and Y of the eye plus the pupil area
            eye1_cols = eye1_data(7:9,:);
            %also define the corresponding type of file
            cam1_type = 'eye';
        else
            %select the columns with the object position, size, and stage
            %status
            eye1_cols = eye1_data(3:end,:);
            %also define the corresponding type of file
            cam1_type = 'stage';
        end
        %do the same for the second camera
        if size(eye2_data,1)<10
           eye2_cols = eye2_data(7:9,:);
           cam2_type = 'eye';
        else
           eye2_cols = eye2_data(3:end,:);
           cam2_type = 'stage';
        end
        %check whether the last frame is missing. If so, copy the second to
        %last frame
        if max(idx_eye1vec)<im_frames
            [~,max_idx] = max(idx_eye1vec);
            idx_eye1vec = cat(2,idx_eye1vec(1:max_idx),im_frames,idx_eye1vec(max_idx+1));
            eye1_cols = cat(2,eye1_cols(:,1:max_idx),eye1_cols(:,max_idx),eye1_cols(:,max_idx+1));
        end
        if max(idx_eye2vec)<im_frames
            [~,max_idx] = max(idx_eye2vec);
            idx_eye2vec = cat(2,idx_eye2vec(1:max_idx),im_frames,idx_eye2vec(max_idx+1));
            eye2_cols = cat(2,eye2_cols(:,1:max_idx),eye2_cols(:,max_idx),eye2_cols(:,max_idx+1));
        end
        %% OFF Plot the stage data
%         if exp_type == 1
%             %Current order (though remember I'm only loading from col 3, see
%             %above)
%             %1: Average loop time
%             %2: Timer (ms)
%             %3: Center X
%             %4: Center Y
%             %5: Orientation
%             %6: Area
%             %7: Stimulus
%             %8: Rep
%             %9: Trial
%             %10: Point
%             %11: Speed
%             %12: Stepper X (from status)
%             %13: Stepper Y (from status)
% 
%             close all
%             %plot the timing signals
%             figure
%             plot(eye2_cols(5,:),'r')
%             hold('on')
%             plot(eye2_cols(6,:),'k')
%             plot(eye2_cols(7,:),'m')
%             plot(eye2_cols(8,:),'c')
%             legend({'Stimulus','Repetition','Trial','Point'})
%             xlabel('Camera Frame')
%         end
        %% Continue with the eye data
        %fill in NaNs at the beginning of the vector and the data to be able to
        %use the splitapply function
        %get the first frame present
        first1 = idx_eye1vec(1);
        first2 = idx_eye2vec(1);

        %remake the vector
        nan_idx_eye1vec = [1:first1-1,idx_eye1vec];
        nan_idx_eye2vec = [1:first2-1,idx_eye2vec];
        %also NaN any values that skip frames
        delta_cam1frame = diff(nan_idx_eye1vec)>1;
        delta_cam2frame = diff(nan_idx_eye2vec)>1;
%         nan_idx_eye1vec([0,delta_cam1frame]==1) = NaN;
%         nan_idx_eye2vec([0,delta_cam2frame]==1) = NaN;
        %extend the data vector by the same amount with NaNs
        nan_eye1_cols = [nan(size(eye1_cols,1),first1-1),eye1_cols];
        nan_eye2_cols = [nan(size(eye2_cols,1),first2-1),eye2_cols];
        
        %find out if there are isolated positions (i.e. not at the
        %beginning or the end) that have NaNs
%         inner_nan1 = find(delta_cam1frame(100:end-100))+100;
        inner_nan1 = find(delta_cam1frame(1:end))+1;

        %consider the entire trace if doing a stage experiment
        if exp_type == 0
            inner_nan2 = find(delta_cam2frame(100:end-100))+100;
        else
            inner_nan2 = find(delta_cam2frame(1:end))+1;
        end
%         inner_nan1 = find(delta_cam1frame);
%         inner_nan2 = find(delta_cam2frame);

       
        %if there are inner nans, insert interpolated frames in the data.
        %If there are more than 10 (and it's not the stage file, cause 
        %lower frame rate), throw and error though
        if length(inner_nan1)>10||(length(inner_nan2)>10 && exp_type==0)
            error('Too many skipped frames in the eye cameras')
        end
        %find where the continuous NaNs begin at the end
        last_nan1 = find(diff(isnan(nan_idx_eye1vec)),1,'last');
        if ~isempty(inner_nan1)

            while ~isempty(inner_nan1)                
                nanc = inner_nan1(1);
                %add an extra value in the index vector and the data vector
                %(via interpolation)
                nan_idx_eye1vec = cat(2,nan_idx_eye1vec(1:nanc-1),...
                    nan_idx_eye1vec(nanc-1)+1,nan_idx_eye1vec(nanc:end));
                %generate the interpolated data point
                interp_point = interp1(1:2,nan_eye1_cols(:,nanc-1:nanc)',1.5)';
                nan_eye1_cols = cat(2,nan_eye1_cols(:,1:nanc-1),...
                    interp_point,nan_eye1_cols(:,nanc:end));
                
                delta_cam1frame = diff(nan_idx_eye1vec)>1;
                
                inner_nan1 = find(delta_cam1frame(1:last_nan1))+1;
                last_nan1 = last_nan1 + 1;
%                 inner_nan1 = find(delta_cam1frame(100:end-100))+100;
%                 inner_nan1 = find(delta_cam1frame);
            end
        end
        
%         %if there were continuous stretches of NaNs that couldn't be
%         %interpolated, linearly interpolate between these frames, but throw
%         %a warning
%         if sum(isnan(nan_idx_eye1vec(1:last_nan1)))> 0
%             warning('Had to interpolate eye1 multiple frames in a row');
%             %find the stretches of NaNs
%             nan_idx = bwlabel(isnan(nan_idx_eye1vec(1:last_nan1)));
%             %get the number of stretches
%             nan_num = length(unique(nan_idx));
%             %for all the stretches
%             for stretch = 1:nan_num
%                 %get the point before the start and after the end
%                 start_stretch = find(nan_idx==stretch,1,'first')-1;
%                 end_stretch = find(nan_idx==stretch,1,'last')+1;
%                 %interpolate between them
%                 nan_idx_eye1vec = cat(2,nan_idx_eye1vec(1:start_stretch),...
%                     nan_idx_eye1vec(start_stretch)+1:nan_idx_eye1vec(end_stretch)-1,...
%                     nan_idx_eye1vec(end_stretch:end));
%                 
%                 
%             end
%             
%         end
        
        %find where the continuous NaNs begin at the end
        last_nan2 = find(diff(isnan(nan_idx_eye2vec)),1,'last');
        if ~isempty(inner_nan2)
            while ~isempty(inner_nan2)
                nanc = inner_nan2(1);
                %add an extra value in the index vector and the data vector
                %(via interpolation)
                nan_idx_eye2vec = cat(2,nan_idx_eye2vec(1:nanc-1),...
                    nan_idx_eye2vec(nanc-1)+1,nan_idx_eye2vec(nanc:end));
                %generate the interpolated data point
                interp_point = interp1(1:2,nan_eye2_cols(:,nanc-1:nanc)',1.5)';
                nan_eye2_cols = cat(2,nan_eye2_cols(:,1:nanc-1),...
                    interp_point,nan_eye2_cols(:,nanc:end));
                
                delta_cam2frame = diff(nan_idx_eye2vec)>1;
                %consider the entire trace if doing a stage experiment
                if exp_type == 0
                    inner_nan2 = find(delta_cam2frame(100:end-100))+100;
                else
                    inner_nan2 = find(delta_cam2frame(1:last_nan2))+1;
                    last_nan2 = last_nan2 + 1;
                end
            end
            
        end
        
        %Now NaN the extrema if they contain deltas larger than 1
        delta_cam1frame = diff(nan_idx_eye1vec)>1;
        delta_cam2frame = diff(nan_idx_eye2vec)>1;
        nan_idx_eye1vec([0,delta_cam1frame]==1) = NaN;
        nan_idx_eye2vec([0,delta_cam2frame]==1) = NaN;
        %trim the traces to contain the same number of scope frames
        min_frames = min([max(nan_idx_eye1vec),max(nan_idx_eye2vec)]);
        nan_idx_eye1vec(nan_idx_eye1vec>min_frames) = NaN;
        nan_idx_eye2vec(nan_idx_eye2vec>min_frames) = NaN;
        
        %define the anonymous function to put the results in a cell
        anon_mode = @(x) {mode(x,2)};
        %take the mode of the eye data to have it one by one with the microscope frames
        eye1_perframe = splitapply(anon_mode,nan_eye1_cols,nan_idx_eye1vec);
        eye2_perframe = splitapply(anon_mode,nan_eye2_cols,nan_idx_eye2vec);
        %concatenate the individual cells across eyes together in a single
        %array
        cam_data = cat(1,cat(2,eye1_perframe{:}),cat(2,eye2_perframe{:}));
        %% If a stage experiment, assemble the trial_info matrix from the just interpolated cam data
        if exp_type == 1
            %trial data columns
            %1: Stimulus
            %2: Repetition
            %3: Trial
            %4: Point
            %extract the relevant info from the cam data
            frame_data = round(cam_data(8:11,:));
            

            %get the number of trials
            trial_num = length(unique(frame_data(3,~isnan(frame_data(3,:)))));
            %get the number of protocols
            stimprot_num = size(stimprot_ids,1);

            %calculate the stim ON trace (scale to match the visual stim
            %code)
            stim_perframe = (frame_data(4,:)>0).*5;

            %fix the stimulus trace
            frame_data(1,frame_data(3,:)>0) =  frame_data(1,frame_data(3,:)>0) + 1;
            %now lift the numbers of all data to start at index 1
            frame_data([2 4],:) = frame_data([2 4],:) + 1;
            %NaN the beginning of the trace before the first stimulus
            frame_data(3,isnan(frame_data(3,:)))= 0;
            start_stim = find(diff(frame_data(3,:)),1,'first');
            frame_data(:,1:start_stim) = NaN;
            %then replace the last trial num (i.e. spont act) by that
            %number instead of 0
            frame_data(3,frame_data(3,:)==0) = trial_num;
            %finally, lift the stim numbers to 1 indexing
            frame_data(1,:) = frame_data(1,:) + 1;
            
            
            %allocate memory for a matrix with the trial information
            trial_info = zeros(trial_num,2);
            
            %for all the trials
            for trials = 1:trial_num
                %fill in the ID of the trial
                trial_info(trials,1) = mode(frame_data(1,frame_data(3,:)==trials));
                %fill in the rep
                trial_info(trials,2) = mode(frame_data(2,frame_data(3,:)==trials));
               
                
            end

%             %plot the data, now in scope frames
%             close all
%             figure
%             plot(frame_data(1,:),'r')
%             hold('on')
%             plot(frame_data(2,:),'k')
%             plot(frame_data(3,:),'m')
%             plot(frame_data(4,:),'c')
% 
%             plot(stim_perframe,'g')
%             legend({'Stimulus','Repetition','Trial','Point','StimON'})
%             xlabel('Scope Frame')
            %check the spont activity. If none present, discard this
            %protocol from the list and the experiment
            if isnan(trial_info(end,1))
                stimprot_ids = stimprot_ids(2:end,:);
                stimprot_num = stimprot_num - 1;
                trial_info = trial_info(1:end-1,:);
                trial_num = trial_num - 1;
            end
        end
        %% Calculate dRoR

        %get the stim starts
        stim_start_ori = find(diff(stim_perframe)>2);
        %and the stim ends
        stim_end_ori = find(diff(stim_perframe)<-2);
        
        %if a stage experiment and spont act, add those limits to the start
        %and end
        if exp_type == 1 && any(trial_info(:,1)==1)
            stim_start_ori = cat(2,stim_start_ori,find(frame_data(1,:)==1,1,'first'));
            stim_end_ori = cat(2,stim_end_ori,find(frame_data(1,:)==1,1,'last'));
        end
        
        %initialize stim_start and stim_end (and do nothing if they are
        %the same length)
        stim_start = stim_start_ori;
        stim_end = stim_end_ori;
        %if the experiment is incomplete (i.e. trials are lost at the
        %beginning or the end)
        if length(stim_start_ori)<length(stim_end_ori)
            %this means there is a starts missing at the beginning. Pad with
            %NaN
            stim_start = [NaN,stim_start_ori];           
            %also kill the matching end
            stim_end(1) = NaN;
        elseif length(stim_start_ori)>length(stim_end_ori)
            %this means there is an end missing at the end. Pad with NaN
            stim_end = [stim_end_ori,NaN];
            %also kill the matching start
            stim_start(end) = NaN;
        end
        
        %determine whether there are trials missing
        if length(stim_start) ~= trial_num
            %based on the cut_flag, deal with truncated files
            switch cut_flag
                case 1 %start is truncated
                    %pad the beginning of the stim_start and stim_end with
                    %NaNs
                    stim_start = cat(2,nan(1,trial_num-length(stim_start)),stim_start);
                    stim_end = cat(2,nan(1,trial_num-length(stim_end)),stim_end);
                case 2 %end is truncated
                    %pad the end of stim_start and stim_end with NaNs
                    stim_start = cat(2,stim_start,nan(1,trial_num-length(stim_start)));
                    stim_end = cat(2,stim_end,nan(1,trial_num-length(stim_end)));
            end
        end

%         %determine the number of frames to take for R0
%         r0_frames = ceil(sample_rate*r0_time);

        %determine the number of cells in the experiment
        cell_num = size(R_detrend,2);

        %initialize the structure
        dRoR = struct();
        
        if exp_type == 0
            %get a vector with the protocol id numbers
            prot_id = cat(1,stimprot_ids{:,2});

            %for all the stimulus protocols
            for protocols = 1:stimprot_num
                
                

                %load the name of the protocol
                dRoR(protocols).name = stimprot_ids{prot_id==sorted_protocols(protocols),1};
                %get the info for the trials in this protocol
                dRoR(protocols).trialInfo = trial_info(trial_info(:,1)==sorted_protocols(protocols),:);
                %get the number of trials in the protocol
                subtrial_num = size(dRoR(protocols).trialInfo,1);
                
%                 %time for R0 calculation (s)
%                 r0_time = 3;
                
                %get the post_stim time for the protocol
                r0_time = dRoR(protocols).trialInfo(1,4);
                %determine the number of frames to take for R0
                r0_frames = ceil(sample_rate*r0_time);

                %get the starts and ends corresponding to these trials
                prot_start = stim_start(trial_info(:,1)==sorted_protocols(protocols));
                prot_end = stim_end(trial_info(:,1)==sorted_protocols(protocols));
                
%                 %if the start or end are empty, skip the iteration
%                 if isempty(prot_start)||isempty(prot_end)
%                     continue
%                 end
                %define the vector of trial durations (including the r0
                %time)
                trial_vec = r0_frames + (prot_end-prot_start);
                switch dRoR(protocols).name
                    case 'LO'
                        %get the max trial duration plus the r0 frames
                        trial_duration = max(trial_vec);
                    case 'RFM'
                        %get the min trial duration plus the r0 frames at
                        %beginning and end
                        trial_duration = min(trial_vec) + r0_frames;
                        trial_vec = ones(1,subtrial_num).*trial_duration;
                    otherwise
                        %get the min trial duration plus the r0 frames
                        trial_duration = min(trial_vec);
                        trial_vec = ones(1,subtrial_num).*trial_duration;
                end

                %load the rep for each trial
                rep_num = length(unique(dRoR(protocols).trialInfo(:,3)));
                %also the number of stimuli
                stim_num = length(unique(dRoR(protocols).trialInfo(:,2)));

                %allocate memory for the trials
                trial_mat = zeros(cell_num,trial_duration,rep_num,stim_num);
                %allocate memory for the cam data (the last dimension is the
                %two cameras, each with x and y eye position and pupil area)
                trial_cam = zeros(trial_duration,rep_num,stim_num,6);

                %for each trial
                for trials = 1:subtrial_num
                    %get the current rep
                    reps = dRoR(protocols).trialInfo(trials,3);
                    %get the current stim
                    stim = dRoR(protocols).trialInfo(trials,2);
                    %determine the trial start
                    trial_start = prot_start(trials)-r0_frames;
                    %if the start is a NaN,
                    if isnan(trial_start)
                        %NaN also the cell in the trial_mat
                        trial_mat(:,:,reps,stim) = NaN;
                        trial_cam(:,reps,stim,:)= NaN;
                    else %transfer the dRoR info
                        %load the trial in the matrix, including the r0 period
                        trial_mat(:,1:trial_vec(trials),reps,stim) = permute(R_detrend(trial_start:trial_start + trial_vec(trials)-1,:),[2 1 3 4]);
                        %get the corresponding cam info
                        trial_cam(1:trial_vec(trials),reps,stim,:)= cam_data(:,trial_start:trial_start + trial_vec(trials)-1)';
                    end

                end

                %calculate dRoR and store the main structure
                R0 = median(trial_mat(:,1:r0_frames,:,:),2);
                dRoR(protocols).data = (trial_mat - R0)./R0;
                %save the camera data
                dRoR(protocols).camData = trial_cam;

        %         %get all the traces from the corresponding protocol
        %         temp_traces = R_detrend();
            end
        else
            %get a vector with the protocol id numbers
            prot_id = cat(1,stimprot_ids{:,2});

            %for all the stimulus protocols
            for protocols = 1:stimprot_num
                %load the name of the protocol
                dRoR(protocols).name = stimprot_ids{protocols,1};
                %get the info for the trials in this protocol
                dRoR(protocols).trialInfo = trial_info(trial_info(:,1)==stimprot_ids{protocols,2},:);
                %get the number of trials in the protocol
                subtrial_num = size(dRoR(protocols).trialInfo,1);

                %get the starts and ends corresponding to these trials
                prot_start = stim_start(trial_info(:,1)==stimprot_ids{protocols,2});
                prot_end = stim_end(trial_info(:,1)==stimprot_ids{protocols,2});
                %if the start or end are empty, skip the iteration
                if isempty(prot_start)||isempty(prot_end)
                    continue
                end

                %get the min trial duration plus the r0 frames
                trial_duration = r0_frames + min(prot_end-prot_start);

                %load the rep for each trial
                rep_num = length(unique(dRoR(protocols).trialInfo(:,2)));
%                 %also the number of stimuli
%                 stim_num = length(unique(dRoR(protocols).trialInfo(:,2)));

                %allocate memory for the trials
                trial_mat = zeros(cell_num,trial_duration,rep_num);
                %allocate memory for the cam data (the last dimension is the
                %two cameras, each with x and y eye position and pupil area)
                trial_cam = zeros(trial_duration,rep_num,size(cam_data,1));

                %for each trial
                for trials = 1:subtrial_num
                    %get the current rep
                    reps = dRoR(protocols).trialInfo(trials,2);
%                     %get the current stim
%                     stim = dRoR(protocols).trialInfo(trials,2);
                    %determine the trial start
                    trial_start = prot_start(trials)-r0_frames;
                    %if the start is a NaN,
                    if isnan(trial_start)
                        %NaN also the cell in the trial_mat
                        trial_mat(:,:,reps) = NaN;
                        trial_cam(:,reps,:)= NaN;
                    else %transfer the dRoR info
                        %load the trial in the matrix, including the r0 period
                        trial_mat(:,:,reps) = permute(R_detrend(trial_start:trial_start + trial_duration-1,:),[2 1 3]);
                        %get the corresponding cam info
                        trial_cam(:,reps,:)= cam_data(:,trial_start:trial_start + trial_duration-1)';
                    end

                end

                %calculate dRoR and store the main structure
                R0 = median(trial_mat(:,1:r0_frames,:),2);
                dRoR(protocols).data = (trial_mat - R0)./R0;
                %save the camera data
                dRoR(protocols).camData = trial_cam;

                %store the number of r0 frames
                dRoR(protocols).r0Frames = r0_frames;
                %store the stim_id (i.e. object)
                dRoR(protocols).stimID = stim_id;
                
            end
        end
        %% Assemble the structure of metadata only with the data from included ROIs

        %calcium/Suite2P related info
        meta_data.caData.ops = ca_data.ops;
        meta_data.caData.stat = ca_data.stat(cat(1,ca_data.stat(:).iscell)==1);
        meta_data.caData.sp = ca_data.sp{exps}(cat(1,ca_data.stat(:).iscell)==1,:);
        meta_data.caData.filename = ca_data.filename;
        %if it's a visual stim
        if exp_type == 0
            %stim info
            meta_data.stimData = stim_data;
        end
        %save the frame averaging constant used
        meta_data.frameAve = frame_averaging;
        %store the type of camera file
        meta_data.camtypes = {cam1_type,cam2_type};
        %% Save the dRoR and the stimulus sequence paired with each frame
        
        %check whether the pre-processing folder already exists in the
        %experiment folder. if not, make it
        folder_list = dir(single_path);
        folder_list = folder_list(3:end);
        folder_list = {folder_list(cat(1,folder_list(:).isdir)==1).name};
        if isempty(folder_list)||~contains(folder_list,'preProcessing')
            mkdir(single_path,'preProcessing')
        end
        %assemble the file name
        save_name = strsplit(single_path,'\');
        save_name = strcat('preProc_',save_name{end-2},'_',save_name{end-1},'_',save_name{end},'.mat');
        %assemble the save path
        save_full = strcat(single_path,'\preProcessing\',save_name);
        
        save(save_full,'stimprot_ids','dRoR','meta_data')
    end
end