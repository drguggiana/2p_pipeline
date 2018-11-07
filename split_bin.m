%% Split the binary file into a series of tif files for processing with Suite2P
%% Clean up
clearvars
close all
%% OFF Load aux file

% %select the files to convert
% tar_path = uipickfiles('FilterSpec','J:\Drago Guggiana Nilo\Data\*.lvd');
% 
% %load the file
% [data, scanrateA, numchannels, timestamp, inputrange, status] = load_lvd(tar_path{1});
%% OFF Plot lvd data
% close all
% %get the number of fields
% field_num = size(data,1);
% figure
% %for all the fields
% for fields = 1:field_num
%     subplot(ceil(sqrt(field_num)),round(sqrt(field_num)),fields)
%     
%     plot(data(fields,:))
% end
% 
% %calculate the number of stimuli
% stim_num = sum(diff(data(1,:))>2.5);
%% Load bin file

%set the number of frames per file
frames_perfile = 200;

%select the experiment folders to convert
tar_path = uipickfiles('FilterSpec','J:\Drago Guggiana Nilo\Data\');
%% Write the tiff files

%get the number of files
exp_num = length(tar_path);

%for all of the files
for experiment = 1:exp_num
    %get the binary file paths from this experiment
    file_list = dir(fullfile(tar_path{experiment},'*.bin'));
    
    
    %define the green and red paths
    g_path = fullfile(file_list(1).folder,file_list(1).name);
    r_path = fullfile(file_list(2).folder,file_list(2).name);
    
    %define the save path
    save_path = tar_path{experiment};
    
    [~,binname] = fileparts(tar_path{experiment});
    reverseStr = '';
    
    %load the green file info (generic)
    finfo=dir(g_path);
    %open the files
    fg = fopen(g_path,'r');
    fr = fopen(r_path,'r');
    %get the x and y resolution of the data
    x_res=fread(fg,1,'int16=>double');
    y_res=fread(fg,1,'int16=>double');
    %get the number of images in the file, times the number of channels
    nbr_images=round(finfo.bytes/x_res/y_res/2);
    %determine the number of files to write
    file_num = ceil(nbr_images/frames_perfile);
    
    %set up a frame counter
    frame_count = 1;
    
    %for all the frames
    for tifffiles = 1:file_num
        %assemble the save file name
%         [~,save_name] = fileparts(tar_path{experiment});
        save_full = strcat(save_path,'\',num2str(tifffiles,'%05u'),'.tif');
        
        %create a tiff object
        tiff_obj = Tiff(save_full,'w8');
        
        %set the object properties
        tiff_props.Photometric = Tiff.Photometric.MinIsBlack;
        tiff_props.Compression = Tiff.Compression.None;
        tiff_props.SampleFormat = Tiff.SampleFormat.Int;
        tiff_props.ImageLength = x_res;
        tiff_props.ImageWidth = y_res;
        tiff_props.BitsPerSample = 16;
        tiff_props.SamplesPerPixel = 1;
        tiff_props.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        
        for frames = 1:frames_perfile% (mod(frame_count,frames_perfile)~=0)&&frame_count< nbr_images+1
            %read the data
            g_data=reshape(fread(fg,y_res*x_res,'int16=>int16'),y_res,x_res)';
            r_data=reshape(fread(fr,y_res*x_res,'int16=>int16'),y_res,x_res)';
            % disp frame loading in command window
            msg = sprintf('Loading %s (x=%1.0fy=%1.0ft=%1.0f): %1.0f\n',binname,x_res,y_res,nbr_images,frame_count);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            
            %pass the properties to the object
            tiff_obj.setTag(tiff_props);
           
            %write the green file
            tiff_obj.write(g_data);
            
            
            %create a new directory
            tiff_obj.writeDirectory;
            
            tiff_obj.setTag(tiff_props);
            %write the red file
            tiff_obj.write(r_data);
            
            
            %if nbr_images is reached, kill the loop
            if frame_count == nbr_images
                break
            end
            
            %if it's not the last frame, add another IFD
            if frame_count ~= nbr_images && frames ~= frames_perfile
                
                %create a new directory
                tiff_obj.writeDirectory;
                
                tiff_obj.setTag(tiff_props);
            end
            
            %update the frame counter
            frame_count = frame_count + 1;
        end
        %close the object
        tiff_obj.close();
    end
    %close the file
    fclose(fg);
    fclose(fr);
    
    
end