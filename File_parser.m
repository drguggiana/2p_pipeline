%% File parser
%% Clean up
clearvars
close all
%% Define whether there is a stim file or not
% stim_file = 1;
%% Set up file paths

%bin path
bin_path = 'J:\Drago Guggiana Nilo\Bin_files';

%aux path
aux_path = 'J:\Drago Guggiana Nilo\Aux Data\';

%stim path
stim_path = 'J:\Drago Guggiana Nilo\Stim data\';
%% Parse aux files

%get the aux file list
aux_folders = strsplit(genpath(aux_path),';')';
separator = cell(length(aux_folders),1);
separator(:) = {'\'};
aux_folders = cellfun(@strsplit,aux_folders,separator,'UniformOutput',0);
%get a vector with the ones going into the actual files
full_paths = cellfun(@size,aux_folders,num2cell(ones(length(aux_folders),1).*2));
aux_folders = aux_folders(full_paths==6);
aux_folders = cat(1,aux_folders{:});
%allocate memory to store the actual files
aux_cell = cell(length(aux_folders),1);
%for all the remaining files
for aux = 1:length(aux_folders)
    
%     %load the contained files in the cell
%     aux_cell{aux,1} = fullfile(aux_folders{aux,:});
%     aux_cell{aux,2} = dir(aux_cell{aux,1});
%     aux_cell{aux,2} = {aux_cell{aux,2}(3:end).name}';
%     aux_cell{aux,3} = ones(size(aux_cell{aux,2})).*aux;

    %load the folder path
    folder_path = fullfile(aux_folders{aux,:});
    %load the files in that folder
    folder_files = dir(strcat(folder_path,'\*.lvd'));
%     folder_files = {folder_files(3:end).name}';
    folder_files = {folder_files.name}';
    
    %extract the different parts of the file names
    [~,aux_name,~] = cellfun(@fileparts,folder_files,'UniformOutput',0);
%     separator = cell(length(aux_name),1);
%     separator(:) = {'_'};
%     aux_parts = cellfun(@strsplit,aux_name,separator,'UniformOutput',0);
%     aux_parts = cat(1,aux_parts{:});
%     aux_experiment = cat(1,aux_parts(:,3));
%     aux_animal = cellfun(@strcat,aux_parts(:,1),separator,aux_parts(:,2),'UniformOutput',0);
    experiment_function = @(x) x(12:15);
    animal_function = @(x) x(1:11);
    aux_experiment = cellfun(experiment_function,aux_name,'UniformOutput',0);
    aux_animal = cellfun(animal_function,aux_name,'UniformOutput',0);
    
    %get the number of experiments
    num_exp = length(aux_experiment);
    
%     %get the coordinates of the individual experiments
%     [~,ia,ic] = unique(aux_experiment);
%     %get the number of experiments
%     num_exp = length(ia);
    %allocate memory to store the per experiment info
    info_perexp = cell(num_exp,7);
    %for all the experiments
    for experiment = 1:num_exp
        %load the cell with the corresponding info
        %animal name
        info_perexp{experiment,1} = aux_animal{experiment};
        %experiment
        info_perexp{experiment,2} = aux_experiment{experiment};
        %folder
        info_perexp{experiment,3} = folder_path;
        %lvd, eye1 and eye2 path
        info_perexp{experiment,4} = folder_files{experiment};
        info_perexp{experiment,5} = strcat(folder_files{experiment}(1:end-3),'eye1');
        info_perexp{experiment,6} = strcat(folder_files{experiment}(1:end-3),'eye2');
        %check if there's a text file. if so, add it
        if exist(fullfile(folder_path,strcat(folder_files{experiment}(1:end-3),'txt')),'file')==2
            info_perexp{experiment,7} = strcat(folder_files{experiment}(1:end-3),'txt');
        end
        
    end
    %load the info into the main cell
    aux_cell{aux} = info_perexp;
end

%concatenate the matrix
aux_cell = cat(1,aux_cell{:});

%create a structure containing the paths and file information
aux_struct = cell2struct(aux_cell,{'animal','experiment','folder','lvd_file','eye1_file','eye2_file','txt_file'},2);

%load the timestamp from each file
%for all the folders
for aux = 1:length(aux_struct)
    %load the timestamp from file
    [~, ~, ~, timestamp, ~, ~] = load_lvd(fullfile(aux_struct(aux).folder,aux_struct(aux).lvd_file));
    timestamp = num2str(timestamp);
    %format the time stamp properly and save
    aux_struct(aux).timestamp = strcat(timestamp(1:4),'_',timestamp(5:6),...
        '_',timestamp(7:8),'_',timestamp(9:10),'_',timestamp(11:12),'_',timestamp(13:14));
    
end
%% TODO: Add movement of the txt files if they exist
%% Matching
% %start with the bin files
% bin_folders = strsplit(genpath(bin_path),';');
% %get rid of the last one, since it's empty
% bin_folders = bin_folders(1:end-1);
% 
% %for each folder
% for folders = 1:length(bin_folders)

%allocate memory for the path structure
path_struct = struct([]);
%get the bin files in each one
%     bin_files = dir(strcat(bin_folders{folders},'\*.ini'));
bin_files = dir(strcat(bin_path,'\*525.ini'));
% %exclude the dot folders
% bin_files = bin_files(3:end);
%for all the files
for files = 1:length(bin_files)
    %read the file
    fileID = fopen(fullfile(bin_files(files).folder,bin_files(files).name),'r');
    file_content = strsplit(fscanf(fileID,'%s'),'"');
    fclose(fileID);
    file_date = file_content{4}(1:10);
    file_date = strrep(file_date,'-','_');
    file_time = file_content{4};
    file_time = strcat(file_time(1:10),'-',file_time(11:18));
    %         file_time = file_time(1:8);
    file_time = datetime(file_time,'InputFormat','yyyy-MM-dd-HH:mm:ss');
    
    %% Aux matching
    %find the coordinate of the time that most resembles the current
    %file
    [aux_times,aux_idx] = sort(datetime({aux_struct(:).timestamp},'InputFormat','yyyy_MM_dd_HH_mm_ss'));
    
    %define an array with times 30s before and after the recorded time
    time_array = [file_time - seconds(30),file_time + seconds(30)];
    
    %determine which files are within the tolerance
    tar_aux_time = aux_idx(isbetween(aux_times,time_array(1),time_array(2)));
%     tar_aux_time = aux_idx(find(file_time>aux_times,1,'last'));
    %% Determine whether there is a stim file and load it if so
    %assemble the folder path for the stim files
    stim_subpath = strcat(stim_path,file_date);

    %get the files in the stim path
    stim_files = dir(stim_subpath);

    %if the directory exists
    if ~isempty(stim_files)
        stim_files = stim_files(cat(1,stim_files(:).isdir)==0);
        %get the dates of the files
        separator = cell(length(stim_files),1);
        separator(:) = {'-'};
        stim_times = cellfun(@strsplit,{stim_files(:).name}',separator,'UniformOutput',0);
        stim_times = cat(1,stim_times{:});
        %         str_extract = @(x) x(12:19);
        %         stim_times = datetime(cellfun(str_extract,stim_times(:,2),'UniformOutput',0),'InputFormat','HH_mm_ss');
        stim_times = datetime(stim_times(:,2),'InputFormat','yyyy_MM_dd_HH_mm_ss');
        
        %define an array with times 30s before and after the recorded time
        time_array = [file_time - seconds(180),file_time + seconds(30)];

        %determine which files are within the tolerance
        tar_stim_time = isbetween(stim_times,time_array(1),time_array(2));
        
        %again, if there is no file, turn off the stim variable
        if ~isempty(tar_stim_time)

    %         %find the first date later than the bin time
    %         tar_stim_time = find(file_time>stim_times,1,'last');
            %         %load the corresponding file
            %         load(fullfile(stim_files(tar_time).folder,stim_files(tar_time).name));

            %set the stim file variable to on
            stim_file = 1;
        else
            stim_file = 0;
        end
    else
        %set the stim file variable to off
        stim_file = 0;
    end
    %% Fill in the structure containing the corresponding paths for every experiment and animal
    %so I know what to put in each solver in the next part
    
    %find the animal corresponding to this file
    path_struct(files).animal = aux_struct(tar_aux_time).animal;
    %load the date of the file
    path_struct(files).date = file_date;
    %load the experiment number
    path_struct(files).expNum = num2str(str2double(aux_struct(tar_aux_time).experiment));
    %load the path to the ini files
    path_struct(files).iniPath525 = fullfile(bin_files(files).folder,bin_files(files).name);
    %create the ini path for the 610 ini file
    path610 = bin_files(files).name;
    path610(end-6:end-4) = '610';
    path_struct(files).iniPath610 = fullfile(bin_files(files).folder,path610);

    %load the path to the imaging files
    %change the ini paths to bin paths
    bin525 = path_struct(files).iniPath525;
    bin525(end-2:end) = 'bin';
    bin610 = path_struct(files).iniPath610;
    bin610(end-2:end) = 'bin';
    path_struct(files).binPath525 = bin525;
    path_struct(files).binPath610 = bin610;
    if stim_file == 1
        %load the path to the stim file
        path_struct(files).stimPath = fullfile(stim_files(tar_stim_time).folder,stim_files(tar_stim_time).name);
    else
        path_struct(files).txtPath = fullfile(aux_struct(tar_aux_time).folder,aux_struct(tar_aux_time).txt_file);
    end
    %load the paths to the aux files
    path_struct(files).lvdPath = fullfile(aux_struct(tar_aux_time).folder,aux_struct(tar_aux_time).lvd_file);
    path_struct(files).eye1Path = fullfile(aux_struct(tar_aux_time).folder,aux_struct(tar_aux_time).eye1_file);
    path_struct(files).eye2Path = fullfile(aux_struct(tar_aux_time).folder,aux_struct(tar_aux_time).eye2_file);
%     %load the stim file time
%     path_struct(files).stimTime = stim_times(tar_stim_time);
%     %load the aux file time
%     path_struct(files).auxTime = aux_times(tar_aux_time);
%     %load the aux experiment number
    
    
end
% end
%return control to the user to check the file structure. if it looks right,
%just run the cell below to actually move the files
return
%% Create folder structure and move files

%define the path to target the files to
data_path = 'J:\Drago Guggiana Nilo\Data\';
%get the folders already there
data_list = dir(data_path);
data_list = data_list(3:end);
%ID the animals present
[animal_list,ia,ib] = unique({path_struct(:).animal});
%get the number of animals
animal_num = length(animal_list);
%get the list of fields
field_list = fieldnames(path_struct);
%remove the first 2 (animal and date)
field_list = field_list(4:end);

%for all the animals
for animals = 1:animal_num
    
    %see whether the animal folder is in there already or not
    if isempty(data_list)||sum(contains({data_list(:).name},animal_list{animals}))==0
        %if it doesn't, create the folder
        mkdir(data_path,animal_list{animals})
    end
    %assemble the animal path
    animal_path = fullfile(data_path,animal_list{animals});
    
    %get a list of the actual files that are in the same animal
    file_list = find(ib==animals);
    %for all the files from this animal on the raw data
    for files = 1:length(file_list)
        
        %get a list of the date folders inside the animal folder
        date_list = dir(animal_path);
        date_list = date_list(3:end);
        %as before, if the folder for the current date is not there, create it
        if isempty(date_list)||any(~contains({date_list(:).name},path_struct(file_list(files)).date))
            %if it doesn't, create the folder
            mkdir(animal_path,path_struct(file_list(files)).date)
        end
        
        %assemble the date path for this file
        date_path = fullfile(animal_path,path_struct(file_list(files)).date);
        %assemble the final path for the files
        destination = fullfile(date_path,path_struct(file_list(files)).expNum);
        
        %create a folder in the date path for the experiment number
        mkdir(date_path,path_struct(file_list(files)).expNum)
        
        %go through the files and transfer them to the final path
        for fields = 1:length(field_list)
            
            %if there's no stim file, skip the iteration
            if isempty(path_struct(file_list(files)).(field_list{fields}))
                continue
            end
            %get the source
            source = path_struct(file_list(files)).(field_list{fields});
            [status,msg] = movefile(source,destination);
            if status ~= 1
                error(msg)
            end
        end
    end
    
end