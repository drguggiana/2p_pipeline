function [imeta_info] = load_eye_monitor_data_UD(file_name)


if nargin < 1
    disp('usage: load_labview_acute(file_name[, file_format])');
    return;
end

meta_info_size=13;

%% try to open file, otherwise quit
[fi, message] = fopen(file_name, 'r', 'ieee-be');

status = 0;
% data = [];
% scanrateA = [];

if fi == -1
    disp('There was a problem reading the following file:');
    disp(file_name);
    disp(message);
    status = -1;
    return;
end

%% get file size

% meta_info = fread(fi, meta_info_size, 'double');
% if meta_info(3)==0 | rem(meta_info(3),1)~=0
%     size_x=meta_info(5);
%     size_y=meta_info(6);
% else
%     size_x=meta_info(5)-meta_info(3);
%     size_y=meta_info(6)-meta_info(4);
% end


fseek(fi, 0, 'eof');
file_size = ftell(fi);
fclose(fi);


%% open the file again - this time to read in the data
fi = fopen(file_name, 'r', 'ieee-be');

%% get all header information


% % double = 8 bytes, uint8 = 1 byte
% % nbr_frames=(file_size)/(size_x*size_y+meta_info_size*8);
% % nbr_frames=round((file_size)/(size_x*size_y+meta_info_size*8));  %change made by Alessandro
% 
% nbr_frames = floor((file_size)/(meta_info_size*8));
% 
% % if load_img
% %     idata=zeros(size_x,size_y,nbr_frames,'uint8');
% % else
% %     idata=[];
% % end
% imeta_info=zeros(meta_info_size,nbr_frames);
% 
% % go_on=1;
% for ind=1:nbr_frames
%     imeta_info(:,ind) = fread(fi, meta_info_size, 'double');
% %     if load_img
% %         idata(:,:,ind) = uint8(reshape(fread(fi, size_x*size_y, 'uint8'),size_x,size_y));
% %     else
% %         fseek(fi,size_x*size_y,'cof');
% %     end
% end
%% Fix any data misalignments (fill in gaps on tracking)

%load the file
imeta_temp = fread(fi, 'double');

%relying on the ms timer value, count the frames
nbr_frames = sum(imeta_temp>10000);
% %allocate memory for the frames
% imeta_info = zeros(meta_info_size,nbr_frames);

%get a vector with the position of said frames
ms_timer = find(imeta_temp>10000);
%find the short frames
short_frames = find(diff(ms_timer)<meta_info_size);

%for all the short frames
for frames = (short_frames-1)'
    %insert a vector with 4 zeros(nothing tracked) at the corresponding
    %position
    imeta_temp = cat(1,imeta_temp(1:frames*meta_info_size+2),[0;0;0;0],imeta_temp(frames*meta_info_size+3:end));
end

%reshape the matrix and output
imeta_info = reshape(imeta_temp,meta_info_size,nbr_frames);

%% read in the data

fclose(fi);


%% EOF
