function [file_size, read_until] = get_file_size(fi)
   
    % go to the end of the file
    fseek(fi, 0, 'eof');
    file_size = ftell(fi);

    % with a sampling rate of 30kHz we set the threshold to 20 minutes
    max_length_tmp = 30000 * 60 * 20;

    % by the default we read in the whole file
    read_until = 0;
    if file_size > max_length_tmp
        read_until = 1;
    end
    
end
%% EOF