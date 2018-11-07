function [iframe_times,shutter_open]=get_iframe_times(imeta_data,shutter_channel,nbr_iframes)

if nargin<3
    nbr_iframes=length(imeta_data);
end

cnt=0;
go_on=1;
shutter_close=1;

while go_on
    cnt=cnt+1;
    if isempty(find(shutter_channel(shutter_close(cnt):end)>2.5,1))
        go_on=0;
    else
        shutter_open(cnt)=find(shutter_channel(shutter_close(cnt):end)>2.5,1)+shutter_close(cnt);
        shutter_close(cnt+1)=find(shutter_channel(shutter_open(cnt):end)<2.5,1)+shutter_open(cnt);
    end
end

if length(shutter_open)~=length(nbr_iframes)
    disp('Number of shutter openings does not match the eye data');
    return
end


for ind=1:length(nbr_iframes)
    iframe_times(sum(nbr_iframes(1:ind-1))+1:sum(nbr_iframes(1:ind)))=...
        imeta_data(2,sum(nbr_iframes(1:ind-1))+1:sum(nbr_iframes(1:ind)))+...
        shutter_open(ind)-imeta_data(2,sum(nbr_iframes(1:ind-1))+1);
end