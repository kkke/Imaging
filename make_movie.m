%% load tiff files
filename = 'test1.tif'
INFO=imfinfo(filename);
j=length(INFO);
x=INFO(1).Width;
y=INFO(1).Height;
data=zeros(y,x,j,'int16');
% data=zeros(y,x,j);
handle=waitbar(0,'Loading image');
for i=1:j
    data(:,:,i)=imread(filename,i);
    waitbar(i/j,handle)
end
close(handle);
disp(filename);
Iminfo.FileName=filename;
%% show original movie
data = double(data);
figure;
for i = 1:size(data,3)
    imagesc(data(:,:,i))
    caxis([500,10000])
    colormap('jet')
    pause(0.1)
end

%% show df/f movie
Im_avg = mean(data(:,:,1:12),3);
Im_dm  = data-repmat(Im_avg,1,1,size(data,3));
Im_dm_m = Im_dm./repmat(Im_avg,1,1,size(data,3));

figure;
for i = 1:size(Im_dm_m,3)
    imshow(Im_dm_m(:,:,i),[-3,3])
    
    pause(0.1)
end