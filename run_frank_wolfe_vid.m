clc
clear all
close all

%% reading video file and save it as a matrix
 FileName = 'Vole1-2.mp4';
 path = fullfile('..', '..', 'DataSet', 'RGB', 'Vole', FileName);
%  x_min = 20;
%  x_max = 212;
%  y_min = 89;
%  y_max = 400;
 vid = VideoReader(path);
 % starting at specific time
 vid.CurrentTime = 200; % in seconds
 i=1;
 currAxes = axes;
 while hasFrame(vid) && i<=500
     vidFrame = readFrame(vid);
%      vidFrame = vidFrame(x_min:x_max, y_min:y_max, :);
     imArray(:,:,i) = rgb2gray(vidFrame);
     i = i+1;
     % playing back the video
%      image(vidFrame, 'Parent', currAxes);
%      pause(1/vid.FrameRate);
 end
 %% extracting the size information of the cropped video and reshaping the video sequence
orgSize = size(imArray);
imNum = orgSize(3);   % number of images in the dataset
imDim = orgSize(1:2); % resolution of each frame
rate = 1;   % down sampling rate
height = imDim(1)/rate;
width = imDim(2)/rate;
dwnSize = height*width/(rate^2);  % dimention of the downsampled frame as a vector
%%

imMat = reshape(imArray,height*width,[]);
nFrm = size(imMat,2);
X = double(imMat(:,1:nFrm));
niter = 40;
tic
[A,E,D]=frank_wolfe_vI(X,niter);
toc
 figure(1);
 for i=1:1:size(A,2)
     im_f = E(:,i);
     im_f = reshape(abs(im_f), height, width);
     im_b = A(:,i);
     im_b = reshape(im_b, height, width);
     im_o = reshape(imMat(:,i), height, width);
     subplot(3,1,1);
     imshow(im_o, []);
     subplot(3,1,2);
     imshow(im_f, []);
     subplot(3,1,3);
     imshow(im_b, []);
     pause(1/10);
 end
%%
F = uint8(reshape(abs(E), height, width, []));
B = uint8(reshape(abs(A), height, width, []));
v = VideoWriter('original.avi', 'Grayscale AVI');
v.FrameRate = 10;
% v.Colormap = 'Grayscale AVI';
open(v)
writer.FrameRate = vid.FrameRate;
writeVideo(v, imArray(:,:,1:nFrm));
close(v);

v = VideoWriter('background.avi', 'Grayscale AVI');
v.FrameRate = 10;
% v.Colormap = 'Grayscale AVI';
open(v)
writer.FrameRate = vid.FrameRate;
writeVideo(v, B);
close(v);

v = VideoWriter('forground.avi', 'Grayscale AVI');
v.FrameRate = 10;
% v.Colormap = 'Grayscale AVI';
open(v)
writer.FrameRate = vid.FrameRate;
writeVideo(v, F);
close(v);
%%
Eb=E;
E(E>0)=0;
cd 'output'
delete *.bmp
cd ..\
for i=1:1:size(E,2)
    im_f = abs(Eb(:,i));
    Mx = max(im_f);
    im_f = reshape(imbinarize(abs(im_f),Mx/3),height, width);
    im_f = imopen(im_f, strel('rectangle', [3,3]));
    im_f = imclose(im_f, strel('rectangle', [5, 5]));
    im_f = imfill(im_f, 'holes');
    im_f = uint8(im_f*255);
    FileName = strcat(num2str(i), '.bmp');
    path = fullfile('output', FileName);
%      im_f = reshape(abs(im_f), height, width);
    imwrite(im_f, path);
 end
