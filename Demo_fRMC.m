 clc
 close all; 
 clear variables;
 %% reading video file and save it as a matrix
 FileNames = {'Video_003'};
 l = length(FileNames);
 tElapsed = nan;
 Disp = 0;
 max_niter = 10;
 gamma1 = 0.8 ;
 gamma2 = 0.8 ;
 savePathSt = fullfile('Data','result_frmc');
%  counter = 1;
 for counter = 1:1:length(FileNames)
     path = fullfile('Data', strcat(FileNames{counter},'.avi'));
     vid = VideoReader(path);
     % starting at specific time
     vid.CurrentTime = 0; % in seconds
     i=1;
     while hasFrame(vid)
         vidFrame = readFrame(vid);
    %      vidFrame = vidFrame(x_min:x_max, y_min:y_max, :);
         imArray(:,:,i) = rgb2gray(vidFrame);
         i = i+1;
     end
    %% extracting the size information of the cropped video and reshaping the video sequence
    orgSize = size(imArray);
    imNum = orgSize(3);   % number of images in the dataset
    imDim = orgSize(1:2); % resolution of each frame
    rate = 1;   % down sampling rate
    height = imDim(1)/rate;
    width = imDim(2)/rate;
    dwnSize = height*width/(rate^2);  % dimention of the downsampled frame as a vector

    %% applying the fRMC on gray scaled videos
    % making the video matrix by concatinating vectorized frames
    imMatG = reshape(imArray(:,:,:),height*width,[]);
    clear imArray
    imMatG = double(imMatG);
    frNum = size(imMatG,2);
    tic
	[A, ~ ] = InfaceExtFrankWolfe(imMatG, gamma1, gamma2, max_niter);
    tElapsed(counter) = toc;
    E = abs(A - imMatG);
    savePath = fullfile(savePathSt);
    save(fullfile(savePath,'Foreground.mat'), 'E');
    save(fullfile(savePath,'Background.mat'), 'A');
    clear A
    %% save the binary mask as video
    power =1;
    coef = 1;
    ForegEn = coef * E .^ power;
%     clear E
    % ForegEn = Foreground;
    Th = (1/5) * max(max(ForegEn));
    % thresholding
    ForegMask = ForegEn > Th;
    ForegMask = reshape(ForegMask,height, width, []);
    % morphologocal processing
    ForegMask = imopen(ForegMask, strel('rectangle', [3,3]));
    ForegMask = imclose(ForegMask, strel('rectangle', [5, 5]));
    ForegMask = imfill(ForegMask, 'holes');
    ForegMask = 255* uint8(ForegMask);
    v = VideoWriter(fullfile(savePath,'forground.avi'), 'Grayscale AVI');
    v.FrameRate = 10;
    % v.Colormap = 'Grayscale AVI';
    open(v)
    writer.FrameRate = vid.FrameRate;
    writeVideo(v, ForegMask);
    close(v);
    for j=1:1:size(ForegMask,3)
    FileName = strcat('fg_',num2str(j, '%.06i'), '.png');
    path = fullfile(savePath, FileName);
    imwrite(ForegMask(:, :, j), path);
    end
    save(fullfile(savePath,'elapse-time.txt'), 'tElapsed','-ascii');
 end