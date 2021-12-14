function [im_bw_out,im_centroid]=...
    SegmentByFiltering(im,opt_seg,im_bp)
% this function takes as an input an image with bright blobs and segments
% it by looking at the eigenvalues of the hessian matrix. Objects can be
% further seperated using a watershed filter based on the object size. An
% optional 3rd input is a smoothed version of the image, which bypasses the
% need to do a bpass filter in this code. 

%% Initialize default parameters, all of these can also be fields in opt_seg
thresh1 =.1; %initial Threshold for binarization
thresh2 = .1; %threshold for last binarization
minObjSize=5; % min object size
filterSize3=[3,3,4]; %bp filter size low f
filterSize=[3,3]; %bp filter size low f
noise=1; % bp filter hi f
prefilter=0; % 0: creates im_bp through band pass filter
gaussFilter=1; 
plot_interm = 0; % plot intermediate result

% parse opt_seg to load fields
if nargin>=2
    Fnames=fieldnames(opt_seg);
    for i=1:length(Fnames)
        eval([Fnames{i} '= opt_seg.' Fnames{i} ';']);
    end
else
    
end

%cuberoot for min obj dimension

imsize=size(im);
nd = length(imsize); %number of dimension for img (2 or 3)

%% if smoothed image is an input, use it, otherwise, do a bpass filter
% tries to separate neurons close by through erosion
im(im<0)=0;
if ~prefilter
    if nd==3
        im_bp=bpass3(im,noise,filterSize3);
    else
        im_bp=bpass(im,noise,filterSize);
    end
else
    if nargin<3
        im_bp=im;
    end
end

if plot_interm
    figure; imagesc(im); title('original');
    figure; imagesc(im_bp); title('band pass filtered');
end
im_bp=normalizeRange(im_bp); % normalize to [0,1]

%% initial threshold
im_bw_uncleaned=im_bp>thresh1; %convert to binary
if plot_interm
    figure; imagesc(im_bw_uncleaned); title('bw');
end
%remove small objects
im_bw=AreaFilter(im_bw_uncleaned,minObjSize,[],6);
if plot_interm
    figure; imagesc(im_bw); title('remove small objects');
end

im_centroid=uint16(normalizeRange(im_bw)); % normalize to [0,1]

% % find connected objects and centroids
% cc=bwconncomp(im_centroid,6);
% blobStats=regionprops(cc,'Area','BoundingBox','Centroid');

%% smooth and binarized
%smooth original img
if nd==3
   	im_sm = smooth3(im,'gaussian',2*gaussFilter+1,gaussFilter);
else
    im_sm = imgaussfilt(im,gaussFilter);
end

im_sm=normalizeRange(im_sm); % normalize to [0,1]

if plot_interm
    figure; imagesc(im_sm); title('smooth');
end

%binarize it
im_bw_out = im_sm>thresh2;



    
end%function
%% 
