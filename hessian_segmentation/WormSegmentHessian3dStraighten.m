function [im_bw_out,im_smooth]=...
    WormSegmentHessian3dStraighten(im,opt_seg,im_smooth)
% this function takes as an input an image with bright blobs and segments
% it by looking at the eigenvalues of the hessian matrix. Objects can be
% further seperated using a watershed filter based on the object size. An
% optional 3rd input is a smoothed version of the image, which bypasses the
% need to do a bpass filter in this code. 

%% Initialize default parameters, all of these can also be fields in opt_seg
thresh1=.1; %initial Threshold
hthresh=-.0; %threshold for trace of hessian.
minObjSize=5; % min object size
maxObjSize=100; % max object size
valleyRatio=.75;
% watershed filter object shapes? is also the value for imhmin
watershedFilter=0; 
filterSize3=[3,3,4]; %bp filter size low f
filterSize=[3,3]; %bp filter size low f
noise=1; % bp filter hi f
pad=3; % pad to take around each sub blob
maxSplit=0; % split objects using regional maxima
minSphericity=.7; % minimum sphericity for splitting.
prefilter=0; % 0: creates im_smooth through band pass filter
gaussFilter=1;
plot_interm = 1; % plot intermediate result

% parse opt_seg to load fields
if nargin>=2
    Fnames=fieldnames(opt_seg);
    for i=1:length(Fnames)
        eval([Fnames{i} '= opt_seg.' Fnames{i} ';']);
    end
else
    
end

%cuberoot for min obj dimension
minObjDim=round(minObjSize^(1/3));
imsize=size(im);
nd = length(imsize); %number of dimension for img (2 or 3)
if nd == 3
    imsize=imsize([2,1,3]);
else
    imsize=imsize([2,1]);
end
%% if smoothed image is an input, use it, otherwise, do a bpass filter
% tries to separate neurons close by through erosion
im(im<0)=0;
if ~prefilter
    if nd==3
        im_smooth=bpass3(im,noise,filterSize3);
    else
        im_smooth=bpass(im,noise,filterSize);
    end
else
    if nargin<3
        im_smooth=im;
    end
end

if plot_interm
    figure; imagesc(im); title('original');
    figure; imagesc(im_smooth); title('smooth');
end
im_smooth=normalizeRange(im_smooth); % normalize to [0,1]

%% initial threshold
im_bw_uncleaned=im_smooth>thresh1; %convert to binary
if plot_interm
    figure; imagesc(im_bw_uncleaned); title('bw');
end
%remove small objects
im_bw=AreaFilter(im_bw_uncleaned,minObjSize,[],6);
if plot_interm
    figure; imagesc(im_bw); title('remove small objects');
end

% find connected objects
cc=bwconncomp(im_bw,6);
blobStats=regionprops(cc,'Area','BoundingBox','Centroid');



%% use hessian to find nuclei in objects
% after initial rough thresholding, loop through each segmented region and
% do a thresholding based on hessian and watersheds. This is faster than
% just shotgun analyzing the entire image.


im_bw_out=zeros(size(im_bw));
    for iblob=1:cc.NumObjects
        %crop out object with pad
        box=floor(blobStats(iblob).BoundingBox);
        if nd==3
            box(1:3)=box(1:3)-[pad,pad,pad];
            box(4:end)=box(4:end)+box(1:3)+[2*pad,2*pad,2*pad];
        else
            box(1:2)=box(1:2)-[pad,pad];
            box(3:end)=box(3:end)+box(1:2)+[2*pad,2*pad];
        end

        %don't overshoot size of image if negative, cast to 0, if larger than
        %imsize, make the edge of the box the edge of the image. 
        box(box<1)=1;
        over_edge=bsxfun(@ge,box((length(box)/2+1):end),imsize);
        if nd==3
            box([false,false,false,over_edge])=imsize(over_edge);
        else
            box([false,false,over_edge])=imsize(over_edge);
        end

        if nd==3
            sub_bw=im_bw(box(2):box(5),box(1):box(4),box(3):box(6));
            sub_im=im_smooth(box(2):box(5),box(1):box(4),box(3):box(6));
        else
            sub_bw=im_bw(box(2):box(4),box(1):box(3));
            sub_im=im_smooth(box(2):box(4),box(1):box(3));
        end
        sub_im=normalizeRange(sub_im);

%         if plot_interm && iblob==1
%             figure; imagesc(sub_im); title('one blob (smooth)');
%         end

        % smooth image and calculate hessian and eigenvalues for segmentation,
        % filter less if prefiltered
        if nd==3 
            if ~prefilter
             sub_im=smooth3(sub_im,'gaussian',2*gaussFilter+1,gaussFilter);
            else
             sub_im=bpass3(sub_im,2,filterSize3);
            end
        else
            if ~prefilter
             sub_im=imgaussfilt(sub_im,gaussFilter); %smooth2a(sub_im,2*gaussFilter+1,gaussFilter); 
             %This function smooths the data in matrixIn using a mean filter over a
             %rectangle of size (2*Nr+1)-by-(2*Nc+1)
            else
             sub_im=bpass(sub_im,2,filterSize);
            end
        end

%         if plot_interm && iblob==1
%             figure; imagesc(sub_im); title('one blob (gaussian smooth)');
%         end
        %clculate hessian matrix for each point
        H=hessianMatrix(sub_im,8);
        %Find the eigenvalues for the hessian at each point which is above
        %threshold. 
        Heig=hessianEig(H,sub_bw);
        Heig(isnan(Heig))=0;
        %find where the trace is below a threshold
        if nd==3
            Htrace=real(Heig(:,:,:,1));
        else
            Htrace=real(Heig(:,:,1));
        end
        hess_bw=Htrace<hthresh ;
%         if plot_interm && iblob==1
%             figure; imagesc(hess_bw(:,:,1)); title('Hessian < 0 (of one blob)');
%         end
        %apply area threshold
        hess_bw=AreaFilter(hess_bw,minObjSize,[],6);

        %% watershed filter shapes

        if watershedFilter
            Jd=-bwdist(~hess_bw);  %make distance map
            %Jd=smooth3(Jd,'gaussian',5,2);
            
            Jd=imhmin(Jd,watershedFilter);
            Jd(~hess_bw)=Inf;
            Jw=watershed(Jd);
            hess_bw=hess_bw.*(Jw>0);
        end

        %%
        %watershed splitting based on local maxima locations
        if maxSplit
            %find regionalmaxima and threshold around that intensity
            subImaxPnts=imregionalmax(sub_im.*hess_bw);

            %make regions around each maxima that was segmented by the hessian.
            %The regions have value of the intensity of that point times the valley
            %ratio
            max_value_im=subImaxPnts.*hess_bw.*sub_im*valleyRatio;
            dilate_kernel=true(minObjDim,minObjDim,minObjDim);
            subImax=imdilate(max_value_im,dilate_kernel);

           subImaxReg=subImaxPnts>subImax & sub_im>(2*thresh1);
            %make labelled mask
            hess_bwlabel=bwlabeln(hess_bw,6);
            %loop through labelled regions
            for iLabel=1:max(hess_bwlabel(:));
                %select specific region and the peaks in that region
                sub_hess_bw=hess_bwlabel==iLabel;
                subsubImax=subImaxReg & sub_hess_bw;

                %if more than one peak is found, watershed split them
                subsubcc=bwconncomp(subsubImax);
                if subsubcc.NumObjects>1
                    maxBW=watershed(bwdist(subsubImax));
                    hess_bw(maxBW==0 & subJm)=0;
                end
            end
        end %maxsplit
        %% after splitting, apply size filters
%         hess_bw2=hess_bw.*sub_bw;
%         hess_bw3=AreaFilter(hess_bw2,minObjSize,[],6);
%         %recurrent function that splits up objects in a BW image if they are not
%         % sufficiently sphereical or have a long axis which is much longer than the
%         % other axes. 
%         opt_rs.minObjSize = minObjSize;
%         opt_rs.maxObjSize = maxObjSize;
%         opt_rs.minSphericity = minSphericity;
%         hess_bw4=regionSplit(hess_bw3,opt_rs);
%         % remove small objects in subImage
%         hess_bw=AreaFilter(hess_bw4,minObjSize,[],6);


        %compile results of all sub images into final segmented mask.
        if nd==3
            im_bw_out(box(2):box(5),box(1):box(4),box(3):box(6))= hess_bw;
        else
            im_bw_out(box(2):box(4),box(1):box(3))= hess_bw;
        end

    end %for each blob
    
end%function
%% 
