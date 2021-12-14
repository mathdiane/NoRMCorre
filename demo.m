clear
gcp;

%simulated video from dNMF
name = '..\dNMF\colab_demo_results\size_64_128_2_1000_30\original.tif';
%size_x_y_z_T_K; K is the number of neurons, value ranges: [0,255]
%original.tif only looks at max along z, so the tif file size is [x y T]

%original default demo data file froom NoRMCorre
%name = 'granule_love2.tif';
if ~exist(name,'file') && strcmpi(name,'granule_love2.tif')% download file if it doesn't exist in the directory
    url = 'https://www.dropbox.com/s/mjmtwn4pdgydkny/granule_love2.tif.zip?dl=1';
    filename = 'granule_love2.tif.zip';
    fprintf('downloading the file...');
    outfilename = websave(filename,url);
    fprintf('...done. Now unzipping...')
    unzip(filename);
    fprintf('done.');
elseif ~exist(name,'file')
    fprintf('[Error]filename is not granule_love2.tif AND file %s does not exists',name);
end

tic; Y = read_file(name); toc; % read the file (optional, you can also pass the path in the function instead of Y)
% Y from 'granule_love2.tif' is 64 x 128 x 4000 (frames)
% value range [38,2773]
Y = single(Y);                 % convert to single precision 
T = size(Y,ndims(Y));
Y = Y - min(Y(:)); %sclae Y s.t. min = 0
%granule_love2.tif's Y now has the value range: [0,2735]

% set parameters 
corr_bidir = 0; %don't correct bi-directional scanning, s.t. shifts_method is still fft; default is 1
% options could be changed in normcorre/normcorre_batch 
% ex: if there's bidirectional scanning and
% the options.shifts_method is 'fft', then
% options.shifts_method would be changed to 'cubic'
% options_XXX_ori: the original options before being changed

%corr_bidir should be set to True for granule_love2.tif
%because its img seems to be from bidirectional scanning when looking at
%the img through eyes
%plus the detected correction is larger (7 pixels)
%corr_bidir can be set to False for dNMF example
%because its img doesn't seem to be from bidirectional scanning when looking at
%the img through eyes
%plus the detected correction is smaller (-2 pixels)


%running time for fft is ~2x of cubic for pw-rigid 
options_rigid_ori = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200, 'correct_bidir', corr_bidir);
options_nonrigid_ori = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200, 'correct_bidir', corr_bidir);
options_nonrigid_ori_plot = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200,'plot_flag',1,'make_avi',1, 'correct_bidir', corr_bidir);
%     'plot_flag          ' % flag for plotting results in real time (default: false)
%     'make_avi           ' % flag for making movie (default: false)
% but they would be very time-consuming (10x running time)
%% perform rigid motion correction (img level translation) (not in parallel)
fprintf('---------- rigid ----------\n')
tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid_ori); toc
%corr_bidir = 1--> cubic shift
%time: 37.949212 seconds for viedo of size 64 x 128 x 4000
%(granule: Offset 7.0e-01 pixels due to bidirectional scanning detected.)
%time:  9.552002 seconds for video of size 64 x 128 x 1000 
%(dNMF: Offset -2.0e-01 pixels due to bidirectional scanning detected. )
%corr_bidir = 0, fft shift
%time:  5.703990 seconds for video of size 64 x 128 x 1000
%% now try non-rigid motion correction (in parallel)
fprintf('---------- pw-rigid(parallel) ----------\n')
% non-rigid motion is approx. by pw-rigid motion through patch level
% translation
tic; [M2,shifts2,template2,options_nonrigid2] = normcorre_batch(Y,options_nonrigid_ori); toc
% the algorithm is implemented in the function ```normcorre.m```. 
% If you have access to the parallel computing toolbox, then the function ```normcorre_batch.m``` can offer speed gains by enabling within mini-batch parallel processing.
%corr_bidir = 1--> cubic shift
%time: 82.998566 seconds for viedo of size 64 x 128 x 4000
%time: 21.795310 seconds for video of size 64 x 128 x 1000
%corr_bidir = 0, fft shift
%time: 51.962705 seconds for video of size 64 x 128 x 1000
%% try non-rigid motion correction (not in parallel)
% not parallel (normcorre) cost around 2x time than parallel(normcorre_batch)
% Diane added
fprintf('---------- pw-rigid ----------\n')
tic; [M3,shifts3,template3,options_nonrigid3] = normcorre(Y,options_nonrigid_ori); toc
%corr_bidir = 1--> cubic shift
%time: 142.769127 seconds for viedo of size 64 x 128 x 4000 
%time: 31.176072 seconds for video of size 64 x 128 x 1000
%corr_bidir = 0, fft shift
%time: 87.228900 seconds for video of size 64 x 128 x 1000
%% plot non-rigid motion correction (not in parallel)
%fprintf('---------- pw-rigid (plot)----------\n')
%tic; [M4,shifts4,template4,options_nonrigid_plot4] = normcorre(Y,options_nonrigid_ori_plot); toc
%corr_bidir = 1--> cubic shift
%time: 336.612847 seconds for video of size 64 x 128 x 1000
%corr_bidir = 0, fft shift
%time: 322.000896 seconds for video of size 64 x 128 x 1000

%% (plot) non-rigid motion correction w/deformable-trapezoidal weights(not in parallel)
fprintf('---------- pw-rigid (deformable-trapezoid function)----------\n')
tic; [M5,shifts5,template5,options_nonrigid5] = d_normcorre(Y,options_nonrigid_ori); toc
%corr_bidir = 0, fft shift
%(not plotting) time: 244.781173 seconds for video of size 64 x 128 x 1000
%(plotting) time:529.707105  seconds for video of size 64 x 128 x 1000
%% compute metrics

nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

corr1 = M3; %M1;
corr2 = M5; %M2;

shifts_c1 = shifts3; %shifts1
shifts_c2 = shifts5; %shifts2

title1 = 'trapezoid';%'rigid';
title2 = 'deformable';%'non-rigid';

[cY,mY,vY] = motion_metrics(Y,10);
[cM1,mM1,vM1] = motion_metrics(corr1,10);
[cM2,mM2,vM2] = motion_metrics(corr2,10);
T = length(cY);
%% plot metrics
figure;
    ax1 = subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,3,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean '+ string(title1)+' corrected','fontsize',14,'fontweight','bold')
    ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean '+ string(title2)+' corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data',string(title1),string(title2)); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel(string(title1)+' corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel(string(title2)+' corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2,ax3],'xy')
%% plot shifts        

shifts_r = squeeze(cat(3,shifts_c1(:).shifts));
shifts_nr = cat(ndims(shifts_c2(1).shifts)+1,shifts_c2(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_nr(:,1,:))';
shifts_y = squeeze(shifts_nr(:,2,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data',title1,title2); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

%% plot a movie with the results
% commented out to save time when running the code
% figure;
% for t = 1:1:T
%     subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     set(gca,'XTick',[],'YTick',[]);
%     drawnow;
%     pause(0.02);
% end