clear
gcp;

%% set file path
output_root_folder = 'outputs\';

%simulated video from Colab's dNMF demo.ipynb
input_subfolder = 'size_64_128_2_1000_30_-80';
input_folder = strcat('dNMF_colab_demo_data\',input_subfolder,'\');
name = strcat(input_folder,'original.tif');
%size_x_y_z_T_K_S; K is the number of neurons. S related to noise std.
%Tif file value ranges: [0,255]
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

%% set parameters 
plot_mean_corr = 0; %plot mean img, correlation comparison btw different methods
plot_shift_t = 0; % plot correlation and displacement over time


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

phase_corr = 0; %default: 0
% In the case of high signal-to-noise ratio (SNR), phase
% correlation can be used by setting phase_flag to 1(logical) of
% options,when performing img registration


%running time for fft is ~2x of cubic for pw-rigid 
options_rigid_ori = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200, 'correct_bidir', corr_bidir, 'phase_flag', phase_corr);
options_nonrigid_ori = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200, 'correct_bidir', corr_bidir, 'phase_flag', phase_corr);
options_nonrigid_ori_plot = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200,'plot_flag',1,'make_avi',1, 'correct_bidir', corr_bidir, 'phase_flag', phase_corr);
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
%fprintf('---------- pw-rigid(parallel) ----------\n')
% non-rigid motion is approx. by pw-rigid motion through patch level
% translation
%tic; [M2,shifts2,template2,options_nonrigid2] = normcorre_batch(Y,options_nonrigid_ori); toc
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

%motion corrected data for comparison
corr0 = M1;
corr1 = M3; %M1;
corr2 = M5; %M2;

shifts_c0 = shifts1;
shifts_c1 = shifts3; %shifts1
shifts_c2 = shifts5; %shifts2

title0 = 'rigid';
short_title1 = 'trapezoid';
title1 = strcat('pw-rigid: ',short_title1);%'rigid';
short_title2 = 'deformable';
title2 = strcat('pw-rigid: ',short_title2);%'non-rigid';


[cY,mY,vY] = motion_metrics(Y,10);
[cM0,mM0,vM0] = motion_metrics(corr0,10);
[cM1,mM1,vM1] = motion_metrics(corr1,10);
[cM2,mM2,vM2] = motion_metrics(corr2,10);
%outputs of motion_metrics:
% cY:           correlation coefficient of each frame with the mean
% mY:           mean image
% ng:           norm of gradient of mean image

output_folder = strcat(output_root_folder, input_subfolder,'\');
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

fprintf('------------------------------\n');
fprintf(strcat(input_subfolder,'\n'));
fprintf('method  corr(mean)  (std)  crisp(mean img)\n');
print_plot_motion_metrics(cY,mY,vY,'original', output_folder, 'original');
print_plot_motion_metrics(cM0,mM0,vM0,title0, output_folder, title0);
print_plot_motion_metrics(cM1,mM1,vM1,title1, output_folder, short_title1);
print_plot_motion_metrics(cM2,mM2,vM2,title2, output_folder, short_title2);
T = length(cY);
%% plot metrics
if plot_mean_corr
    figure;
        ax1 = subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
        ax2 = subplot(2,3,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean '+ string(short_title1)+' corrected','fontsize',14,'fontweight','bold')
        ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean '+ string(short_title2)+' corrected','fontsize',14,'fontweight','bold')
        subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data',string(short_title1),string(short_title2)); title('correlation coefficients','fontsize',14,'fontweight','bold')
        subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
            xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel(string(short_title1)+' corrected','fontsize',14,'fontweight','bold');
        subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
            xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel(string(short_title2)+' corrected','fontsize',14,'fontweight','bold');
        linkaxes([ax1,ax2,ax3],'xy')
end
%% plot shifts        
if plot_shift_t
    shifts_r = squeeze(cat(3,shifts_c1(:).shifts));
    shifts_nr = cat(ndims(shifts_c2(1).shifts)+1,shifts_c2(:).shifts);
    shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
    shifts_x = squeeze(shifts_nr(:,1,:))';
    shifts_y = squeeze(shifts_nr(:,2,:))';

    patch_id = 1:size(shifts_x,2);
    str = strtrim(cellstr(int2str(patch_id.')));
    str = cellfun(@(x) ['patch # ',x],str,'un',0);

    % below are both roughly -0.7 for dNMF
    % means the correction improve more (cM1-cY or cM2-cY is larger) 
    % when img correlation is low (cY is small)
    % i.e. when there are more motion artifacts
    % corr(cY, cM1-cY)
    % corr(cY, cM2-cY)

    % below is roughly -0.3 for dNMF; 
    % didn't observe similar pattern as the chunk above
    % corr(cY, cM2 - cM1)
    figure;
        ax0 = subplot(411); plot(1:T,cY-max(cY),1:T,(cM2 - cM1)/max(cM2 - cM1)); legend('shifted raw data (max = 0)',strcat('(',short_title2,'-',short_title1,')/diff_{max}')); title(strcat('correlation diff.( ',short_title2,' - ',short_title1,')'),'fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])     
        ax1 = subplot(412); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data',short_title1,short_title2); title('correlation coefficients','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])       
        ax2 = subplot(413); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax3 = subplot(414); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
                xlabel('timestep','fontsize',14,'fontweight','bold')
        linkaxes([ax0,ax2,ax3],'x')
end
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

