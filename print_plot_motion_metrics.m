function print_plot_motion_metrics(cY,mY,ng,method_name,output_folder,output_file_tag)
    fprintf(method_name);
    fprintf('  '+string(mean(cY))+'  '+string(std(cY))+'  '+string(ng) +'\n');
    
    fig = figure;
    imagesc(mY); axis equal; axis tight;
    title(sprintf('Mean Image (%s)',method_name),'fontweight','bold','fontsize',14);
    saveas(fig, strcat(output_folder,'mean_img_',output_file_tag,'.png'));

end

