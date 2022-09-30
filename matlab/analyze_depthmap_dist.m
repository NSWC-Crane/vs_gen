format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

commandwindow;

%% select input file

startpath = fileparts(startpath);
file_filter = {'*.txt','Text Files';'*.*','All Files' };

[data_file, data_path] = uigetfile(file_filter, 'Select Scene Generation Parameters File', startpath);
if(data_path == 0)
    return;
end

 
%% load the dll/so file

lib_path = strcat(startpath,'\build\Release\');
lib_name = 'vs_gen';
header_file = 'vs_gen_lib.h';

% p = parpool('AttachedFiles',{strcat(startpath,'\include\',header_file)}); 

% if(~libisloaded(lib_name))
% %     spmd
%         [notfound, warnings] = loadlibrary(strcat(lib_path,lib_name,'.dll'), strcat(startpath,'\include\',header_file));
% %     end
% end
% 
% if(~libisloaded(lib_name))
%     fprintf('\nThe %s library did not load correctly!\n',  strcat(lib_path,lib_name,'.dll'));    
%     return;
% end
% 
% % initialize the generator using the file
% calllib(lib_name,'init_vs_gen_from_file',fullfile(data_path, data_file));

% image size
img_w = 64;
img_h = 64;

% img_f1 = uint8(zeros(img_h * img_w * 3, 1));
% img_f2 = uint8(zeros(img_h * img_w * 3, 1));
% dm = uint8(zeros(img_h * img_w, 1));

% % create the correct matlab pointers to pass into the function
% img_f1_t = libpointer('uint8Ptr', img_f1);
% img_f2_t = libpointer('uint8Ptr', img_f2);
% dm_t = libpointer('uint8Ptr', dm);
% 
% 
% % void get_vs_minmax(unsigned short* min_dm_value, unsigned short* max_dm_value);
% min_dm_value_t = libpointer('uint16Ptr', 0);
% max_dm_value_t = libpointer('uint16Ptr', 0);
% calllib(lib_name,'get_vs_minmax', min_dm_value_t, max_dm_value_t);

min_dm_value = 0; %double(min_dm_value_t.Value);
max_dm_value = 22; %double(max_dm_value_t.Value);


%% grab some data
% number of images
N = 90000*64;

dm_hist = zeros(N, max_dm_value + 1);
dm_hist2 = zeros(N, max_dm_value + 1);
dm_hist3 = zeros(N, max_dm_value + 1);


% y = 0.00000000261125x3 - 0.00000239054362x2 + 0.00085286458333x + 0.01452380952383
%  64 x  64: shape_scale = 0.060; good
% 128 x 128: shape_scale = 0.090; good
% 256 x 256: shape_scale = 0.120; good 
% 512 x 512: shape_scale = 0.175; good
% shape_scale = 0.00000000261125*(img_w*img_w*img_w) - 0.00000239054362*(img_w*img_w) + 0.00085286458333*img_w + 0.01452380952383;
shape_scale = 0.06;

num_bins = 23;

tic;
fprintf('Starting Scene Generation ...\n');

for idx=1:N
    
    fprintf('Iteration: %04d\n', idx);
    
    if(~libisloaded(lib_name))
    %     spmd
%             [notfound, warnings] = loadlibrary(strcat(lib_path,lib_name,'.dll'), strcat(startpath,'\include\',header_file), 'mfilename', strcat(startpath,'\matlab\','vs_gen_prototype_file'));
            [notfound, warnings] = loadlibrary(strcat(lib_path,lib_name,'.dll'), @vs_gen_prototype_file);
    %     end
    end

    if(~libisloaded(lib_name))
        fprintf('\nThe %s library did not load correctly!\n',  strcat(lib_path,lib_name,'.dll'));    
        %return;
    end

    % initialize the generator using the file
    calllib(lib_name,'init_vs_gen_from_file',fullfile(data_path, data_file));

    img_f1 = uint8(zeros(img_h * img_w * 3, 1));
    img_f2 = uint8(zeros(img_h * img_w * 3, 1));
    dm = uint8(zeros(img_h * img_w, 1));

    % create the correct matlab pointers to pass into the function
    img_f1_t = libpointer('uint8Ptr', img_f1);
    img_f2_t = libpointer('uint8Ptr', img_f2);
    dm_t = libpointer('uint8Ptr', dm);


%     % void get_vs_minmax(unsigned short* min_dm_value, unsigned short* max_dm_value);
%     min_dm_value_t = libpointer('uint16Ptr', 0);
%     max_dm_value_t = libpointer('uint16Ptr', 0);
%     calllib(lib_name,'get_vs_minmax', min_dm_value_t, max_dm_value_t);


    %if(~libisloaded(lib_name))
    %    loadlibrary(strcat(lib_path,lib_name,'.dll'), strcat(startpath,'\include\',header_file));
    %    calllib(lib_name,'init_vs_gen_from_file',fullfile(data_path, data_file));
    % generate the scene
    %calllib(lib_name,'generate_vs_scene', 0.1, shape_scale, img_w, img_h, img_f1_t, img_f2_t, dm_t);        
    %else
    % generate the scene
    calllib(lib_name,'generate_vs_scene', 0.1, shape_scale, img_w, img_h, img_f1_t, img_f2_t, dm_t);
    %end   
        
    % get the depthmap images
%     dm = double(dm_t.Value);
    dm = reshape(dm_t.Value, [img_h, img_w])';

    dm = dm(16:47, 16:47);
    dm = dm(:);
    
    dm_tmp = zeros(1, max_dm_value + 1);
    dm_tmp2 = zeros(1, max_dm_value + 1);
    dm_tmp3 = zeros(1, max_dm_value + 1);

    for jdx=1:(max_dm_value + 1)
        dm_tmp(1, jdx) = dm_tmp(1, jdx) + sum(dm==(jdx-1));
    end
    
    dm_hist(idx, :) = dm_tmp;
%     
%     dm_hist2(idx, :) = histcounts(dm, num_bins);
%     
%     
%     for jdx=1:numel(dm)
%         dm_tmp3(1, dm(jdx)+1) = dm_tmp3(1, dm(jdx)+1) + 1;      
%     end
%     dm_hist3(idx, :) = dm_tmp3;
    
%     fprintf('.');
%     if(mod(idx, 100) == 0)
%         fprintf('\n');
%     end
        
end

toc
%%
dm_hist_sum = 8*sum(dm_hist, 1);
dm_hist2_sum = 8*sum(dm_hist2, 1);
dm_hist3_sum = 8*sum(dm_hist3, 1);


if(false)
    % deinterleave the pointers and stack to create the images
    img_f1 = cat(3, reshape(img_f1_t.Value(3:3:end), [img_h, img_w])', reshape(img_f1_t.Value(2:3:end), [img_h, img_w])', reshape(img_f1_t.Value(1:3:end), [img_h, img_w])');
    img_f2 = cat(3, reshape(img_f2_t.Value(3:3:end), [img_h, img_w])', reshape(img_f2_t.Value(2:3:end), [img_h, img_w])', reshape(img_f2_t.Value(1:3:end), [img_h, img_w])');
    
    dm = reshape(dm_t.Value, [img_h, img_w])';
    
    figure(plot_num); 
    image(img_f1);
    axis off;
    plot_num = plot_num + 1;
    
    figure(plot_num); 
    image(img_f2);
    axis off;
    plot_num = plot_num + 1;
    
    figure(plot_num); 
    imagesc(dm); 
    colormap(gray((max_dm_value + 1)  - min_dm_value));
    axis off;
    plot_num = plot_num + 1;
end

fprintf('\nComplete!\n\n');
%%
% plot the hist of the dm values

hist_bins = min_dm_value:1:max_dm_value;

figure(plot_num)
set(gcf,'position',([100,100,1200,600]),'color','w')
hold on
box on
grid on

% b1 = bar(hist_bins(1:end), [dm_hist_sum; dm_hist2_sum; dm_hist3_sum]);
b1 = bar(hist_bins(1:end), dm_hist_sum);

set(gca,'fontweight','bold','FontSize',13);

% X-Axis
xlim([hist_bins(1)-1, hist_bins(end)+1]);
xticks([hist_bins(1):1:hist_bins(end)]);
%xtickangle(90);
xlabel(strcat('Depth Map Value'),'fontweight','bold')

% Y-Axis
ytickformat('%2.1f');
ylabel('Depth Map Count','fontweight','bold');

b1(1).FaceColor = 'b';
b1(2).FaceColor = 'r';

title(strcat('Ground Truth Distribution - Shape Scale =',32,num2str(shape_scale,'%1.3f')), 'fontweight','bold','FontSize',15);
% lgd = legend([b1(1), b1(2)],'Training Data', 'hist2', 'location','southoutside', 'orientation', 'horizontal');
lgd = legend([b1(1)],'Training Data', 'location','southoutside', 'orientation', 'horizontal');

ax = gca;
ax.Position = [0.07 0.16 0.90 0.79];

%print(plot_num, '-dpng', fullfile(save_dir,'dm_combined.png'));

plot_num = plot_num + 1;

return;
%% unload the library from memory after we're done with it
fprintf('Unloading %s\n', lib_name);
unloadlibrary(lib_name);

