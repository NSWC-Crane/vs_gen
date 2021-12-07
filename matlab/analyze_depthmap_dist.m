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

if(~libisloaded(lib_name))
    [notfound, warnings] = loadlibrary(strcat(lib_path,lib_name,'.dll'), strcat(startpath,'\include\',header_file));
end

if(~libisloaded(lib_name))
    fprintf('\nThe %s library did not load correctly!',  strcat(lib_path,lib_name,'.dll'));    
    return;
end

% initialize the generator using the file
calllib(lib_name,'init_vs_gen_from_file',fullfile(data_path, data_file));

% number of images
N = 2000;

% image size
img_w = 64;
img_h = 64;

img_f1 = uint8(zeros(img_h * img_w * 3, 1));
img_f2 = uint8(zeros(img_h * img_w * 3, 1));
dm = uint8(zeros(img_h * img_w, 1));

% create the correct matlab pointers to pass into the function
img_f1_t = libpointer('uint8Ptr', img_f1);
img_f2_t = libpointer('uint8Ptr', img_f2);
dm_t = libpointer('uint8Ptr', dm);


% void get_vs_minmax(unsigned short* min_dm_value, unsigned short* max_dm_value);
min_dm_value_t = libpointer('uint16Ptr', 0);
max_dm_value_t = libpointer('uint16Ptr', 0);
calllib(lib_name,'get_vs_minmax', min_dm_value_t, max_dm_value_t);

min_dm_value = double(min_dm_value_t.Value);
max_dm_value = double(max_dm_value_t.Value);


%% grab some data
dm_hist = zeros(max_dm_value + 1,1);

shape_scale = 0.057;

fprintf('Starting Scene Generation ...\n');

for idx=1:N
    
    % generate the scene
    calllib(lib_name,'generate_vs_scene', 0.1, shape_scale, img_w, img_h, img_f1_t, img_f2_t, dm_t);

    % get the depthmap images
    dm = double(dm_t.Value);
    
    for jdx=1:(max_dm_value + 1)
        dm_hist(jdx,1) = dm_hist(jdx,1) + sum(dm==(jdx-1));
    end 
    
    fprintf('.');
    if(mod(idx, 100) == 0)
        fprintf('\n');
    end
        
end

% % deinterleave the pointers and stack to create the images
% img_f1 = cat(3, reshape(img_f1_t.Value(3:3:end), [img_h, img_w])', reshape(img_f1_t.Value(2:3:end), [img_h, img_w])', reshape(img_f1_t.Value(1:3:end), [img_h, img_w])');
% img_f2 = cat(3, reshape(img_f2_t.Value(3:3:end), [img_h, img_w])', reshape(img_f2_t.Value(2:3:end), [img_h, img_w])', reshape(img_f2_t.Value(1:3:end), [img_h, img_w])');
% 
% dm = reshape(dm_t.Value, [img_h, img_w])';
% 
% figure(plot_num); 
% image(img_f1);
% axis off;
% plot_num = plot_num + 1;
% 
% figure(plot_num); 
% image(img_f2);
% axis off;
% plot_num = plot_num + 1;
% 
% figure(plot_num); 
% imagesc(dm); 
% colormap(gray((max_dm_value + 1)  - min_dm_value));
% axis off;
% plot_num = plot_num + 1;

fprintf('\nComplete!\n\n');

% plot the hist of the dm values

hist_bins = min_dm_value:1:max_dm_value;

figure(plot_num)
set(gcf,'position',([100,100,1300,600]),'color','w')
hold on
box on
grid on

b1 = bar(hist_bins(1:end), dm_hist);
set(gca,'fontweight','bold','FontSize',13);

% X-Axis
xlim([hist_bins(1)-1, hist_bins(end)+1]);
xticks([hist_bins(1):1:hist_bins(end)]);
%xtickangle(90);
xlabel(strcat('Depth Map Value'),'fontweight','bold')

% Y-Axis
ytickformat('%2.1f');
ylabel('Depth Map Ratio','fontweight','bold');

b1(1).FaceColor = 'b';

title(strcat('Ground Truth Distribution:',32,num2str(shape_scale,'%1.3f')), 'fontweight','bold','FontSize',16);
lgd = legend([b1(1)],'Training Data', 'location','southoutside', 'orientation', 'horizontal');

ax = gca;
ax.Position = [0.05 0.16 0.90 0.79];

%print(plot_num, '-dpng', fullfile(save_dir,'dm_combined.png'));

plot_num = plot_num + 1;

return;
%%
unloadlibrary(lib_name);