format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[scriptpath,  filename, ext] = fileparts(full_path);
plot_count = 1;
line_width = 1.0;

cd(scriptpath);
commandwindow;

%% load the dll/so file

lib_path = '../build/Release/';
lib_name = 'vs_gen';
hfile = 'vs_gen_lib.h';

fprintf('loading library: %s\n', strcat(lib_path, lib_name, '.dll'));
fprintf('header file: %s\n', strcat('../include/', hfile));

if(~libisloaded(lib_name))
    [notfound, warnings] = loadlibrary(strcat(lib_path, lib_name, '.dll'), strcat('../include/', hfile));
end

if(~libisloaded(lib_name))
   fprintf('\nThe %s library did not load correctly!',  strcat(lib_path, lib_name,'.dll'));    
end

% show all the available functions and 
% libfunctionsview(lib_name);
% pause(1);


% initialize the generator using the file
fprintf('calling init_vs_gen_from_file function with: %s\n', '../blur_params_v23a.txt');
calllib(lib_name,'init_vs_gen_from_file', '../blur_params_v23a.txt');

% set the seed for the random generator
%calllib(lib_name,'set_vs_seed', 3851);

% void generate_vs_scene(double scale, unsigned int img_w, unsigned int img_h, unsigned char* img_f1_t,  unsigned char* img_f2_t, unsigned char* dm_t);
img_w = uint32(64);
img_h = uint32(64);
img_f1 = uint8(zeros(img_h * img_w * 3, 1));
img_f2 = uint8(zeros(img_h * img_w * 3, 1));
dm = uint8(zeros(img_h * img_w, 1));

% create the correct matlab pointers to pass into the function
img_f1_t = libpointer('uint8Ptr', img_f1);
img_f2_t = libpointer('uint8Ptr', img_f2);
dm_t = libpointer('uint8Ptr', dm);


%%
num_crops = 64;
num_iterations = 90;
num_bins = 23;
vs_scale = 0.1;
shape_scale = 0.06;

dm_hist = zeros(1, num_bins);

for idx=1:num_iterations
    for jdx=1:num_crops
        % generate the scene
        calllib(lib_name,'generate_vs_scene', vs_scale, shape_scale, img_w, img_h, img_f1_t, img_f2_t, dm_t);

        dm = reshape(dm_t.Value, [img_h, img_w])';
        dm = dm(16:47, 16:47);
        
        % get the counts of each depthmap value
        dm_hist = dm_hist + 8*histcounts(dm, num_bins);
    end
end
% deinterleave the pointers and stack to create the images
% img_f1 = cat(3, reshape(img_f1_t.Value(3:3:end), [img_h, img_w])', reshape(img_f1_t.Value(2:3:end), [img_h, img_w])', reshape(img_f1_t.Value(1:3:end), [img_h, img_w])');
% img_f2 = cat(3, reshape(img_f2_t.Value(3:3:end), [img_h, img_w])', reshape(img_f2_t.Value(2:3:end), [img_h, img_w])', reshape(img_f2_t.Value(1:3:end), [img_h, img_w])');

dm = reshape(dm_t.Value, [img_h, img_w])';

% figure; image(img_f1);
% figure; image(img_f2);
figure; imagesc(dm); colormap(gray(23));








bp = 1;

%%
fprintf('unloading library...\n');
unloadlibrary(lib_name);

fprintf('complete!\n');
