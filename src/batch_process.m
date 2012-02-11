clear;

%% Files to process
data_dir = ['..' filesep '3DPalm_zonly'];
file_mask = [data_dir filesep '*.zonly'];
file_list = dir(file_mask);
file_list = char({file_list.name});

%% Allocate memory to store results
num_of_horizontal_layers = 8;
num_of_samples = length(file_list);
max_depth = zeros(num_of_samples,1,'single');
hc_area = zeros(num_of_samples,num_of_horizontal_layers,'int32');

%% Prepare for parallel computing
matlabpool;

%% Extract features from each sample
tic;
parfor i = 1:num_of_samples    
    [max_depth(i,1),hc_area(i,:)] = extract_global_feature([data_dir filesep file_list(i,:)], num_of_horizontal_layers);
end
toc;

%% Close all parallel workers
matlabpool close;

%% Save result to file
save_file = 'features.mat';
save(save_file,'file_list','max_depth','hc_area');