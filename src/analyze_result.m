save_file = 'features.mat';
load(save_file);

x = 2:1:20;
hist(max_depth(1:200),x);