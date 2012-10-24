% main function to call restaore_image()

input = 'noisy.jpg';
src = imread(input, 'jpg');

%% parameters
max_diff = 200;
weight_diff = 0.02;
iterations = 10;
covar = 100;
%%
dst = restore_image(src, covar, max_diff, weight_diff, iterations);