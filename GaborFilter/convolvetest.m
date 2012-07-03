I=rgb2gray(imread('lena.jpg','jpg'));

mask = [-1 -1 -1;-1 8 -1;-1 -1 -1;];

figure();
imshow(I);

J = conv2(I, mask, 'valid');

figure();
imshow(J);