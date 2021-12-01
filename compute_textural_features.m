function [im_entropy, im_contrast, im_correlation, im_energy, im_homogeneity, im_excited] = compute_textural_features(im_stack)
%COMPUTE_TEXTURAL_FEATURES Summary of this function goes here
%   Detailed explanation goes here


s = [size(im_stack,1), size(im_stack,2)];
NUM = size(im_stack,3);

im_entropy = zeros(1,NUM);
im_contrast = zeros(1,NUM);
im_correlation = zeros(1,NUM);
im_energy = zeros(1,NUM);
im_homogeneity = zeros(1,NUM);
im_excited = zeros(1,NUM);

im_mean = mean(im_stack(:));
im_std = std(im_stack(:));

for j = 1:NUM
    
    im = im_stack(:,:,j);
    im = (im-im_mean)/im_std;
    
    im_gc = graycoprops(graycomatrix(im));
    im_entropy(j) = entropy(im);
    im_contrast(j) = im_gc.Contrast;
    im_correlation(j) = im_gc.Correlation;
    im_energy(j) = im_gc.Energy;
    im_homogeneity(j) = im_gc.Homogeneity;
    im_excited(j) = sum(im(:) > 0)/sum(im(:) > -1000);
    
end

end

