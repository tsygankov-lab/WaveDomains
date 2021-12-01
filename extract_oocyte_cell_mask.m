clear; clc;

Im = imread('Z:\Siarhei Hladyshau\Integrative_Model_of_Cell_Morphodynamics\actin_cortical_waves\starfish_oocyte\ncb3251-sv16\629.tif');
Im = double(Im>0);
Im = clean_mask(Im, 10);

windowSize = 51;
kernel = ones(windowSize) / windowSize ^ 2;
Im = conv2(single(Im), kernel, 'same');
Im = double(Im > 0.5); 

Im = imresize(Im, [512, 512], 'nearest');
Im = double(Im>0);

figure();
imagesc(Im);

imwrite(Im, 'Z:\Siarhei Hladyshau\Integrative_Model_of_Cell_Morphodynamics\actin_cortical_waves\selforganized_spiral_waves\cell_mask_oocyte\oocyte_mask.tif');

function BW_corr = clean_mask(BW, min_area)
    BW_corr = double(BW);
    
    %CC = bwconncomp(BW);
    CC = bwconncomp(BW,4);
    for i = 1:CC.NumObjects
        idx = CC.PixelIdxList{i};
        if length(idx)< min_area
            BW_corr(idx) = 0;
        end
    end
        
    CC = bwconncomp(~BW);
    for i = 1:CC.NumObjects
        idx = CC.PixelIdxList{i};
        if length(idx)< min_area
            BW_corr(idx) = 1;
        end
    end
    BW_corr = double(boolean(BW_corr));
end
