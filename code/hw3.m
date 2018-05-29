imgR = imread('../data/part1/uttower/right.jpg');
imgL = imread('../data/part1/uttower/left.jpg');

%convert the images to grayscale
%if(size(imgR,3) == 1)
    imgR = rgb2gray(imgR);
%end

%if(size(imgL,3) == 1)
    imgL = rgb2gray(imgL);
%end

%convert to double
imgR = im2double(imgR);
imgL = im2double(imgL);

radius = 1;

[featuresLeft, rowLeft, colLeft] = harris(imgL,2,0.003,radius,1);
[featuresRight, rowRight, colRight] = harris(imgR,2,0.003,radius,1);

[rl, cl] = size(imgL);
[rr, cr] = size(imgR);

%center_pL = zeros(size(colLeft));
radius = 20; %it can't be less than 3 because for homography we are using 2*numMatches x 9, and less than 3 will result in exceeds matrix error

featureDescriptionL = Describe_Features2(imgL, rowLeft, colLeft, radius);
featureDescriptionR = Describe_Features2(imgR, rowRight, colRight, radius);

%Finding the matching Features
descriptor_Distance = dist2(featureDescriptionL, featureDescriptionR);
[~,distance_id] = sort(descriptor_Distance(:), 'ascend');
number_matches = 150;
bestMatches = distance_id(1:number_matches);
[rowIdx_inDistMatrix, colIdx_inDistMatrix] = ind2sub(size(descriptor_Distance), bestMatches);
imgLFeature_idx = rowIdx_inDistMatrix;
imgRFeature_idx = colIdx_inDistMatrix;
 
match_rL = rowLeft(imgLFeature_idx);
match_cL = colLeft(imgLFeature_idx);
match_rR = rowRight(imgRFeature_idx);
match_cR = colRight(imgRFeature_idx);

% Display lines connecting the matched features
plot_r = [match_rL, match_rR];
plot_c = [match_cL, match_cR + cl];
figure; imshow([imgL imgR]); hold on; title('Mapping of top matched features');
hold on; 
plot(match_cL, match_rL,'ys');           
plot(match_cR + cl, match_rR, 'ys'); 
for i = 1:number_matches             
    plot(plot_c(i,:), plot_r(i,:));
end

imgLMatchFeatPts = [match_cL, match_rL, ones(150,1)];
imgRMatchFeatPts = [match_cR, match_rR, ones(150,1)];

[H, inlierIndices] = homography(imgLMatchFeatPts,imgRMatchFeatPts);

% Warp image
homographyTransform = maketform('projective', H);
imgLTransformed = imtransform(imgL, homographyTransform);

[h1, w1, numChannelL] = size(imgL);
[h2, w2, numChannelR] = size(imgR);

corners = [ 1 1 1;
           w1 1 1;
           w1 h1 1;
           1 h1 1];

    hCoord = (corners * H);
    dimensionM = size(hCoord, 2) - 1;
    normCoordM = bsxfun(@rdivide,hCoord,hCoord(:,end));
    warpCorners = normCoordM(:,1:dimensionM);

     minX = min( min(warpCorners(:,1)), 1);
    maxX = max( max(warpCorners(:,1)), w2);
    minY = min( min(warpCorners(:,2)), 1);
    maxY = max( max(warpCorners(:,2)), h2);
   
    xResRange = minX : maxX; 
    yResRange = minY : maxY; 

    [x,y] = meshgrid(xResRange,yResRange) ;
    Hinv = inv(H);

    warpedScaleFactor = Hinv(1,3) * x + Hinv(2,3) * y + Hinv(3,3);
    warpX = (Hinv(1,1) * x + Hinv(2,1) * y + Hinv(3,1)) ./ warpedScaleFactor ;
    warpY = (Hinv(1,2) * x + Hinv(2,2) * y + Hinv(3,2)) ./ warpedScaleFactor ;


    if numChannelL == 1
        blendLeft = interp2( im2double(imgL), warpX, warpY, 'cubic') ;
        blendRight = interp2( im2double(imgR), x, y, 'cubic') ;
    end
  
    blendWeight = ~isnan(blendLeft) + ~isnan(blendRight) ;
    blendLeft(isnan(blendLeft)) = 0 ;
    blendRight(isnan(blendRight)) = 0 ;
    stitchedImg = (blendLeft + blendRight) ./ blendWeight ;
    
figure, imshow(stitchedImg);
title('Alignment by homography');