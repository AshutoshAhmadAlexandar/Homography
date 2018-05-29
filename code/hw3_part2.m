DrawDebug = true;
NormalizePts = true;
TrueMatches = true; 

% load images 
ImgL = imread('../data/part2/house1.jpg');
ImgR = imread('../data/part2/house2.jpg');

matches = load('../data/part2/house_matches.txt'); 

numberMatches = size(matches,1);

if(DrawDebug)
    
    figure; imshow([ImgL ImgR]); hold on; title('Overlay detected features (corners)');
    hold on; plot(matches(:,1),matches(:,2),'ys'); plot(matches(:,3) + size(ImgL,2), matches(:,4), 'ys'); 
    
    figure; imshow([ImgL ImgR]); hold on;
    plot(matches(:,1), matches(:,2), '+r');
    plot(matches(:,3)+size(ImgL,2), matches(:,4), '+r');
    line([matches(:,1) matches(:,3) + size(ImgL,2)]', matches(:,[2 4])', 'Color', 'r');
end

if (TrueMatches)  
    display('All matches are fitting to all');
    
    cartCoordL = matches(:,1:2);
    [numCoordinates, dimension] = size(cartCoordL);
    homoCoordL = ones(numCoordinates, dimension+1);
    homoCoordL(:,1 : dimension) = cartCoordL(:,1:dimension);
    x1 = homoCoordL;
   
    cartCoordR = matches(:,3:4);
    [numCoordinates, dimension] = size(cartCoordR);
    homoCoordR = ones(numCoordinates, dimension+1);
    homoCoordR(:,1 : dimension) = cartCoordR(:,1:dimension);
    x2 = homoCoordR;
    
    if NormalizePts        
        center = mean(x1(:,1:2)); 
        offset = eye(3);
        offset(1,3) = -center(1); 
        offset(2,3) = -center(2); 
        sX= max(abs(x1(:,1)));
        sY= max(abs(x1(:,2)));
        scale = eye(3);
        scale(1,1)=1/sX;
        scale(2,2)=1/sY;                  
        transform_1 = scale * offset;
        x1_norm = (transform_1 * x1')';
        
        center = mean(x2(:,1:2)); 
        offset = eye(3);
        offset(1,3) = -center(1); 
        offset(2,3) = -center(2); 
        sX= max(abs(x2(:,1)));
        sY= max(abs(x2(:,2)));
        scale = eye(3);
        scale(1,1)=1/sX;
        scale(2,2)=1/sY;                       
        transform_2 = scale * offset;
        x2_norm = (transform_2 * x2')';
        
        x1 = x1_norm;
        x2 = x2_norm;
    end
    
    u1 = x1(:,1);
    v1 = x1(:,2);
    u2 = x2(:,1);
    v2 = x2(:,2);
   
    %group at least 8 known matches together in a useful form
    temp = [ u2.*u1, u2.*v1, u2, v2.*u1, v2.*v1, v2, u1, v1, ones(size(matches,1), 1)];
  
    [~,~,V] = svd(temp);
    f_vec = V(:,9);
    
    F = reshape(f_vec, 3,3); %reshape the 9x1 vec into the 3x3 fund matrix
    [U, S, V] = svd(F);
    S(end) = 0;
    F = U*S*V';
    
    if NormalizePts
        F = transform_2' * F * transform_1;
    end

    
end
residuals = residual_error2(F,matches);
display(['Mean residual is: ' , num2str(mean(residuals))]);

L = (F * [matches(:,1:2) ones(numberMatches,1)]')'; 

L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); 
point_line_dist = sum(L .* [matches(:,3:4) ones(numberMatches,1)],2);
closest_point = matches(:,3:4) - L(:,1:2) .* repmat(point_line_dist, 1, 2);

pt1 = closest_point - [L(:,2) -L(:,1)] * 10; 
pt2 = closest_point + [L(:,2) -L(:,1)] * 10;

figure;
imshow(ImgR); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_point(:,1)]', [matches(:,4) closest_point(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

cam_MatrixL = load('../data/part2/house1_camera.txt');
[~, ~, V] = svd(cam_MatrixL);
cameraCenterL = V(:,end);
CoordL = cameraCenterL';
dimensionL = size(CoordL, 2) - 1;
normCoordL = bsxfun(@rdivide,CoordL,CoordL(:,end));
cartCoordL = normCoordL(:,1:dimensionL);
cameraCenterL = cartCoordL;

cam_MatrixR = load('../data/part2/house2_camera.txt');
[~, ~, V] = svd(cam_MatrixR);
cameraCenterR = V(:,end);
CoordR = cameraCenterR';
dimensionR = size(CoordR, 2) - 1;
normCoordR = bsxfun(@rdivide,CoordR,CoordR(:,end));
cartCoordR = normCoordR(:,1:dimensionR);
cameraCenterR = cartCoordR;

cartCoordL = matches(:,1:2);
[numCoordinates, dimension] = size(cartCoordL);
hCoordL = ones(numCoordinates, dimension+1);
hCoordL(:,1 : dimension) = cartCoordL(:,1:dimension);
x1 = hCoordL;
    
    
cartCoordR = matches(:,3:4);
[numCoordinates, dimension] = size(cartCoordR);
hCoordR = ones(numCoordinates, dimension+1);
hCoordR(:,1 : dimension) = cartCoordR(:,1:dimension);
x2 = hCoordR;

numberMatches = size(x1,1);
triangPoints = zeros(numberMatches, 3);
projPointsImgL = zeros(numberMatches, 2);
projPointsImgR = zeros(numberMatches, 2);

for i = 1:numberMatches
    pt1 = x1(i,:);
    pt2 = x2(i,:);
    crossProductMat1 = [  0   -pt1(3)  pt1(2); pt1(3)   0   -pt1(1); -pt1(2)  pt1(1)   0  ];
    crossProductMat2 = [  0   -pt2(3)  pt2(2); pt2(3)   0   -pt2(1); -pt2(2)  pt2(1)   0  ];    
    Eqns = [ crossProductMat1*cam_MatrixL; crossProductMat2*cam_MatrixR ];
    
    [~,~,V] = svd(Eqns);
    triangPointHomo = V(:,end)'; 
    dimensionT = size(triangPointHomo, 2) - 1;
    normCoordT = bsxfun(@rdivide,triangPointHomo,triangPointHomo(:,end));
    triangPoints(i,:) = normCoordT(:,1:dimensionT);
    
    %residual calculations
    hCoord1 = (cam_MatrixL * triangPointHomo')';
    dimension1 = size(hCoord1, 2) - 1;
    normCoord1 = bsxfun(@rdivide,hCoord1,hCoord1(:,end));
    projPointsImgL(i,:) = normCoord1(:,1:dimension1);
    
    hCoord2 = (cam_MatrixR * triangPointHomo')';
    dimension2 = size(hCoord2, 2) - 1;
    normCoord2 = bsxfun(@rdivide,hCoord2,hCoord2(:,end));
    projPointsImgL(i,:) = normCoord2(:,1:dimension2);
    
end

% plot the triangulated points and the camera centers
triangulationF(triangPoints, cameraCenterL, cameraCenterR);

distances1 = diag(dist2(matches(:,1:2), projPointsImgL));
distances2 = diag(dist2(matches(:,3:4), projPointsImgR));
display(['Mean Residual for Image 1(left side): ', num2str(mean(distances1))]);
display(['Mean Residual for Image 2(right side): ', num2str(mean(distances2))]);
