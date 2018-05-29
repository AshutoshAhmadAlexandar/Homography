function [ H, inlierIndex ] = homography( imgLMatchFeatPts, imgRMatchFeatPts )

    parameters.numIterations = 150;      %the number of iterations to run
    parameters.subsetSize = 4;          %number of matches to use each iteration
    parameters.inlierDistThreshold = 10;   %the minimum distance for an inlier
    parameters.minInlierRatio = .3;     %minimum inlier ratio required to store a fitted model

    [H, inlierIndex] = ransac1(parameters, imgLMatchFeatPts, imgRMatchFeatPts, @Homography_fit, @residual_error);
    
    display('Number of inliers:');
    display(length(inlierIndex));
    display('Average residual for the inliers:')
    display(mean(residual_error(H, imgLMatchFeatPts(inlierIndex,:), imgRMatchFeatPts(inlierIndex,:))));
end