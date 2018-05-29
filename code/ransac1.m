function [ bestFit, inlierIndex ] = ransac1( parameters, x, y, Homography_fit, residual_Error )

    [number_matches, ~] = size(x);
    numInliersEachIteration = zeros(parameters.numIterations,1);
    storedModels = {};
    
    for i = 1 : parameters.numIterations;
        subsetIndices = randsample(number_matches, parameters.subsetSize);
        x_subset = x(subsetIndices, :);
        y_subset = y(subsetIndices, :);
            
        %fit a model to that subset
        model = Homography_fit(x_subset, y_subset);
        
        residualErrors = residual_Error(model, x, y);
        
        inlierIndex = find(residualErrors < parameters.inlierDistThreshold);      

        %record the number of inliers
        numInliersEachIteration(i) = length(inlierIndex);
        
        currentInlierRatio = numInliersEachIteration(i)/number_matches;
        if currentInlierRatio >=  parameters.minInlierRatio
            x_inliers = x(inlierIndex, :);
            y_inliers = y(inlierIndex, :);
            storedModels{i} = Homography_fit(x_inliers, y_inliers);
        end
    end
    
    bestIteration = find(numInliersEachIteration == max(numInliersEachIteration));
    bestIteration = bestIteration(1); 
    bestFit = storedModels{bestIteration};
    
    residualErrors = residual_Error(bestFit, x, y);
    inlierIndex = find(residualErrors < parameters.inlierDistThreshold);
end