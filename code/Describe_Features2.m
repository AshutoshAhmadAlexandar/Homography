function [ featureDescription ] = Describe_Features2( img, row, col, radius )

    numberFeatures = length(row); %number of features
    featureDescription = zeros(numberFeatures, (2 * radius + 1)^2);

    padHelper = zeros(2 * radius + 1); 
    padHelper(radius + 1, radius + 1) = 1;
    paddedImg = imfilter(img, padHelper, 'replicate', 'full');

    %Extract the neighborhoods around the found features
    for i = 1 : numberFeatures
        rowRange = row(i) : row(i) + 2 * radius;
        colRange = col(i) : col(i) + 2 * radius;
        neighborhood = paddedImg(rowRange, colRange);
        flattenedFeatureVec = neighborhood(:);
        featureDescription(i,:) = flattenedFeatureVec;
    end
    
    %Normalize all descriptors to have zero mean and unit standard deviation
    featureDescription = zscore(featureDescription')';
end
