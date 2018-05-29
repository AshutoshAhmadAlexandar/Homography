function H = Homography_fit(pts1_homogenous, pts2_homogenous)

    if size(pts1_homogenous) ~= size(pts2_homogenous)
        error('Number of matched features in the subset supplied to fit_homography does not match for both images')
    end 
    
    [number_matches, ~] = size(pts1_homogenous);
    
    
    A = []; % will be 2*numMatches x 9
    for i = 1:number_matches
        p1 = pts1_homogenous(i,:);
        p2 = pts2_homogenous(i,:);
        
        % 2x9 matrix to append onto A. 
        A_i = [ zeros(1,3)  ,   -p1     ,   p2(2)*p1;
                    p1      , zeros(1,3),   -p2(1)*p1];
        A = [A; A_i];        
    end
    
    %solve for A*h = 0
    [~,~,eigenVecs] = svd(A); 
    h = eigenVecs(:,9);     
    H = reshape(h, 3, 3);   
    H = H ./ H(3,3);        
    
end