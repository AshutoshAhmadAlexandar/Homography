function residuals = residual_error(H, homoCoord1, homoCoord2)
    transformedPoints = homoCoord1 * H;
    
    lambda_t =  transformedPoints(:,3); 
    lambda_2 = homoCoord2(:,3);    
    cartDistX = transformedPoints(:,1) ./ lambda_t - homoCoord2(:,1) ./ lambda_2;
    cartDistY = transformedPoints(:,2) ./ lambda_t - homoCoord2(:,2) ./ lambda_2;
    residuals = cartDistX .* cartDistX + cartDistY .* cartDistY;
end