function [  ] = triangulationF( triangPoints, cameraCenterL, cameraCenterR )

    figure; axis equal;  hold on; 
    plot3(-triangPoints(:,1), triangPoints(:,2), triangPoints(:,3), '.r');
    plot3(-cameraCenterL(1), cameraCenterL(2), cameraCenterL(3),'*g');
    plot3(-cameraCenterR(1), cameraCenterR(2), cameraCenterR(3),'*b');
    grid on; xlabel('x'); ylabel('y'); zlabel('z'); axis equal;
    
end