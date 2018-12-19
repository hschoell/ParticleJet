clear all; close all; clc;

i = 1;
M = zeros(10000,5);
dt = 0.01;
St = 1;

while i <= 10000
    
    y_pos = rand*0.1 + 0.45;
    
    M(i,:) = [1,0,y_pos,exp(-(y_pos - 0.5).^2/(0.05^2)),0];
    
    for j = 1:i
        if M(j,1) == 1
            M(j,2) = M(j,2) + M(j,4)*dt;
            M(j,3) = M(j,3) + M(j,5)*dt;
            M(j,4) = rand;
            k1 = (1/St)*((rand - 0.5) + M(j,5));
            middle = M(j,5) + k1*dt;
            k2 = (1/St)*((rand - 0.5) + middle);
            M(j,5) = middle + k2*dt;
        end
        
        if M(j,1) == 1 && (M(j,2) < 0 || M(j,2) > 2 || M(j,3) < 0 || M(j,3) > 1)
            M(j,1) = 0;
            M(j,2) = -10;
            M(j,3) = -10;
        end
    end
    
   
    
    plot(M(:,2),M(:,3),'.');
    xlim([0,2]);
    ylim([0,1]);
    
    i = i + 1;
    
    pause(0.015);
    
end
