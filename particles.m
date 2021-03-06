function M = particles(M, u, v, i, dx, dy, dt, St)
    %Add one new point into the flow per timestep
    y_pos = rand*0.1 + 0.45;
    M(i,:) = [1,0.0001,y_pos,exp(-(y_pos - 0.5).^2/(0.05^2)),0];
    %vector is: [in use?, x_pos, y-_pos, x_vel, y_vel]
    
    %update every particle currently in flow
    for j = 1:i
        if M(j,1) == 1
            %find u and v at the particle position            
            x_pos = ceil(M(j,2)/dx);
            y_pos = ceil(M(j,3)/dy);            
            
            %use linear interpolation to find velocities
            if x_pos > 128
                u_val = 0.5*(u(x_pos, y_pos) + u(x_pos - 1, y_pos));
            else
                u_val = 0.5*(u(x_pos, y_pos) + u(x_pos + 1, y_pos));
            end
            
            if y_pos > 64
                v_val = 0.5*(v(x_pos, y_pos) + v(x_pos, y_pos - 1));
            else
                v_val = 0.5*(v(x_pos, y_pos) + v(x_pos, y_pos + 1));
            end
            
            %update
            %x velocity
            k1 = (1/St)*(u_val + M(j,4));
            middle = M(j,4) + k1*dt/2;
            k2 = (1/St)*(u_val + middle);
            M(j,4) = M(j,4) + k2*dt;
            
            %y velocity
            k1 = (1/St)*(v_val + M(j,5));
            middle = M(j,5) + k1*dt/2;
            k2 = (1/St)*(v_val + middle);
            M(j,5) = M(j,5) + k2*dt;
            
            %x position
            k1 = M(j,4);
            middle = M(j,2) + k1*dt/2;
            k2 = M(j,4);
            M(j,2) = M(j,2) + k2*dt;
            
            %y position
            k1 = M(j,5);
            middle = M(j,3) + k1*dt/2;
            k2 = M(j,5);
            M(j,3) = M(j,3) + k2*dt;            
        end
        
        %discard any particles outside the domain
        if M(j,1) == 1 && (M(j,2) < 0 || M(j,2) > 2 || M(j,3) < 0 || M(j,3) > 1)
            M(j,1) = 0;
            M(j,2) = -10;
            M(j,3) = -10;
        end
    end

end