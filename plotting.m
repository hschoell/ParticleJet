 %plotting
    %particles
    figure(1);
    hold on
    hs=streamslice(x(2:end),y(2:end),u(1:nx,1:ny)',v(1:nx,1:ny)');
    set(hs,'color','k','linewidth',1)
    plot(M(:,2),M(:,3),'.');
    xlim([0,2]);
    ylim([0,1]);
    hold off
    %u velocity
    figure(2)
    surf(u')
    %v velocity
    figure(3)
    surf(v')