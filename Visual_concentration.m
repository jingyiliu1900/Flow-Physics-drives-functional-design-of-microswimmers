function Visual_concentration(r, theta, c_ses, c_mot, Sh_ses, Sh_mot)
    % Transfer into cartesian coordinates
    z = (r)*cos(theta);
    x = (r)*sin(theta);
    % Bottom half domain 
    x_extend = [x,-x];
    z_extend = [z,z];

    c_ses_extend = [c_ses, c_ses]; % sessile
    c_mot_extend = [c_mot, c_mot]; % sessile

    theta2 = 0:pi/100:2*pi;
    figure('color','w')
    colormap(hot)
    subplot(1,2,1)
    pcolor(x_extend, z_extend, c_ses_extend)
    shading interp; colorbar;
    xlim([-5,5]); ylim([-5,5]); axis square; set(gca,'xtick',[]);set(gca,'ytick',[]);
    Circle1 = cos(theta2);   Circle2 = sin(theta2);
    hold on; plot(Circle1,Circle2, 'linewidth',1,'color','k')
    title(['sessile, Sh=', num2str(Sh_ses)]); hold off

    subplot(1,2,2)
    pcolor(x_extend, z_extend, c_mot_extend)
    colormap(hot)
    shading interp; colorbar;
    xlim([-5,5]); ylim([-5,5]); axis square; set(gca,'xtick',[]);set(gca,'ytick',[]);
    Circle1 = cos(theta2);   Circle2 = sin(theta2);
    hold on; plot(Circle1,Circle2, 'linewidth',1,'color','k')
    title(['motile, Sh=', num2str(Sh_mot)]); hold off
end
