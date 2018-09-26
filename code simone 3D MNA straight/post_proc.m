maxiter=1000;
close all
for k=0:maxiter
    load(['convergence_',num2str(k,'%03d')])
    load(['configuration_',num2str(k,'%03d')])
%     load(['VMstress_',num2str(k,'%03')])
    scrsz = get(groot,'ScreenSize');
    f1=figure(2);
    set(f1,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
    hold on
    plot(outeriter,perfo,'bo','MarkerFaceColor','b')
    plot(outeriter,Vol*100,'ro','MarkerFaceColor','r')
    plot(outeriter,(GKSl)*100,'ko','MarkerFaceColor','k')
%     plot(outeriter,C,'mo','MarkerFaceColor','m')
    title(['volfrac = ',num2str(Vol*100),' %, \Delta TSFC =',num2str(perfo),'% ,G_{KS} \times 100 = ',num2str((GKSl)*100),' % , iter = ', num2str(k)])
    grid on
    xlabel('iter')
    legend('\Delta TSFC %','Volume Fraction %','G_{KS} \times 100')
%     axis([0 k+1 -1 10000])
    print(['convergence_',num2str(k,'%04d')],'-dpng')
    h = figure(1);clf
    set(h,'Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
    set(h,'Color',[1 1 1]);
    axes1 = axes('Parent',h);
    hold on; h24_patch = patch('Parent',axes1,'Vertices',coordo,'Faces',FACES_shown,'FaceVertexCData',[0.4 0.8 1],'FaceColor','flat'); axis equal; axis off;
    quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
    text(1000,0,0,'x')
    text(0,1000,0,'y')
    text(0,0,1000,'z')
    view([27.6 18]);
    % Set the remaining axes properties
    set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[434 342.3 684.6]);
    drawnow;
    hold off
     print(['configuration_',num2str(k,'%04d')],'-dpng')
%       h = figure(6);clf
%     set(h,'Color',[1 1 1]);
%     set(h,'Position',[1 1 scrsz(3)/2 scrsz(4)/2])
%     axes1 = axes('Parent',h);
%     hold on;  patch('Parent',axes1,'Vertices',coordo,'Faces',FACES_shown,'FaceVertexCData',ColorOdsR,'FaceColor','interp','EdgeColor','none'); axis equal; axis off; drawnow;
%     colormap jet;
%     axis equal
%     axis off
%     s1=quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k');
%     s1=quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k');
%     s1=quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k');
%     text(1000,0,0,'x')
%     text(0,1000,0,'y')
%     text(0,0,1000,'z')
%     view([27.6 18]);
%     % Set the remaining axes properties
%     set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
%         'PlotBoxAspectRatio',[434 342.3 684.6]);
%     C=colorbar;
%     caxis([-1 0.5])
%     C.Limits=[-1 0.5];
%     C.Label.String = 'x(\sigma_{VM}/\sigma_{lim}-1) ';
%     C.FontWeight='bold';
%     C.Location='southoutside';
%     drawnow;
%     title(['Von Mises stress at iteration ',num2str(k)])
%     print(['VMstress_',num2str(k,'%03d')],'-dpng')
 
    
end
%  load(['VMstress_',num2str(k,'%03d')])
%   h = figure(6);clf
%     set(h,'Color',[1 1 1]);
%     set(h,'Position',[1 1 scrsz(3)/2 scrsz(4)/2])
%     axes1 = axes('Parent',h);
%     hold on;  patch('Parent',axes1,'Vertices',coordo,'Faces',FACES_shown,'FaceVertexCData',ColorOdsR,'FaceColor','interp','EdgeColor','none'); axis equal; axis off; drawnow;
%     colormap jet;
%     axis equal
%     axis off
%     s1=quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k');
%     s1=quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k');
%     s1=quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k');
%     text(1000,0,0,'x')
%     text(0,1000,0,'y')
%     text(0,0,1000,'z')
%     view([27.6 18]);
%     % Set the remaining axes properties
%     set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
%         'PlotBoxAspectRatio',[434 342.3 684.6]);
%     C=colorbar;
%     caxis([-1 0.5])
%     C.Limits=[-1 0.5];
%     C.Label.String = 'x(\sigma_{VM}/\sigma_{lim}-1) ';
%     C.FontWeight='bold';
%     C.Location='southoutside';
%     drawnow;
%     title(['Von Mises stress at iteration ',num2str(k)])
%     print(['VMstress_',num2str(k,'%03d')],'-dpng')
for k=0:360
     h = figure(1);
%     set(h,'Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
%     set(h,'Color',[1 1 1]);
%     axes1 = axes('Parent',h);
%     hold on; h24_patch = patch('Parent',axes1,'Vertices',coordo,'Faces',FACES_shown,'FaceVertexCData','c','FaceColor','flat'); axis equal; axis off;
%     quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
%     quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
%     quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
%     text(1000,0,0,'x')
%     text(0,1000,0,'y')
%     text(0,0,1000,'z')
    view([27.6+k 18]);
    % Set the remaining axes properties
%     set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
%         'PlotBoxAspectRatio',[434 342.3 684.6]);
%     drawnow;
%     hold off
     print(['angle/configuration_',num2str(maxiter,'%03d'),'_angle',num2str(k,'%03d')],'-dpng')
end