% Enhancing the display                    
  DisplayGrid = 1;
  Size1 = 14;
  Size2 = 16;
  
 % a = axis; axis([Fmin Fmax a(3) a(4)]);
  set(gca,'XGrid','off')
  set(gca,'YGrid','off')
  set(gca,'FontName','Arial')    
  set(gca,'LineWidth',2)    
  set(gca,'FontSize',Size1)    
  set(get(gca,'XLabel'),'fontsize',Size2,'FontWeight','bold','FontName','Arial')
  %a = get(get(gca,'XLabel'),'Position');
  %a(2) = a(2) - abs(a(2))*0.01;
  %set(get(gca,'XLabel'),'Position',a);  
  set(get(gca,'YLabel'),'fontsize',Size2,'FontWeight','bold','FontName','Arial')
  set(get(gca,'Title'),'fontsize',Size2,'FontWeight','bold','FontName','Arial')
  set(gca,'Box','On')  
  set(gcf,'WindowStyle','normal')
  set(gcf,'Position',[300 240 1000 600])
  %set(findobj(gca,'Type','line'),'LineWidth',[LineWidth])
  grid on