function plotmat2(i1,i2,i3,dh,dv,nh,nv,mat,str)
% This function plots 2D image
h= subplot(i1,i2,i3); imagesc(0:(nh-1),0:(nv-1),mat);     
axis ij;
if verLessThan('matlab','8.4')
    %Valid in Matlab 2013
    xtk = dh*get(h,'XTick');        set(h,'XTickLabel',xtk)
    ytk = dv*get(h,'YTick');        set(h,'YTickLabel',ytk)
    ax = gca;                       set(ax,'UserData',[dh,dv]);
    h=gcf;                          set(h,'ResizeFcn','axTickUpdate');

    set(gca,'XAxisLocation','top')
    colorbar;       colormap(flipud(gray));
    title(str);
    xlabel('X (m)');    ylabel('Z (m)')

else
%%% Valid in Matlab 2016
    h.XTickLabel=dh.*h.XTick;
    h.YTickLabel = dv*h.YTick;

    h=gcf;
    ax = gca;
    ax.UserData.dh = dh;
    ax.UserData.dv = dv;

    h.SizeChangedFcn = 'axTickUpdate';

    set(gca,'XAxisLocation','top')
    hcb=colorbar;       
    set(get(hcb,'Title'),'String',str)
    colormap(flipud(gray));
    %title(str);
    xlabel('X (m)');    ylabel('Z (m)')
end
end