close all;
% left1 =0.0;
% left2 = 0.35;
% left3= 0.7;
% bottom= 0.3;
% width = 0.33;
% height = 0.5;
tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
subplot(1,3,1);
h=pcolor(X2,Z2,squeeze(w(st:en,end/2,:)));
set(h,'edgecolor','none');

%pbaspect([1 1.5 1]);
xlabel('x(m)');
ylabel('d(m)');
colorbar;
caxis([0 2]);
colormap(othercolor('BuDOr_18'))
 title('w(m/s)');
 %set(gca,'OuterPosition',[left1,bottom, width, height])
subplot(1,3,2);
h=pcolor(X2,Z2,squeeze(s(st:en,end/2,:)));
set(h,'edgecolor','none');
caxis([33 35]);
 %pbaspect([1 1.4 1]); 
 title('S (psu)');
 xlabel('x(m)');
ylabel('d(m)');
shading interp
colorbar;
%set(gca,'OuterPosition',[left2,bottom, width, height])
colormap(othercolor('BuDRd_18'))
subplot(1,3,3);
h=pcolor(X2,Z2,squeeze(velmag(st:en,end/2,:)));
set(h,'edgecolor','none');
%caxis([0 .3]);
%caxis([1005 1020]);
 %pbaspect([1 1.5 1]);
 %title('t(deg C)');
 title('vel. mag (m/s)');
 xlabel('x(m)');
ylabel('d(m)');
colorbar;
colormap(othercolor('BuDRd_18'))
%set(gca,'OuterPosition',[left3,bottom, width, height])
sgtitle(titl); 
%set(gcf,'position',get(0,'ScreenSize'))
x0=10;
y0=10;
width=1200;
height=300;
set(gcf,'position',[x0,y0,width,height])

file=sprintf('%sX05m_%ds.png',outdir,timeframe);
    
saveas(gca,file);