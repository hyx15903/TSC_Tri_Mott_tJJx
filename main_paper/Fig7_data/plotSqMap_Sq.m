function plotSqMap_Sq()
clear;
clc;
FONTSIZE = 16;
POS_X = 0.26;
POS_INIT_X = 0.05;
POS_Y = 0.39;
POS_INIT_Y = 0.06;

PLOT_RANGE = 0.8;
COLOR_TYPE = 'jet';
LINE_COLOR = 'black';
LINE_WIDTH = 0.4;
MARKERSIZE = 7;
%// surf plot SqMap
figure()

%// 1
pos = [POS_INIT_X (POS_INIT_Y*2.+POS_Y) POS_X POS_Y];
subplot('position',pos)
Test = load(['./Kspace_tJ_Spin_Stru_Nx36Ny6j20.01jk0.05SU2.txt']);
kx = Test(:,1)./(2.*pi);
ky = Test(:,2)./(2.*pi);
Sq = Test(:,3);
[X,Y] = meshgrid(-PLOT_RANGE:0.01:PLOT_RANGE,-PLOT_RANGE:0.01:PLOT_RANGE); % 0.01 for sampling step length
Z=griddata(kx,ky,Sq,X,Y,'v4');
h = surf(X,Y,Z); 
hold on;
set(h,'edgecolor','none');
axis([-PLOT_RANGE,PLOT_RANGE,-PLOT_RANGE,PLOT_RANGE]);
caxis([0, 3.5]);
pbaspect([1 1 1]);
colormap(COLOR_TYPE);
set(get(gca, 'XLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
xlabel('$k_{x}/2\pi$');
set(get(gca, 'YLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
ylabel('$k_{y}/2\pi$');
box on;
text(0.67,0.12,5, '$\mathbf{K}$', 'Interpreter','latex','FontSize',14)
text(0.55,0.35,5, '$\mathbf{M}$', 'Interpreter','latex','FontSize',14)

line([-0.333,0.333],[1./sqrt(3.),1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.333,0.667],[1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.667,0.333],[0,-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,0.333],[-1./sqrt(3.),-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,-0.667],[-1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.667,-0.333],[0,1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
set(gca, 'XTick', [-0.5 0 0.5])
set(gca, 'YTick', [-0.5 0 0.5])
xticklabels(gca,{'$-0.5$','$0$','$0.5$'})
set(gca, 'fontsize', FONTSIZE,'TickLabelInterpreter','latex')
set(gca, 'linewidth', 2)
set(get(gca, 'Children'), 'linewidth', 2)
title('$(a) J_{2}=0.01, J_{\chi}=0.05$','FontSize',FONTSIZE,...
    'Interpreter','latex');
view(2);
ax=gca;
plot3(0.667,0,5,'.k','MarkerSize',25)
plot3(0.5,(1/sqrt(3.)/2.),5,'.k','MarkerSize',25)
for x=-16:1:16
    for y=-6:1:6
        plot3((y/12.)+(2./12.*x/2.),(y/12.*sqrt(3.))-(2./sqrt(3.)/12.*x/2.),5,'.k','MarkerSize',MARKERSIZE)
    end
end
surf(ax)

%// 2
pos = [(POS_INIT_X*2.+POS_X) (POS_INIT_Y*2.+POS_Y) POS_X POS_Y];
subplot('position',pos)
Test = load(['./Kspace_tJ_Spin_Stru_Nx36Ny6j20.05jk0.05SU2.txt']);
kx = Test(:,1)./(2.*pi);
ky = Test(:,2)./(2.*pi);
Sq = Test(:,3);
[X,Y] = meshgrid(-1:0.01:1,-1:0.01:1); % 0.01 for sampling step length
Z=griddata(kx,ky,Sq,X,Y,'v4');
hold on;
h = surf(X,Y,Z); 
set(h,'edgecolor','none');
axis([-PLOT_RANGE,PLOT_RANGE,-PLOT_RANGE,PLOT_RANGE]);
caxis([0, 3.5]);
pbaspect([1 1 1]);
colormap(COLOR_TYPE);
box on;

line([-0.333,0.333],[1./sqrt(3.),1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.333,0.667],[1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.667,0.333],[0,-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,0.333],[-1./sqrt(3.),-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,-0.667],[-1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.667,-0.333],[0,1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
set(get(gca, 'XLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
xlabel('$k_{x}/2\pi$');
set(get(gca, 'YLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
ylabel('$k_{y}/2\pi$');
set(gca, 'XTick', [-0.5 0 0.5])
set(gca, 'YTick', [-0.5 0 0.5])
xticklabels(gca,{'$-0.5$','$0$','$0.5$'})
set(gca, 'fontsize', FONTSIZE,'TickLabelInterpreter','latex')
set(gca, 'linewidth', 2);
set(get(gca, 'Children'), 'linewidth', 2);
title('$(b) J_{2}=0.05, J_{\chi}=0.05$','FontSize',FONTSIZE,...
    'Interpreter','latex');
view(2);
ax=gca;
for x=-16:1:16
    for y=-6:1:6
        plot3((y/12.)+(2./12.*x/2.),(y/12.*sqrt(3.))-(2./sqrt(3.)/12.*x/2.),5,'.k','MarkerSize',MARKERSIZE)
    end
end
surf(ax)


%// 3
pos = [(POS_INIT_X*3.+POS_X*2.) (POS_INIT_Y*2.+POS_Y) POS_X POS_Y];
subplot('position',pos)
Test = load(['./Kspace_tJ_Spin_Stru_Nx36Ny6j20.05jk0.1SU2.txt']);
kx = Test(:,1)./(2.*pi);
ky = Test(:,2)./(2.*pi);
Sq = Test(:,3);
[X,Y] = meshgrid(-1:0.01:1,-1:0.01:1); % 0.01 for sampling step length
Z=griddata(kx,ky,Sq,X,Y,'v4');
hold on;
h = surf(X,Y,Z); 
set(h,'edgecolor','none');
axis([-PLOT_RANGE,PLOT_RANGE,-PLOT_RANGE,PLOT_RANGE]);
caxis([0, 3.5]);
pbaspect([1 1 1]);
colormap(COLOR_TYPE);
colorbar;
set(get(gca, 'XLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
xlabel('$k_{x}/2\pi$');
set(get(gca, 'YLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
ylabel('$k_{y}/2\pi$');
box on;

line([-0.333,0.333],[1./sqrt(3.),1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.333,0.667],[1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.667,0.333],[0,-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,0.333],[-1./sqrt(3.),-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,-0.667],[-1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.667,-0.333],[0,1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
set(gca, 'XTick', [-0.5 0 0.5])
set(gca, 'YTick', [-0.5 0 0.5])
xticklabels(gca,{'$-0.5$','$0$','$0.5$'})
set(gca, 'fontsize', FONTSIZE,'TickLabelInterpreter','latex')
set(gca, 'linewidth', 2);
set(get(gca, 'Children'), 'linewidth', 2);
title('$(c) J_{2}=0.05, J_{\chi}=0.1$','FontSize',FONTSIZE,...
    'Interpreter','latex');
view(2);
ax=gca;
set(ax, 'Position', pos);
for x=-16:1:16
    for y=-6:1:6
        plot3((y/12.)+(2./12.*x/2.),(y/12.*sqrt(3.))-(2./sqrt(3.)/12.*x/2.),5,'.k','MarkerSize',MARKERSIZE)
    end
end
surf(ax)

%// 4
pos = [POS_INIT_X POS_INIT_Y POS_X POS_Y];
subplot('position',pos)
Test = load(['./Kspace_tJ_Spin_Stru_Nx36Ny6j20.1jk0.05SU2.txt']);
kx = Test(:,1)./(2.*pi);
ky = Test(:,2)./(2.*pi);
Sq = Test(:,3);
[X,Y] = meshgrid(-1:0.01:1,-1:0.01:1); % 0.01 for sampling step length
Z=griddata(kx,ky,Sq,X,Y,'v4');
hold on;
h = surf(X,Y,Z); 
set(h,'edgecolor','none');
axis([-PLOT_RANGE,PLOT_RANGE,-PLOT_RANGE,PLOT_RANGE]);
caxis([0, 3.5]);
pbaspect([1 1 1]);
colormap(COLOR_TYPE);
box on;

line([-0.333,0.333],[1./sqrt(3.),1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.333,0.667],[1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.667,0.333],[0,-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,0.333],[-1./sqrt(3.),-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,-0.667],[-1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.667,-0.333],[0,1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
set(get(gca, 'XLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
xlabel('$k_{x}/2\pi$');
set(get(gca, 'YLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
ylabel('$k_{y}/2\pi$');
set(gca, 'XTick', [-0.5 0 0.5])
set(gca, 'YTick', [-0.5 0 0.5])
xticklabels(gca,{'$-0.5$','$0$','$0.5$'})
set(gca, 'fontsize', FONTSIZE,'TickLabelInterpreter','latex')
set(gca, 'linewidth', 2);
set(get(gca, 'Children'), 'linewidth', 2);
title('$(d) J_{2}=0.1, J_{\chi}=0.05$','FontSize',FONTSIZE,...
    'Interpreter','latex');
view(2);
ax=gca;
for x=-16:1:16
    for y=-6:1:6
        plot3((y/12.)+(2./12.*x/2.),(y/12.*sqrt(3.))-(2./sqrt(3.)/12.*x/2.),5,'.k','MarkerSize',MARKERSIZE)
    end
end
surf(ax)

%// 5
pos = [(POS_INIT_X*2.+POS_X) POS_INIT_Y POS_X POS_Y];
subplot('position',pos)
Test = load(['./Kspace_tJ_Spin_Stru_Nx36Ny6j20.1jk0.1SU2.txt']);
kx = Test(:,1)./(2.*pi);
ky = Test(:,2)./(2.*pi);
Sq = Test(:,3);
[X,Y] = meshgrid(-1:0.01:1,-1:0.01:1); % 0.01 for sampling step length
Z=griddata(kx,ky,Sq,X,Y,'v4');
hold on;
h = surf(X,Y,Z); 
set(h,'edgecolor','none');
axis([-PLOT_RANGE,PLOT_RANGE,-PLOT_RANGE,PLOT_RANGE]);
caxis([0, 3.5]);
pbaspect([1 1 1]);
colormap(COLOR_TYPE);
box on;

line([-0.333,0.333],[1./sqrt(3.),1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.333,0.667],[1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.667,0.333],[0,-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,0.333],[-1./sqrt(3.),-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,-0.667],[-1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.667,-0.333],[0,1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
set(get(gca, 'XLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
xlabel('$k_{x}/2\pi$');
set(get(gca, 'YLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
ylabel('$k_{y}/2\pi$');
set(gca, 'XTick', [-0.5 0 0.5])
set(gca, 'YTick', [-0.5 0 0.5])
xticklabels(gca,{'$-0.5$','$0$','$0.5$'})
set(gca, 'fontsize', FONTSIZE,'TickLabelInterpreter','latex')
set(gca, 'linewidth', 2);
set(get(gca, 'Children'), 'linewidth', 2);
title('$(e) J_{2}=0.1, J_{\chi}=0.1$','FontSize',FONTSIZE,...
    'Interpreter','latex');
view(2);
ax=gca;
for x=-16:1:16
    for y=-6:1:6
        plot3((y/12.)+(2./12.*x/2.),(y/12.*sqrt(3.))-(2./sqrt(3.)/12.*x/2.),5,'.k','MarkerSize',MARKERSIZE)
    end
end
surf(ax)

%// 6
pos = [(POS_INIT_X*3.+POS_X*2.) POS_INIT_Y POS_X POS_Y];
subplot('position',pos)
Test = load(['./Kspace_tJ_Spin_Stru_Nx36Ny6j20.2jk0.2SU2.txt']);
kx = Test(:,1)./(2.*pi);
ky = Test(:,2)./(2.*pi);
Sq = Test(:,3);
[X,Y] = meshgrid(-1:0.01:1,-1:0.01:1); % 0.01 for sampling step length
Z=griddata(kx,ky,Sq,X,Y,'v4');
hold on;
h = surf(X,Y,Z); 
set(h,'edgecolor','none');
axis([-PLOT_RANGE,PLOT_RANGE,-PLOT_RANGE,PLOT_RANGE]);
caxis([0, 3.5]);
pbaspect([1 1 1]);
colormap(COLOR_TYPE);
colorbar;
set(get(gca, 'XLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
xlabel('$k_{x}/2\pi$');
set(get(gca, 'YLabel'), 'FontSize', FONTSIZE,'Interpreter','latex');
ylabel('$k_{y}/2\pi$');
box on;

line([-0.333,0.333],[1./sqrt(3.),1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.333,0.667],[1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([0.667,0.333],[0,-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,0.333],[-1./sqrt(3.),-1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.333,-0.667],[-1./sqrt(3.),0],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
line([-0.667,-0.333],[0,1./sqrt(3.)],[4,4],'Color',LINE_COLOR,'LineWidth', LINE_WIDTH)
set(gca, 'XTick', [-0.5 0 0.5])
set(gca, 'YTick', [-0.5 0 0.5])
xticklabels(gca,{'$-0.5$','$0$','$0.5$'})
set(gca, 'fontsize', FONTSIZE,'TickLabelInterpreter','latex')
set(gca, 'linewidth', 2);
set(get(gca, 'Children'), 'linewidth', 2);
title('$(f) J_{2}=0.2, J_{\chi}=0.2$','FontSize',FONTSIZE,...
    'Interpreter','latex');
view(2);
ax=gca;
set(ax, 'Position', pos);
for x=-16:1:16
    for y=-6:1:6
        plot3((y/12.)+(2./12.*x/2.),(y/12.*sqrt(3.))-(2./sqrt(3.)/12.*x/2.),5,'.k','MarkerSize',MARKERSIZE)
    end
end
surf(ax)

exportgraphics(ax,'test.png','Resolution',1200)
end
