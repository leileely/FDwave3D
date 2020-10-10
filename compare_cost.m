% comparison of the computational cost
% reproducing Figure 6 of the related manuscript

% model size
x=[121*121*121 171*171*171 215*215*215 401*401*94 341*341*291];

% processing time for pure-MATLAB and C-MEX functions
time=[1635.4 4757.1 9865.8 4112.32 12071.54]./[801 801 801 201 241];
ctime=[72.12 245.44 511.4 824.22 2171.54]./[101 101 101 101 101];


fig = figure;
left_color = [1 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

%yyaxis left

loglog(x,time,'ro-','linewidth',2,'markersize',10);
hold on
loglog(x,ctime,'bo-','linewidth',2,'markersize',10);
ylabel('Normalized processing time (second)');
axis([10^6 10^8 0 10^2])
%set(gca,'xtick',[10^6 121*121*121 441*441*134 341*341*291],'xticklabel',[{} {'Homogeneous'} {'overthrust'}   {'Layered'} {}],'fontsize',16);

xlabel('Total number of nodes')
legend('pure-MATLAB','C-MEX');
set(gcf,'pos',[100 100 800 600]);
set(gca,'fontsize',14);
%title({'Homogeneous' 'Overthrust'  'Layered'},'fontsize',12);





