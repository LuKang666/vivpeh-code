%% fig.9_12 (共振时电压和功率的矩随beta和eta的变化曲面图,sigma=0.01)
clc; clear;
load('D:\Program Files\MATLAB\work\Script\VIVPEHs\data\Fig9-12.mat');
swEPSfigure;

beta = linspace(0.1, 3, 20);
EV = zeros(20, 20);
EP = zeros(20, 20);
for i = 1:20
    EV(i, :) = eval(['V', mat2str(i)]);
    EP(i, :) = eval(['P', mat2str(i)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
mesh(eta, beta, EV);
colormap('cool');
shading interp;
set(gca, 'xlim', [0, 3.1], ...
    'ylim', [0, 3.1], ...
    'LineWidth', 3.5, ...
    'view', [67.5, 30]);
xlabel('$\eta$'); ylabel('$\overline{\beta}$'); zlabel('$E[\overline{V}^2]$');
xticks([1, 2, 3]); yticks([0, 1, 2, 3]);
set(gcf, 'position', [0, 0, 500, 510]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
mesh(eta, beta, EP);
colormap('cool');
shading interp;
set(gca, 'xlim', [0, 3.1], ...
    'ylim', [0, 3.1], ...
    'LineWidth', 3.5, ...
    'view', [67.5, 30]);
xlabel('$\eta$'); ylabel('$\overline{\beta}$'); zlabel('$E[\overline{P}]$');
xticks([0, 1, 2, 3]); yticks([0, 1, 2, 3]);
set(gcf, 'position', [0, 0, 500, 510]);
