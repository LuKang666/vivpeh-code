% swEPSfigure.m
%
%   Set the default font names and sizes for the eps figures
%       prepared for Scientific Word
%   In SW, a 65-50% reduction of the figures is normally done
%   Full LaTeX commands can be used in the labels, legends, etc.
%
%


set(groot,'DefaultAxesFontName','Times New Roman');
set(groot,'DefaultTextFontName','Times New Roman');
set(groot,'DefaultAxesFontSize',22);
set(groot,'DefaultTextFontSize',22);
set(groot,'DefaultLineLineWidth',3.0);
set(groot,'DefaultLineMarkerSize',12);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultlegendinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultFigurePosition',[50,100,1200,600]);    %改变Figure宽度为515,原值为560
%%%%%%%%%%%%%%%%%%%%%%%
ONEFIG='set(gcf,''unit'',''normalized'',''position'',[0.1,0.2,0.33,0.45]); set(gca,''Position'',[0.11,0.115,0.85,0.86]);';
%%%%%%%%%%%%%%%%%%%%%%%
TWOFIGS='set(gcf,''unit'',''normalized'',''position'',[0.1,0.2,0.56,0.4]);';
FIG1='set(gca,''Position'',[0.068,0.125,0.425,0.85])';
FIG2='set(gca,''Position'',[0.562,0.125,0.425,0.85])';
%%%%%%%%%%%%%%%%%%%%%%%
FOURFIGS='set(gcf,''unit'',''normalized'',''position'',[0.1,0.2,0.45,0.65]);';
FIG11='set(gca,''Position'',[0.082,0.57,0.41,0.41]);';
FIG12='set(gca,''Position'',[0.57,0.57,0.41,0.41]);';
FIG21='set(gca,''Position'',[0.082,0.08,0.41,0.41]);';
FIG22='set(gca,''Position'',[0.57,0.08,0.41,0.41]);';

% 以下是原文件设置
%set(0,'DefaultAxesFontName','Times New Roman');

%set(0,'DefaultTextFontName','Times New Roman');

%set(0,'DefaultAxesFontSize',24);
%set(0,'DefaultTextFontSize',24);

%set(0,'defaulttextinterpreter','latex');


% disp('  ');
% disp('  Changing Default Font to Times New Roman');
% disp('  Changing Default Font Size to 24');
% disp('  ');
