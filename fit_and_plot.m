function fit_and_plot(Ax,Ay,Az,...
    write_storex,write_storey,write_storez)
% This   

x = zeros(size(write_storez,1),1);
for i = 1:size(write_storez,1)
    x(i,1) = ...   
    [Ax*mean(write_storex{i,1},1)+...
        Ay*mean(write_storey{i,1},1)+...
        Az*mean(write_storez{i,1},1)]/3;   
end

[NMR_file,NMR_path] = uigetfile('*.csv',...
    ['Select file with '...
    'chemical shifts of the aromatic - non-aromatic ',...
    'listed in the same order as in the lists of atoms']);
if isequal(NMR_file,0)
   disp('No NMR file provided');
else
   disp(['User selected ', fullfile(NMR_path,NMR_file)]);
   y = csvread(fullfile(NMR_path,NMR_file));
end


text_font_size = 8;
SS = 50;
 
format long
b1 = x\y; % b1 is the slope  

yCalc1 = b1*x;
Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);

temp_handle = get(0,'ScreenSize');
fig = figure('units','centimeters','position',...
    [temp_handle(3:4)/2-4.5,9,9]);
movegui('center');
fig.Position(1) = fig.Position(1)+9.5;


g = scatter(x,y,SS,'filled');
g.MarkerFaceColor = [1 1 1]; 
g.MarkerEdgeColor = [1 1 1]; 
hold on


if size(x,1) > 1
    marg_x = 0.15;
    marg_y = 0.15;
xlim([(min(x)-(max(x)-min(x))*marg_x),...
    (max(x)+(max(x)-min(x))*marg_x)]);

ylim([(min(y)-(max(y)-min(y))*marg_y),...
    (max(y)+(max(y)-min(y))*marg_y)]);
else
xlim([-0.2 0.2])
ylim([-2 2])
disp('x and y limits set to smallest [2 2] value')
end


p1 = plot([min(x)-1;max(x)+1],...
    [min(x)-1;max(x)+1]*b1,'b-');
p1.LineWidth = 1;
hold on 

p2 = plot([min(x)-10;max(x)+10],...
    [0;0],'k-');
p2.LineWidth = 0.5;
hold on 

p3 = plot([0;0],[2*min(y)-10;2*max(y)+10],...
    'k-');
p3.LineWidth = 0.5;
hold on 

color_list = [... % coloring for spectator atoms 
    0.1255,0.4667,0.7098;
    0.9608,0.4941,0.1255;
    0.1765,0.6275,0.2824;
    0.8392,0.1608,0.1569;
    0.5725,0.4078,0.6745;
    0.5490,0.3412,0.2980;
    0.8471,0.4784,0.6941;
    0.4980,0.4980,0.4980;
    0.7373,0.7451,0.1961;
    0.1137,0.7490,0.8118];

for i = 1:size(x,1)
    g1 = scatter(x(i,1),y(i,1),SS,'filled');
    g1.MarkerFaceColor = color_list(i,:); 
    g1.MarkerEdgeColor = [0 0 0]; 
    g1.MarkerFaceAlpha = 0.5;
    hold on
end


b1fak = num2str(b1);

if b1 > 0
MyString2 = ['y = ',b1fak(1:4),' x'];
else
MyString2 = ['y = ',char(8211),b1fak(2:5),' x'];
end

posit = [0.5 0.85 0 0];
if isnan(Rsq1)
MyString = ['             '];
else
Rfak = num2str(Rsq1);
MyString = ['R^{2} = ',Rfak(1:4),'  '];
end

title_2 = text(mean(xlim),...
    max(ylim)-0.12*(max(ylim)-min(ylim)),...
    MyString2, ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom',...
     'FontSize',text_font_size,...
     'BackgroundColor', 'white',...
     'EdgeColor','none');
 
title_3 = text(mean(xlim),...
    max(ylim)-0.18*(max(ylim)-min(ylim)),...
    MyString, ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom',...
     'FontSize',text_font_size,...
     'EdgeColor','none');

xlabel(['effect of the ring current / ',...
    'ppm·nA^{',char(8211),'1} T'],...
    'fontname', 'Arial');
ylabel(['{\Delta}','{\delta} NMR chemical shift / ppm'])

x_list = xticks;
y_list = yticks;

idx_zero = find(xticks < 0.001 & xticks > -0.001);
x_list(idx_zero) = 0;

x_list_str = string(x_list);
y_list_str = string(y_list);


% Replace ugly hyphen with n-dash for
% negative values. 
for i = 1:size(x_list,2)
if x_list(1,i) < 0
    x_list_str(1,i) = ...
        [char(8211),num2str(abs(x_list(1,i)))];
end
end

for i = 1:size(y_list,2)
if y_list(1,i) < 0
    y_list_str(1,i) = ...
        [char(8211),num2str(abs(y_list(1,i)))];
end
end
   
xticklabels(x_list_str);
yticklabels(y_list_str);

set(gca,... 
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8,...
    'FontName','Ariel')
set(gca, 'Position', get(gca, 'OuterPosition') - ...
    get(gca,'TightInset')*[-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

hold off



end



