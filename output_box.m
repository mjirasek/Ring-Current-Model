function continue_fitting = output_box(Ax,Ay,Az,...
    write_storex,write_storey,write_storez)
% This function fits a linear function to the
% experimental NMR data and RCGF. The data and
% fit are then ploted.

x = zeros(size(write_storez,1),1);
for i = 1:size(write_storez,1)
    x(i,1) = ...
        [Ax*mean(write_storex{i,1},1)+...
        Ay*mean(write_storey{i,1},1)+...
        Az*mean(write_storez{i,1},1)]/3;
end


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


answer_list = "The averaged RCGF for the atoms are:";
answer_list_console = answer_list;

for i = 1:size(x,1)
    answer_list = [answer_list;
        strcat("\color[rgb]{",num2str(color_list(i,:)),"}atoms no. ",num2str(i),...
        " \color{black} = ", string(x(i,1))," ppm T / nA")];
end

for i = 1:size(x,1)
    answer_list_console = [answer_list_console;
        strcat("atoms no. ",num2str(i),...
        " = ", string(x(i,1))," ppm T / nA")];
end

opts.Interpreter = 'tex';
opts.Default = 'OK and exit';

disp(answer_list_console);
answer = questdlg(answer_list, ...
    'RCGF for selected atoms', ...
    'OK and exit','continue to fit to NMR data',opts);
% Handle response



switch answer
    case 'OK and exit'
        continue_fitting = 0;
    case 'continue to fit to NMR data'
        continue_fitting = 1;
end
opts.Interpreter = 'Default';

end
