% This script executes Ring Current Model (RCM)
% on the molecule with provided xyz.
% To execute corretly, MATLAB version 2019a or newer is required.
%
% If used in publication please cite as: doi:10.xxxx/xxxxxxx
%
% GitHub page:
% https://github.com/mjirasek
%
% Questions regarding the script (how to run / found bugs)
% adress to m.jirasek@seznam.cz
%
% Alternatively, Prof. Harry Anderson can be reached on
% harry.anderson@chem.ox.ac.uk
%
% Michael Jirasek, PhD.

clear
close all

msg1 = 'Select xyz file';
disp(msg1);
[xyz_file,xyz_path] = uigetfile('*.xyz',msg1);
if isequal(xyz_file,0)
    disp('No xyz file provided');
    return
else
    disp(['User selected ', fullfile(xyz_path,xyz_file)]);
    xyz_crude = readcell( (fullfile(xyz_path,xyz_file)),...
        'FileType', 'text');
    geomTXT = string(xyz_crude(3:end,1));
    geomNO = cell2mat(xyz_crude(3:end,2:4));
    xyz = geomNO;
end

msg2 = ['Select file with numbering of atoms ',...
    'for which the induced magnetic filed ',...
    'should be calculated'];
disp(msg2);
[idx_file,idx_path] = uigetfile('*.csv',msg2);
if isequal(idx_file,0)
    disp('No connectivity file provided');
else
    disp(['User selected ', fullfile(idx_path,idx_file)]);
    idx = csvread(fullfile(idx_path,idx_file));
end

msg3 = 'Select connectivity file';
disp(msg3);
[conn_file,conn_path] = uigetfile('*.csv',msg3);
if isequal(conn_file,0)
    disp('No connectivity file provided');
else
    disp(['User selected ', fullfile(conn_path,conn_file)]);
    conn = csvread(fullfile(conn_path,conn_file));
end

separation = 0.7; % distance of the current paths
%  from the molecular plane. Use 0 for the single loop model.
show_ring_path = 0; % highlights the RCM in molecule
show_spectator_atoms = 1; % highlights examined atoms


[write_storex,write_storey,write_storez] = ...
    RCM(xyz,conn,idx,separation); % returns RCGF(x,y,z)

[Ax,Ay,Az] = area_finder(xyz,conn,idx); % returns cross-section areas

geometry_plot(geomTXT,geomNO,xyz_file,idx,...
    conn,show_ring_path,show_spectator_atoms); % geometry visualisation

continue_fitting = output_box(Ax,Ay,Az,...
    write_storex,write_storey,write_storez);

if continue_fitting
    fit_and_plot(Ax,Ay,Az,...
        write_storex,write_storey,write_storez);
else
end




    
