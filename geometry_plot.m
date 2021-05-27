function geometry_plot(...
    geomTXT,geomNO,xyz_file,idx,conn,...
    show_ring_path,show_spectator_atoms)
% This function produces simple visualisation
% of the provided geometry. The arguments
% show_ring_path and show_spectator_atoms
% serve as a swich to turn off (= 0) or
% on (= 1) visualisations of the RCM path
% and coloring of the spectator atoms,
% respectively, according to the provided
% csv files.


BT = 0.9; % bond thickness
BT2 = 0.7; % C-H bond thickness
SS = 17; % sizes of atom
SS2 = 100; % sizes of spectator atom
ringT = 6; % ring current thickness
trp = 0.5; % transparency
spec_trp = 0.5; % spectator atoms transparency
text_font_size = 8; % font size for the text

% defining the xyz coordinates of individual atoms
c = geomNO(find(geomTXT== "C"),:);
h = geomNO(find(geomTXT== "H"),:);
n = geomNO(find(geomTXT== "N"),:);
s = geomNO(find(geomTXT== "S"),:);
f = geomNO(find(geomTXT== "F"),:);
o = geomNO(find(geomTXT== "O"),:);
zn = geomNO(find(geomTXT== "Zn"),:);
at_rest = geomNO(find(geomTXT ~= "C" & geomTXT ~= "H"...
    & geomTXT ~= "N" & geomTXT ~= "S" ...
    & geomTXT ~= "F" & geomTXT ~= "O" ...
    & geomTXT ~= "Zn"),:);

% this adjusts pop-up window in the centre of the screen
temp_handle = get(0,'ScreenSize');
fig = figure('units','centimeters','position',...
    [temp_handle(3:4)/2-4.5,9,9]);
movegui('center');
view([0 90]);

distMat = squareform(pdist(geomNO)); % distance matrix
adjMat = distMat < 1.75; % cutoff for drawing bonds
for i=1:size(adjMat,1)
    for j=1:size(adjMat,2)
        if i ~= j && adjMat(i,j) == 1
            hold on
            g = plot3([geomNO(i,1);geomNO(j,1)],...
                [geomNO(i,2);geomNO(j,2)],...
                [geomNO(i,3);geomNO(j,3)],'k');
            set(g,'LineStyle','-');
            if geomTXT(i,1) == "H" & ...
                    geomTXT(j,1) == "C" || ...
                    geomTXT(i,1) == "C" & ...
                    geomTXT(j,1) == "H"
                set(g,'LineWidth',BT2);
                hold on
            else
                set(g,'LineWidth',BT);
                hold on
            end
        else
            % do nothing
        end
    end
end

% drawing atoms
s1 = scatter3(c(:,1),c(:,2),c(:,3),SS,'filled');
s1.MarkerFaceColor = [0 0 0];
s1.MarkerEdgeColor = [0 0 0];
hold on
s2 = scatter3(s(:,1),s(:,2),s(:,3),SS,'filled');
s2.MarkerFaceColor = [1 165/255 0];
s2.MarkerEdgeColor = [0 0 0];
hold on
s3 = scatter3(h(:,1),h(:,2),h(:,3),SS,'filled');
s3.MarkerFaceColor = 0.8*[1 1 1];
s3.MarkerEdgeColor = [0 0 0];
hold on
s4 = scatter3(n(:,1),n(:,2),n(:,3),SS,'filled');
s4.MarkerFaceColor = [0 0 1];
s4.MarkerEdgeColor = [0 0 0];
hold on
s5 = scatter3(f(:,1),f(:,2),f(:,3),SS,'filled');
s5.MarkerFaceColor = [178 102 255]/255;
s5.MarkerEdgeColor = [0 0 0];
hold on
s6 = scatter3(o(:,1),o(:,2),o(:,3),SS,'filled');
s6.MarkerFaceColor = [255 0 0]/255;
s6.MarkerEdgeColor = [0 0 0];
hold on
s7 = scatter3(zn(:,1),zn(:,2),zn(:,3),SS,'filled');
s7.MarkerFaceColor = [0 155 155]/255;
s7.MarkerEdgeColor = [0 0 0];
hold on
s8 = scatter3(at_rest(:,1),at_rest(:,2),at_rest(:,3),SS,'filled');
s8.MarkerFaceColor = [0.5 0.5 0.5];
s8.MarkerEdgeColor = [0 0 0];
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

if show_spectator_atoms % show spectator atoms
    for nuclei = 1:size(idx,2)
        a = idx(:,nuclei);
        a = a(find(a~=0));
        
        t1 =  scatter3(geomNO(a,1),...
            geomNO(a,2),geomNO(a,3),SS2,'filled');
        t1.MarkerFaceColor = color_list(nuclei,:);
        t1.MarkerFaceAlpha = spec_trp;
    end
else
end

if show_ring_path % show path (the path is only single loop)
    for i = 1:size(conn,1)
        r = plot3(...
            [geomNO(conn(i,2),1);geomNO(conn(i,3),1)],...
            [geomNO(conn(i,2),2);geomNO(conn(i,3),2)],...
            [geomNO(conn(i,2),3);geomNO(conn(i,3),3)]);
        r.Color = 'red';
        r.LineStyle = '-';
        r.LineWidth = abs(conn(i,1))*ringT;
        r.Color(4)=trp;
        hold on
    end
else
end



newChr = strrep(xyz_file,'_',' ');
newChr = newChr(1:end-4);

% simple annotation of the figure with the name of the xyz file
posit = [0.02 0.98 0 0];
h1 = annotation('textbox',...
    posit, 'String',...
    newChr, 'FitBoxToText',....
    true,...
    'FontSize',8,...
    'EdgeColor','none');

hold off
axis 'equal'
axis tight
set(gca,'xtick',[],'ytick',[])
xlim([-15 15])
ylim([-15 15])
axis off
set(gca, 'Position', ...
    get(gca, 'OuterPosition') - ...
    get(gca,'TightInset')...
    *[-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

end



