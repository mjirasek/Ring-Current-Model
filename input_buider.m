clear
close all

% setup of the script
show_connectivity = 0;
show_heads = 0;
show_spectator_atoms = 0;
show_numbering = 0;
show_hydrogens = 1;
print_figure = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    max_special = ceil(max(max(geomNO)))+0.5;
    min_special = floor(min(min(geomNO)))-0.5;
    max_lims = ceil(max(geomNO))+0.5;
    min_lims = floor(min(geomNO))-0.5;
end

if show_connectivity
    msg3 = 'Select connectivity file';
    disp(msg3);
    [conn_file,conn_path] = uigetfile('*.csv',msg3);
    if isequal(conn_file,0)
        disp('No connectivity file provided');
    else
        disp(['User selected ', fullfile(conn_path,conn_file)]);
        conn = csvread(fullfile(conn_path,conn_file));
    end
else
end

figsize = 2;
size_faktor = 15/max_special;

BT = 0.4*figsize*size_faktor; % bond thickness
BT2 = 0.3*figsize*size_faktor; % C-H bond thickness
SS = 17*figsize*size_faktor; % sizes of atom
SS2 = 100*figsize*size_faktor; % sizes of spectator atom
ringT = 4*figsize*size_faktor; % ring current thickness
trp = 0.5; % transparency
spec_trp = 0.5; % spectator atoms transparency
text_font_size = 8*figsize*size_faktor; % font size for the text

% defining the xyz coordinates of individual atoms
c = geomNO(find(geomTXT== "C"),:);
non_h = geomNO(find(geomTXT~= "H"),:);
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
    [temp_handle(3:4)/2-9*figsize/2,9*figsize,9*figsize]);
ax1 = axes;
movegui('center');
view([83 21]);

if show_hydrogens
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
else
    
    distMat = squareform(pdist(non_h)); % distance matrix
    adjMat = distMat < 1.75; % cutoff for drawing bonds
    for i=1:size(adjMat,1)
        for j=1:size(adjMat,2)
            if i ~= j && adjMat(i,j) == 1
                hold on
                g = plot3([non_h(i,1);non_h(j,1)],...
                    [non_h(i,2);non_h(j,2)],...
                    [non_h(i,3);non_h(j,3)],'k');
                set(g,'LineStyle','-');
                
                set(g,'LineWidth',BT);
                hold on
            else
                % do nothing
            end
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
if show_hydrogens
    s3 = scatter3(h(:,1),h(:,2),h(:,3),SS,'filled');
    s3.MarkerFaceColor = 0.8*[1 1 1];
    s3.MarkerEdgeColor = [0 0 0];
    hold on
else
end

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

if show_numbering
    offset = 0.07;
    for i = 1:size(geomNO,1)
        if show_hydrogens
            txt(i) = text(geomNO(i,1)+offset,geomNO(i,2)+offset,...
                geomNO(i,3)+offset,strcat(geomTXT(i,1),",",num2str(i)));
            txt(i).Color = 'black';
            txt(i).BackgroundColor = 'white';
            txt(i).EdgeColor = 'black';
            txt(i).FontSize = 9;
            txt(i).Margin = 1;
            txt(i).VerticalAlignment = 'baseline';
            hold on
        else
            if geomTXT(i) ~= "H"
                txt(i) = text(geomNO(i,1)+offset,geomNO(i,2)+offset,...
                    geomNO(i,3)+offset,strcat(geomTXT(i,1),",",num2str(i)));
                txt(i).Color = 'black';
                txt(i).BackgroundColor = 'white';
                txt(i).EdgeColor = 'black';
                txt(i).FontSize = 9;
                txt(i).Margin = 1;
                txt(i).VerticalAlignment = 'baseline';
                hold on
            else
            end
        end
    end
else
end


if show_spectator_atoms % show spectator atoms
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
        for nuclei = 1:size(idx,2)
            a = idx(:,nuclei);
            a = a(find(a~=0));
            
            t1 =  scatter3(geomNO(a,1),...
                geomNO(a,2),geomNO(a,3),SS2,'filled');
            t1.MarkerFaceColor = color_list(nuclei,:);
            t1.MarkerFaceAlpha = spec_trp;
        end
    end
else
end

if show_connectivity % show path (the path is only single loop)
    conn_color_map = [...
        252,12,12;...
        0,17,255;...
        204,0,255;...
        95,252,12;...
        102,102,102]/255;
    ring_current_weights = flipud(unique(conn(:,1)));
    for i = 1:size(conn,1)
        r(i) = plot3(...
            [geomNO(conn(i,2),1);geomNO(conn(i,3),1)],...
            [geomNO(conn(i,2),2);geomNO(conn(i,3),2)],...
            [geomNO(conn(i,2),3);geomNO(conn(i,3),3)]);
        r(i).Color = conn_color_map(...
            find(ring_current_weights ==  conn(i,1)),:);
        r(i).LineStyle = '-';
        r(i).LineWidth = abs(conn(i,1))*ringT;
        r(i).Color(4)=trp;
        hold on
    end
    if ~show_heads
        leg = legend([r(1:size(ring_current_weights,1))]',...
            num2str(ring_current_weights,5));
        leg.Position(1:2) = leg.Position(1:2) - 0.05;
    else
    end
    
    if show_heads
        start = geomNO(conn(:,2),1:3);
        stop = geomNO(conn(:,3),1:3);
        head_clr = zeros(size(r,2),3);
        for i = 1:size(r,2)
            head_clr(i,1:3) = r(i).Color;
        end
        dvec=stop-start;
        dis=sqrt(sum(dvec.^2,2));
        hv=min(dis)*0.35;
        cosrang=acos(dvec(:,3)./dis)*180/pi;
        nvec=[-dvec(:,2) dvec(:,1) zeros(size(dis))];
        hheads=[];
        hhgrd=[];
        pv=dis-hv;
        for i=1:length(dis)
            [xi,yi,zi] = cylinder([tan(30/180*pi),0],10);
            xi=xi*hv;yi=yi*hv;zi=zi*hv+pv(i);
            [rx,ry,rz] = rotatedata(xi,yi,zi,nvec(i,:),cosrang(i),[0,0,0]);
            cx=start(i,1)+rx;cy=start(i,2)+ry;cz=start(i,3)+rz;
            hheads(i)=surf(cx,cy,cz,'edgecolor','none','facecolor',head_clr(i,:)*0.6);
        end
    else
    end
else
end

set(gca,'xtick',[],'ytick',[],'ztick',[])
axis equal
view([-8 24])
box off
axis tight
set(gca,'visible','off')

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

ax = gca;
ax.XLim = [min_lims(1) max_lims(1)];
ax.YLim = [min_lims(2) max_lims(2)];
ax.ZLim = [min_lims(3) max_lims(3)];

axis vis3d
axis off
ax.Position = ax.OuterPosition - ax.TightInset...
    *[-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];

if print_figure
    print(['porphyrin',...
        num2str(show_connectivity),...
        num2str(show_heads),...
        num2str(show_spectator_atoms),...
        num2str(show_numbering),...
        num2str(show_hydrogens),...
        '.png'],'-dpng', '-r300')
end



