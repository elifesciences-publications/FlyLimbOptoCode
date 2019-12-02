% plotUMAPFigures.m: This script reads in and plots the UMAP embedding data

% Read in the .csv file
% data_path = '/Users/bdeangelis/Desktop/Datasets/OptoMethodDatasets/csv/UMAP_embedding_limb-data_1000-neighbors_3-dims.csv';
data_path = '/Users/bdeangelis/Desktop/Datasets/OptoMethodDatasets/csv/UMAP_embedding_body-data_1000-neighbors_3-dims.csv';
% data_path = '/Users/bdeangelis/Desktop/Datasets/OptoMethodDatasets/csv/UMAP_embedding_all-data_1000-neighbors_3-dims.csv';
data = readtable(data_path);

% Create new columns with the appropriate units for the post-stimulus
% NOTE: Samples were taken at 150 fps. To get the cumulative displacements
% we need to multiply by 1/150.
data.post_v_perp_cum_mm = data.post_v_perp * (1/150);
data.post_v_rot_cum_deg = data.post_v_rot * (1/150) * (180/pi); %Converts from radians to degrees

% Define some colormaps
[ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningPaperColormaps();

% Define some plot limits
centroidLimits = [-8 8];

% TODO: Set the final view
% Set the viewing angle for all figures
viewport = [-30 20];
viewport_side = [26 10];

% Filter the centroid data based on a variable of interest
 data = data(data.pre_v_par > 0.5,:);

% Pull out the centroid data
centroidUMAP = [data.x, data.y, data.z];

%% Create 3D plots of the data - Body Velocity Plots

% Figure 4A
% Colored by prestimulus mean forward walking speed
MakeFigure;
scattern(centroidUMAP, 10, data.pre_v_par, 'filled');
axis('equal');
caxis([0 30]); % Update
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'v_{||}(t<0) (mm/s)');
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);
zlim(centroidLimits);
view(viewport);

% Figure 4B
% Colored by post-stimulus mean forward walking speed
MakeFigure;
scattern(centroidUMAP, 10, data.post_v_par, 'filled');
axis('equal');
caxis([0 30]); % Update
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'v_{||}(t>0) (mm/s)');
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);
zlim(centroidLimits);
view(viewport);

% Figure 4C
% Colored by post-stimulus cumulative lateral displacement (View 1) 
MakeFigure;
scattern(centroidUMAP, 10, data.post_v_perp_cum_mm, 'filled');
axis('equal');
caxis([-.75 .75]);
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'Cumulative Lateral Displacement (t>0) (mm)');
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);
zlim(centroidLimits);
view(viewport);


% Figure 4D
% Colored by post-stimulus cumulative lateral displacement (View 2) 
MakeFigure;
scattern(centroidUMAP, 10, data.post_v_perp_cum_mm, 'filled');
axis('equal');
caxis([-.75 .75]);
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'Cumulative Lateral Displacement (t>0) (mm)');
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);
zlim(centroidLimits);
view(viewport_side);

% Figure 4E
% Colored by post-stimulus cumulative rotation (View 1)
MakeFigure;
scattern(centroidUMAP, 10, data.post_v_rot_cum_deg, 'filled');
axis('equal');
caxis([-40 40]);
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'Cumulative Rotation (t>0) (Degrees)');
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);
zlim(centroidLimits);
view(viewport); 

% Figure 4F
% Colored by post-stimulus cumulative rotation (View 2)
MakeFigure;
scattern(centroidUMAP, 10, data.post_v_rot_cum_deg, 'filled');
axis('equal');
caxis([-40 40]);
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'Cumulative Rotation (t>0) (Degrees)');
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);
zlim(centroidLimits);
view(viewport_side); 

%% Comparing experimental conditions

% Define all the hit types
conds = unique(data.Condition);

% Figure 4G
MakeFigure;
for i = 1:length(conds)
    hold on;
    if strcmp(conds{i}, 'Bristle-Chrimson')
        curUMAP = centroidUMAP(strcmp(data.Condition, conds{i}),:);
        h = scattern(curUMAP, 100, 'b', 'marker', '.');
        set(h, 'MarkerEdgeAlpha', .2, 'MarkerFaceAlpha', .2)
        
    elseif strcmp(conds{i}, 'Gr5a-Chrimson')
        curUMAP = centroidUMAP(strcmp(data.Condition, conds{i}),:);
        h = scattern(curUMAP, 100, 'r', 'marker', '.');
        set(h, 'MarkerEdgeAlpha', .2, 'MarkerFaceAlpha', .2)
                
    else
        curUMAP = centroidUMAP(strcmp(data.Condition,conds{i}),:);
        h = scattern(curUMAP, 100, [.5,.5,.5], 'marker', '.');
        set(h, 'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha', .1)

    end
end
axis('equal');    
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);
zlim(centroidLimits);
view(viewport);
grid on;
legend(conds)

% Figure 4H
MakeFigure;
for i = 1:length(conds)
    hold on;
    if strcmp(conds{i}, 'Bristle-Chrimson')
        curUMAP = centroidUMAP(strcmp(data.Condition, conds{i}),:);
        h = scattern(curUMAP, 100, 'b', 'marker', '.');
        set(h, 'MarkerEdgeAlpha', .2, 'MarkerFaceAlpha', .2)
        
    elseif strcmp(conds{i}, 'Gr5a-Chrimson')
        curUMAP = centroidUMAP(strcmp(data.Condition, conds{i}),:);
        h = scattern(curUMAP, 100, 'r', 'marker', '.');
        set(h, 'MarkerEdgeAlpha', .2, 'MarkerFaceAlpha', .2)
                
    else
        curUMAP = centroidUMAP(strcmp(data.Condition,conds{i}),:);
        h = scattern(curUMAP, 100, [.5,.5,.5], 'marker', '.');
        set(h, 'MarkerEdgeAlpha', .1, 'MarkerFaceAlpha', .1)

    end
end
axis('equal');    
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);
zlim(centroidLimits);
view(viewport_side);
grid on;
legend(conds)