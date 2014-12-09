%% Load exteral functions
path_backup = path();
addpath('..');
addpath('../src');

%% Load Template
templateFilename = 'template.stl';
[mesh_template.vertices, mesh_template.faces] = READ_stl(templateFilename);
[mesh_template.vertices, mesh_template.faces] = patchslim(mesh_template.vertices, mesh_template.faces);

%% Load Target
targetFilename = 'target.stl';
[mesh_target.vertices, mesh_target.faces] = READ_stl(targetFilename);
[mesh_target.vertices, mesh_target.faces] = patchslim(mesh_target.vertices, mesh_target.faces);

%% Roughly align template and target (if not already done)
% here you have to put in some form of manual alignment (matrix transformation)
% or algorithm that roughly alignes the two datasets
mesh_template.vertices = bsxfun(@minus, mesh_template.vertices, mean(mesh_template.vertices));
mesh_target.vertices = bsxfun(@minus, mesh_target.vertices, mean(mesh_target.vertices));

%% Plot template and target (to test if they are roughly aligned)
%myPlot(mesh_template.vertices, mesh_template.faces, mesh_target.vertices, mesh_target.faces, 0.5, 1.0);

%% Reconstruction 
[reconstruction_vertices, weights] = nonRigidICP(mesh_template, mesh_target);

%% Write reconstruction to STL-file
stlwrite('nonRigidICP_Reconstruction.stl', mesh_template.faces, reconstruction_vertices);

%% Redo changes to Path
path(path_backup);
    
