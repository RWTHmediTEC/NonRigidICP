function nonRigidICPDemo()
    
    %% Load external functions

    path_backup = path();
    prefix = fileparts(mfilename('fullpath')); % path to current m-file
    addpath([prefix '/..']);
    addpath([prefix '/../src']);


    %% Load Template

    templateFilename = 'template.stl';
    templateTriangulation = stlread(templateFilename);
    template_mesh.vertices = templateTriangulation.Points;
    template_mesh.faces = templateTriangulation.ConnectivityList;


    %% Load Target

    targetFilename = 'target.stl';
    targetTriangulation = stlread(targetFilename);
    target_mesh.vertices = targetTriangulation.Points;
    target_mesh.faces = targetTriangulation.ConnectivityList;


    %% Roughly align template and target (if not already done)

    % Here you have to put in some form of manual alignment (matrix
    % transformation) or algorithm that roughly alignes the two datasets.

    template_mesh.vertices = bsxfun(@minus, template_mesh.vertices, mean(template_mesh.vertices));
    target_mean = mean(target_mesh.vertices);
    target_mesh.vertices = bsxfun(@minus, target_mesh.vertices, target_mean);


    %% Plot template and target (to test if they are roughly aligned)

    % myPlot(template_mesh.vertices, template_mesh.faces, target_mesh.vertices, target_mesh.faces, 0.5, 1.0);
    %
    % or
    %
    % plotPolygonMesh(template_mesh, 'Lighting', true)
    % plotPolygonMesh(target_mesh)
    % axis ('equal', 'off')


    %% Reconstruction 

    % Define some stiffness parameters if the default parameters do not fit your needs.
    % alpha = [1e9 1e7 1e5 1000 100 10 1]';
    alpha = 1e8 * exp(-(0:21)');

    [reconstruction, inliers] = nonRigidICP(template_mesh, target_mesh, 'alpha', alpha);

    % Shift the reconstruction in such a manner that it is aligned with the
    % original target data set.
    reconstruction.vertices = bsxfun(@plus, reconstruction.vertices, target_mean);


    %% Write reconstruction to STL-file
    reconstructionTriangulation = triangulation(template_mesh.faces, reconstruction.vertices);
    stlwrite(reconstructionTriangulation, 'nonRigidICP_Reconstruction.stl');


    %% Rervert path changes
    path(path_backup);
end
