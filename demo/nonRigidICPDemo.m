function nonRigidICPDemo()
    %% Load external functions

    path_backup = path();
    prefix = fileparts(mfilename('fullpath')); % path to current m-file
    addpath([prefix '/..']);
    addpath([prefix '/../src']);


    %% Load Template

    templateFilename = 'template.stl';
    [template_mesh.vertices, template_mesh.faces] = READ_stl(templateFilename);
    [template_mesh.vertices, template_mesh.faces] = patchslim(template_mesh.vertices, template_mesh.faces);


    %% Load Target

    targetFilename = 'target.stl';
    [target_mesh.vertices, target_mesh.faces] = READ_stl(targetFilename);
    [target_mesh.vertices, target_mesh.faces] = patchslim(target_mesh.vertices, target_mesh.faces);


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
    t = 0:15;
    alpha = 1e8 * exp(-t');

    [reconstruction, weights] = nonRigidICP(template_mesh, target_mesh, 'alpha', alpha);

    % Shift the reconstruction in such a manner that it is aligned with the
    % original target data set.
    reconstruction.vertices = bsxfun(@plus, reconstruction.vertices, target_mean);


    %% Write reconstruction to STL-file

    stlwrite('nonRigidICP_Reconstruction.stl', template_mesh.faces, reconstruction.vertices);


    %% Rervert path changes
    path(path_backup);
end

