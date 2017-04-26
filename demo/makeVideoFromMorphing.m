function varargout = makeVideoFromMorphing(output_filename, source_mesh, target_mesh, varargin)

% MAKEVIDEOFROMMORPHING   Create a video file from a morphing process
%
%    makeVideoFromMorphing(output_filename, source_mesh, target_mesh, ...)
%    frames = makeVideoFromMorphing(__)
%
%    This function creates a video file from the process of non-rigidly
%    aligning two triangulated 3D meshes using the nonRigidICP algorithm.
%    Optional parameters can be passed, all unknown parameters are
%    forwarded to nonRigidICP.
%
%    Additional optional arguments can be passed in the form of pairs with
%    the following meanings:
%
%       'SourceMeshOpacity'        - Opacity of the source mesh's faces
%                                    when displayed. Must be a value
%                                    between 0 and 1. Defaults to 1.
%
%       'TargetMeshOpacity'        - Opacity of the target mesh's faces
%                                    when displayed. Must be a value
%                                    between 0 and 1. Defaults to 0.5.
%
%       'View'                     - Viewing angles [az, el] for a
%                                    three-dimensional plot. Defaults to
%                                    a value of [10 30].
%
%    Example
%
%       importer = LibraryImporter();
%       importer.import('MeshProcessing')
%       importer.import('StatisticalShapeModel')
%       template = shiftMeshToOrigin(loadSTL('template.stl'));
%       target = shiftMeshToOrigin(loadSTL('target.stl'));
%       frames = makeVideoFromMorphing('test.avi', template, target);
%       figure;
%       movie(fig, frames, 1, 2)

% Copyright 2015, 2017 Chair of Medical Engineering, RWTH Aachen University
% Written by Christoph Hänisch (haenisch@hia.rwth-aachen.de)
% Version 1.1
% Last changed on 2017-04-26.
% License: Modified BSD License (BSD license with non-military-use clause)

    %% Import external libraries

    importer = LibraryImporter();
    importer.addFolder('..');
    importer.addFolder('../src');

    %% Parse the input parameters

    if mod(length(varargin), 2)
        error('Parameter list must have an even length.')
    end

    known_keywords = {'SourceMeshOpacity', 'TargetMeshOpacity', 'View'};

    input_args = {};
    unknown_args = {};

    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case lower(known_keywords)
                input_args = [input_args, varargin(i), varargin(i+1)]; %#ok<AGROW>
            otherwise
                unknown_args = [unknown_args, varargin(i), varargin(i+1)]; %#ok<AGROW>
        end
    end

    parser = inputParser;
    addParameter(parser, 'SourceMeshOpacity', 1.0, @(x)validateattributes(x,{'numeric'},{'>=',0},{'<=',1}, 'makeVideoFromMorphing'));
    addParameter(parser, 'TargetMeshOpacity', 0.5, @(x)validateattributes(x,{'numeric'},{'>=',0},{'<=',1}, 'makeVideoFromMorphing'));
    addParameter(parser, 'View', [10 30]);
    parse(parser, input_args{:});

    %% Perform the morphing and save store the obtained geometries

    figure
    hold on
    patch_handle_target_mesh = plotPolygonMesh(target_mesh, 'Lighting', true, 'ColorScheme', 'Transparent');
    patch_handle_target_mesh.FaceAlpha = parser.Results.TargetMeshOpacity;
    patch_handle_source_mesh = plotPolygonMesh(source_mesh, 'Color', 'g', 'ColorScheme', 'Transparent');
    patch_handle_source_mesh.FaceAlpha = parser.Results.SourceMeshOpacity;
    axis off equal
    view(parser.Results.View)
    title('Initial situtation')
    drawnow
    frames(1) = getframe(gcf);

    nonRigidICP(source_mesh, target_mesh, unknown_args{:}, 'callback', @callback);
    
    % fig = figure;
    % movie(fig, frames, 4)

    if nargout == 1
        varargout{1} = frames;
    end

    return


    %% Helper functions

    function callback(data)
        % Field names of data:
        % - alpha
        % - computation_time
        % - gradual_change
        % - inner_iteration
        % - iteration
        % - number_of_inliers
        % - overall_computation_time
        % - transformed_mesh

        delete(patch_handle_source_mesh)
        patch_handle_source_mesh = plotPolygonMesh(data.transformed_mesh, 'Color', 'g', 'ColorScheme', 'Transparent');
        patch_handle_source_mesh.FaceAlpha = parser.Results.SourceMeshOpacity;
        title(sprintf('Iteration %d, inner iteration %d, alpha %f, diff %f, inliers %d', ...
                      data.iteration, ...
                      data.inner_iteration, ...
                      data.alpha, ...
                      data.gradual_change, ...
                      data.number_of_inliers))
        drawnow
        frames(end+1) = getframe(gcf);
    end
    
end