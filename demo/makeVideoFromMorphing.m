function varargout = makeVideoFromMorphing(output_filename, source_mesh, target_mesh, varargin)

% MAKEVIDEOFROMMORPHING   Create a video file from a morphing process
%
%    makeVideoFromMorphing(output_filename, source_mesh, target_mesh, ...)
%    frames = makeVideoFromMorphing(__)
%
%    This function creates a video (file) from the process of non-rigidly
%    aligning two triangulated 3D meshes using the nonRigidICP algorithm.
%    Optional parameters can be passed, all unknown parameters are
%    forwarded to nonRigidICP.
%
%    If 'output_filename' is set to the empty array [] no file is written.
%
%    Additional optional arguments can be passed in the form of pairs with
%    the following meanings:
%
%       'FrameRate'                - Specifies the frame rate of the
%                                    resulting video file. Must be a
%                                    positive scalar. Defaults to 5.
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
%       frames = makeVideoFromMorphing('test.avi', template, target, ...
%                           'alpha', [1e9 1e7 1e5 1e3 100 10 1 0.1 0.01]');
%       save frames.mat frames
%       figure;
%       movie(frames, 1, 2)
%       exit

% Copyright 2015, 2017 Chair of Medical Engineering, RWTH Aachen University
% Written by Christoph Hänisch (haenisch@hia.rwth-aachen.de)
% Version 1.1
% Last changed on 2017-04-27.
% License: Modified BSD License (BSD license with non-military-use clause)

    %% Import external libraries

    importer = LibraryImporter();
    importer.import('Tools')
    importer.addFolder('..');
    importer.addFolder('../src');

    %% Parse the input parameters
    
    assert(ischar(output_filename) || isempty(output_filename), 'Output filename must be a character string or empty.')

    if mod(length(varargin), 2)
        error('Parameter list must have an even length.')
    end

    known_keywords = {'FrameRate', 'SourceMeshOpacity', 'TargetMeshOpacity', 'View'};

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
    addParameter(parser, 'FrameRate', 5, @(x) isscalar(x) && isnumeric(x) && x > 0);
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
    legend({'Target mesh', 'Source mesh'})
    fullscreen
    drawnow
    frames(1) = getframe(gcf);

    nonRigidICP(source_mesh, target_mesh, unknown_args{:}, 'callback', @callback);
    
    %% Write out the movie

    if ~isempty(output_filename)
        video_obj = VideoWriter(output_filename);
        video_obj.FrameRate = parser.Results.FrameRate;
        video_obj.open()
        for i = 1:length(frames)
            video_obj.writeVideo(frames(i).cdata);
        end
        video_obj.close()
    end

    %% Set output arguments accordingly

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
        title({sprintf('Iteration %d, inner iteration %d', data.iteration, data.inner_iteration), ...
               sprintf('alpha %.4f, diff %.4f, inliers %d', data.alpha, data.gradual_change, data.number_of_inliers)})
        drawnow
        frames(end+1) = getframe(gcf);
    end
    
end