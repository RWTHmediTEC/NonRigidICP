function makeVideoFromMorphing(output_filename, source_mesh, target_mesh, varargin)

% MAKEVIDEOFROMMORPHING   Create a video file from a morphing process
%
%    MAKEVIDEOFROMMORPHING(output_filename, source_mesh, target_mesh, ...)
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

% Copyright 2015 Chair of Medical Engineering, RWTH Aachen University
% Written by Christoph Hänisch (haenisch@hia.rwth-aachen.de)
% Version 1.0
% Last changed on 2015-08-06.
% License: Modified BSD License (BSD license with non-military-use clause)

    %% Parse the input parameters

    if mod(length(varargin), 2)
        error('Parameter list must have an even legth.')
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


    %% Load exteral functions
    path_backup = path();
    prefix = fileparts(mfilename('fullpath')); % path to current m-file
    addpath([prefix '/..']);
    addpath([prefix '/../src']);


    %% Perform the morphing and save store the obtained geometries

    figure
    patch_handle_target_shape = plotPolygonMesh(target_mesh, 'Lighting', true);
    patch_handle_target_shape.FaceAlpha = parser.Results.TargetMeshOpacity;
    axis off equal
    view(parser.Results.View)
    ax = gca;
    ax.NextPlot = 'replaceChildren';

    frames = [];

    nonRigidICP(source_mesh, target_mesh, unknown_args, 'callback', @callback);
    
    fig = figure;
    movie(fig, frames, 4)


    %% Rervert path changes
    path(path_backup);

    return


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

        patch_handle_source_shape = plotPolygonMesh(target_shape, 'color', 'g');
        patch_handle_source_shape.FaceAlpha = parser.Results.SourceMeshOpacity;
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