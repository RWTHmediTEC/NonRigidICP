function [transformed_mesh, inliers, accumulated_transformations] = nonRigidICP(template_mesh, target_mesh, varargin)

% NONRIGIDICP   Non-rigid alignment of two meshes (e.g. from STL files)
%
%    NONRIGIDICP(template_mesh, target_mesh)
%    NONRIGIDICP(template_mesh, target_mesh, ...)
%
%    This function non-rigidly aligns two triangulated 3D meshes and is
%    based on the paper: "Optimal Step Nonrigid ICP Algorithm for Surface
%    Registration" by B. Amberg, S. Romdhani and T. Vetter from 2007.
%
%    The positional input arguments 'template_mesh' and 'target_mesh' must
%    be structs having the non-empty fields 'vertices' and 'faces'. Face
%    normals can be optionally given in the field 'normals'. The structure
%    must obey the following form:
%
%       'mesh.vertices'        - 3D points in space [n x 3 double]
%       'mesh.faces'           - indeces to the resp. vertices [m x 3 int]
%       'mesh.normals'         - 3D vectors [m x 3 double]
%
%    Note: The vertices must be unique such that faces sharing the same
%    vertex must hold the same index. Use a function like patchslim to
%    adapt your input data accordingly.
%
%    Additional optional arguments can be passed in the form of pairs with
%    the following meanings:
%
%       'alpha'                    - Array of stiffness parameters.
%                                    Defaults to an exponentially
%                                    descending curve (see the internal
%                                    function getAlpha() for more details).
%                                    The values are highly application
%                                    specific and may range from intervals
%                                    like [1; 1e8] to [1e-5; 100].
%
%       'callback'                 - Function handle to callback function.
%                                    This function is called after each
%                                    iteration and the following parameters
%                                    are passed: deformed template mesh,
%                                    current iteration, current alpha
%                                    value.
%
%       'epsilon'                  - Convergence parameter. Has to be
%                                    positive and defaults to a value of
%                                    0.0002.
%
%       'initial_transformation'   - 3-by-4 transformation matrix that
%                                    roughly aligns the two meshes (the
%                                    template vertices are transformed).
%                                    Not needed, if the two input meshes
%                                    are already roughly aligned. Defaults
%                                    to the matrix [eye(3), [0; 0; 0]].
%
%       'max_dist'                 - Maximum vertex-to-surface-distance for
%                                    a nearest neighbor search to classify
%                                    a point as an inlier. Has to be
%                                    positive and defaults to a value of 20
%                                    units.
%
%       'max_iter'                 - Maximal number of iterations per alpha
%                                    step. Has to be positive and defaults
%                                    to a value of 100 units.
%
%       'max_normal_diff'          - Maximum normal difference in degrees
%                                    between vertex and surface to classify
%                                    a point as an inlier. Has to be
%                                    positive and defaults to a value of 22
%                                    degrees.
%
%       'verbosity'                - Denotes the amount of informational
%                                    text put out by this function. The
%                                    higher the value the more information
%                                    is put out. Currently, only values of
%                                    zero and one are supported. The
%                                    default value is one.
%
%       'weights'                  - A column vector filled with weights
%                                    corresponding to the n vertices of the
%                                    mesh. The higher the weight, the more
%                                    the particular vertex contributes to
%                                    the final morphing result. A value of
%                                    zero effectively regards a vertex as
%                                    outlier. Defaults to ones(n, 1).
%
%
%    Outputs:
%
%       'transformed_mesh'             - aligned template mesh
%
%       'inliers'                      - template vertex weights of the
%                                        last iteration. A value of zero
%                                        denotes an outlier, a value of one
%                                        denotes an inlier.
%
%       'accumulated_transformations'  - transformation which transforms
%                                        the original template_mesh such
%                                        that it is aligned with the
%                                        target_mesh
%
%
%    Example:
%
%       % Load in the data.
%       [template.vertices, template.faces] = READ_stl('template.stl');
%       [target.vertices, target.faces] = READ_stl('target.stl');
%
%       % Remove the duplicates (mandatory).
%       [template.vertices, template.faces] = ...
%           patchslim(template.vertices, template.faces);
%       [target.vertices, target.faces] =  ...
%           patchslim(target.vertices, target.faces);
%
%       % Perform the fitting.
%       transformed_target = nonRigidICP(template, target);
%
%
%    See also patchslim (located on Matlab File Exchange)

% Copyright 2014, 2015 Chair of Medical Engineering, RWTH Aachen University
% Written by Erik Noorman and Christoph Hänisch (haenisch@hia.rwth-aachen.de)
% Version 1.4.1
% Last changed on 2015-07-20.
% License: Modified BSD License (BSD license with non-military-use clause)

    %% Parse the input parameters
    parser = inputParser;
    addParameter(parser, 'alpha', getAlpha(), @(x) validateattributes(x,{'numeric'},{'>',0},'nonRigidICP'));
    addParameter(parser, 'callback', [], @(x) isempty(x) || isa(x, 'function_handle'));
    addParameter(parser, 'epsilon', 0.0002, @(x) validateattributes(x,{'numeric'},{'>',0},'nonRigidICP'));
    addParameter(parser, 'initial_transformation', [eye(3), [0;0;0]], @(x) validateattributes(x,{'numeric'},{'size',[3,4]}));
    addParameter(parser, 'max_dist', 20, @(x) validateattributes(x,{'numeric'},{'>',0},'nonRigidICP'));
    addParameter(parser, 'max_iter', 100, @(x) validateattributes(x,{'numeric'},{'>',0},'nonRigidICP'));
    addParameter(parser, 'max_normal_diff', 22, @(x) validateattributes(x,{'numeric'},{'>',0},'nonRigidICP'));
    addParameter(parser, 'weights', [], @(x)  validateattributes(x,{'numeric'},{'size',[NaN,1]}));
    addParameter(parser, 'verbosity', 1, @(x) validateattributes(x,{'numeric'},{'>=',0},'nonRigidICP'));
    parse(parser, varargin{:});

    %% Load exteral functions
    path_backup = path();
    prefix = fileparts(mfilename('fullpath')); % path to current m-file
    addpath([prefix '/src']);

    %% Initialise variables
    n = size(template_mesh.vertices, 1); % # template vertices

    alpha = parser.Results.alpha; % Stiffness parameter list
    callback = parser.Results.callback;
    epsilon = parser.Results.epsilon;
    initial_transformation = parser.Results.initial_transformation;
    max_dist = parser.Results.max_dist; % distance in mm
    max_iter = parser.Results.max_iter; % maximal # iteration for fixed stiffness
    max_normal_diff = parser.Results.max_normal_diff; % angle in degrees
    verbosity_level = parser.Results.verbosity;
    weights = parser.Results.weights;
    
    if isempty(weights)
        weights = ones(n,1);
    end

    % Initialize accumulated transformations and template_mesh.vertices
    if ~isequal(initial_transformation, [eye(3), [0;0;0]])
        % Transform the template_mesh.vertices with the initial transformation matrix.
        template_vertices_homogeneous = [template_mesh.vertices, ones(n, 1)];
        transformed_vertices_homogeneous = template_vertices_homogeneous * initial_transformation';
        template_mesh.vertices = transformed_vertices_homogeneous(:, 1:3);

        % If the initial transformation implies a reflection flip the face
        % normals such that they point to the former direction (e.g. outwards).
        if sign(det(initial_transformation(1:3,1:3))) == -1
            template_mesh.faces = [template_mesh.faces(:,1), template_mesh.faces(:,3), template_mesh.faces(:,2)];
        end
    end
    X = repmat(initial_transformation', n, 1);
    accumulated_transformations = X;

    %% TESTING: meshes are roughly aligned? (uncomment following line for testing)
    %myPlot(template_mesh.vertices, template_mesh.faces, target_mesh.vertices, target_mesh.faces, 0.5, 1.0); pause();

    % Create structures for efficient point2planeDist calculation
    kd_tree_target = KDTreeSearcher(target_mesh.vertices, 'distance', 'euclidean'); % vertex distance KD-Tree
    v2f_target = buildVerticesToFacesList(target_mesh.vertices, target_mesh.faces); % vertex-face mapping
    if isfield(target_mesh, 'normals')
        target_mesh.vertex_normals = STLVertexNormals(target_mesh.faces, target_mesh.vertices, target_mesh.normals);
    else
        target_mesh.vertex_normals = STLVertexNormals(target_mesh.faces, target_mesh.vertices);% vertex normals
    end

    v2f_template = buildVerticesToFacesList(template_mesh.vertices, template_mesh.faces); % vertex-face mapping - needed to calculate (area) vertex weights
    e = buildEdgeList(template_mesh.faces); % the cost function needs a template edge list
    m = size(e, 1); % # template edges
    i = reshape([1:m; 1:m], 2*m, 1);
    j = reshape(e', 2*m, 1);
    s = repmat([-1; 1], m, 1);
    M_init = sparse(i, j, s); % [m x n] (M_init with (-1|1) for every edge)

    %% Loop over stiffness alpha^i element {a^1,...,a^k}, a^i > a^(i+1)
    if verbosity_level > 0
        fprintf('\niter\talpha\tdiff\tinlier\ttime\n');
    end

    stop_the_outermost_loop = false; % used to stop the outer loop if there is virtually no morphing in any of the computing steps

    for a = alpha'
    
        if stop_the_outermost_loop
            break;
        end

        %% TESTING: uncomment following line for ploting the current meshes
        %myPlot(template_mesh.vertices, template_mesh.faces, target_mesh.vertices, target_mesh.faces, 0.5, 1.0); pause();

        % do while |X^(i) - X^(i-1)| < epsilon
        for iter = 1:max_iter
            t1 = tic;

            X_old = X;
            % Find preliminary correspondences
            template_mesh.vertex_normals = STLVertexNormals(template_mesh.faces, template_mesh.vertices); % compute vertex normals (face area weighted)

            [U, inliers, ~] = closestPointsOnSurface(...
                template_mesh, target_mesh, ...
                v2f_target, kd_tree_target, ...
                max_dist, max_normal_diff);

            %% Build Error Matrix and find solution X = A/B and transform model
            % Construct sparse matrices: M, G, W, D, A, B
            M = a * M_init; % [m x n] ((-1|1) for every edge else 0)

            gamma = 1; % Perhaps set gamma to meaningful value e.g. (1/)radius of data (tested: seems to make no difference)
            G = diag([1,1,1,gamma]);

            vWeight = weights.*inliers.*computeVertexWeight(template_mesh, v2f_template);
            W = sparse(1:n, 1:n, vWeight); % [n x n]
            %W = sparse(1:n, 1:n, w); % [n x n]

            % From the vertices (x_i, y_i, z_i) create the matrix
            %
            %     | x1 y1 y2 1  0  0  0  0  0  0  0  0  ... |
            % D = | 0  0  0  0  x2 y2 y3 1  0  0  0  0  ... |
            %     | 0  0  0  0  0  0  0  0  x3 y3 z3 1  ... |
            %     | ...                                     |
            template_vertices_homogeneous = [template_mesh.vertices, ones(n,1)];
            vec = @(matrix) reshape(matrix, [], 1); % stacks the columns of 'matrix' on top of each other resulting in a column vector
            i = kron((1:n), ones(1, 4)); % i = [1 1 1 1 2 2 2 2 3 3 3 3 ...];
            j = 1:4*n; % j = [1 2 3 4 5 ... 4*n];
            D = sparse(i, j, vec(template_vertices_homogeneous')); % n x 4n

            % kron(M,G) [4m   x 4n]
            % W*D       [n    x 4n]
            % A         [4m+n x 4n]
            A = [kron(M,G); W*D];

            % 0    [4m   x 3]
            % W*U  [n    x 3]
            % B    [4m+n x 3]
            B = [sparse(4*m,3); W*U];

            % Solve linear system
            X = A\B;

            % Transform template points
            template_mesh.vertices = full(D*X);

            %%
            % Track point alignments

            XHom = [X, repmat([0;0;0;1], n, 1)];
            sX = reshape(XHom', 16*n, 1);
            accuXHom = [accumulated_transformations, repmat([0;0;0;1], n, 1)];
            sXaccu = reshape(accuXHom', 16*n, 1);

            i = reshape(repmat(reshape(1:4*n, 4, n), 4, 1), 16*n, 1); % 123412341234123456785678...
            j = reshape(repmat(1:4*n, 4, 1), 16*n, 1); % 111122223333...
            XBigMat = sparse(i, j, sX); % n x 4n
            accuXBigMat = sparse(i, j, sXaccu); % n x 4n

            accumulated_transformations = repmat(eye(4), 1, n) * (XBigMat * accuXBigMat);
            accumulated_transformations = accumulated_transformations(1:3,:)'; % throw away every homogenous stuff [0 0 0 1]

            % Compute new do-while condition
            normXdiff = norm(full(X_old-X))/n;

            if verbosity_level > 0
                fprintf('%d:\t%.0f\t%.4f\t%d\t%.1f\n', iter, a, normXdiff, sum(inliers), toc(t1)); % print iter # and time
            end
            
            % Invoke callback function

            if ~isempty(callback)
                transformed_mesh.vertices = template_mesh.vertices;
                transformed_mesh.faces = template_mesh.faces;
                callback(transformed_mesh, iter, a);
            end

            % do-while condition still holds?
            if normXdiff < epsilon/10000
                % The change rate is so marginal that we can break out from *all* iterations
                stop_the_outermost_loop = true;
                break;
            elseif normXdiff < epsilon
                break;
            end
        end
    end

    transformed_mesh.vertices = template_mesh.vertices;
    transformed_mesh.faces = template_mesh.faces;

    %% Revert path changes
    path(path_backup);
end


function alpha = getAlpha()
    % Stiffness parameter alpha choosen as in:
    % N. Hasler, C.Stoll, M. Sunkel, B. Rosenhahn, and H.-P. Seidel
    % 'A Statistical Model of Human Pose and Body Shape' 2009
    tMax = 15;
    k0 = 3000;
    kInf = 0.01;
    lambda = log(k0/kInf)/tMax;
    alpha = zeros(tMax,1);
    for t = 1:tMax
        alpha(tMax-t+1) = k0 * exp(lambda * t);
    end
end


function [closest_points, inliers, distances] = closestPointsOnSurface(source_mesh, target_mesh, v2f_target, kd_tree_target, max_dist, max_normal_diff)
    % This function finds for every vertex in the source mesh the closest
    % point on the surface of target mesh and sets the values of the
    % inliers vector to 'false' if one of the following condition holds:
    % 1. The normal differenz is higher than the value 'max_normal_diff'
    % 2. The point-to-plane distance is higher than 'max_dist'
    % 3. The closest target vertex is a border vertex

    % Initialize closest_points with closest vertices
    [closest_points_idx, distances] = knnsearch(kd_tree_target, source_mesh.vertices); % nearest neighbor search
    closest_points = target_mesh.vertices(closest_points_idx, :);

    % Set weights to 'false' if point normal differences are too high
    normal_diff = acos(dot(source_mesh.vertex_normals, target_mesh.vertex_normals(closest_points_idx,:), 2)) / pi * 180; % in degrees
    inliers = normal_diff < max_normal_diff;

    % For each vertex in the source mesh find the closest surface point (on
    % those triangles that touch the closest vertex of the target mesh) 
    for i = 1:size(source_mesh.vertices, 1)
        if inliers(i)
            vertex = source_mesh.vertices(i,:);
            closest_faces_idx = v2f_target{closest_points_idx(i)}; % all triangles that include the closest vertex

            % Test if vertex lies on border
            number_faces = size(closest_faces_idx, 2);
            number_vertices = size(unique(reshape(target_mesh.faces(closest_faces_idx,:), number_faces*3, 1)), 1);
            if number_faces+1 ~= number_vertices
                inliers(i) = false;
                continue;
            end

            for closest_face_idx = closest_faces_idx
                triangle_idx = target_mesh.faces(closest_face_idx,:);
                tri = [target_mesh.vertices(triangle_idx(1),:); ...
                       target_mesh.vertices(triangle_idx(2),:); ...
                       target_mesh.vertices(triangle_idx(3),:)];
                % pointTriangleDistance: yields similar values compared to CloudCompare but slightly bigger
                [current_distance, closest_point] = pointTriangleDistance(tri, vertex); % calc NN to triangle distance
                if distances(i) > current_distance % keep the smallest distance
                    distances(i) = current_distance;
                    closest_points(i, :) = closest_point';
                end
            end
        end
    end

    inliers = inliers & distances < max_dist; % check if closest points are not to far away
end


function v2f = buildVerticesToFacesList(v, f)
    % Builds a mapping for every vertex to all faces that incorporate that vertex
    v2f = cell(size(v,1),1);
    for i=1:size(f,1) % every face
        v2f{f(i,1)} = [v2f{f(i,1)}, i];
        v2f{f(i,2)} = [v2f{f(i,2)}, i];
        v2f{f(i,3)} = [v2f{f(i,3)}, i];
    end
end


function e = buildEdgeList(f)
    % Build edge list from faces
    m = size(f,1);
    e = zeros(m*3,2);

    idx=1;
    for i=1:m
        face = sort(f(i,:));
        e(idx,:) = face(1,[1,2]);
        e(idx+1,:) = face(1,[1,3]);
        e(idx+2,:) = face(1,[2,3]);
        idx = idx + 3;
    end
    e = unique(e,'rows');
end


function vWeight = computeVertexWeight(mesh, v2f)
    n = size(mesh.vertices, 1); % #vertices
    vWeight = zeros(n, 1);

    m = size(mesh.faces, 1); % #faces
    A = zeros(m, 3); % 1. vertices of faces
    B = zeros(m, 3); % 2. vertices of faces
    C = zeros(m, 3); % 3. vertices of faces

    for i=1:m
        A(i,:) = mesh.vertices(mesh.faces(i,1),:);
        B(i,:) = mesh.vertices(mesh.faces(i,2),:);
        C(i,:) = mesh.vertices(mesh.faces(i,3),:);
    end
    % Norm(cross( a-c, b-c )) / 2 % triangle surface in 3D (a,b,c vert of tria)
    fWeight = sqrt(sum(abs(cross( A-C, B-C )).^2,2))./2; % surface area per face

    for i=1:n
        vWeight(i) = sum(fWeight(v2f{i}, 1));
    end
end
