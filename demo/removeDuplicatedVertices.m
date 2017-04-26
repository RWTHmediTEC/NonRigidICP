function simplified_mesh = removeDuplicatedVertices(mesh, order)

% REMOVEDUPLICATEDVERTICES   Removes duplicated vertices in surface meshes
% 
%    simplified_mesh = removeDuplicatedVertices(mesh)
%    simplified_mesh = removeDuplicatedVertices(mesh, order)
%
%    This function finds and removes duplicated vertices in a surface mesh.
%    The passed mesh must be a struct of the following form:
%
%       mesh.vertices     The vertices are given as vertically stacked rows
%                         (i.e. an n-by-3 matrix) where each row denotes a
%                         vertex location [x y z] in space.
%
%       mesh.faces        The faces are also given in the form of
%                         vertically stacked rows. Each row contains three
%                         integer indeces [id1 id2 id3] referencing the
%                         respective vertices in 'target.vertices'.
%
%       mesh.normals      Face normals are given as an m-by-3 matrix of
%                         type double. This field is optionally.
%
%    The second parameter is optional and specifies the order of the
%    vertices in the reduced set. A value of 'ascending' means that the
%    rows are sorted ascending. A value of 'stable' means that the original
%    order is retained. If not specified, order defaults to 'stable'.

% Written by Christoph Hänisch (haenisch@hia.rwth-aachen.de)
% Version 1.1
% Last changed on 2015-08-07
% Licence: Modified BSD License (BSD license with non-military-use clause)

    %% Check input parameter

    if ~isstruct(mesh)
        error('Input mesh has a wrong format.')
    end

    if ~isfield(mesh, 'vertices')
        error('Input mesh does not provide field variable ''vertices''.');
    end

    if ~isfield(mesh, 'faces')
        error('Input mesh does not provide field variable ''faces''.');
    end

    if nargin < 2
        order = 'stable';
    else
        if ~any(strcmp(order, {'ascending', 'stable'}))
            error('Order must be ''ascending'' or ''stable''.')
        end
    end

    %% Simplify mesh

    switch order
        case 'ascending'
            [simplified_mesh.vertices, ~, index_mapping] = unique(mesh.vertices, 'rows', 'sorted');
        case 'stable'
            [simplified_mesh.vertices, ~, index_mapping] = unique(mesh.vertices, 'rows', 'stable');
    end
    simplified_mesh.faces = index_mapping(mesh.faces);
    
    if isfield(mesh, 'normals')
        simplified_mesh.normals = mesh.normals;
    end

end