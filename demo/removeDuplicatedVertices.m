function simplified_mesh = removeDuplicatedVertices(mesh)

% REMOVEDUPLICATEDVERTICES   Removes duplicated vertices in surface meshes
% 
%    simplified_mesh = removeDuplicatedVertices(mesh)
%
%    This function finds and removes duplicated vertices in a surface mesh.
%    The passed mesh must be and the returned mesh is a struct of the
%    following form:
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
%                         type double.

% Written by Christoph Hänisch (haenisch@hia.rwth-aachen.de)
% Version 1.0
% Last changed on 2015-04-27
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

    %% Simplify mesh

    [simplified_mesh.vertices, ~, index_mapping] = unique(mesh.vertices, 'rows');
    simplified_mesh.faces = index_mapping(mesh.faces);
    
    if isfield(mesh, 'normals')
        simplified_mesh.normals = mesh.normals;
    end

end