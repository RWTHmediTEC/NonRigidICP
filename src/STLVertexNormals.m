function [vn]=STLVertexNormals(f, v, fn)

% STLVertexNormals
%
% The function STLVertexNormals takes in a tesselated mesh model in the
% faces, vertices format, and calculates the per-vertex normals. The
% default method is to calculate face-weighted per-vertex normals by
% averaging the non-normalized cross-product for each triangle that touches
% a vertex. The slightly faster (but less pretty at sharp edges) method is to average the
% precomputed per-face normals passed as an argument.
%
% Usage:
%   [vn]=STLVertexNormals(f, v, fn);
%
% Inputs:
%   f, the triangle faces in the model
%   v, the vertices for each triangle in the model
%   fn, the face-normals for each triangle in the model (optional input)
%
% If fn is specified, the routine should be slightly faster, but less
% accurate.
%
% Outputs:
%   vn, the vertex-normals
%
% Francis Esmonde-White, 11 October, 2011, francis@esmonde-white.com

% mostly based on:
% http://www.devmaster.net/forums/showthread.php?t=414

if ~exist('fn','var')
    a = v(f(:,1),:);
    b = v(f(:,2),:);
    c = v(f(:,3),:);
    fn = cross((b-a),(c-a)); % un-normalized face normals
%     fn = fn./repmat(sqrt(sum(fn.^2,2)),[1,3]); % normalized face normals
end

vn = zeros(size(v)); % preallocate for vertex normals
vn(f(:,1),:) = vn(f(:,1),:) + fn; % add the normals for the first point in each triangle
vn(f(:,2),:) = vn(f(:,2),:) + fn; % add the normals for the second point in each triangle
vn(f(:,3),:) = vn(f(:,3),:) + fn; % add the normals for the third point in each triangle

vn = vn./repmat(sqrt(sum(vn.^2,2)),[1,3]); % normalized vertex normals
% Based on:
% http://stackoverflow.com/questions/1061276/how-to-normalize-a-vector-in-matlab-efficiently-any-related-built-in-function
