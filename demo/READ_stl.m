function [v,f,n] = READ_stl(stlFILENAME)
% READ_stlascii  Read mesh data in the form of an <*.stl> file
%==========================================================================
% FILENAME:          READ_stl.m
% AUTHOR:            Adam H. Aitkenhead
% INSTITUTION:       The Christie NHS Foundation Trust
% CONTACT:           adam.aitkenhead@physics.cr.man.ac.uk
% DATE:              29th March 2010
% PURPOSE:           Read mesh data in the form of an <*.stl> file.
%
% USAGE:
%
%     [v,f,n] = READ_stl(stlFILENAME)
%
% INPUT PARAMETERS:
%
%     stlFILENAME       (String)    filename of the STL file
%
%
% OUTPUT PARAMETERS:
%
%     v = vertices      (N x 3)     x-y-z positions
%     f = faces         (N/3 x 3)   v1-v2-v3 indices
%     n = normals       (N/3 x 3)   per face
%
%==========================================================================

%==========================================================================
% STL ascii file format
%======================
% ASCII STL files have the following structure.  Technically each facet
% could be any 2D shape, but in practice only triangular facets tend to be
% used.  The present code ONLY works for meshes composed of triangular
% facets.
%
% solid object_name
% facet normal x y z
%   outer loop
%     vertex x y z
%     vertex x y z
%     vertex x y z
%   endloop
% endfacet
%
% <Repeat for all facets...>
%
% endsolid object_name
%==========================================================================

%==========================================================================
% STL binary file format
%=======================
% Binary STL files have an 84 byte header followed by 50-byte records, each
% describing a single facet of the mesh.  Technically each facet could be
% any 2D shape, but that would screw up the 50-byte-per-facet structure, so
% in practice only triangular facets are used.  The present code ONLY works
% for meshes composed of triangular facets.
%
% HEADER:
% 80 bytes:  Header text
% 4 bytes:   (int) The number of facets in the STL mesh
%
% DATA:
% 4 bytes:  (float) normal x
% 4 bytes:  (float) normal y
% 4 bytes:  (float) normal z
% 4 bytes:  (float) vertex1 x
% 4 bytes:  (float) vertex1 y
% 4 bytes:  (float) vertex1 z
% 4 bytes:  (float) vertex2 x
% 4 bytes:  (float) vertex2 y
% 4 bytes:  (float) vertex2 z
% 4 bytes:  (float) vertex3 x
% 4 bytes:  (float) vertex3 y
% 4 bytes:  (float) vertex3 z
% 2 bytes:  Padding to make the data for each facet 50-bytes in length
%   ...and repeat for next facet... 
%==========================================================================

stlFORMAT = IDENTIFY_stl_format(stlFILENAME);

%Load the STL file:
if strcmp(stlFORMAT,'ascii')
  [v,f,n] = READ_stlascii(stlFILENAME);
elseif strcmp(stlFORMAT,'binary')
  [v,f,n] = READ_stlbinary(stlFILENAME);
end %if

end %function
%==========================================================================


%==========================================================================
function [stlFORMAT] = IDENTIFY_stl_format(stlFILENAME)
% IDENTIFY_stl_format  Test whether an stl file is ascii or binary

% Open the file:
fidIN = fopen(stlFILENAME);

% Check the file size first, since binary files MUST have a size of 84+(50*n)
fseek(fidIN,0,1);         % Go to the end of the file
fidSIZE = ftell(fidIN);   % Check the size of the file

if rem(fidSIZE-84,50) > 0
    
  stlFORMAT = 'ascii';

else

  % Files with a size of 84+(50*n), might be either ascii or binary...
    
  % Read first 80 characters of the file.
  % For an ASCII file, the data should begin immediately (give or take a few
  % blank lines or spaces) and the first word must be 'solid'.
  % For a binary file, the first 80 characters contains the header.
  % It is bad practice to begin the header of a binary file with the word
  % 'solid', so it can be used to identify whether the file is ASCII or
  % binary.
  fseek(fidIN,0,-1);        % Go to the start of the file
  firsteighty = char(fread(fidIN,80,'uchar')');

  % Trim leading and trailing spaces:
  firsteighty = strtrim(firsteighty);

  % Take the first five remaining characters, and check if these are 'solid':
  firstfive = firsteighty(1:min(5,length(firsteighty)));

  % Double check by reading the last 80 characters of the file.
  % For an ASCII file, the data should end (give or take a few
  % blank lines or spaces) with 'endsolid <object_name>'.
  % If the last 80 characters contains the word 'endsolid' then this
  % confirms that the file is indeed ASCII.
  if strcmp(firstfive,'solid')
  
    fseek(fidIN,-80,1);     % Go to the end of the file minus 80 characters
    lasteighty = char(fread(fidIN,80,'uchar')');
  
    if findstr(lasteighty,'endsolid')
      stlFORMAT = 'ascii';
    else
      stlFORMAT = 'binary';
    end
  
  else
    stlFORMAT = 'binary';
  end
  
end

% Close the file
fclose(fidIN);

end %function
%==========================================================================


%==========================================================================
function [v,f,n] = READ_stlascii(stlFILENAME)
% READ_stlascii  Read mesh data in the form of an ascii <*.stl> file

% Read the ascii STL file
fidIN = fopen(stlFILENAME);
fidCONTENTcell = textscan(fidIN,'%s','delimiter','\n');                  %Read all the file
fidCONTENT = fidCONTENTcell{:}(logical(~strcmp(fidCONTENTcell{:},'')));  %Remove all blank lines
fclose(fidIN);

% Read the vector normals
if nargout >= 2
  stringNORMALS = char(fidCONTENT(logical(strncmp(fidCONTENT,'facet normal',12))));
  n  = str2num(stringNORMALS(:,13:end));
end

% Read the vertex coordinates
facetTOTAL       = sum(strcmp(fidCONTENT,'endfacet'));
stringVERTICES   = char(fidCONTENT(logical(strncmp(fidCONTENT,'vertex',6))));
coordVERTICESall = str2num(stringVERTICES(:,7:end));
cotemp           = zeros(3,facetTOTAL,3);
cotemp(:)        = coordVERTICESall;
coordVERTICES    = shiftdim(cotemp,1);

v = reshape([coordVERTICES(:,:,1),coordVERTICES(:,:,2),coordVERTICES(:,:,3)]',3,[])';
f = reshape(1:size(v,1),3,[])';

end %function
%==========================================================================


%==========================================================================
function [v, f, n] = READ_stlbinary(stlFILENAME)
% This function reads an STL file in binary format into vertex and face
% matrices v and f.
%
% USAGE: [v, f, n,] = stlread(stlFILENAME);
%
% v contains the vertices for all triangles [3*n x 3].
% f contains the vertex lists defining each triangle face [n x 3].
% n contains the normals for each triangle face [n x 3].
% c is optional and contains color rgb data in 5 bits [n x 3].
%
% To see plot the 3D surface use:
%   patch('Faces',f,'Vertices',v,'FaceVertexCData',c);
% or
%   plot3(v(:,1),v(:,2),v(:,3),'.');
%
% Duplicate vertices can be removed using:
%   [v, f]=patchslim(v, f);
%
% For more information see:
%  http://www.esmonde-white.com/home/diversions/matlab-program-for-loading-stl-files
%
% Based on code originally written by:
%    Doron Harlev
% and combined with some code by:
%    Eric C. Johnson, 11-Dec-2008
%    Copyright 1999-2008 The MathWorks, Inc.
% 
% Re-written and optimized by Francis Esmonde-White, May 2010.

fid=fopen(stlFILENAME, 'r'); %Open the file, assumes STL Binary format.

fread(fid,80,'uchar=>schar'); % Read file title
numFaces=fread(fid,1,'int32'); % Read number of Faces

T = fread(fid,inf,'uint8=>uint8'); % read the remaining values
fclose(fid);

% Each facet is 50 bytes
%  - Three single precision values specifying the face normal vector
%  - Three single precision values specifying the first vertex (XYZ)
%  - Three single precision values specifying the second vertex (XYZ)
%  - Three single precision values specifying the third vertex (XYZ)
%  - Two color bytes (possibly zeroed)

% 3 dimensions x 4 bytes x 4 vertices = 48 bytes for triangle vertices
% 2 bytes = color (if color is specified)

trilist = 1:48;

ind = reshape(repmat(50*(0:(numFaces-1)),[48,1]),[1,48*numFaces])+repmat(trilist,[1,numFaces]);
Tri = reshape(typecast(T(ind),'single'),[3,4,numFaces]);

n=squeeze(Tri(:,1,:))';
n=double(n);

v=Tri(:,2:4,:);
v = reshape(v,[3,3*numFaces]);
v = double(v)';

f = reshape(1:3*numFaces,[3,numFaces])';

end
%==========================================================================

