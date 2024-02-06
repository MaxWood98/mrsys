%Build Full Faces Function 

%Version = 0.1
%Updated = 06-12-23

%Max Wood 2023
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [mesh] = build_full_faces(mesh)
    mesh.faces_full = zeros(mesh.nface,mesh.max_vtx_inface);
    for ii=1:mesh.nface
        mesh.faces_full(ii,1:mesh.faces{ii}.nvertex) = mesh.faces{ii}.vertices(:);
    end
    mesh.faces_full(mesh.faces_full == 0) = nan;
end