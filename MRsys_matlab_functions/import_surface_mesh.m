%Load Surface Mesh Function 

%Version = 0.2
%Updated = 06-12-23

%Max Wood 2023
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [mesh] = import_surface_mesh(filename)

    %Open file
    fid = fopen(filename,'r');
    
    %Read quantities 
    itemdat = textscan(fid,'%d %d %d',1);
    lenreadface = textscan(fid,'%d',1);
    lenreadface = lenreadface{1};
    
    %Initialise mesh
    mesh.nvertex = itemdat{1};
    mesh.nface = itemdat{2};
    mesh.ndim = itemdat{3};
    mesh.faces = cell(mesh.nface,1);
    mesh.vertices = zeros(mesh.nvertex,mesh.ndim);
    
    %Read vertices 
    if mesh.ndim == 3
        vtxdat = textscan(fid,'%f %f %f',mesh.nvertex);
        mesh.vertices(:,1) = vtxdat{1,1};
        mesh.vertices(:,2) = vtxdat{1,2};
        mesh.vertices(:,3) = vtxdat{1,3};
    elseif mesh.ndim == 2
        vtxdat = textscan(fid,'%f %f',mesh.nvertex);
        mesh.vertices(:,1) = vtxdat{1,1};
        mesh.vertices(:,2) = vtxdat{1,2};
    end 
    
    %Read faces 
    face_data = textscan(fid,'%d',lenreadface);
    face_data = face_data{1,1};

    %Extract faces
    row = 1;
    maxfaceNvtx = 0;
    for ii=1:mesh.nface 
        nvface = face_data(row,1);
        if nvface > maxfaceNvtx
            maxfaceNvtx = nvface;
        end
        mesh.faces{ii}.nvertex = nvface;
        row = row + 1;
        mesh.faces{ii}.vertices = face_data(row:row+nvface-1,1);
        row = row + nvface;
    end
    mesh.max_vtx_inface = maxfaceNvtx;
    
    %Close file
    fclose(fid);
end