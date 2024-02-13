%MRsys interface script 
clearvars

%% Perform subdivision 

%Options
surfname = 'test_geometry\cube_open';
levels = 3;

%Call 
sdcommand = ['mrsys sdm ',surfname,' ',num2str(levels)];
system(sdcommand);

%Import surface
[mesh] = import_surface_mesh('subd_data\mrsys_surface');
[mesh] = build_full_faces(mesh);

%Plot
cla reset
patch('vertices',mesh.vertices,'faces',mesh.faces_full,'facecolor',[0.7,0.7,0.7])
axis equal
axis tight
xlabel('x')
ylabel('y')
zlabel('z')