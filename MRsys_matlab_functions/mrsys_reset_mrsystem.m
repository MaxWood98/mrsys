%MRsys reset multiresolution system control nets function 

%Version = 0.1
%Updated = 14-12-23

%Max Wood 2023
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [] = mrsys_reset_mrsystem(system_fname,surface_fname,display_toggle)
    cmd_str = ['mrsys rmrs ',system_fname,' ',surface_fname];
    if display_toggle == 1
        system(cmd_str);
    else
        [~,~] = system(cmd_str);
    end
end