%MRsys deform multiresolution system function 

%Version = 0.1
%Updated = 14-12-23

%Max Wood 2023
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [] = mrsys_deform_mrsystem(system_fname,level2def,deformation_fname,display_toggle)
    cmd_str = ['mrsys pmrs ',system_fname,' ',num2str(level2def),' ',deformation_fname];
    if display_toggle == 1
        system(cmd_str);
    else
        [~,~] = system(cmd_str);
    end
end