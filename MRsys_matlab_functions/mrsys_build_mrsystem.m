%MRsys build multiresolution system function 

%Version = 0.1
%Updated = 14-12-23

%Max Wood 2023
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [] = mrsys_build_mrsystem(surf_name,ncrs_max,display_toggle)
    cmd_str = ['mrsys bmrs ',surf_name,' ',num2str(ncrs_max)];
    if display_toggle == 1
        system(cmd_str);
    else
        [~,~] = system(cmd_str);
    end
end