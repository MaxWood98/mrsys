%MRsys multiresolution system surface gradient evaluation function 

%Version = 0.1
%Updated = 14-12-23

%Max Wood 2023
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [surfgrad] = mrsys_evaluate_mrsystem_surfgrad(grad_fname,system_fname,level2eval,display_toggle)
    
    %Evaluate
    cmd_str = ['mrsys gradients ',system_fname,' ',num2str(level2eval)];
    if display_toggle == 1
        system(cmd_str);
    else
        [~,~] = system(cmd_str);
    end
    
    %Import
    surfgrad = load(grad_fname);
end