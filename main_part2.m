%% anotate the profiles by marking species
[fname, pname] = uigetfile([pwd '/*.mat'], 'Select mat file');
sigma_integrate_band = 1.0;
data = load([pname fname]); % load data for completion chek

is_done = false;
if isfield(data.profileData, 'aggregateSelectedArea')
    is_done = ~strcmp(questdlg('Pocket and monomers have already been selected. Redo?','Redo?' ,'No','Yes', 'Yes'),'Yes');
end

if ~is_done
    step2_anotate_profiles(fname, pname, sigma_integrate_band);
else
    disp('Nothing to do');
end
