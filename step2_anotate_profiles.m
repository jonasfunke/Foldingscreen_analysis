%% anotate the profiles by marking species
[fname, pname] = uigetfile([pwd '/*.mat'], 'Select mat file');

data = load([pname fname]); % load data


is_done = false;
if isfield(data.profileData, 'aggregateSelectedArea')
    is_done = ~strcmp(questdlg('Pocket and monomers have already been selected. Redo?','Redo?' ,'No','Yes', 'Yes'),'Yes');
end

if ~is_done
    step2(data, fname, pname)
else
    disp('Nothing to do');
end

function step2(data, fname, pname) 

    % select aggregates and monomer bands
    data.profileData = select_species(data.profileData, data.gelData, 1.0); 

    % integrate aggregates, smear, and monomer bands
    data.profileData = integrate_species(data.profileData, 1.0);
    
    % save
    disp(['Saving to ' pname fname]);
    save([pname fname], '-struct', 'data');
    
    fname_fractions = [fname(1:end-4) '_fractions.out'];
    disp(['Saving fractions to ' pname fname_fractions])
    fid = fopen([pname fname_fractions],'w'); 
    fprintf(fid,'%s\n','#Monomer Smear Pocket MonomerFraction SmearFraction PocketFraction');
    fclose(fid);
    tmp = [data.profileData.monomerTotal data.profileData.smearTotal data.profileData.pocketTotal];
    fraction_tmp = tmp./(tmp(:,1)+tmp(:,2)+tmp(:,3));
    dlmwrite([pname fname_fractions], [tmp fraction_tmp], 'delimiter', '\t', '-append');
    disp('Done')

end
