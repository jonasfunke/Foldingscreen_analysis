function step2_anotate_profiles(fname, pname, sigma_integrate)
% step2: IFS analysis pipeline
% anotate the profiles by marking species

    data = load([pname filesep fname]); % load data
    % select aggregates and monomer bands
    data.profileData.sigma_integrate = sigma_integrate;
    data.profileData = select_species(data.profileData, data.gelData, data.gelInfo); 

    % integrate aggregates, smear, and monomer bands
    data.profileData = integrate_species(data.profileData);
    
    % save
    disp(['Saving to ' pname filesep fname]);
    save([pname filesep fname], '-struct', 'data');
    
    fname_fractions = [fname(1:end-4) '_fractions.out'];
    disp(['Saving fractions to ' pname fname_fractions])
    fid = fopen([pname fname_fractions],'w'); 
    fprintf(fid,'%s\n','#Monomer Smear Pocket MonomerFraction SmearFraction PocketFraction');
    fclose(fid);
    tmp = [data.profileData.monomerTotal data.profileData.smearTotal data.profileData.pocketTotal];
    fraction_tmp = tmp./(tmp(:,1)+tmp(:,2)+tmp(:,3));
    dlmwrite([pname filesep fname_fractions], [tmp fraction_tmp], 'delimiter', '\t', '-append');
    disp('Done')

end
