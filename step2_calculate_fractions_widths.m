%% analyze the profiles

close all, clear all, clc
[fname pname] = uigetfile('/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/*.mat', 'Select mat file');

data = load([pname fname]); % load data


if ~isfield(data.profileData, 'aggregateSelectedArea')
    % select aggregates and monomer bands
    data.profileData = select_species(data.profileData, data.gelData, 1.5); 

    % integrate aggregates, smear, and monomer bands
    data.profileData = integrate_species(data.profileData, 1.5);
    %get_best_folding(data.profileData, data.gelInfo, data.gelData, true);
    
    save([pname fname], '-struct','data')
    
    % write intergrated intensities to file
    tmp = [data.profileData.monomerTotal data.profileData.smearTotal data.profileData.pocketTotal];
    fraction_tmp = tmp./(tmp(:,1)+tmp(:,2)+tmp(:,3));
    
    fid = fopen([pname fname(1:end-4) '_data.out'],'w'); 
    fprintf(fid,'%s\n','#Monomer Smear Pocket MonomerFraction SmearFraction PocketFraction');
    fclose(fid);
    dlmwrite([pname fname(1:end-4) '_data.out'], [tmp fraction_tmp], 'delimiter', '\t', '-append');
    
    disp(['Data saved to ' pname fname])

else
    display('Pocket and monomers have already been selected.')
end

%% dertime the best folding conditions
[data_tmp, cur_fig] = get_best_folding(data.profileData, data.gelInfo, data.gelData, true);
tmp = strsplit(pname, filesep);
print(cur_fig, '-dpdf', [pname tmp{end-1} '_analysis.pdf']); %save figure
data.foldingAnalysis = data_tmp;
save([pname fname], '-struct','data')
disp(['Data saved to ' pname fname])
disp('Done')


