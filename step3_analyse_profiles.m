%% analyse the profiles
[fname, pname] = uigetfile([pwd '/*.mat'], 'Select mat file');
data = load([pname fname]);

step3(data, fname, pname)

function step3(data, fname, pname)
    [data_tmp, cur_fig] = get_best_folding(data.profileData, data.gelInfo, data.gelData, true);

    tmp = strsplit(pname, filesep);
    print(cur_fig, '-dpdf', [pname tmp{end-1} '_analysis.pdf']); %save figure
    data.foldingAnalysis = data_tmp;

    disp(['Saving to ' pname fname])
    save([pname fname], '-struct', 'data')
    disp('Done')
end

