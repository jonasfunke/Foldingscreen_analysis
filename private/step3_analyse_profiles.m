function step3_analyse_profiles(fname, pname)
% step3: IFS analysis pipeline
% analyse the profiles
    data = load([pname filesep fname]);
    [data_tmp, cur_fig] = get_best_folding(data.profileData, data.gelInfo, data.gelData, true);

    tmp = strsplit(pname, filesep);
    print(cur_fig, '-dpdf', [pname filesep tmp{end-1} '_analysis.pdf']); %save figure
    data.foldingAnalysis = data_tmp;

    disp(['Saving to ' pname filesep fname])
    save([pname filesep fname], '-struct', 'data')
    disp('Done')
end

