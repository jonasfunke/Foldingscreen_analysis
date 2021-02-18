%% analyse the profiles
[fname, pname] = uigetfile([pwd '/*.mat'], 'Select mat file');
name = fname(1:end-9);
step3_analyse_profiles(name, pname);