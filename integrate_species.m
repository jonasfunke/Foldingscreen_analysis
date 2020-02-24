function [profileData] = integrate_species(profileData, sigma_integrate)
% compute fractions
%   Detailed explanation goes here

    %% integrate aggregates
    pocket_sums = zeros(length(profileData.profiles), 1);
    pocket_boundaries = zeros(length(profileData.profiles),2);

    for i=1:length(profileData.profiles)
        pocket_boundaries(i,:) = [max(1,round(profileData.aggregateFit.b1-sigma_integrate*profileData.aggregateFit.c1)) ...
            round(profileData.aggregateFit.b1+sigma_integrate*profileData.aggregateFit.c1)  ];
        pocket_sums(i) = sum(profileData.fullProfiles{i}(pocket_boundaries(i,1):pocket_boundaries(i,2) ));
    end
    %% integrate monomers

    monomer_sums = zeros(length(profileData.profiles),1);
    monomer_boundaries = zeros(length(profileData.profiles),2);
    for i=1:length(profileData.profiles)
        monomer_boundaries(i,:) = [max(1,round(profileData.monomerFits{i}.b1-sigma_integrate*profileData.monomerFits{i}.c1)) ...
            round(profileData.monomerFits{i}.b1+sigma_integrate*profileData.monomerFits{i}.c1)  ];
        monomer_sums(i) = sum(profileData.fullProfiles{i}(monomer_boundaries(i,1):monomer_boundaries(i,2)));
    end  

     %% intergrate smear
        
    smear_sums = zeros(length(profileData.profiles),1);
    smear_boundaries = zeros(length(profileData.profiles),2);
    for i=1:length(profileData.profiles)
        smear_boundaries(i,:) = [pocket_boundaries(i,2) ...
            monomer_boundaries(i,1)];
        smear_sums(i) = sum(profileData.fullProfiles{i}(smear_boundaries(i,1):smear_boundaries(i,2)));
    end

    profileData.monomerBoundaries = monomer_boundaries;
    profileData.pocketBoundaries = pocket_boundaries;
    profileData.monomerTotal = monomer_sums;
    profileData.pocketTotal = pocket_sums;
    profileData.smearTotal = smear_sums;
    profileData.smearBoundaries = smear_boundaries;
    profileData.sigma_integrate = sigma_integrate;

    
    %subplot(5,1,5)
        %tmp = zeros(length(profileData.profiles),1);
        %for i=1:length(profileData.profiles)
        %    tmp(i) = fits{i}.b1;
        %end
        %plot(tmp, '.-')
        
    %% Find best lane
    
    
    
end

