function [profileData] = integrate_species(profileData)
% @ step2
% compute fractions of selected species per lane

    %% integrate aggregates
    n = length(profileData.profiles);
    pocket_sums = zeros(n, 1);
    pocket_boundaries = zeros(n,2);
    mu = profileData.aggregateFit(2);
    sig = profileData.aggregateFit(3) * profileData.sigma_integrate;
    for i=1:n
        if profileData.monomerFits(i,1) == 0.0 %ignored lane
            pocket_sums(i) = 1.0;
        else 
            pocket_boundaries(i,:) = [max(1,round(mu-sig)) round(mu+sig)];
            pocket_sums(i) = sum(profileData.fullProfiles{i}(pocket_boundaries(i,1):pocket_boundaries(i,2) ));
        end
    end
    
    %% integrate monomers
    monomer_sums = zeros(n,1);
    monomer_boundaries = zeros(n,2);
    for i=1:n
        if profileData.monomerFits(i,1) == 0.0 
            monomer_sums(i) = 0.0;
        else    
            mu = profileData.monomerFits(i,2);
            sig = profileData.monomerFits(i,3) * profileData.sigma_integrate;
            monomer_boundaries(i,:) = [max(1,round(mu-sig)) round(mu+sig)];
            monomer_sums(i) = sum(profileData.fullProfiles{i}(monomer_boundaries(i,1):monomer_boundaries(i,2)));
        end
    end

    %% intergrate smear
    smear_sums = zeros(n,1);
    smear_boundaries = zeros(n,2);
    for i=1:n
        if profileData.monomerFits(i,1) == 0.0 
            smear_sums(i) = 0.0;
        else 
            smear_boundaries(i,:) = [pocket_boundaries(i,2) monomer_boundaries(i,1)];
            smear_sums(i) = sum(profileData.fullProfiles{i}(smear_boundaries(i,1):smear_boundaries(i,2)));
        end
    end

    profileData.monomerBoundaries = monomer_boundaries;
    profileData.pocketBoundaries = pocket_boundaries;
    profileData.monomerTotal = monomer_sums;
    profileData.pocketTotal = pocket_sums;
    profileData.smearTotal = smear_sums;
    profileData.smearBoundaries = smear_boundaries;
    
end

