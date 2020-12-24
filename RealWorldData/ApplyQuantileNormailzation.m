%% Process labels to find which are cancer tumors and which normal tissue

load('KIRCPANCANLabels.mat');
KIRC_CancerOrNotLabels = zeros(1, length(KIRCPANCANLabels));
for i=1:length(KIRCPANCANLabels)
    a = char(KIRCPANCANLabels(i));
    b = a(end-1:end); %Extract last two characters
    if strcmp(b,'11')==1
        KIRC_CancerOrNotLabels(i) = 1; %1's are normal
    elseif strcmp(b,'01')==1
        KIRC_CancerOrNotLabels(i) = 2; %2's are cancer
    end
end
normal_idx = find(KIRC_CancerOrNotLabels==1);
cancer_idx = find(KIRC_CancerOrNotLabels==2);

%% Apply Quantile Normalization:
load('KIRCPANCANFeatures');
KIRCPANCANFeatures_QN = (quantilenorm(KIRCPANCANFeatures'))';
KIRC_p_vals_QN = zeros(1,size(KIRCPANCANFeatures_QN,2));
for i=1:size(KIRCPANCANFeatures_QN,2)
    KIRC_p_vals_QN(i)=ranksum(KIRCPANCANFeatures_QN(cancer_idx,i),KIRCPANCANFeatures_QN(normal_idx,i));
end

%% Save values
save KIRC_p_vals_QN.mat KIRC_p_vals_QN; 

figure
hist(KIRC_p_vals_QN,30)
title('Distribution of Quantile Normalized P-values')
load('KIRC_p_vals.mat')
figure
hist(KIRC_p_vals,30)
title('Original Distribution of P-values')