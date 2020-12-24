%% Remove weird character from genes
for i=1:length(KIRCPANCANGenes)
    a = char(KIRCPANCANGenes(i));
    a(regexp(a,'[?,|]'))=[];
    KIRCPANCANGenes(i) = a;
end

%% Check if gene names match
ThereOrNot = zeros(length(KIRCPANCANGenes),1);
for i=1:length(KIRCPANCANGenes)
    ThereOrNot(i) = sum(strcmp(PANCANGenes,KIRCPANCANGenes(i)));
end

ThereOrNot2 = zeros(length(PANCANGenes),1);
for i=1:length(PANCANGenes)
    ThereOrNot2(i) = sum(strcmp(KIRCPANCANGenes,PANCANGenes(i)));
end

length(find(ThereOrNot==0))
length(find(ThereOrNot2==0))

%% Reorder the cancer specific gene names to match global data set

reorder_KIRC = zeros(length(KIRCPANCANGenes),1);
for i=1:length(KIRCPANCANGenes)
    reorder_KIRC(i) = find(PANCANGenes(i)==KIRCPANCANGenes);
end

OrderedKIRCPANCANGenes = KIRCPANCANGenes(reorder_KIRC);
sum(strcmp(PANCANGenes,OrderedKIRCPANCANGenes)) %Check they match

OrderedKIRCPANCANFeatures = KIRCPANCANFeatures(:,reorder_KIRC);

%% If everything checks out, just save the reordered everything

KIRCPANPANFeatures = OrderedKIRCPANCANFeatures;
save KIRCPANCANFeatures.mat KIRCPANCANFeatures;
KIRCPANPANGenes = OrderedKIRCPANCANGenes;
save KIRCPANCANGenes.mat KIRCPANCANGenes;

save KIRC_p_vals.mat KIRC_p_vals; %These were re-ordered on creation
save KIRCPANCANLabels.mat KIRCPANCANLabels; %No need to re-order


