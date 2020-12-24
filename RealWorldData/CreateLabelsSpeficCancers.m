%% Process labels to find which are cancer tumors and which normal tissue

THCA_CancerOrNotLabels = zeros(1, length(THCAPANCANLabels));
for i=1:length(THCAPANCANLabels)
    a = char(THCAPANCANLabels(i));
    b = a(end-1:end); %Extract last two characters
    if strcmp(b,'11')==1
        THCA_CancerOrNotLabels(i) = 1; %1's are normal
    elseif strcmp(b,'01')==1
        THCA_CancerOrNotLabels(i) = 2; %2's are cancer
    end
end
normal_idx = find(THCA_CancerOrNotLabels==1);
cancer_idx = find(THCA_CancerOrNotLabels==2);

%% Remove weird character from genes

for i=1:length(THCAPANCANGenes)
    a = char(THCAPANCANGenes(i));
    a(regexp(a,'[?,|]'))=[];
    THCAPANCANGenes(i) = a;
end

%% Check if gene names match

load('PANCANGenes.mat');
ThereOrNot = zeros(length(THCAPANCANGenes),1);
for i=1:length(THCAPANCANGenes)
    ThereOrNot(i) = sum(strcmp(PANCANGenes,THCAPANCANGenes(i)));
end

ThereOrNot2 = zeros(length(PANCANGenes),1);
for i=1:length(PANCANGenes)
    ThereOrNot2(i) = sum(strcmp(THCAPANCANGenes,PANCANGenes(i)));
end

length(find(ThereOrNot==0))
length(find(ThereOrNot2==0))

%% Reorder the cancer specific gene names to match global data set

reorder_THCA = zeros(length(THCAPANCANGenes),1);
for i=1:length(THCAPANCANGenes)
    reorder_THCA(i) = find(PANCANGenes(i)==THCAPANCANGenes);
end

OrderedTHCAPANCANGenes = THCAPANCANGenes(reorder_THCA);
sum(strcmp(PANCANGenes,OrderedTHCAPANCANGenes)) %Check they match

OrderedTHCAPANCANFeatures = THCAPANCANFeatures(:,reorder_THCA);

%% Make sure order of genes matches full data set before calculating this

THCA_p_vals = zeros(1,size(OrderedTHCAPANCANFeatures,2));
for i=1:size(OrderedTHCAPANCANFeatures,2)
    THCA_p_vals(i)=ranksum(OrderedTHCAPANCANFeatures(cancer_idx,i),OrderedTHCAPANCANFeatures(normal_idx,i));
end

%% If everything checks out, just save the reordered everything

THCAPANPANFeatures = OrderedTHCAPANCANFeatures;
save THCAPANCANFeatures.mat THCAPANCANFeatures;
THCAPANPANGenes = OrderedTHCAPANCANGenes;
save THCAPANCANGenes.mat THCAPANCANGenes;

save THCA_p_vals.mat THCA_p_vals; %These were re-ordered on creation
save THCAPANCANLabels.mat THCAPANCANLabels; %No need to re-order
