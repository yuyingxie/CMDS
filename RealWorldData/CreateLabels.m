PANCANProcesedLabels = strings(11069,1);
counter = zeros(11069,1); %to check labels are unique
for i=1:11069
    % Extract TSS code:
    [token,remain] = strtok(PANCANLabels{1,i},'-');
    token2 = strtok(remain,'-');
    % Compare to TSS codes associated with each cancer type:
    for j=1:830
        if strcmp(token2,LabelLookUpTable.TSSCode{j})==1
           PANCANProcesedLabels(i) = LabelLookUpTable.CancerType{j};
           counter(i) = counter(i)+1;
        end
    end
end

%% Create Digit Labels from the text labels from PANCANDigitCategoryLookupTable:

PANCANDigitLabels = zeros(11069,1);
for i=1:11069
    for j=1:33
        if strcmp(PANCANProcesedLabels(i),PANCANDigitCategoryLookupTable.CancerType{j})==1
            PANCANDigitLabels(i) = PANCANDigitCategoryLookupTable.CancerTypeDigitLabel(j);
        end
    end
end