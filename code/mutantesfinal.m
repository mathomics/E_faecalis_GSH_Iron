% This function performs the analysis described in the article "Reduced
% glutathione levels in Enterococcus faecalis trigger metabolic and
% transcriptional compensatory adjustments during iron exposure" (Víctor
% Aliaga-Tobar, Jorge Torres, Sebastián Mendoza, , Gabriel Gálvez, Jaime
% Ortega, Sebastián Gómez, Valentina Parra, Felipe Arenas, Alejandro Maass,
% Anne Siegel, Mauricio González, and Mauricio Latorre). In this analysis,
% in silico mutant strains are created and flux distributions are computed
% and compared with the wildtype strain. 
%
% INPUTS:
%    iEF671.mat:             Enterococcus faecalis model from the paper by
%                            Veith et al, 2015
%                            (https://doi.org/10.1128/AEM.03279-14)
%    tablaID.xlsx:           List of exchange reaction names for substrates
%                            initially in the media together with its
%                            correspoding lower bounds
%
% OUTPUTS:
%    flux_distributions.xlsx:   Excel File with the flux distributions of
%                               the different strains
%    specific_growth_rate_comparison.pdf: PDF file with a plot comparing
%                                         specific growth rates for the
%                                         different strains
%    gsh_mutant.mat:             Model with the strain ∆gsh
%    arc_mutant.mat:             Model with the strain ∆arc
%    gsh_arc_mutant.mat:         Model with the strain ∆gsh∆arc
%
% .. Authors: - Jorge Torres, Sebastián Mendoza 01/03/2025

%% Initialization of COBRA TOOLBOX, GUROBI and our toolbox

% initCobraToolbox
% changeCobraSolver('gurobi')
% cd('..')
% init_e_faecalis_gsh_iron_toolbox;
baseFolder = getBaseFolder_e_faecalis_ghs_iron_toolbox;
cd([baseFolder filesep 'results'])

%% Perform analysis

% load model
load([baseFolder filesep 'data' filesep 'iEF671_original.mat']);

% Display some information about the model
disp(model);

% change objective function to biomass prodduction
model.c(find(contains(lower(model.rxns),'bio'))) = 1;

% find positions of exchange reations 
posEX = find(startsWith(model.rxns, 'Ex_'));
media = [model.rxns(posEX) model.rxnNames(posEX) num2cell(model.lb(posEX))];
media2 = [model.rxns model.rxnNames num2cell(model.lb)];

% reset lower bounds
model.lb(posEX) = zeros(length(posEX),1);

% load lower bounds from excel file and set lower bounds in the model
tablaID = readtable([baseFolder filesep 'data' filesep 'tablaID.xlsx']);
for i = 1:height(tablaID)
    reactionID = tablaID.ID{i};  % reaction ID
    value = tablaID.Value(i);    % reaction lower bound
    rxnIndex = find(strcmp(model.rxns, reactionID));
    if ~isempty(rxnIndex)
        model.lb(rxnIndex) = value;  % change lower bound
    else
        fprintf('Reaction %s not found in the model.\n', reactionID);
    end
end

% check the new lower bounds
medianew = [model.rxns(posEX) model.rxnNames(posEX) num2cell(model.lb(posEX))];

% perfom Flux Balance Analysis for wild-type strain
fba_wt_aux = optimizeCbModel(model);

% compute scaled reduced costs for wild-type strain
pos_scaled_reduced_cost = find(abs(fba_wt_aux.w.*fba_wt_aux.x)>1e-6);
scale_reduced_costs = fba_wt_aux.w.*fba_wt_aux.x;
scale_reduced_cost_rxns = model.rxns(pos_scaled_reduced_cost);
scale_reduced_cost_values = scale_reduced_costs(pos_scaled_reduced_cost);pos_scaled_reduced_cost = find(abs(fba_wt_aux.w.*fba_wt_aux.x)>1e-6);

% perfom parsimonious Flux Balance Analysis for wild-type strain
[~, ~, ~, fba_wildtype] = pFBA_own(model, 'skipclass', 1,'tol',10^-7);

display(fba_wildtype)
for i = 1:length(model.grRules)
    if isempty(model.grRules{i}) 
        model.grRules{i}= ''; 
    end
end


% find the gene associated to gene EF3089 (glutathione synthetase)
pos_gsh = find(contains(model.grRules,'EF3089'));%GSH 3089

% create mutant ∆gsh
gsh_mutant = removeRxns(model,model.rxns(pos_gsh));

% save the model for the mutant strain ∆gsh
save([baseFolder filesep 'results' filesep 'models' filesep 'gsh_mutant'],'gsh_mutant');

% perform parsimonious Flux Balance Analysis for strain ∆gsh
[~, ~, ~, fba_mutante_gsh] = pFBA_own(gsh_mutant, 'skipclass', 1,'tol',10^-7);

% find the genes associated to arginine catabolism
pos_arg1 = find(contains(model.grRules,'EF0104'));%arginine deiminasa | catalyzes the degradation of arginine to citruline
pos_arg2 = find(contains(model.grRules,'EF0105'));%ornithine carbamoyltransferase|catalyzes the formation of ornithine and carbamylphosphate from citrulline in the arginine catabolic pathway
pos_arg3 = find(contains(model.grRules,'EF0106'));%carbamate kinase|catalyzes the reversible synthesis of carbamate and ATP from carbamoyl phosphate and ADP
pos_lip = find(contains(model.rxns ,'GLYCDH'));
pos_1 = find(contains(model.grRules,'EF0255'));
pos2 = find(contains(model.grRules,'EF0641'));

% perform parsimonious Flux Balance Analysis for strain ∆arc
arc_mutant = removeRxns(model,unique([model.rxns(pos_arg1); model.rxns(pos_arg3); model.rxns(pos_arg2(2))]));

% save the model for the mutant strain ∆arc
save([baseFolder filesep 'results' filesep 'models' filesep 'arc_mutant'],'arc_mutant');

% perform parsimonious Flux Balance Analysis for strain ∆arc
[~, ~, ~, fba_mutante_arginina] = pFBA_own(arc_mutant, 'skipclass', 1,'tol',10^-7);

% create mutant ∆gsh∆arc
pos_dko = unique([pos_gsh; pos_arg1; pos_arg2(2); pos_arg3]);
gsh_arc_mutant = removeRxns(model,model.rxns(pos_dko));
% save the model for the mutant strain for the double mutant strain ∆gsh∆arc
save([baseFolder filesep 'results' filesep 'models' filesep 'gsh_arc_mutant'],'gsh_arc_mutant');

% perform parsimonious Flux Balance Analysis for the double mutant strain ∆gsh∆arc
[~, ~, ~, fba_DKO_mutante] = pFBA_own(gsh_arc_mutant, 'skipclass', 1,'tol',10^-7);

%visualize fluxes for different strains
colors = [
    0.1059  0.6196  0.4667; % Color 1
    0.8510  0.3725  0.0078; % Color 2
    0.4588  0.4392  0.7020; % Color 3
    0.9059  0.1608  0.5412  % Color 4
];
labels = {'WT', '\Delta \it{gsh}', '\Delta \it{arg}', 'DKO'};
% Create single figure
figure;

biomass_fluxes = [
    fba_wildtype(find(strcmp(model.rxns, 'BIOMASS'))), 
    fba_mutante_gsh(find(strcmp(gsh_mutant.rxns, 'BIOMASS'))), 
    fba_mutante_arginina(find(strcmp(arc_mutant.rxns, 'BIOMASS'))), 
    fba_DKO_mutante(find(strcmp(gsh_arc_mutant.rxns, 'BIOMASS')))
];

hold on; % Hold on to overlay bars
for i = 1:length(biomass_fluxes)
    b = bar(i, biomass_fluxes(i));
    b.FaceColor = colors(i, :);
    b.FaceAlpha = 0.85; % Transparency
end
hold off;

% Set the x-axis labels
set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels);

% Label the axes
xlabel('Strains');
ylabel('Specific growth rate [1/h]');

% Set the title
title('Specific growth rate for different strains');

% save figure
saveas(gcf, [baseFolder filesep 'results' filesep 'plots' filesep 'specific_growth_rate_comparison.pdf']);    

%ATPrxn 
atpIndex = find(strcmp(model.mets, 'atp_c'));  
fluxes = fba_wildtype;
atpStoichiometry = model.S(atpIndex, :);
atpConsumingReactions = {};
atpProducingReactions = {};
for i = 1:length(model.rxns)
    if atpStoichiometry(i) < 0 && fluxes(i) > 0
        atpConsumingReactions{end+1} = model.rxns{i};
    elseif atpStoichiometry(i) > 0 && fluxes(i) < 0
        atpConsumingReactions{end+1} = model.rxns{i};
    elseif atpStoichiometry(i) > 0 && fluxes(i) > 0
        atpProducingReactions{end+1} = model.rxns{i};
    elseif atpStoichiometry(i) < 0 && fluxes(i) < 0
        atpProducingReactions{end+1} = model.rxns{i};
    end
end
disp('ATP Consuming Reactions:');
disp(atpConsumingReactions);
disp('ATP Producing Reactions:');
disp(atpProducingReactions);
%ADP rxn
adpIndex = find(strcmp(model.mets, 'adp_c'));  
fluxes = fba_wildtype;
adpStoichiometry = model.S(adpIndex, :);
adpConsumingReactions = {};
adpProducingReactions = {};
for i = 1:length(model.rxns)
    if adpStoichiometry(i) < 0 && fluxes(i) > 0
        adpConsumingReactions{end+1} = model.rxns{i};
    elseif adpStoichiometry(i) > 0 && fluxes(i) < 0
        adpConsumingReactions{end+1} = model.rxns{i};
    elseif adpStoichiometry(i) > 0 && fluxes(i) > 0
        adpProducingReactions{end+1} = model.rxns{i};
    elseif adpStoichiometry(i) < 0 && fluxes(i) < 0
        adpProducingReactions{end+1} = model.rxns{i};
    end
end
disp('ADP Consuming Reactions:');
disp(adpConsumingReactions);
disp('ADP Producing Reactions:');
disp(adpProducingReactions);
%NADP rxn
nadpIndex = find(strcmp(model.mets, 'nadp_c'));  
fluxes = fba_wildtype;
nadpStoichiometry = model.S(nadpIndex, :);
nadpConsumingReactions = {};
nadpProducingReactions = {};
for i = 1:length(model.rxns)
    if nadpStoichiometry(i) < 0 && fluxes(i) > 0
        nadpConsumingReactions{end+1} = model.rxns{i};
    elseif nadpStoichiometry(i) > 0 && fluxes(i) < 0
        nadpConsumingReactions{end+1} = model.rxns{i};
    elseif nadpStoichiometry(i) > 0 && fluxes(i) > 0
        nadpProducingReactions{end+1} = model.rxns{i};
    elseif nadpStoichiometry(i) < 0 && fluxes(i) < 0
        nadpProducingReactions{end+1} = model.rxns{i};
    end
end
disp('NADP Consuming Reactions:');
disp(nadpConsumingReactions);
disp('NADP Producing Reactions:');
disp(nadpProducingReactions);
%NAD rxn
nadIndex = find(strcmp(model.mets, 'nad_c'));  
fluxes = fba_wildtype;
nadStoichiometry = model.S(nadIndex, :);
nadConsumingReactions = {};
nadProducingReactions = {};
for i = 1:length(model.rxns)
    if nadStoichiometry(i) < 0 && fluxes(i) > 0
        nadConsumingReactions{end+1} = model.rxns{i};
    elseif nadStoichiometry(i) > 0 && fluxes(i) < 0
        nadConsumingReactions{end+1} = model.rxns{i};
    elseif nadStoichiometry(i) > 0 && fluxes(i) > 0
        nadProducingReactions{end+1} = model.rxns{i};
    elseif nadStoichiometry(i) < 0 && fluxes(i) < 0
        nadProducingReactions{end+1} = model.rxns{i};
    end
end
disp('NAD Consuming Reactions:');
disp(nadConsumingReactions);
disp('NAD Producing Reactions:');
disp(nadProducingReactions);
%NADH rxn
nadhIndex = find(strcmp(model.mets, 'nadh_c'));  
fluxes = fba_wildtype;
nadhStoichiometry = model.S(nadhIndex, :);
nadhConsumingReactions = {};
nadhProducingReactions = {};
for i = 1:length(model.rxns)
    if nadhStoichiometry(i) < 0 && fluxes(i) > 0
        nadhConsumingReactions{end+1} = model.rxns{i};
    elseif nadhStoichiometry(i) > 0 && fluxes(i) < 0
        nadhConsumingReactions{end+1} = model.rxns{i};
    elseif nadhStoichiometry(i) > 0 && fluxes(i) > 0
        nadhProducingReactions{end+1} = model.rxns{i};
    elseif nadhStoichiometry(i) < 0 && fluxes(i) < 0
        nadhProducingReactions{end+1} = model.rxns{i};
    end
end
disp('NADH Consuming Reactions:');
disp(nadhConsumingReactions);
disp('NADH Producing Reactions:');
disp(nadhProducingReactions);
%NADPH rxn
nadphIndex = find(strcmp(model.mets, 'nadph_c'));  
fluxes = fba_wildtype;
nadphStoichiometry = model.S(nadphIndex, :);
nadphConsumingReactions = {};
nadphProducingReactions = {};
for i = 1:length(model.rxns)
    if nadphStoichiometry(i) < 0 && fluxes(i) > 0
        nadphConsumingReactions{end+1} = model.rxns{i};
    elseif nadphStoichiometry(i) > 0 && fluxes(i) < 0
        nadphConsumingReactions{end+1} = model.rxns{i};
    elseif nadphStoichiometry(i) > 0 && fluxes(i) > 0
        nadphProducingReactions{end+1} = model.rxns{i};
    elseif nadphStoichiometry(i) < 0 && fluxes(i) < 0
        nadphProducingReactions{end+1} = model.rxns{i};
    end
end
disp('NADPH Consuming Reactions:');
disp(nadphConsumingReactions);
disp('NADPH Producing Reactions:');
disp(nadphProducingReactions);



% Wild-type
rxns_wt = model.rxns;
fluxes_wt = fba_wildtype;
wt_table = table(rxns_wt, fluxes_wt, 'VariableNames', {'Reactions', 'Fluxes'});

% Mutant GSH
rxns_gsh = gsh_mutant.rxns;
fluxes_gsh = fba_mutante_gsh;
gsh_table = table(rxns_gsh, fluxes_gsh, 'VariableNames', {'Reactions', 'Fluxes'});

% Mutant Arginine
rxns_arg = arc_mutant.rxns;
fluxes_arg = fba_mutante_arginina;
arg_table = table(rxns_arg, fluxes_arg, 'VariableNames', {'Reactions', 'Fluxes'});

% Double Knockout Mutant
rxns_dko = gsh_arc_mutant.rxns;
fluxes_dko = fba_DKO_mutante;
dko_table = table(rxns_dko, fluxes_dko, 'VariableNames', {'Reactions', 'Fluxes'});

%Write to Excel file
filename = 'flux_distributions.xlsx';
writetable(wt_table, [baseFolder filesep 'results' filesep 'tables' filesep filename], 'Sheet', 'Wild_type');
writetable(gsh_table, [baseFolder filesep 'results' filesep 'tables' filesep filename], 'Sheet', 'Mutant_GSH');
writetable(arg_table, [baseFolder filesep 'results' filesep 'tables' filesep filename], 'Sheet', 'Mutant_Arginine');
writetable(dko_table, [baseFolder filesep 'results' filesep 'tables' filesep filename], 'Sheet', 'Double_Knockout');

%%%%%

fba_wt_gsh = optimizeCbModel(gsh_mutant);
fba_wt_arg = optimizeCbModel(arc_mutant);
fba_wt_dko = optimizeCbModel(gsh_arc_mutant);

%compute reduced costs and scaled reduced costs for the strains
pos_scaled_reduced_cost1 = find(abs(fba_wt_gsh.w.*fba_wt_gsh.x)>1e-6);
scale_reduced_costs1 = fba_wt_gsh.w.*fba_wt_gsh.x;
scale_reduced_cost_rxns1 = gsh_mutant.rxns(pos_scaled_reduced_cost1);
scale_reduced_cost_values1 = scale_reduced_costs1(pos_scaled_reduced_cost1);

pos_scaled_reduced_cost2 = find(abs(fba_wt_arg.w.*fba_wt_arg.x)>1e-6);
scale_reduced_costs2 = fba_wt_arg.w.*fba_wt_arg.x;
scale_reduced_cost_rxns2 = arc_mutant.rxns(pos_scaled_reduced_cost2);
scale_reduced_cost_values2 = scale_reduced_costs2(pos_scaled_reduced_cost2);

pos_scaled_reduced_cost3 = find(abs(fba_wt_dko.w.*fba_wt_dko.x)>1e-6);
scale_reduced_costs3 = fba_wt_dko.w.*fba_wt_dko.x;
scale_reduced_cost_rxns3 = gsh_arc_mutant.rxns(pos_scaled_reduced_cost3);
scale_reduced_cost_values3 = scale_reduced_costs3(pos_scaled_reduced_cost3);
disp('')

%polyamines pathway
% List of metabolites
mets = {'orn_L_c', 'ptrc_c', 'spmd_c', 'arg_L', 'agma_c', 'cbp_c', 'amet_c', 'p15da_c'};

% Initialize a cell array to hold reactions for each metabolite
reactions = cell(length(mets), 2);

% Loop through each metabolite
for i = 1:length(mets)
    % Find indexes of reactions that involve the current metabolite
    rxnIndexes = find(contains(model.mets, mets{i}));
    
    % Get the IDs of these reactions
    reactions{i, 1} = mets{i};
    reactions{i, 2} = model.rxns(rxnIndexes);
end

% Display reactions for each metabolite
reactions

% Extract reaction IDs and names
rxnIDs = model.rxns;
rxnNames = model.rxnNames;

% Extract reactions and their corresponding gene reaction rules
reactions = model.rxns;
geneReactionRules = model.grRules;

% Create a table with the reactions and their gene reaction rules
reactionTab = table(reactions, geneReactionRules, 'VariableNames', {'Reactions', 'GeneReactionRules'});

% Display the table
disp(reactionTab);

% Crear una tabla con las reacciones y sus límites
reactionIDs = model.rxns;
lowerBounds = model.lb;
upperBounds = model.ub;

% Create table with reactions and bounds
reactionlimit = table(reactionIDs, lowerBounds, upperBounds, ...
                      'VariableNames', {'ReactionID', 'LowerBound', 'UpperBound'});

% save table as excel file
writetable(reactionlimit, [baseFolder filesep 'results' filesep 'tables' filesep 'reactions_bounds.xlsx']);

disp('')
