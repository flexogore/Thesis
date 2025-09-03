%model = createModel;
mod = readSBML('iCre1355_auto.xml', 1000);
% Define the reaction list
rxn_list = {'CO2t', 'CO2th', 'RBPCh', 'BFBPh', 'RPIh', 'PRUK', 'GLPThi', 'GAPDH_nadp_hi',...
    'ATPSh', 'CBFC', 'CEF', 'FNORh', 'PSIblue', 'PSIIblue', 'CS', 'ACONT', 'SUCDHh',...
    'FUMm', 'ENO', 'PGM', 'PYK', 'PYRtm', 'PDHam1mi', 'PDHam2mi', 'PDHe2r','PDHe3mr',...
    'CSm', 'ACONTm', 'ICDHm', 'AKGDHmi', 'AKGDHe2r', 'SUCDHm', 'FUMm',...
    'NADHOR_2m', 'CYOR_q8_m', 'CYOO6m', 'ATPSm', 'PPCKm', 'PYKm'};

%candidate rxns , 'ACS', 'CS', 'ACONT, 'ICL', 'SUCDHh', 'FUMm'

% Extract the subnetwork
selected_rxn_indices = ismember(mod.rxns, rxn_list);
selected_rxns = mod.rxns(selected_rxn_indices);
model = extractSubNetwork(mod, selected_rxns);

% Add reactions
model = addReaction(model,'EX_co2_e','metaboliteList',{'co2[e]'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'CO2tm','metaboliteList',{'co2[m]', 'co2[c]'},'stoichCoeffList',[-1, 1], 'reversible',true);
%model = addReaction(model,'ACCOAtm','metaboliteList',{'accoa[m]', 'coa[c]', 'accoa[c]', 'coa[m]'},'stoichCoeffList',[-1, -1, 1, 1], 'reversible',true);
model = addReaction(model,'ACCOAtm','metaboliteList',{'accoa[m]', 'accoa[c]'},'stoichCoeffList',[-1, 1], 'reversible',true);
model = addReaction(model,'COAtm','metaboliteList',{'coa[c]', 'coa[m]'},'stoichCoeffList',[-1, 1], 'reversible',true);


model = addReaction(model,'ICL','metaboliteList',{'icit[c]', 'glx[c]', 'succ[c]'},'stoichCoeffList',[-1, 1, 1], 'reversible',false);

model = addReaction(model,'T_T3P','metaboliteList',{'3pg[h]','3pg[c]'},'stoichCoeffList',[-1, 1], 'reversible',false); %rev
model = addReaction(model,'T_FBAp','metaboliteList',{'3pg[c]','fdp_B[c]'},'stoichCoeffList',[-2, 1], 'reversible',false); %rev

model = addReaction(model,'PGKh','metaboliteList',{'3pg[h]', '13dpg[h]'},'stoichCoeffList',[-1, 1], 'reversible',true); %rev
model = addReaction(model,'TPIh','metaboliteList',{'g3p[h]','dhap[h]'},'stoichCoeffList',[-1, 1], 'reversible',true); %rev
model = addReaction(model,'FBAh','metaboliteList',{'g3p[h]','dhap[h]','fdp_B[h]'},'stoichCoeffList',[-1, -1, 1], 'reversible',true); %rev


model = addReaction(model,'TKT1h','metaboliteList',{'g3p[h]', 's7p[h]', 'xu5p_D[h]', 'r5p[h]'},'stoichCoeffList',[-1, -1, 1, 1], 'reversible',true); %rev
model = addReaction(model,'TKT2h','metaboliteList',{'g3p[h]', 'f6p_B[h]', 'xu5p_D[h]', 'e4p[h]'},'stoichCoeffList',[-1, -1, 1, 1], 'reversible',true); %rev
model = addReaction(model,'RPEh','metaboliteList',{'xu5p_D[h]', 'ru5p_D[h]'},'stoichCoeffList',[-1, 1], 'reversible',true); %rev
model = addReaction(model,'TAh','metaboliteList',{'e4p[h]','f6p_B[h]','g3p[h]','s7p[h]'},'stoichCoeffList',[-1,-1, 1, 1], 'reversible',true); %rev

model = addReaction(model,'PGIAh','metaboliteList',{'f6p_B[h]','g6p_A[h]'},'stoichCoeffList',[-1, 1], 'reversible',true); %rev
model = addReaction(model,'PGMTh','metaboliteList',{'g6p_A[h]','g1p[h]'},'stoichCoeffList',[-1, 1], 'reversible',true); %rev
model = addReaction(model,'Biomass','metaboliteList',{'adpglc[h]'},'stoichCoeffList' ,[-1], 'reversible',false); %rev
model = addReaction(model,'MDH_nadp_h','metaboliteList',{'mal_L[h]', 'nadp[h]', 'h[h]', 'nadph[h]', 'oaa[h]'},'stoichCoeffList' ,[-1, -1, 1, 1, 1], 'reversible',false); %rev

model = addReaction(model,'PGIp_c','metaboliteList',{'fdp_B[c]'},'stoichCoeffList',[-1], 'reversible',false); %rev

model = addReaction(model,'EX_h2o','metaboliteList',{'h2o[h]'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'EX_h2o_u','metaboliteList',{'h2o[u]'},'stoichCoeffList',[1], 'reversible',false);


model = addReaction(model,'T_h','metaboliteList',{'h[h]', 'h[c]'},'stoichCoeffList',[-1, 1], 'reversible',false);
model = addReaction(model,'EX_pi','metaboliteList',{'pi[h]'},'stoichCoeffList',[-1], 'reversible',false);
model = addReaction(model,'EX_nadph','metaboliteList',{'nadph[h]'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'EX_nadp','metaboliteList',{'nadp[h]'},'stoichCoeffList',[-1], 'reversible',false);
model = addReaction(model,'EX_ppi','metaboliteList',{'ppi[h]'},'stoichCoeffList',[-1], 'reversible',false);

model = addReaction(model,'EX_atp','metaboliteList',{'atp[h]'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'EX_adp','metaboliteList',{'adp[h]'},'stoichCoeffList',[-1], 'reversible',false);

model = addReaction(model,'EX_photon437','metaboliteList',{'photon437[u]'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'EX_photon438','metaboliteList',{'photon438[u]'},'stoichCoeffList',[1], 'reversible',false);

model = addReaction(model,'EX_o2D','metaboliteList',{'o2D[u]'},'stoichCoeffList',[-1], 'reversible',false);
%model = addReaction(model,'T_atp_h','metaboliteList',{'atp[h]', 'atp[c]'},'stoichCoeffList',[-1, 1], 'reversible',false);
model = addReaction(model,'T_atp_m','metaboliteList',{'atp[m]', 'atp[c]'},'stoichCoeffList',[-1, 1], 'reversible',false);

model = addReaction(model,'MALSc','metaboliteList',{'accoa[c]', 'glx[c]', 'h2o[c]', 'coa[c]', 'h[c]', 'mal_L[c]'},'stoichCoeffList',[-1, -1, -1, 1, 1, 1], 'reversible',false);
model = addReaction(model,'MDHc','metaboliteList',{'mal_L[c]', 'nad[c]', 'h[c]', 'nadh[c]', 'oaa[c]'},'stoichCoeffList',[-1, -1, 1, 1, 1], 'reversible',false);
model = addReaction(model,'MDHm','metaboliteList',{'mal_L[m]', 'nad[m]', 'h[m]', 'nadh[m]', 'oaa[m]'},'stoichCoeffList',[-1, -1, 1, 1, 1], 'reversible',false);
model = addReaction(model,'SUCLm','metaboliteList',{'adp[m]', 'pi[m]', 'succoa[m]', 'atp[m]', 'coa[m]', 'h[m]', 'succ[m]'},'stoichCoeffList',[-1, -1, -1, 1, 1, 1, 1], 'reversible',false);


model = addReaction(model,'EX_nad','metaboliteList',{'nad[c]'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'EX_nadh','metaboliteList',{'nadh[c]'},'stoichCoeffList',[-1], 'reversible',false);

model = addReaction(model,'T_succ_h','metaboliteList',{'succ[c]', 'succ[h]'},'stoichCoeffList',[-1, 1], 'reversible',false);
%model = addReaction(model,'T_succ_m','metaboliteList',{'succ[c]', 'succ[m]'},'stoichCoeffList',[-1, 1], 'reversible',false);
model = addReaction(model,'EX_fad','metaboliteList',{'fad[h]'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'EX_fadh2','metaboliteList',{'fadh2[h]'},'stoichCoeffList',[-1], 'reversible',false);
model = addReaction(model,'T_fum','metaboliteList',{'fum[h]', 'fum[m]'},'stoichCoeffList',[-1, 1], 'reversible',false);
%model = addReaction(model,'T_fum','metaboliteList',{'fum[h]'},'stoichCoeffList',[-1], 'reversible',false);
model = addReaction(model,'EX_h2o_m','metaboliteList',{'h2o[m]'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'T_mal','metaboliteList',{'mal_L[m]', 'mal_L[h]'},'stoichCoeffList',[-1, 1], 'reversible',false);

model = addReaction(model,'EX_oaa','metaboliteList',{'oaa[h]'},'stoichCoeffList',[-1], 'reversible',false);
%model = addReaction(model,'EX_succ_m','metaboliteList',{'succ[m]'},'stoichCoeffList',[-1], 'reversible',false);
%model = addReaction(model,'EX_mal_m','metaboliteList',{'mal_L[m]'},'stoichCoeffList',[-1], 'reversible',false);


model = addExchangeRxn(model, 'adp[c]', -1000, 1000);
model = addExchangeRxn(model, 'atp[c]', -1000, 1000);
model = addExchangeRxn(model, 'h[m]', -1000, 1000);
%model = addExchangeRxn(model, 'coa[m]', -1000, 1000);

model = addExchangeRxn(model, 'h2o[m]', -1000, 1000);
%model = addExchangeRxn(model, 'accoa[m]', -1000, 1000);
model = addExchangeRxn(model, 'pi[m]', -1000, 1000);
model = addExchangeRxn(model, 'adp[m]', -1000, 1000);
model = addExchangeRxn(model, 'fad[m]', -1000, 1000);
model = addExchangeRxn(model, 'fadh2[m]', -1000, 1000);
model = addExchangeRxn(model, 'o2[m]', -1000, 1000);
model = addExchangeRxn(model, 'nad[m]', -1000, 1000);
model = addExchangeRxn(model, 'nadh[m]', -1000, 1000);

%model = changeRxnBounds(model, 'MDHm', 0, 'l');
%model = changeRxnBounds(model, 'SUCLm', 0, 'l');

%model = addReaction(model, 'EX_succ', ...
%    'metaboliteList', {'succ[c]'}, ...
%    'stoichCoeffList', [1], ...
%    'reversible', false);

% Drain h[c] if it's unused
model = addReaction(model, 'EX_h_c', ...
    'metaboliteList', {'h[c]'}, ...
    'stoichCoeffList', [-1], ...
    'reversible', false);
% Set objective
% model = changeObjective(model, 'Biomass');

%model = changeObjective(model, 'CS');
%model.c(strcmp(model.rxns, 'Biomass')) = 1;

%solus = optimizeCbModel(model, 'max', 'zero');

[mini,maxi] = fluxVariability(model, 0);

view = table(model.rxns, mini, maxi);

writeSBML(model, 'working_model.xml');

options.nStepsPerPoint = 30;
options.nPointsReturned = 20;
options.nWarmupPoints = 10;

[samples,s] = sampleCbModel(model, '', 'ACHR', options);