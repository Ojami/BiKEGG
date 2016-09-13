clc
clear
changeCobraSolver('glpk');
model=readCbModel('iAF1260');
model=changeObjective(model,'BIOMASS_Ec_iAF1260_core_59p81M');
model = changeRxnBounds(model,'EX_glc__D(e)', -8.7,'l');
model = changeRxnBounds(model,'EX_o2(e)', -11.9,'l');
sAr=optimizeCbModel(model,'max','one');
% printFluxVector(model,sAr.x,true,true)
inflxAr=sAr.x;

model = changeRxnBounds(model,'EX_glc__D(e)', -14.9,'l');
model = changeRxnBounds(model,'EX_o2(e)', 0,'b');
sAnr=optimizeCbModel(model,'max','one');
% printFluxVector(model,sAnr.x,true,true)
inflxAnr=sAnr.x;
Inflx = [inflxAr,inflxAnr];
Bigg = model.rxns;
ModelName = 'UniModel';