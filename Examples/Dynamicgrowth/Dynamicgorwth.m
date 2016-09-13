clc
clear
changeCobraSolver('glpk');
model=readCbModel('iAF1260');
model=changeObjective(model,'BIOMASS_Ec_iAF1260_core_59p81M');
% glucose (g/l)
Cglcexp=[1.91,1.883,1.8488,1.779,1.744,1.627,1.569,1.291,0.988,0.5698]/180.15*1000;
% time (h)
t=[0.053,1,1.5,2,2.5,3,3.5,4,4.5,5];
% biomass(gdcw/h)
Xexp=[0.01667,0.0294,0.03957,0.0596,0.094,0.1178,0.1775,0.2554,0.359,0.54113];

for i=1:length(t)-1
    dt(i)=t(i+1)-t(i);
    muexp(i)=(log(Xexp(i+1)/Xexp(i)))/dt(i);
    xglc=((Cglcexp(i+1)-Cglcexp(i))*muexp(i))/((exp(dt(i)*muexp(i))-1)*Xexp(i));
    model=changeRxnBounds(model,'EX_glc__D(e)',xglc,'l');
    s=optimizeCbModel(model,'max','one');
    Inflx(:,i)=s.x;
end
Bigg = model.rxns;
ModelName = 'UniModel';