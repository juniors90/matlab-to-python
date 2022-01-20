%% Function of initializing population
function p = Initialazation()

global No_pop VgMin VgMax No_generator NTrans NTransStep TransTap NQComp QCompMin QCompMax

V = VgMin + rand(No_pop, No_generator) * (VgMax - VgMin);
TT = ceil(NTransStep * rand(No_pop, NTrans));
T = TransTap(TT);
tic
for ii = 1:No_pop
    for jj = 1:NQComp
        QComp(ii, jj) = QCompMin(jj) + rand * (QCompMax(jj) - QCompMin(jj));
    end
end

QCompValue = ceil(1 + rand(No_pop, NQComp) * (length(QComp) - 1));
p = [V TT QComp];