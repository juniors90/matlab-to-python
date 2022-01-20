%% Main mfile should be run
clc 
clear all
definParameters();
global No_Cap_Type NBus No_pop Iter Cap_Price Ke Loaddata Strdata T_OffPeak T_Medium T_Peak NLoadLevel Kp tic
PLoss = zeros(No_pop, 1);
f = zeros(No_pop, 1);
LoadDataBase = Loaddata(:, 3);
LoadOffPeak = 0.3 * LoadDataBase;
LoadMedium = 0.6 * LoadDataBase;
LoadPeak = LoadDataBase;
Loaddata(:, 3) = LoadOffPeak;

%%%% Evaluating Initial conditions

[PLossOutOffPeak0, VbusOutOffPeak0, IsecOut0] = DLF(Strdata, Loaddata);
Loaddata(:, 3) = LoadMedium;
[PLossOutMedium0, VbusOutMedium0, IsecOut0] = DLF(Strdata, Loaddata);
Loaddata(:, 3) = LoadPeak;
[PLossOutPeak0, VbusOutPeak0, Isec0] = DLF(Strdata, Loaddata);
EnergyLossIni = T_OffPeak * PLossOutOffPeak0 +
T_Medium * PLossOutMedium0 + ...
    T_Peak * PLossOutPeak0; Loaddata(:, 3) = LoadDataBase;
%%%
p = ceil(rand(No_pop, NBus - 1) * No_Cap_Type); % % % Initial popoulation
pop = Cap_Mvar_determine(p); % % % Allocation MVAr to the generated population

for i = 1:size(p, 1)
    pop(i, :) = Cap_Mvar_determine(p(i, :));
    Load(:, 1) = LoadOffPeak - (pop(i, :))'; Load(:, 2) = LoadMedium - (pop(i, :))';
    Load(:, 3) = LoadPeak -
    (pop(i, :))'; Total_Cap_Price
    = sum(Cap_Price((p(i, :))));

    for il = 1:NLoadLevel
        Loaddata(:, 3) = Load(:, il);
        [PLoss(i, il), Vbus, Isec(i, il, :)] = DLF(Strdata, Loaddata); % % % Running load flow
        PenaltyVoltageL(i, il) = PenV(Vbus); % % % Calculating amount of penalties
    end

    PenaltyVoltage(i) = sum(PenaltyVoltageL(i, :), 2);
    f(i) = Ke * (T_OffPeak * PLoss(i, 1) + T_Medium * PLoss(i, 2) + T_Peak * PLoss(i, 3)) + Kp * PLoss(i, 1) + Total_Cap_Price; % % % Calculating objective function
    f(i) = f(i) + PenaltyVoltage(i);
end

PBest = p; PBestValue = f;
[GTeacherValue, index] = min(f);
GTeacher = PBest(index, :); % % % The best solution
Xmean = mean(p);

for k = 1:Iter
    k
    [f, p, GTeacher, GTeacherValue, Xmean, PenaltyVoltage, PenaltyVoltageBest] = UpdateSolutions(GTeacher, p, Xmean, f, PenaltyVoltage, LoadOffPeak, LoadMedium, LoadPeak);
    %%% Generating new solutions
    fff(k) = GTeacherValue;
end

toc
ij = 1:Iter;
hold on
plot(ij, fff, 'r')
%% Function of defining input parameters
function definParameters()
    global No_Cap_Type Cap_MVar NBus No_pop Iter Cap_Price
    VLoadMax VLoadMin PF Loaddata Strdata pMax pMin Ke Kp T
    Kl T_OffPeak T_Medium T_Peak NLoadLevel
    No_Cap_Type = 7; % % % Number of capacitor types
    Cap_MVar = 4 * [0 150 300 450 600 900 1200]; % % % MVar of
    capacitors

    Cap_Price = 4 * [0 750 975 1140 1320 1650 2040]; % % % Price
    of capacitors
    No_pop = 100; % % % Number of population
    Iter = 200; % % % Iteration number
    VLoadMax = 1.1; % % % Upper voltage bound
    VLoadMin = 0.9; % % % Lower voltage bound
    PF = 5000; % % % Penalty factor
    %%% Bus P Q
    Loaddata = [
            2  1840 460
            3  980 340
            4  1790 446
            5  1598 1840
            6  1610 600
            7  780 110
            8  1150 60
            9  980 130
            10 1640 200
            ];
    %%% From Bus To Bus Length R X Imax Cap
    Strdata = [1 2 1 0.1233 0.4126 0 0
            2 3 1 0.014 0.6051 0 0
            3 4 1 0.7463 1.205 0 0
            4 5 1 0.6984 0.6084 0 0
            5 6 1 1.9831 1.7276 0 0
            6 7 1 0.9053 0.7886 0 0
            7 8 1 2.0552 1.164 0 0
            8 9 1 4.7953 2.716 0 0
            9 10 1 5.3434 3.0264 0 0
            ];
    NBus = size(Loaddata, 1) + 1; % % % Number of buses
    pMax = No_Cap_Type; % % % Maxiumum bound of populations
    pMin = 1; % % % Minimum bound of populations
    Ke = 0.06; % % % Coefficent of energy loss
    Kp = 300; % % % coefficent of power loss
    T = 8760; % % % time period
    Kl = 168; % % %
    T_OffPeak = 3000; % % % Off peak hours
    T_Medium = 5300; % % % Medium load hours
    T_Peak = 460; % % % Peak hours
    NLoadLevel = 3; % % % Number of load levels
    %% Function of updating solutions
    function [f, p, GTeacher, GTeacherValue, Xmean,
        PenaltyVoltage, PenaltyVoltageBest] =
    UpdateSolutions(GTeacher, p, Xmean, f,
    PenaltyVoltage, LoadOffPeak, LoadMedium, LoadPeak)
    global Cap_Price Ke Kp Strdata Loaddata T_OffPeak
    T_Medium T_Peak NLoadLevel No_Cap_Type

    for i = 1:size(p, 1) % % % % % % % % % % % % % % % Teacher phase % % % % % % % % % % % % % % %
        TF = round(1 + rand); pnew(i, :) = p(i, :) +
        rand(1, size(p, 2)) .* (GTeacher - TF * Xmean);
        pnew(i, :) = round(pnew(i, :));

        for k = 1:size(p, 2)

            if pnew(i, k) > No_Cap_Type
                pnew(i, k) = No_Cap_Type;
            elseif pnew(i, k) < 1
                pnew(i, k) = 1;
            end

        end

        pop(i, :) = Cap_Mvar_determine(pnew(i, :));
        Load(:, 1) = LoadOffPeak - (pop(i, :))';
        Load(:, 2) = LoadMedium - (pop(i, :))';
        Load(:, 3) = LoadPeak - (pop(i, :))';
        Total_Cap_Price = sum(Cap_Price((pnew(i, :))));

        for il = 1:NLoadLevel
            Loaddata(:, 3) = Load(:, il);
            [PLoss(i, il), Vbus, Isec(i, il, :)] = DLF(Strdata, Loaddata);
            PenaltyVoltageL(i, il) = PenV(Vbus);
        end

        PenaltyVoltageNew(i) = sum(PenaltyVoltageL(i, :), 2);
        fnew(i) = Ke * (T_OffPeak * PLoss(i, 1) + T_Medium * PLoss(i, 2) + T_Peak * PLoss(i, 3)) + Kp * PLoss(i, 1) + Total_Cap_Price;
        fnew(i) = fnew(i) + PenaltyVoltageNew(i);

        if fnew(i) < f(i)
            p(i, :) = pnew(i, :); f(i) = fnew(i);
            PenaltyVoltage(i) =
            PenaltyVoltageNew(i);
        end

        j = round(1 + rand * (i - 1)); % % % % % % % % % % % Student phase % % % % % % % % % % % % % % % % % % % % %

        if j ~= i

            if f(i) < f(j)
                pnew(i, :) = p(i, :) + rand(1, size(p, 2)) .* (p(i, :) - p(j, :));
            else
                pnew(i, :) = p(i, :) + rand(1, size(p, 2)) .* (p(j, :) - p(i, :));
            end

            pnew(i, :) = round(pnew(i, :));

            for k = 1:size(p, 2)

                if pnew(i, k) > No_Cap_Type
                    pnew(i, k) = No_Cap_Type;
                elseif pnew(i, k) < 1
                    pnew(i, k) = 1;
                end

            end

            pop(i, :) = Cap_Mvar_determine(pnew(i, :));
            Load(:, 1) = LoadOffPeak - (pop(i, :))';
            Load(:, 2) = LoadMedium - (pop(i, :))';
            Load(:, 3) = LoadPeak - (pop(i, :))';
            Total_Cap_Price = sum(Cap_Price((pnew(i, :))));

            for il = 1:NLoadLevel
                Loaddata(:, 3) = Load(:, il);
                [PLoss(i, il), Vbus, Isec(i, il, :)] = DLF(Strdata, Loaddata);
                PenaltyVoltageL(i, il) = PenV(Vbus);
            end

            PenaltyVoltageNew(i) = sum(PenaltyVoltageL(i, :), 2);
            fnew(i) = Ke * (T_OffPeak * PLoss(i, 1) + T_Medium * PLoss(i, 2) + T_Peak * PLoss(i, 3)) + Kp * PLoss(i, 1) + Total_Cap_Price;
            fnew(i) = fnew(i) + PenaltyVoltageNew(i);

            if fnew(i) < f(i)
                p(i, :) = pnew(i, :);
                f(i) = fnew(i);
                PenaltyVoltage(i) =
                PenaltyVoltageNew(i);
            end

        end

    end

[GTeacherValue, index] = min(f);
GTeacher = p(index, :);
PenaltyVoltageBest = PenaltyVoltage(index(1));
Xmean = mean(p);
%% Function of backward forward load flow
function [PLoss, Vbus, Isec] = DLF(Strdata, Loaddata)
    % Strdata->> 1-from/2-to/3-Length(km)/4-R(ohm/km)/5-
    X(ohm / km) / 6 - Imax(Amp) / 7 - Capacitor (kvar)
    % Loaddata->> 1-bus/2-P(kw)/3-Q(kw)
    PLoss = [];
    Nsec = length(Strdata(:, 1)); %Number of sections (or to buses)
    Vbase = 23000; %V base of the system (v)
    Isec = zeros(Nsec, 1);
    Vbus = Vbase * ones(Nsec, 1);
    Cbus = zeros(Nsec, 1);
    Sbus = zeros(Nsec, 1);
    Rsec = Strdata(:, 4) .* Strdata(:, 3);
    Xsec = Strdata(:, 5) .* Strdata(:, 3);
    Zsec = Rsec + i * Xsec;
%=============================Algorithm================ ===
BI = zeros(Nsec, Nsec + 1);
BI(1, 1) = 1;
BV = BI;
for k = 1:Nsec
    BI(:, Strdata(k, 2)) = BI(:, Strdata(k, 1));
    BI(k, Strdata(k, 2)) = 1;
    BV(:, Strdata(k, 2)) = BV(:, Strdata(k, 1));
    BV(k, Strdata(k, 2)) = Zsec(k);
end

        BI(:, 1) = [];
        BV(:, 1) = [];
        BV = BV.';
        Cbus(Strdata(:, 2)) = Strdata(:, 7) * 1000;
        Sbus(Loaddata(:, 1)) = (Loaddata(:, 2) + i * Loaddata(:, 3)) * 100
        0; Cbus(1, :) = [];
        Sbus(1, :) = [];
        Iter = 0;

        NERROR = 1;
        S_bus = Sbus - i * (Cbus .* (Vbus / Vbase).^2); %for P
        constant Ibus = conj(S_bus ./ (sqrt(3) * Vbus));

        while ((Iter < 100) && (NERROR > 1e-5))
            Iter = Iter + 1;
            OldIbus = Ibus;
            VD = sqrt(3) * (BV * BI) * Ibus;
            Isec = BI * Ibus;
            Vbus = Vbase - VD;
            S_bus = Sbus - i * (Cbus .* (Vbus / Vbase).^2); %for P
            constant Ibus = conj(S_bus ./ (sqrt(3) * Vbus));
            NERROR = max(max(abs(Ibus - OldIbus)));
        end

%======================================================== =
LossSec = 3 * abs(Isec).^2 .* (Rsec) / 1000;
PLoss = sum(LossSec);
Vbus = abs(Vbus) / Vbase; % voltage of to buses return
%% Function of allocating MVAr to the generated population
function pop = Cap_Mvar_determine(p)
    global Cap_MVar NBus
    for i = 1:size(p, 1)
        pop_row = p(i, :);
        pop_row_MVar = zeros(1, NBus - 1);
        for j = 1:NBus - 1
            pop_row_MVar(j) = Cap_MVar(pop_row(j));
        end
        pop(i, :) = pop_row_MVar;
    end

%% Function of applying upper and lower bounds of population
function p = ApplyingConstraint(p)
    global No_Cap_Type
    for i = 1:size(p, 1)
        for j = 1:size(p, 2)
            if p(i, j) > No_Cap_Type
                p(i, j) = No_Cap_Type;
            elseif p(i, j) < 1
                p(i, j) = 1;
            end
        end
    end

%% Function of penalizing infeasible solutions
function PenaltyVoltage = PenV(Vbus)
    global VLoadMax VLoadMin PF
    for i = 1:size(Vbus, 1)
        if (Vbus(i) > VLoadMax) || (Vbus(i) < VLoadMin)
            Penalty(i) = PF;
        else
            Penalty(i) = 0;
        end
    end

PenaltyVoltage = sum(Penalty);
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
