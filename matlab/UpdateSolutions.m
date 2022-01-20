%% Function of updating solutions
function [f, p, GTeacher, GTeacherValue, Xmean, PenaltyVoltage, PenaltyVoltageBest] = UpdateSolutions(GTeacher, p, Xmean, f, PenaltyVoltage, LoadOffPeak, LoadMedium, LoadPeak)
global Cap_Price Ke Kp Strdata Loaddata T_OffPeak T_Medium T_Peak NLoadLevel No_Cap_Type

for i = 1:size(p, 1) % % % % % % % % % % % % % % % Teacher phase % % % % % % % % % % % % % % %
    TF = round(1 + rand); pnew(i, :) = p(i, :) + rand(1, size(p, 2)) .* (GTeacher - TF * Xmean);
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
        PenaltyVoltage(i) = PenaltyVoltageNew(i);
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
            PenaltyVoltage(i) = PenaltyVoltageNew(i);
        end

    end

end

[GTeacherValue, index] = min(f);
GTeacher = p(index, :);
PenaltyVoltageBest = PenaltyVoltage(index(1));
Xmean = mean(p);