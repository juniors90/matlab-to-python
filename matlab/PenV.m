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