%% Function of backward forward load flow
function [PLoss, Vbus, Isec] = DLF(Strdata, Loaddata)
% Strdata->> 1-from/2-to/3-Length(km)/4-R(ohm/km)/5-X(ohm / km) / 6 - Imax(Amp) / 7 - Capacitor (kvar)
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
S_bus = Sbus - i * (Cbus .* (Vbus / Vbase).^2); %for P constant 
Ibus = conj(S_bus ./ (sqrt(3) * Vbus));

while ((Iter < 100) && (NERROR > 1e-5))
    Iter = Iter + 1;
    OldIbus = Ibus;
    VD = sqrt(3) * (BV * BI) * Ibus;
    Isec = BI * Ibus;
    Vbus = Vbase - VD;
    S_bus = Sbus - i * (Cbus .* (Vbus / Vbase).^2); %for P constant 
    Ibus = conj(S_bus ./ (sqrt(3) * Vbus));
    NERROR = max(max(abs(Ibus - OldIbus)));
end

%======================================================== =
LossSec = 3 * abs(Isec).^2 .* (Rsec) / 1000;
PLoss = sum(LossSec);
Vbus = abs(Vbus) / Vbase; % voltage of to buses 
return