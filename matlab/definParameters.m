%% Function of defining input parameters
function definParameters()

global No_Cap_Type Cap_MVar NBus No_pop Iter Cap_Price VLoadMax VLoadMin PF Loaddata Strdata pMax pMin Ke Kp T Kl T_OffPeak T_Medium T_Peak NLoadLevel

No_Cap_Type = 7; % % % Number of capacitor types
Cap_MVar = 4 * [0 150 300 450 600 900 1200]; % % % MVar of capacitors
Cap_Price = 4 * [0 750 975 1140 1320 1650 2040]; % % % Price of capacitors
No_pop = 100; % % % Number of population
Iter = 200; % % % Iteration number
VLoadMax = 1.1; % % % Upper voltage bound
VLoadMin = 0.9; % % % Lower voltage bound
PF = 5000; % % % Penalty factor
%%% Bus P Q
Loaddata = [2 1840 460
            3 980 340
            4 1790 446
            5 1598 1840
            6 1610 600
            7 780 110
            8 1150 60
            9 980 130
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