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