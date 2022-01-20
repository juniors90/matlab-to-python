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