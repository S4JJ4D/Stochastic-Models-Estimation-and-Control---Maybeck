res = cell(1,10);

for k=1:10
    res{k} = CEP(k);
end

data = vertcat(res{:});
dataTable = array2table(data, 'VariableNames', ...
    {'k (sx/sy)', 'CEP/sx', 'CEP1/sx', 'CEP2/sx', 'CEP1/sx (% error)', 'CEP2/sx (% error)'});
home;
dataTable

function res = CEP(k)

% cep1 = @(sx,sy) .588*(sx + sy);
% cep2 = @(sx,sy) .563*sx + .614*sy;

cep1_ratio = @(k) .588*(1 + 1/k);
cep2_ratio = @(k) .563 + .614/k;


gamma0 = cep1_ratio(k);

fprintf('\nInitial Guess: %.4f\n', gamma0);

options = optimoptions('fsolve','Display','iter');

[cep,fval] = fsolve(@(x) F(x,k), gamma0, options);


fprintf('\nCEP_1 = %.5f\nCEP_2 = %.5f\nR =     %.5f\n', ...
    cep1_ratio(k),cep2_ratio(k),cep);

fprintf('\nCEP_1 ERROR %% = %.5f\nCEP_2 ERROR %% = %.5f\n', ...
    abs((cep1_ratio(k)-cep)/cep)*100,abs((cep2_ratio(k)-cep)/cep)*100);


res = [k, cep, cep1_ratio(k), cep2_ratio(k), ...
    abs((cep1_ratio(k)-cep)/cep)*100,abs((cep2_ratio(k)-cep)/cep)*100];

    function f = F(gamma,k)
        f = integral(@(theta) exp(-1/2 * (gamma^2) ./ (cos(theta).^2 + 1/k^2 * sin(theta).^2)), 0, 2*pi) - pi;
    end
end