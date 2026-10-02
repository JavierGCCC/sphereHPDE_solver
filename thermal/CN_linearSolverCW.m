function [T_unit, T_0, info]=CN_linearSolverCW(F, Q_unit)
t_ref=tic;

Q_unit_a = Q_unit; Q_unit_a(end) = 0;
Q_unit_b = zeros(size(Q_unit)); Q_unit_b(end) = Q_unit(end);

X = F\[Q_unit_a, Q_unit_b];

T_unit = X(:,1);
T_0    = X(:,2);


elapsed_time=toc(t_ref);
info.elapsedTime=elapsed_time;

end