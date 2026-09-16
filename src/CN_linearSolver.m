function [T, info]=CN_linearSolver(F_d, F_a, Q, Tamb, M, L, U)
t_ref=tic;
s=size(Q);
%Initial condition.
T=zeros(s(1),s(2));
T(:,1)=Tamb*ones(s(1),1);
switch nargin
    case 4
        h = waitbar(0, 'Solving heat transfer...');
        for i=1:s(2)
            T(:,i+1)=F_d\(F_a*T(:,i)+Q(:,i));
            waitbar(i / s(2), h);
        end
        close(h);
        elapsed_time=toc(t_ref);
    case 5
        F_d_cond=M*F_d;
        F_a_cond=M*F_a;
        Q_cond=M*Q;
        h = waitbar(0, 'Solving heat transfer...');
        for i=1:s(2)
            T(:,i+1)=F_d_cond\(F_a_cond*T(:,i)+Q_cond(:,i));
            waitbar(i / s(2), h);
        end
        close(h);
        elapsed_time=toc(t_ref);
    case 7
        U_cond=M*U;
        h = waitbar(0, 'Solving heat transfer...');
        for i=1:s(2)
            S=L\(F_a*T(:,i)+Q(:,i));
            T(:,i+1)=U_cond\(M*S);
            waitbar(i / s(2), h);
        end
        close(h);
        elapsed_time=toc(t_ref);
    otherwise
        disp('Incorrect number of inputs');
        elapsed_time=NaN;
end
info.elapsedTime=elapsed_time;

end