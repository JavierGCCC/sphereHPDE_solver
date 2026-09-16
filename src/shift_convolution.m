function dt = shift_convolution(tau1, tau2, gamma, display)
b = 1 / tau1;
a = pi / tau2^2;
 f_raw = @(t) exp(-b * t) .*(erfc(sqrt(a) * (-t + b / (2 * a))));
 % encontrar t_max (donde f(t) es máximo)
    t_grid = linspace(0, 5*tau2, 1e6); %garantiza precision 1e-5tau2
    [~, idx] = max(f_raw(t_grid));
    t_max = t_grid(idx);
    f_max = f_raw(t_max);

    %reescalado. 
    f=@(t) f_raw(t)/f_max;
    

     % ecuación para el desplazamiento
    target = gamma;
    eq = @(dt) f(-dt) - target;

    % búsqueda de raíz (dt > 0)
    dt = fzero(eq, [t_max,10*tau2]);
    


if display 
    % ======= GRAFICA =======
    figure;
    hold on;
    plot(t_grid, f_raw(t_grid)/f_max, 'b', 'LineWidth', 1.3); % función normalizada
    yline(target, '--k', sprintf('γ = %.1e', gamma));
    xline(dt, '--r', sprintf('Δt = %.2e s', dt), 'LabelVerticalAlignment','bottom');
    plot(t_grid, f_raw(t_grid - dt)/f_max, 'r--', 'LineWidth', 1.3); % función desplazada
    xline(t_max, ':g', sprintf('t_{max} = %.2e s', t_max));
    xlabel('t [s]');
    ylabel('f(t)/f_{max}');
    title('Desplazamiento temporal Δt');
    legend('f(t)/f_{max}','γ','Δt','f(t+Δt)','t_{max}','Location','best');
    grid on; box on;
    hold off;
end
end


% gamma=1e-3;
% tau1=1.7e-12;
% tau2=30*tau1;
% 
% dt=shift_convolution(tau1,tau2,gamma,true);
