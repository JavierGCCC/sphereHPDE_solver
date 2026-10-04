function dt = shift_convolution(tau1, tau2, gamma, display)
if nargin < 4 || isempty(display)
    display = false;
end
b = 1 / tau1;
a = pi / tau2^2;
f_raw = @(t) exp(-b * t) .*(erfc(sqrt(a) * (-t + b / (2 * a))));
 % find t_max (where f(t) is maximum)
    t_grid = linspace(0, 5*tau2, 1e6); %guarantees precision: 1e-5tau2
    [~, idx] = max(f_raw(t_grid));
    t_max = t_grid(idx);
    f_max = f_raw(t_max);

    %rescaling. 
    f=@(t) f_raw(t)/f_max;

    %equation for displacement.
    target = gamma;
    eq = @(dt) f(-dt) - target;

    % rooots (dt > 0)
    dt = fzero(eq, [t_max,10*tau2]);

 % ======= Visualization =======
if display 
    figure;
    hold on;
    plot(t_grid, f_raw(t_grid)/f_max, 'b', 'LineWidth', 1.3); 
    yline(target, '--k', sprintf('γ = %.1e', gamma));
    xline(dt, '--r', sprintf('Δt = %.2e s', dt),...
        'LabelVerticalAlignment','bottom');
    plot(t_grid, f_raw(t_grid - dt)/f_max, 'r--', 'LineWidth', 1.3);
    xline(t_max, ':g', sprintf('t_{max} = %.2e s', t_max));
    xlabel('t [s]');
    ylabel('f(t)/f_{max}');
    title('Desplazamiento temporal Δt');
    legend('f(t)/f_{max}','γ','Δt','f(t+Δt)','t_{max}','Location','best');
    grid on; box on;
    hold off;
end
end

