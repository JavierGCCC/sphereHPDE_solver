
% Primer plot: data_peq
x_original = data_peq(:,1);
x_interp = linspace(min(x_original), max(x_original), 1000);  % Eje interpolado

figure;
hold on;
for i = 1:3
    y_original = data_peq(:,i+1);
    y_interp = interp1(x_original, y_original, x_interp, 'spline');
    plot(x_interp, y_interp, 'LineWidth', 2);
end
%title('Interpolated Curves - data\_peq');
xlabel('\lambda (nm)'); 
ylabel('\sigma_{abs} (m^2)');
legend('\delta_{Au}=5 nm','\delta_{Au}=10 nm','\delta_{Au}=15 nm');
legend boxoff;
hold off;
set(gca,'FontSize',18); set(gca,'LineWidth',2);box on;

% Segundo plot: data_gran
x_original = data_gran(:,1);
x_interp = linspace(min(x_original), max(x_original), 1000);  % Eje interpolado

figure;
hold on;
for i = 1:3
    y_original = data_gran(:,i+1);
    y_interp = interp1(x_original, y_original, x_interp, 'spline');
    plot(x_interp, y_interp, 'LineWidth', 2);
end
%title('Interpolated Curves - data\_gran');
xlabel('\lambda (nm)'); 
ylabel('\sigma_{abs} (m^2)');
hold off;
set(gca,'FontSize',18); set(gca,'LineWidth',2); box on;

