clear;close all;clc

tic
c0 = 340;
rho0 = 1.225;

nx = [0 0 1 1 0 1 0 2 2 1]';
ny = [0 1 0 1 2 2 3 0 1 3]';

lx = 0.1; % Cavity height [m]
ly = 0.2;  % Cavity width  [m]
lz = 0.15; % Cavity depth  [m]

fn = c0/2 .* sqrt((nx./lx).^2 + (ny./ly).^2); % Modal frequencies of infinite long rect duct


theta = 0;
% f = [100, 103, 106, 109, 112, 115, 118, 122, 125, 128, 132, 136, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 206, 212, 218, 224, 230, 236, 243, 250, 258, 265, 272, 280, 290, 300, 307, 315, 325, 335, 345, 355, 365, 375, 387, 400, 412, 425, 437, 450, 462, 475, 487, 500, 515, 530, 545, 560, 580, 600, 615, 630, 650, 670, 690, 710, 730, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1e3, 1.03e3, 1.06e3, 1.09e3, 1.12e3, 1.15e3, 1.18e3, 1.22e3, 1.25e3, 1.28e3, 1.32e3, 1.36e3, 1.4e3, 1.45e3, 1.5e3, 1.55e3, 1.6e3, 1.65e3, 1.7e3, 1.75e3, 1.8e3, 1.85e3, 1.9e3, 1.95e3, 2e3, 2.06e3, 2.12e3, 2.18e3, 2.24e3, 2.3e3, 2.36e3, 2.43e3, 2.5e3, 2.58e3, 2.65e3, 2.72e3, 2.8e3, 2.9e3, 3e3, 3.07e3, 3.15e3, 3.25e3, 3.35e3, 3.45e3, 3.55e3, 3.65e3, 3.75e3, 3.87e3, 4e3, 4.12e3, 4.25e3, 4.37e3, 4.5e3, 4.62e3, 4.75e3, 4.87e3, 5e3];
f = [50, 51.5, 53, 54.5, 56, 58, 60, 61.5, 63, 65, 67, 69, 71, 73, 75, 77.5, 80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5, 100, 103, 106, 109, 112, 115, 118, 122, 125, 128, 132, 136, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 206, 212, 218, 224, 230, 236, 243, 250, 258, 265, 272, 280, 290, 300, 307, 315, 325, 335, 345, 355, 365, 375, 387, 400, 412, 425, 437, 450, 462, 475, 487, 500, 515, 530, 545, 560, 580, 600, 615, 630, 650, 670, 690, 710, 730, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1e3, 1.03e3, 1.06e3, 1.09e3, 1.12e3, 1.15e3, 1.18e3, 1.22e3, 1.25e3, 1.28e3, 1.32e3, 1.36e3, 1.4e3, 1.45e3, 1.5e3, 1.55e3, 1.6e3, 1.65e3, 1.7e3, 1.75e3, 1.8e3, 1.85e3, 1.9e3, 1.95e3, 2e3, 2.06e3, 2.12e3, 2.18e3, 2.24e3, 2.3e3, 2.36e3, 2.43e3, 2.5e3, 2.58e3, 2.65e3, 2.72e3, 2.8e3, 2.9e3, 3e3, 3.07e3, 3.15e3, 3.25e3, 3.35e3, 3.45e3, 3.55e3, 3.65e3, 3.75e3, 3.87e3, 4e3, 4.12e3, 4.25e3, 4.37e3, 4.5e3, 4.62e3, 4.75e3, 4.87e3, 5e3, 5.15e3, 5.3e3, 5.45e3, 5.6e3, 5.8e3, 6e3, 6.15e3, 6.3e3, 6.5e3, 6.7e3, 6.9e3, 7.1e3, 7.3e3, 7.5e3, 7.75e3, 8e3];
w = 2*pi*f;
k = w./c0;
S = lx*ly;

kz = sqrt(k.^2 - (pi/lx)^2 - (pi/ly)^2);

%% Simplified TL taking only first mode into consideration
%%RIGID GROUND

TLg = zeros(length(k), 1);

parfor i=1:length(k)

% % V1 FUNC
% Normalized Radiation Impedance of a baffled rigid piston on a rigid ground
   eps = 1e-6; % Small regularization factor to avoid singularity

   ki = k(i);

f4g = @(x1, y1, x2, y2) ...
    (exp(-1i*ki*sqrt((x1-x2).^2 + (y1-y2).^2)) ./ sqrt((x1-x2).^2 + (y1-y2).^2 + eps) + ...
     exp(-1i*ki*sqrt((x1+x2).^2 + (y1-y2).^2)) ./ sqrt((x1+x2).^2 + (y1-y2).^2 + eps));
  
   q4g = integralN(f4g,0, lx, 0, ly, 0, lx, 0, ly,'AbsTol',1e-5,'RelTol',1e-3);
   
   Zg = (1i*ki)/(2*pi*S)* q4g;

% Transmission Loss
TLg(i) = 10*log10((cos(theta)*abs((Zg+1)^2 * exp(1i*ki*lz) - (Zg-1)^2 * exp(-1i*ki*lz))^2) / (64*real(Zg)));

end




%% Simplified TL taking only first mode into consideration
%%100% ABSORPTIVE GROUND

TLf = zeros(length(k), 1);

parfor i=1:length(k)

% % V1 FUNC
% Normalized Radiation Impedance of a baffled rigid piston on a absorptive ground
   eps = 1e-6; % Small regularization factor to avoid singularity

   ki = k(i);

f4f = @(x1, y1, x2, y2) ...
    (exp(-1i*ki*sqrt((x1-x2).^2 + (y1-y2).^2)) ./ sqrt((x1-x2).^2 + (y1-y2).^2 + eps));
 
   q4f = integralN(f4f,0, lx, 0, ly, 0, lx, 0, ly,'AbsTol',1e-5,'RelTol',1e-3);
   
   Zf = (1i*ki)/(2*pi*S)* q4f;

% Transmission Loss
TLf(i) = 10*log10((cos(theta)*abs((Zf+1)^2 * exp(1i*ki*lz) - (Zf-1)^2 * exp(-1i*ki*lz))^2) / (16*real(Zf)));

end



toc
%%


% Plot the full curve without markers
figure;
semilogx(f, TLf, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
hold on;
semilogx(f, TLg, '-.', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);



% Add markers at specific points
marker_indices = 1:10:length(f); % Plot markers every 10th point
semilogx(f(marker_indices), TLg(marker_indices), '^', ...
    'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'MarkerFaceColor', [0.8500 0.3250 0.0980]);

semilogx(f(marker_indices), TLf(marker_indices), '^', ...
    'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'MarkerFaceColor', [0 0.4470 0.7410]);

grid on
ylim([-15 10])
ylabel('TL [dB]','Interpreter','latex')
xlabel('Frequency [Hz]','Interpreter','latex')
set(gca, 'FontSize', 14, 'Linewidth', 1,'TickLabelInterpreter','latex');
legend('Eq. (25) + Eq. (24), absorptive ground','Eq. (23) + Eq. (22), rigid ground','Interpreter', 'latex', 'Location', 'southwest')
