clear;
clc;
n = 1;

% Parameter Settings
for alpha = 4e-3:-1e-4:-4e-3
    G = 1000;
    P_cell = 200;       % Intracellular hydrostatic pressure
    miu = 0.4;
    E_ecm = 1300;
    can = 0.3;
    lambda = 1e-2;      % Contractility of apical actin band
    % alpha = 1e-3;     % Lateral membrane tension (Loop variable)
    gamma = -5e-4;      % Basal adhesion tension energy density
    
    V0 = 2*sqrt(3)*1e-15;   % Cell volume with parameters
    A = 5*10^-23;           % Gaussian non-adhesive matrix modulus

    x = linspace(5e-6, 30e-6, 200);
    y = linspace(5e-6, 30e-6, 200);
    % [r_1, r_2] = meshgrid(x, y);

    for j = 1:length(x)
        for i = j:length(x)
            h(i,j) = V0 ./ (x(i)^2 + x(i).*y(j) + y(j).^2); % Cell height
            c(i,j) = (x(i) - y(j)) ./ (y(j).*h(i,j));       % Apical curvature
            r(i,j) = 1./c(i,j) + h(i,j);
            c_1(i,j) = 1./r(i,j);                           % Basal curvature
            
            % Pressure term
            K_press(i,j) = (P_cell + can.*sqrt(abs(c(i,j)))).^2 / (24*G);
            K_press(i,i) = 1;

            % Total Energy Function
            f(i,j) = 0.5*sqrt(3)*lambda.*y(j).^2 ...
                   + sqrt(3)/2*gamma.*x(i).^2 ...
                   - sqrt(3)*alpha.*(x(i)+y(j)).*h(i,j) ...
                   + 2*A./(x(i).*y(j)) ...
                   + A./h(i,j).^2 ...
                   + K_press(i,j)*r(i,j).*(x(i).^6)./(y(j).^4);
        end
    end

    % Mask upper triangle to ensure physical constraints (e.g., R_b >= R_a or similar logic)
    for i = 1:length(x)
        for j = i:length(x)
            f(i,j) = NaN;
        end
    end

    [val, liner] = min(f(:));
    [a, b] = ind2sub(size(f), liner); % Convert linear index to subscripts
    
    a1(n) = a;
    b1(n) = b;
    R_1(n) = x(a);
    R_2(n) = y(b);
    rk(n) = r(a,b);
    H(n) = V0 / (R_1(n).^2 + R_1(n).*R_2(n) + R_2(n).^2); % Cell height
    cri(n) = (R_1(n) - R_2(n)) ./ (R_2(n).*H(n));         % Apical curvature
    f1(n) = f(a,b);
    k(n) = alpha + 4e-3;
    u(n) = 1/cri(n);                                      % Lumen radius
    q(n) = floor(2*pi*u(n)/R_2(n));                       % Number of cells in 2D ring
    q_t(n) = floor(4/3*pi*(rk(n)^3 - u(n)^3)/10^-15);     % Total cell count in 3D shell
    m(n) = 4*pi*rk(n)^2 / (q_t(n)*sqrt(3)/2*R_1(n)^2);    % Area incompatibility metric
    
    n = n + 1;
end

% Assume we have the following data:
innerRadii = u; 
outerRadii = rk;     
numDivisions = q; 

% ========== Publication-Quality Visualization ==========
figure('Color', [0.98 0.98 0.98], 'Position', [100, 100, 800, 800]); % Light gray background
videoFileName = 'cell_grow.avi';
fps = 20;

vidOut = VideoWriter(videoFileName, 'MPEG-4');
vidOut.Quality = 100;
vidOut.FrameRate = fps;
open(vidOut);

% Pre-calculate global display range (ensure absolute stability)
max_display_radius = 3.2 * 1.2 * max([outerRadii, innerRadii]);

% Alpha range (consistent with original loop)
alpha_min = -4e-3;  % Note: Original loop is 4e-3:-1e-4:-4e-3
alpha_max = 4e-3;
width_min = 1.0;    % Red line width when alpha is minimum
width_max = 8.0;    % Red line width when alpha is maximum

for i = length(q) :-1:1
    clf('reset');
    set(gcf, 'Color', [0.98 0.98 0.98]); % Light gray background
    
    % Increase sampling points for extreme smoothness
    theta = linspace(0, 2*pi, 1000); 
    sita = linspace(0, 2*pi, numDivisions(i)+1);
    
    % === Core Structure (Minimalist Design) ===
    hold on;
    
    % Calculate linewidth based on alpha/k value
    % Note: k was defined as alpha + 4e-3, so k ranges from 0 to 8e-3.
    % Adjust scaling factor as needed for visual preference.
    linewidth = 1 + 1 * k(length(q)+1-i) / 1e-3;
    
    % Use anti-aliasing smooth drawing
    % Inner Ring (Apical/Lumen side) - Red
    plot(2*innerRadii(i)*cos(theta), 2*innerRadii(i)*sin(theta), ...
        'Color', [0.8 0 0], 'LineWidth', linewidth);
    
    % Outer Ring (Basal side) - Blue
    plot(2*outerRadii(i)*cos(theta), 2*outerRadii(i)*sin(theta), ...
        'Color', [0.0 0.4470 0.7410], 'LineWidth', 2.5, 'HandleVisibility', 'off');
    
    % Lateral Connections - Green radial lines
    for j = 1:numDivisions(i)
        plot([2*innerRadii(i)*cos(sita(j)), 2*outerRadii(i)*cos(sita(j))], ...
             [2*innerRadii(i)*sin(sita(j)), 2*outerRadii(i)*sin(sita(j))], ...
            'Color', [0.25 0.75 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    
    % === Key Fix: Remove bottom text and add border ===
    % 1. Set axis limits
    axis equal;
    xlim([-max_display_radius, max_display_radius]);
    ylim([-max_display_radius, max_display_radius]);
    
    % 2. Get current axes object
    ax = gca;
    
    % 3. Correctly remove axis labels
    delete(ax.XLabel);   % Delete X-axis label
    delete(ax.YLabel);   % Delete Y-axis label
    set(ax, 'XTick', [], 'YTick', []); % Remove ticks
    
    % 4. Set border color to light gray (professional publication standard)
    set(ax, 'XColor', [0.3, 0.3, 0.3], 'YColor', [0.3, 0.3, 0.3]);
    
    % 5. Extra protection: Remove any residual text objects except title
    textObjs = findobj(ax, 'Type', 'text');
    for j = 1:length(textObjs)
        % Keep only the title if necessary, or delete all and recreate title below
        delete(textObjs(j));
    end
    
    % Title (Fixed format, showing only cell number)
    title(sprintf('Cell Number = %d', q_t(i)), ...
        'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold', ...
        'Color', [0.1 0.1 0.1], 'Position', [0, 0.92*max_display_radius, 0], ...
        'HandleVisibility', 'off');
    
    drawnow;
    writeVideo(vidOut, getframe(gcf));
end

close(vidOut);
clearvars -except videoFileName;

% --- Commented out plotting commands for reference ---
% scatter(-1*k, R_1, 'm');
% % axis([-8e-3, 8e-3, 6e-6, 16e-6])
% xlabel('Lateral Tension');
% ylabel("Basal Edge Length");

% scatter(-1*k, R_1/R_2, 'm');
% % axis([-8e-3, 8e-3, 6e-6, 16e-6])
% xlabel('Lateral Tension');
% ylabel("Basal/Apical Edge Ratio");

% scatter(-1*k, R_2, "g");
% xlabel('Lateral Tension');
% ylabel("Apical Edge Length");

% scatter(-k, H, "cyan");
% xlabel('Lateral Tension');
% ylabel("Cell Height");

% scatter(-k, cri, "red");
% xlabel('Lateral Tension');
% ylabel('Apical Curvature');

% scatter(k, rk, "y");
% xlabel('Lateral Tension');
% ylabel('Organoid Outer Radius');

% scatter(k, u, "black");
% xlabel('Lateral Tension');
% ylabel('Organoid Lumen Radius');