classdef HATASIM_Tool_EN < matlab.apps.AppBase
    % Properties definition
    properties (Access = public)
        UIFigure               matlab.ui.Figure
        TopPanel               matlab.ui.container.Panel
        BottomPanel            matlab.ui.container.Panel
        ParamPanel             matlab.ui.container.Panel
        ResultPanel            matlab.ui.container.Panel
        
        % Parameter setup UI controls
        BalloonParamsPanel     matlab.ui.container.Panel
        CableParamsPanel       matlab.ui.container.Panel
        PrecisionParamsPanel   matlab.ui.container.Panel
        WindParamsPanel        matlab.ui.container.Panel
        DiagramAxes            matlab.ui.control.UIAxes
        
        % Result display UI controls
        PathAxes         matlab.ui.control.UIAxes
        TensionAxes            matlab.ui.control.UIAxes
        OutputPanel            matlab.ui.container.Panel
        OutputTextArea         matlab.ui.control.TextArea
        QueryResultPanel       matlab.ui.container.Panel
        QueryResultTextArea    matlab.ui.control.TextArea
        
        % Common controls
        SolveButton            matlab.ui.control.Button
        ParamNavButton         matlab.ui.control.Button
        ResultNavButton        matlab.ui.control.Button
        QueryLabel             matlab.ui.control.Label
        QueryEditField         matlab.ui.control.NumericEditField
        QueryButton            matlab.ui.control.Button
        
        % Parameter input fields
        D_balloonEdit          matlab.ui.control.NumericEditField
        m_envEdit              matlab.ui.control.NumericEditField
        m_loadEdit             matlab.ui.control.NumericEditField
        d_cableEdit            matlab.ui.control.NumericEditField
        rho_cableEdit          matlab.ui.control.NumericEditField
        KEdit                  matlab.ui.control.NumericEditField
        EEdit                  matlab.ui.control.NumericEditField
        h_step_initEdit        matlab.ui.control.NumericEditField
        h_step_factorEdit      matlab.ui.control.NumericEditField
        toleranceEdit          matlab.ui.control.NumericEditField
        v_refEdit              matlab.ui.control.NumericEditField
        h_refEdit              matlab.ui.control.NumericEditField
        h_groundEdit           matlab.ui.control.NumericEditField
        alphaEdit              matlab.ui.control.NumericEditField
        
        % Data storage
        CalculationData        struct
        MaxLOriginal           double
        MaxLWinded             double
        QueryLines             = {}
    end
    
    methods (Access = private)
        % Calculation function
        function calculate(app)
            % Get parameter values
            D_balloon = app.D_balloonEdit.Value;
            m_env = app.m_envEdit.Value;
            m_load = app.m_loadEdit.Value;
            d_cable = app.d_cableEdit.Value;
            rho_cable = app.rho_cableEdit.Value;
            K = app.KEdit.Value;
            E = app.EEdit.Value;
            h_step_init = app.h_step_initEdit.Value;
            h_step_factor = app.h_step_factorEdit.Value;
            tolerance = app.toleranceEdit.Value;
            v_ref = app.v_refEdit.Value;
            h_ref = app.h_refEdit.Value;
            h_ground = app.h_groundEdit.Value;
            alpha = app.alphaEdit.Value;
            
            % Constants
            g = 9.8;
            R_air = 287;
            R_He = 2077;
            P0 = 101325;
            T0 = 288.15;
            gamma = 0.0065;
            mu0 = 1.79e-5;
            S = 110.4;
            rho0_air = 1.225;
            
            % Part 1: Calculate max tether length without wind
            fun = @(h) rho0_air*((T0 - gamma*(h - h_ground))/T0).^((g/(R_air*gamma))-1) * (1/6)*pi*D_balloon^3 * g ...
                - ( (P0*((T0 - gamma*(h - h_ground))/T0).^(g/(R_air*gamma)))/(R_He*(T0 - gamma*(h - h_ground))) * (1/6)*pi*D_balloon^3 ...
                + m_env + m_load + rho_cable*(pi/4)*d_cable^2*(h - h_ground - D_balloon/2) ) * g;
            
            options = optimset('TolX', 1e-3, 'Display', 'off');
            max_L_no_wind = fzero(fun, [h_ground + D_balloon/2, 10000], options);
            app.MaxLOriginal = max_L_no_wind;
            
            % Part 2: Wind-affected adaptive step calculation
            all_rope_data = struct('H',[],'X_total',[],'L_original',[],'L_deformed',[],'rope_X',[],'rope_Z',[],'F_segments',[]);
            H_results = [];
            X_results = [];
            L_original_results = [];
            L_deformed_results = [];
            F_values = [];
            
            current_h = h_ground + D_balloon/2;
            h_step = h_step_init;
            last_valid_h = current_h;
            max_iterations = 999999999;
            
            % Progress dialog
            progress = uiprogressdlg(app.UIFigure, 'Title','Calculating...', 'Message','Computing ascent path, tether shape, tension...');
            
            for iter = 1:max_iterations
                if h_step < 0.01
                    break;
                end
                
                attempted_h = current_h + h_step;
                
                [valid, X_total, L_original, L_deformed, rope_X, rope_Z, F_segments] = ...
                    app.calculateAtHeight(attempted_h, h_ground, D_balloon, m_env, m_load, d_cable, ...
                    rho_cable, K, v_ref, h_ref, alpha, g, R_air, R_He, P0, T0, gamma, mu0, S, rho0_air, tolerance, E);
                
                if valid
                    new_data = struct('H', attempted_h, 'X_total', X_total, 'L_original', L_original, ...
                        'L_deformed', L_deformed, 'rope_X', rope_X, 'rope_Z', rope_Z, 'F_segments', F_segments);
                    all_rope_data(end+1) = new_data;
                    
                    H_results(end+1) = attempted_h;
                    X_results(end+1) = X_total;
                    L_original_results(end+1) = L_original;
                    L_deformed_results(end+1) = L_deformed;
                    F_values(end+1) = F_segments(1);
                    
                    last_valid_h = attempted_h;
                    current_h = attempted_h;
                else
                    current_h = last_valid_h;
                    h_step = h_step * h_step_factor;
                end
            end
            
            % Store results
            app.CalculationData = struct(...
                'all_rope_data', all_rope_data, ...
                'H_results', H_results, ...
                'X_results', X_results, ...
                'L_original_results', L_original_results, ...
                'L_deformed_results', L_deformed_results, ...
                'F_values', F_values);
            
            % Save max tethered length with wind
            if ~isempty(L_original_results)
                app.MaxLWinded = L_original_results(end);
            else
                app.MaxLWinded = 0;
            end
            
            % Part 3: Visualization
            % Path and tether shape plot
            cla(app.PathAxes);
            plot(app.PathAxes, X_results, H_results - h_ground, 'k:', 'LineWidth', 2.5);
            hold(app.PathAxes, 'on');
            
            if ~isempty(all_rope_data)
                final_data = all_rope_data(end);
                corrected_rope_X = final_data.X_total - final_data.rope_X;
                plot(app.PathAxes, corrected_rope_X, final_data.rope_Z - h_ground, 'r-', 'LineWidth', 3);
            end
            
            title(app.PathAxes, 'Ascent Path & Tether Shape');
            xlabel(app.PathAxes, 'Drift Distance (m)');
            ylabel(app.PathAxes, 'Height Above Ground (m)');
            grid(app.PathAxes, 'on');
            box(app.PathAxes, 'on');
            axis(app.PathAxes, 'equal');
            legend(app.PathAxes, {'Path', 'Final Tether Shape'}, 'Location', 'best');
            hold(app.PathAxes, 'off');
            
            % Tension distribution plot
            cla(app.TensionAxes);
            plot(app.TensionAxes, L_original_results, F_values, 'k:', 'LineWidth', 2.5);
            hold(app.TensionAxes, 'on');
            
            if ~isempty(all_rope_data)
                final_data = all_rope_data(end);
                original_length_positions = (0:K:(length(final_data.F_segments)-1)*K);
                plot(app.TensionAxes, original_length_positions, flip(final_data.F_segments), 'r-', 'LineWidth', 3);
            end
            
            title(app.TensionAxes, 'Tension Distribution');
            xlabel(app.TensionAxes, 'Tether Length from Ground (m)');
            ylabel(app.TensionAxes, 'Tension (N)');
            grid(app.TensionAxes, 'on');
            box(app.TensionAxes, 'on');
            legend(app.TensionAxes, {'Tension', 'Final Tension Distribution'});
            hold(app.TensionAxes, 'off');
            
            % ===== Part 4: Output results =====
            cable_area = pi/4 * d_cable^2;
            
            outputStr = sprintf('\n【Ground Altitude: %.2f m】\n', h_ground);
            outputStr = [outputStr sprintf('·Max Height (No Wind): %.2f m (AGL: %.2f m)\n', max_L_no_wind, max_L_no_wind - h_ground)];
            
            if ~isempty(L_original_results)
                outputStr = [outputStr sprintf('·Max Original Length (Wind): %.2f m\n', L_original_results(end))];
                outputStr = [outputStr sprintf('·Max Deformed Length (Wind): %.2f m\n', L_deformed_results(end))];
            end
            
            if ~isempty(L_original_results) && ~isempty(L_deformed_results)
                total_elongation = L_deformed_results(end) - L_original_results(end);
                elongation_ratio = total_elongation / L_original_results(end) * 100;
                outputStr = [outputStr sprintf('·Total Elongation: %.4f m\n', total_elongation)];
                outputStr = [outputStr sprintf('·Elongation Ratio: %.4f%%\n', elongation_ratio)];
            end
            
            if ~isempty(H_results)
                outputStr = [outputStr sprintf('·Max Altitude: %.2f m (AGL: %.2f m)\n', H_results(end), H_results(end) - h_ground)];
            end
            
            if ~isempty(X_results)
                outputStr = [outputStr sprintf('·Max Drift : %.2f m\n', X_results(end))];
            end
            
            if ~isempty(F_values)
                [max_F_value, max_idx] = max(F_values);
                max_stress_MPa = max_F_value / cable_area / 1e6;
                outputStr = [outputStr sprintf('·Max Tension During Ascent: %.2f N (%.2fMPa) (at %.2f m MSL, %.2f m AGL, original length %.2f m)\n',...
                    max_F_value, max_stress_MPa, H_results(max_idx), H_results(max_idx) - h_ground, L_original_results(max_idx))];
            end
            
            app.OutputTextArea.Value = outputStr;
            
            % Initialize query result area
            app.QueryResultTextArea.Value = 'Enter tether length and click Query to see results...';
            
            % Close progress dialog
            close(progress);
            
            % Switch to results panel
            app.switchToResultPanel();
        end
        
        % Helper calculation function
        function [valid, X_total, L_original, L_deformed, rope_X, rope_Z, F_segments] = ...
                calculateAtHeight(~, h, h_ground, D_balloon, m_env, m_load, d_cable, rho_cable, K, ...
                v_ref, h_ref, alpha, g, R_air, R_He, P0, T0, gamma, mu0, S, rho0_air, tolerance, E)
            
            valid = false;
            X_total = 0;
            L_original = 0;
            L_deformed = 0;
            rope_X = 0;
            rope_Z = h - D_balloon/2;
            F_segments = 0;
            
            T = T0 - gamma*h;
            if T <= 0
                return;
            end
            
            rho_air = rho0_air*(T/T0).^((g/(R_air*gamma))-1);
            mu = mu0*(T/T0).^1.5*(T0 + S)./(T + S);
            rho_He = (P0*(T/T0).^(g/(R_air*gamma)))/(R_He*T);
            
            V_He = (1/6)*pi*D_balloon^3;
            m_He = rho_He*V_He;
            W_balloon = (m_He + m_env + m_load)*g;
            F_b = rho_air*V_He*g;
            
            if F_b <= W_balloon
                return;
            end
            
            v = v_ref*(h/h_ref)^alpha;
            Re_balloon = rho_air*v*D_balloon/mu;
            if Re_balloon <= 1
                Cd_balloon = 24/Re_balloon;
            elseif Re_balloon <= 1e3
                Cd_balloon = 24/Re_balloon*(1 + 0.15*Re_balloon^0.687);
            elseif Re_balloon <= 2e5
                Cd_balloon = 0.47;
            else
                Cd_balloon = 0.2;
            end
            F_d = 0.5*rho_air*v^2*(pi/4)*D_balloon^2*Cd_balloon;
            
            X_total = 0;
            L_deformed = 0;
            h_accumulated = 0;
            W_rope = 0;
            n_segments = 0;
            
            rope_X = 0;
            rope_Z = h - D_balloon/2;
            F_segments = (F_b - W_balloon)/cos(atan(F_d/(F_b - W_balloon)));
            theta = atan(F_d/(F_b - W_balloon));
            F = F_segments(1);
            
            while (rope_Z(1) - h_accumulated) > h_ground + tolerance
                n_segments = n_segments + 1;
                
                delta_L = (F * K) / (E * (pi * d_cable^2 /4));
                L_segment = K + delta_L;
                
                h_seg = rope_Z(1) - h_accumulated - 0.5 * L_segment * cos(theta);
                v_seg = v_ref * (max(h_seg, 0.1)/h_ref)^alpha;
                
                d_new = d_cable * sqrt(K / (K + delta_L));
                Re_cable = rho_air * v_seg * cos(theta) * d_new / mu;
                
                if Re_cable <= 1
                    Cd_n = 8 * pi / (Re_cable * (2 - log(Re_cable)));
                elseif Re_cable <= 1e2
                    Cd_n = 1 + 10 / (Re_cable^0.667);
                elseif Re_cable <= 2e5
                    Cd_n = 1.20;
                else
                    Cd_n = 0.3;
                end
                
                D_n = 0.5 * rho_air * v_seg^2 * cos(theta)^2 * d_new * L_segment * Cd_n;
                segment_weight = rho_cable * (pi/4) * d_cable^2 * K * g;
                W_rope = W_rope + segment_weight;
                
                F_horiz = F * sin(theta) + D_n * cos(theta);
                F_vert = F * cos(theta) - segment_weight - D_n * sin(theta);
                
                if F_vert <= 0
                    return;
                end
                
                theta_new = atan(F_horiz / F_vert);
                F_new = F_horiz / sin(theta_new);
                
                X_total = X_total + L_segment * sin(theta);
                L_deformed = L_deformed + L_segment;
                h_accumulated = h_accumulated + L_segment * cos(theta);
                
                rope_X(end+1) = X_total;
                rope_Z(end+1) = rope_Z(1) - h_accumulated;
                F_segments(end+1) = F_new;
                
                theta = theta_new;
                F = F_new;
            end
            
            L_original = n_segments * K;
            valid = (rope_Z(1) - h_accumulated) <= h_ground + tolerance;
            if (F_b - (W_balloon + W_rope)) <= 0
                valid = false;
            end
        end
        
        % Panel switching functions
        function switchToParamPanel(app)
            app.ParamPanel.Visible = 'on';
            app.ResultPanel.Visible = 'off';
            app.ParamNavButton.BackgroundColor = [0.4 0.8 0.4]; % Green
            app.ResultNavButton.BackgroundColor = [0.9 0.9 0.9]; % Gray
            
            % Hide query section
            app.QueryLabel.Visible = 'off';
            app.QueryEditField.Visible = 'off';
            app.QueryButton.Visible = 'off';
        end
        
        function switchToResultPanel(app)
            app.ParamPanel.Visible = 'off';
            app.ResultPanel.Visible = 'on';
            app.ParamNavButton.BackgroundColor = [0.9 0.9 0.9]; % Gray
            app.ResultNavButton.BackgroundColor = [0.4 0.8 0.4]; % Green
            
            % Show query section
            app.QueryLabel.Visible = 'on';
            app.QueryEditField.Visible = 'on';
            app.QueryButton.Visible = 'on';
        end
        
        % Reset function
        function resetParameters(app)
            % Reset all parameters to defaults
            app.D_balloonEdit.Value = 5;
            app.m_envEdit.Value = 5;
            app.m_loadEdit.Value = 5;
            app.d_cableEdit.Value = 0.006;
            app.rho_cableEdit.Value = 970;
            app.KEdit.Value = 0.1;
            app.EEdit.Value = 95e9;
            app.h_step_initEdit.Value = 1;
            app.h_step_factorEdit.Value = 0.01;
            app.toleranceEdit.Value = 0.001;
            app.v_refEdit.Value = 3;
            app.h_refEdit.Value = 10;
            app.h_groundEdit.Value = 0;
            app.alphaEdit.Value = 0.2;
            
            % Clear results
            app.CalculationData = struct();
            app.OutputTextArea.Value = '';
            app.QueryResultTextArea.Value = 'Enter tether length and click Query to see results...';
            cla(app.PathAxes);
            cla(app.TensionAxes);
            
            % Enable parameter editing
            app.setParametersEditable(true);
            
            % Clear query lines
            for i = 1:length(app.QueryLines)
                if isvalid(app.QueryLines{i})
                    delete(app.QueryLines{i});
                end
            end
            app.QueryLines = {};
            
            % Reset button state
            app.SolveButton.Text = 'Calculate';
            app.SolveButton.BackgroundColor = [1 0.6 0.2]; % Orange
        end
        
        % Set parameter editability
        function setParametersEditable(app, editable)
            state = 'on';
            if ~editable
                state = 'off';
            end
            
            app.D_balloonEdit.Enable = state;
            app.m_envEdit.Enable = state;
            app.m_loadEdit.Enable = state;
            app.d_cableEdit.Enable = state;
            app.rho_cableEdit.Enable = state;
            app.KEdit.Enable = state;
            app.EEdit.Enable = state;
            app.h_step_initEdit.Enable = state;
            app.h_step_factorEdit.Enable = state;
            app.toleranceEdit.Enable = state;
            app.v_refEdit.Enable = state;
            app.h_refEdit.Enable = state;
            app.h_groundEdit.Enable = state;
            app.alphaEdit.Enable = state;
        end
    end
    
    % App initialization
    methods (Access = private)
        function createComponents(app)

            % Get screen size
            screenSize = get(0, 'ScreenSize');
            screenWidth = screenSize(3);
            screenHeight = screenSize(4);

            % Create main window
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [(screenWidth - 1200) / 2, (screenHeight - 800) / 2, 1200, 800];
            app.UIFigure.Name = 'HATASIM Tool';
            app.UIFigure.Resize = 'on';
            
            % Top navigation panel
            app.TopPanel = uipanel(app.UIFigure);
            app.TopPanel.BackgroundColor = [0.25 0.5 0.75];
            app.TopPanel.Position = [1, 761, 1200, 40];
            
            % Navigation buttons
            app.ParamNavButton = uibutton(app.TopPanel, 'push', ...
                'Text', 'Parameters', ...
                'Position', [24, 2, 144, 36], ...
                'FontSize', 18, 'FontWeight', 'bold', ...
                'BackgroundColor', [0.4 0.8 0.4], ... % Green
                'ButtonPushedFcn', createCallbackFcn(app, @ParamNavButtonPushed, true));
            
            app.ResultNavButton = uibutton(app.TopPanel, 'push', ...
                'Text', 'Results', ...
                'Position', [180, 2, 144, 36], ...
                'FontSize', 18, 'FontWeight', 'bold', ...
                'BackgroundColor', [0.9 0.9 0.9], ... % Gray
                'ButtonPushedFcn', createCallbackFcn(app, @ResultNavButtonPushed, true));
            
            % Query section
            app.QueryLabel = uilabel(app.TopPanel, ...
                'Text', 'Original Tether Length (m):', ...
                'Position', [360, 4, 200, 32], ...
                'FontSize', 16, 'FontColor', [1 1 1], ...
                'Visible', 'off');
            
            app.QueryEditField = uieditfield(app.TopPanel, 'numeric', ...
                'Position', [560, 4, 96, 32], ...
                'Limits', [0 10000], ...
                'Value', 0, ...
                'FontSize', 16, ...
                'Visible', 'off');
            
            app.QueryButton = uibutton(app.TopPanel, 'push', ...
                'Text', 'Query', ...
                'Position', [666, 4, 96, 32], ...
                'FontSize', 16, 'FontWeight', 'bold', ...
                'BackgroundColor', [1, 0.65, 0], ... % Orange
                'Visible', 'off', ... 
                'ButtonPushedFcn', createCallbackFcn(app, @QueryButtonPushed, true));
            
            % Bottom panel
            app.BottomPanel = uipanel(app.UIFigure);
            app.BottomPanel.BackgroundColor = [0.25 0.5 0.75]; % Blue
            app.BottomPanel.Position = [1, 1, 1200, 40];
            
            % Footer text
            uilabel(app.BottomPanel, 'Text', 'Supported by An-Group    /    Author: Zhongwei Ni    /    Version: 1.1', ...
                'Position', [200, 2, 800, 36], ...
                'FontSize', 18, 'FontColor', [1 1 1], ...
                'HorizontalAlignment', 'center');
            
            % ========== Parameters Panel ==========
            app.ParamPanel = uipanel(app.UIFigure);
            app.ParamPanel.Position = [1, 41, 1200, 720];
            app.ParamPanel.BackgroundColor = [0.96, 0.98, 1]; % Light blue
            
            % Balloon parameters panel
            app.BalloonParamsPanel = uipanel(app.ParamPanel);
            app.BalloonParamsPanel.Title = 'Balloon Parameters';
            app.BalloonParamsPanel.FontWeight = 'bold';
            app.BalloonParamsPanel.TitleAlignment = 'center';
            app.BalloonParamsPanel.FontSize = 16;
            app.BalloonParamsPanel.BackgroundColor = [0.8 0.9 1]; % Light blue
            app.BalloonParamsPanel.Position = [48, 400, 324, 280];
            
            % Balloon parameter inputs
            uilabel(app.BalloonParamsPanel, 'Text', 'Diameter (m):', 'Position', [40, 190, 132, 28], 'FontSize', 14);
            app.D_balloonEdit = uieditfield(app.BalloonParamsPanel, 'numeric', ...
                'Value', 5, 'Position', [180, 190, 96, 28], 'FontSize', 14);
            
            uilabel(app.BalloonParamsPanel, 'Text', 'Envelope Mass (kg):', 'Position', [40, 130, 132, 28], 'FontSize', 14);
            app.m_envEdit = uieditfield(app.BalloonParamsPanel, 'numeric', ...
                'Value', 5, 'Position', [180, 130, 96, 28], 'FontSize', 14);
            
            uilabel(app.BalloonParamsPanel, 'Text', 'Payload Mass (kg):', 'Position', [40, 70, 132, 28], 'FontSize', 14);
            app.m_loadEdit = uieditfield(app.BalloonParamsPanel, 'numeric', ...
                'Value', 5, 'Position', [180, 70, 96, 28], 'FontSize', 14);
            
            % Tether parameters panel
            app.CableParamsPanel = uipanel(app.ParamPanel);
            app.CableParamsPanel.Title = 'Tether Parameters';
            app.CableParamsPanel.FontWeight = 'bold';
            app.CableParamsPanel.TitleAlignment = 'center';
            app.CableParamsPanel.FontSize = 16;
            app.CableParamsPanel.BackgroundColor = [0.8 0.9 1];
            app.CableParamsPanel.Position = [48, 40, 324, 320];
            
            % Tether parameter inputs
            uilabel(app.CableParamsPanel, 'Text', 'Diameter (m):', 'Position', [36, 240, 132, 28], 'FontSize', 14);
            app.d_cableEdit = uieditfield(app.CableParamsPanel, 'numeric', ...
                'Value', 0.006, 'Position', [180, 240, 96, 28], 'FontSize', 14);
            
            uilabel(app.CableParamsPanel, 'Text', 'Density (kg/m³):', 'Position', [36, 180, 132, 28], 'FontSize', 14);
            app.rho_cableEdit = uieditfield(app.CableParamsPanel, 'numeric', ...
                'Value', 970, 'Position', [180, 180, 96, 28], 'FontSize', 14);
            
            uilabel(app.CableParamsPanel, 'Text', 'Segment Length (m):', 'Position', [36, 120, 132, 28], 'FontSize', 14);
            app.KEdit = uieditfield(app.CableParamsPanel, 'numeric', ...
                'Value', 0.1, 'Position', [180, 120, 96, 28], 'FontSize', 14);
            
            uilabel(app.CableParamsPanel, 'Text', 'Elastic Modulus (Pa):', 'Position', [36, 60, 140, 28], 'FontSize', 14);
            app.EEdit = uieditfield(app.CableParamsPanel, 'numeric', ...
                'Value', 95e9, 'Position', [180, 60, 96, 28], 'FontSize', 14);
            
            % Precision parameters panel
            app.PrecisionParamsPanel = uipanel(app.ParamPanel);
            app.PrecisionParamsPanel.Title = 'Precision Parameters';
            app.PrecisionParamsPanel.FontWeight = 'bold';
            app.PrecisionParamsPanel.TitleAlignment = 'center';
            app.PrecisionParamsPanel.FontSize = 16;
            app.PrecisionParamsPanel.BackgroundColor = [0.8 0.9 1];
            app.PrecisionParamsPanel.Position = [412, 400, 324, 280];
            
            % Precision parameter inputs
            uilabel(app.PrecisionParamsPanel, 'Text', 'Initial Step (m):', 'Position', [40, 190, 120, 28], 'FontSize', 14);
            app.h_step_initEdit = uieditfield(app.PrecisionParamsPanel, 'numeric', ...
                'Value', 1, 'Position', [176, 190, 100, 28], 'FontSize', 14);
            
            uilabel(app.PrecisionParamsPanel, 'Text', 'Step Factor:', 'Position', [40, 130, 120, 28], 'FontSize', 14);
            app.h_step_factorEdit = uieditfield(app.PrecisionParamsPanel, 'numeric', ...
                'Value', 0.01, 'Position', [176, 130, 100, 28], 'FontSize', 14);
            
            uilabel(app.PrecisionParamsPanel, 'Text', 'Tolerance (m):', 'Position', [40, 70, 120, 28], 'FontSize', 14);
            app.toleranceEdit = uieditfield(app.PrecisionParamsPanel, 'numeric', ...
                'Value', 0.001, 'Position', [176, 70, 100, 28], 'FontSize', 14);
            
            % Wind parameters panel
            app.WindParamsPanel = uipanel(app.ParamPanel);
            app.WindParamsPanel.Title = 'Wind Parameters';
            app.WindParamsPanel.FontWeight = 'bold';
            app.WindParamsPanel.TitleAlignment = 'center';
            app.WindParamsPanel.FontSize = 16;
            app.WindParamsPanel.BackgroundColor = [0.8 0.9 1];
            app.WindParamsPanel.Position = [412, 40, 324, 320];
            
            % Wind parameter inputs
            uilabel(app.WindParamsPanel, 'Text', 'Ref. Speed (m/s):', 'Position', [40, 240, 132, 28], 'FontSize', 14);
            app.v_refEdit = uieditfield(app.WindParamsPanel, 'numeric', ...
                'Value', 3, 'Position', [176, 240, 100, 28], 'FontSize', 14);
            
            uilabel(app.WindParamsPanel, 'Text', 'Ref. Height (m):', 'Position', [40, 180, 132, 28], 'FontSize', 14);
            app.h_refEdit = uieditfield(app.WindParamsPanel, 'numeric', ...
                'Value', 10, 'Position', [176, 180, 100, 28], 'FontSize', 14);
            
            uilabel(app.WindParamsPanel, 'Text', 'Ground Altitude (m):', 'Position', [40, 120, 132, 28], 'FontSize', 14);
            app.h_groundEdit = uieditfield(app.WindParamsPanel, 'numeric', ...
                'Value', 0, 'Position', [176, 120, 100, 28], 'FontSize', 14);
            
            uilabel(app.WindParamsPanel, 'Text', 'Shear Exponent:', 'Position', [40, 60, 132, 28], 'FontSize', 14);
            app.alphaEdit = uieditfield(app.WindParamsPanel, 'numeric', ...
                'Value', 0.2, 'Position', [176, 60, 100, 28], 'FontSize', 14);
            
            % Diagram area
            app.DiagramAxes = uiaxes(app.ParamPanel);
            title(app.DiagramAxes, 'System Diagram', 'FontSize', 18, 'FontWeight', 'bold');
            app.DiagramAxes.Position = [784, 180, 360, 500];
            
            % Draw diagram
            plot(app.DiagramAxes, [-10, 8, 23, 35, 44, 50], [0, 9, 22.5, 40.5, 63, 90], 'k-', 'LineWidth', 2); % Tether
            hold(app.DiagramAxes, 'on');
            rectangle(app.DiagramAxes, 'Position', [43.75 83.75 12.5 12.5], 'Curvature', [1 1], 'FaceColor', [0.8 0.9 1]); % Balloon
            plot(app.DiagramAxes, [-10 60], [0 0], 'k-', 'LineWidth', 4); % Ground
            text(app.DiagramAxes, 30, 90, 'Balloon', 'FontSize', 14);
            text(app.DiagramAxes, 8, 22.5, 'Tether', 'FontSize', 14);
            text(app.DiagramAxes, 44.75, 5, 'Ground', 'FontSize', 14);
            text(app.DiagramAxes, 21, 60, 'V_{ref}', 'FontSize', 14);
            quiver(app.DiagramAxes, 20, 55, 10, 0, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 1);
            hold(app.DiagramAxes, 'off');
            axis(app.DiagramAxes, 'equal');
            app.DiagramAxes.XLim = [-20 70];
            app.DiagramAxes.YLim = [-10 105];
            app.DiagramAxes.XTick = [];
            app.DiagramAxes.YTick = [];
            app.DiagramAxes.Box = 'on';
            
            % Solve button
            app.SolveButton = uibutton(app.ParamPanel, 'push', ...
                'Text', 'Calculate', ...
                'Position', [876, 80, 180, 64], ...
                'FontSize', 28, 'FontWeight', 'bold', ...
                'BackgroundColor', [0.6 0.8 1], ...  % Blue
                'ButtonPushedFcn', createCallbackFcn(app, @SolveButtonPushed, true));
            
            % ========== Results Panel ==========
            app.ResultPanel = uipanel(app.UIFigure);
            app.ResultPanel.Position = [1, 41, 1200, 720];
            app.ResultPanel.BackgroundColor = [0.96, 0.98, 1]; % Light blue
            app.ResultPanel.Visible = 'off';
            
            % Path plot
            app.PathAxes = uiaxes(app.ResultPanel);
            app.PathAxes.Position = [48, 300, 528, 400];
            title(app.PathAxes, 'Ascent Path & Tether Shape', 'FontSize', 14, 'FontWeight', 'bold');
            
            % Tension plot
            app.TensionAxes = uiaxes(app.ResultPanel);
            app.TensionAxes.Position = [624, 300, 528, 400];
            title(app.TensionAxes, 'Tension Distribution', 'FontSize', 14, 'FontWeight', 'bold');
                        
            % Output panel
            app.OutputPanel = uipanel(app.ResultPanel, 'Title', 'Final Results', 'FontSize', 14, ...
                'FontWeight', 'bold', 'TitleAlignment', 'center');
            app.OutputPanel.Position = [80, 30, 480, 240];
            app.OutputPanel.BackgroundColor = [0.5 0.75 1]; % Light blue
            
            app.OutputTextArea = uitextarea(app.OutputPanel);
            app.OutputTextArea.Position = [0, 0, 480, 218];
            app.OutputTextArea.FontName = 'Consolas';
            app.OutputTextArea.FontSize = 13;
            app.OutputTextArea.BackgroundColor = [1 1 1]; % White
            app.OutputTextArea.Value = 'Results will appear here...';
            
            % Query result panel
            app.QueryResultPanel = uipanel(app.ResultPanel, 'Title', 'Query Results', 'FontSize', 14, ...
                'FontWeight', 'bold', 'TitleAlignment', 'center');
            app.QueryResultPanel.Position = [640, 30, 480, 240];
            app.QueryResultPanel.BackgroundColor = [0.5 0.75 1]; % Light blue
            
            app.QueryResultTextArea = uitextarea(app.QueryResultPanel);
            app.QueryResultTextArea.Position = [0, 0, 480, 218];
            app.QueryResultTextArea.FontName = 'Consolas';
            app.QueryResultTextArea.FontSize = 13;
            app.QueryResultTextArea.BackgroundColor = [1 1 1]; % White
            app.QueryResultTextArea.Value = 'Enter tether length and click Query to see results...';
            
            % Show main window
            app.UIFigure.Visible = 'on';
        end
    end
    
    % Callback functions
    methods (Access = private)
        function ParamNavButtonPushed(app, ~)
            app.switchToParamPanel();
        end
        
        function ResultNavButtonPushed(app, ~)
            app.switchToResultPanel();
        end
        
        function SolveButtonPushed(app, ~)
            if strcmp(app.SolveButton.Text, 'Calculate')
                % Disable parameter editing
                app.setParametersEditable(false);
                
                % Start calculation
                app.calculate();
                
                % Change button to reset
                app.SolveButton.Text = 'RESET';
                app.SolveButton.BackgroundColor = [1 0.65 0]; % Orange
            else
                % Reset parameters
                app.resetParameters();
            end
        end
        
        function QueryButtonPushed(app, ~)
            if isempty(app.CalculationData) || isempty(app.CalculationData.all_rope_data)
                uialert(app.UIFigure, 'Please complete calculation first', 'Error');
                return;
            end

            % Use max wind-affected tether length as upper limit
            max_rope_length = app.MaxLWinded;
            L_input = app.QueryEditField.Value;
            if L_input < 0 || L_input > max_rope_length
                uialert(app.UIFigure, sprintf('Length out of range (0-%.2f m)', max_rope_length), 'Error');
                return;
            end

            % Find closest calculated length
            [~, idx] = min(abs(app.CalculationData.L_original_results - L_input));
            closest_L = app.CalculationData.L_original_results(idx);
            closest_data = app.CalculationData.all_rope_data(idx+1);

            % Calculate stress
            cable_area = pi/4 * app.d_cableEdit.Value^2;
            stress_MPa = closest_data.F_segments(1) / cable_area / 1e6;

            % Display results
            resultStr = sprintf('\n【Query Length %.2f m】\n', L_input);
            resultStr = [resultStr sprintf('·Matched Length: %.2f m\n', closest_L)];
            resultStr = [resultStr sprintf('·Altitude: %.2f m (MSL) | %.2f m AGL\n', ...
                closest_data.H, closest_data.H - app.h_groundEdit.Value)];
            resultStr = [resultStr sprintf('·Drift Distance: %.2f m\n', closest_data.X_total)];
            resultStr = [resultStr sprintf('·Tension at the balloon–tether interface: %.2f N (%.2f MPa)\n', ...
                closest_data.F_segments(1), stress_MPa)];

            app.QueryResultTextArea.Value = resultStr;

            % Clear previous query lines
            for i = 1:length(app.QueryLines)
                if isvalid(app.QueryLines{i})
                    delete(app.QueryLines{i});
                end
            end
            app.QueryLines = {};

            % Visualize shape comparison
            hold(app.PathAxes, 'on');
            corrected_X = closest_data.X_total - closest_data.rope_X;
            p1 = plot(app.PathAxes, corrected_X, closest_data.rope_Z - app.h_groundEdit.Value, 'b-', 'LineWidth', 2.5);

            % Add tension distribution for query length
            hold(app.TensionAxes, 'on');
            K = app.KEdit.Value;
            original_length_positions = (0:K:(length(closest_data.F_segments)-1)*K);
            p2 = plot(app.TensionAxes, original_length_positions, flip(closest_data.F_segments), 'b-', 'LineWidth', 2.5);
            hold(app.TensionAxes, 'off');

            % Store line handles
            app.QueryLines{end+1} = p1;
            app.QueryLines{end+1} = p2;

            % Update legends
            legend(app.PathAxes, 'off');
            legend(app.PathAxes, {'Path', 'Final Shape', 'Query Shape'}, 'Location', 'best');

            legend(app.TensionAxes, 'off');
            legend(app.TensionAxes, {'Tension at the balloon–tether interface', 'Final Distribution', 'Query Distribution'}, 'Location', 'best');
        end
    end
    
    % App entry point
    methods (Access = public)
        function app = HATASIM_Tool_EN
            createComponents(app)
            registerApp(app, app.UIFigure);
            
            % Initialize default view
            app.switchToParamPanel();
            
            if nargout == 0
                clear app
            end
        end
    end
end