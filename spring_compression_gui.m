function spring_compression_gui
    % initialization
    a_1 = 50; % mm
    b_1 = 80; % mm
    x0_1 = 20;
    k_1 = (250/1E3*9.8)/(80); % N/mm
    num_springs = 2; % default - 2 spring

    fig = figure('Position', [100, 100, 800, 600]);

    uicontrol('Style', 'text', 'Position', [20, 550, 100, 20], 'String', 'Number of Springs:', 'FontSize', 12);
    num_springs_popup = uicontrol('Style', 'popupmenu', 'Position', [130, 550, 100, 20], 'String', {'2', '3'}, 'Callback', @update_num_springs);

    create_input_panel(20, 500, 'x0_1', x0_1, @update_x0_1);
    create_input_panel(20, 450, 'k_1', k_1, @update_k_1);

    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);

    update_plot();

    function create_input_panel(x, y, label, init_value, callback)
        uicontrol('Style', 'text', 'Position', [x, y, 100, 20], 'String', label, 'FontSize', 12);
        uicontrol('Style', 'edit', 'Position', [x+100, y, 100, 20], 'String', num2str(init_value), 'Callback', callback);
        uicontrol('Style', 'slider', 'Position', [x+200, y, 200, 20], 'Min', 0, 'Max', 100, 'Value', init_value, 'Callback', callback);
    end

    function update_num_springs(src, ~)
        num_springs = str2double(src.String{src.Value});
        switch num_springs
            case 2
                delete(findobj('Tag', 'x0_2'));
                delete(findobj('Tag', 'k_2'));
                create_input_panel(20, 400, 'x0_2', x0_1, @update_x0_2);
                create_input_panel(20, 350, 'k_2', k_1, @update_k_2);
            case 3
                delete(findobj('Tag', 'x0_2'));
                delete(findobj('Tag', 'k_2'));
                delete(findobj('Tag', 'x0_3'));
                delete(findobj('Tag', 'k_3'));
                create_input_panel(20, 400, 'x0_2', x0_1, @update_x0_2);
                create_input_panel(20, 350, 'k_2', k_1, @update_k_2);
                create_input_panel(20, 300, 'x0_3', x0_1, @update_x0_3);
                create_input_panel(20, 250, 'k_3', k_1, @update_k_3);
        end
        update_plot();
    end

    function update_x0_1(src, ~)
        x0_1 = str2double(src.String);
        update_plot();
    end

    function update_k_1(src, ~)
        k_1 = str2double(src.String);
        update_plot();
    end

    function update_x0_2(src, ~)
        x0_2 = str2double(src.String);
        update_plot();
    end

    function update_k_2(src, ~)
        k_2 = str2double(src.String);
        update_plot();
    end

    function update_x0_3(src, ~)
        x0_3 = str2double(src.String);
        update_plot();
    end

    function update_k_3(src, ~)
        k_3 = str2double(src.String);
        update_plot();
    end

    function update_plot()
        x1_1 = [0:1:60]; % Displacement
        F_y_1 = zeros(size(x1_1));

        switch num_springs
            case 2
                x0 = [x0_1, x0_2];
                k = [k_1, k_2];
            case 3
                x0 = [x0_1, x0_2, x0_3];
                k = [k_1, k_2, k_3];
        end

        for i = 1:num_springs
            EF_0_i = sqrt(3)*sqrt(b_1^2 - (b_1-x0(i)/2)^2 )*ones(size(x1_1));
            EF_i_i = sqrt(3)*sqrt(b_1^2 - ((2*b_1-x0(i)-x1_1)/2).^2  );

            delta_i = EF_i_i - EF_0_i;
            F_EF_i = k(i) .* delta_i;
            F_EB_i = F_EF_i .* cos(pi/6);

            psi_i = asin((2*b_1-x0(i)-x1_1)/(2*b_1));
            F_EB1_i = F_EB_i ./ (2*cos(psi_i));
            F_B1_i = F_EB1_i .* sin(psi_i);

            F_y_1 = F_y_1 + 6*F_B1_i; % 3 springs, 6 supports
        end

        subplot(ax1);
        cla;
        plot(x1_1, F_y_1, 'LineWidth', 1.5);
        grid on;
        box on;
        xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
        legend(arrayfun(@(i) ['Spring ', num2str(i)], 1:num_springs, 'UniformOutput', false));

%         subplot(ax2);
%         cla;
%         plot(x1_1, cumsum(F_y_1), 'LineWidth', 1.5);
%         grid on;
%         box on;
%         xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
%         ylabel('Cumulative Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
    end
end