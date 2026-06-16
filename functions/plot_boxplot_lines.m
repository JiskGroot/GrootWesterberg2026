function plot_boxplot_lines(data_up, data_p, compartment_name, y_lim)
    % Color selection based on compartment name
    if strcmp(compartment_name, "upper")
        color_to_use = [51 0 255]/255;
    elseif strcmp(compartment_name, "middle")
        color_to_use = [255 128 0]/255;
    elseif strcmp(compartment_name, "lower")
        color_to_use = [60 255 0]/255;
    else
        color_to_use = [0 0 0];
    end
    
    positions = [1, 2];
    data_all = [data_up, data_p];
    
    % Boxplot in current axes
    boxplot(data_all, 'Positions', positions, 'Widths', 0.5);
    hold on; % Allow overlay of points and lines
    
    % Set box colors
    boxes = findobj(gca, 'Tag', 'Box');
    set(boxes, 'Color', color_to_use);
    
    % Plot individual points and connecting lines
    for i = 1:length(data_up)
        plot(positions(1), data_up(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
        plot(positions(2), data_p(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
        plot(positions, [data_up(i), data_p(i)], 'Color', [0.7 0.7 0.7], 'LineStyle', '-', 'LineWidth', 1);
    end
    
    % Axes formatting
    set(gca, 'XTick', positions, 'XTickLabel', {'Unpredictable', 'Predictable'});
    ylabel('Value');
    ylim([0, y_lim]);
    set(gca, "TickDir", "out")
    % title(['Laminar Compartment: ', compartment_name]);
    
    % Statistical test: Wilcoxon signed-rank test
    [p,~,stats] = signrank(data_up, data_p);
    
    % Display p-value and significance
    if p < 0.05
        median_diff = median(data_p - data_up);
        if median_diff > 0
            sig_msg = 'Significant increase';
        elseif median_diff < 0
            sig_msg = 'Significant decrease';
        else
            sig_msg = 'No median change';
        end
        text(1.5, y_lim*0.85, sprintf('p = %.3f\nzval = %.2f\n%s', p, stats.zval, sig_msg), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    else
        text(1.5, y_lim*0.95, sprintf('p = %.3f', p), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
end
