function data_out = level_of_analysis(data_in, metadata_in, level, conditions, statistic)

if nargin < 5
    statistic = @nanmean;
end

sessions = unique(metadata_in.session)';
n_time   = size(data_in, 2);
n_sess   = numel(sessions);

switch level

    case "session"
        data_out = nan(max(sessions), n_time);
        for ii = sessions
            data_out(ii,:) = statistic(data_in(conditions & metadata_in.session == ii, :));
        end

    case "session2d"
        data_out = nan(15, n_time, max(sessions));
        for ii = sessions
            for jj = 1:15
                data_out(jj,:,ii) = statistic(data_in(conditions & ...
                    metadata_in.session == ii & metadata_in.channel == jj, :));
            end
        end

    case "cat_trials"
        n_trials_c = sum(conditions & metadata_in.channel == 1);
        data_out = nan(15, n_time, n_trials_c);
        for jj = 1:15
            data_out(jj,:,:) = data_in(conditions & metadata_in.channel == jj, :)';
        end

    case {"upper", "middle", "lower"}
        switch level
            case "upper";  ch_range = 1:5;
            case "middle"; ch_range = 6:10;
            case "lower";  ch_range = 11:15;
        end
        n_chans  = numel(ch_range);
        data_out = nan(n_sess * n_chans, n_time);
        row = 1;
        for ii = sessions
            for jj = ch_range
                mask = conditions & metadata_in.session == ii & metadata_in.channel == jj;
                if sum(mask) > 2
                    data_out(row,:) = statistic(data_in(mask, :));
                end
                row = row + 1;
            end
        end

    case "channel"
        data_out = nan(n_sess * 15, n_time);
        row = 1;
        for ii = sessions
            for jj = 1:15
                data_out(row,:) = statistic(data_in(conditions & ...
                    metadata_in.session == ii & metadata_in.channel == jj, :));
                row = row + 1;
            end
        end

    case "channel_2"
        channels = unique(metadata_in.channel)';
        data_out = nan(numel(channels), n_time);
        for k = 1:numel(channels)
            data_temp = data_in(conditions & metadata_in.channel == channels(k), :);
            if any(size(data_temp) == 1)
                data_out(k,:) = data_temp;
            elseif size(data_temp, 1) > 1
                data_out(k,:) = statistic(data_temp);
            end
        end
end
