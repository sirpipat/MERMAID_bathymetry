function stats = pp2026_figure2(obs_struct, obsmasterdir, synmasterdir, ...
    presmasterdir, handpick, bath, stats)
% stats = PP2026_FIGURE2(obs_struct, obsmasterdir, synmasterdir, ...
%     presmasterdir, handpick, handpick, bath, stats_input)
%
% Makes figure 2A and 2B for Pipatprathanporn+2026 paper.
%
% INPUT:
% obs_struct        a struct containing
% - fcorners            [lower upper] corner frequencies
% - CCmaxs              maximum correlation coefficients for
%                       [flat bath] cases
% - t_shifts            optimal time shifts for [flat bath] cases
% - metadata            SAC Headers associated to the obsfile
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% presmasterdir     the master directory for the synthetic pressure
%                   seismograms will be created and sorted into IRIS eent
%                   ID folders
% handpick          whether to activate the handpick mode for optimal
%                   timeshift. [Default: false]
%     false -- picks the one that give the maximum correlation coefficient
%     true  -- you are prompts to accept the automaic pick [Y/N]
%              if you reject (N), you will be given options to choose from
%              a few candidates. You may accept (Y), move to the next (N)
%              or previous (P) to find a better fit. The seismograms will
%              be time shifted to give you a visual guide. You may quit (Q)
%              the plotting in case you have made a mistake and want to
%              start over. You may activate the keyword statement (K) to
%              debug.
% bath              source for bathymetry plots
%     'GEBCO' [default] -- interpolated bathymetry grid from GEBCO
%     'SPECFEM'         -- tapered bathymetry grid used for SPECFEM3D run
% stats_input       Input correlation travel-time measurement see STATS
%
% OUTPUT:
% stats             struct with a following fields
% - is_low              whether this entry is a LOWCC run
% - id                  which entry from obs_struct (input) to refer to
% - mean                "mean" elevation of [1D 2D] bathymetry 
%                       If bath == 'SPECFEM', the second column is the edge
%                       elevation of the interface, used for FK injection.
% - std                 standard deviation elevation of [1D 2D] bathymetry
% - t_shifts            time shifts applied to the [2D 3D] synthetics. Note
%                       that 2D synthetics use 1D bathymetry, and 3D
%                       synthetics use 2D bathymetry.
% - ccmaxs              A 5-column matrix with a following columns
%                           1. 2D to observed over [-5 5] window
%                           2. 3D to observed over [-5 5] window
%                           3. 2D to 3D over [-10 30] window
%                           4. 2D to observed over [-10 30] window
%                           5. 3D to observed over [-10 30] window
%
% Last modified by sirawich-at-princeton.edu, 07/16/2025

defval('handpick', false)
defval('bath', 'GEBCO')
defval('stats', [])

keywords = {'ascend', 'descend'};
keyfolders = {'LOWCC', 'HIGHCC'};
letterlabel = {'a', 'b', 'c', 'd', 'e'};

% record pick stats
if isempty(stats)
    have_stats_input = false;

    stats.is_low = nan(60, 1);
    stats.id = nan(60, 1);
    stats.mean = nan(60, 2);
    stats.std = nan(60, 2);
    stats.t_shifts = nan(60, 2);
    stats.ccmaxs = nan(60, 5);
else
    have_stats_input = true;
end

for pp = 1:1
    ii_cases_low_cc = [1 2 13 20 28] + 0 *5 * (pp - 1);%[1 2 14 9 10];
    ii_cases_high_cc = [1 2 12 5 29] + 0 * 5 * (pp - 1);%[1 2 8 11 14];
    ii_cases_master = [ii_cases_low_cc; ii_cases_high_cc];
    
    min_snr = 15;
    min_gcarc = 20;
    min_depth = 33;
    wh_valid = and(and(obs_struct.snr(:,2) > min_snr, ...
        obs_struct.metadata.GCARC > min_gcarc), ...
        obs_struct.metadata.EVDP > min_depth);
    ic_list = indeks((1:length(obs_struct.CCmaxs(:,2)))', wh_valid);
    for kk = 1:2
        [~, ic] = sort(obs_struct.CCmaxs(wh_valid,2), keywords{kk});
        ic = ic_list(ic);
    
        ii_cases = ii_cases_master(kk,:);
        
        figure(kk);
        set(gcf, 'Units', 'inches', 'Position', [0 1 12 7])
        clf
        
        xposshift = [0 0.2 0.4 0.6 0.8];
        
        for ii = 1:5
            ii_stats = (kk-1) * 30 + ii_cases(ii);

            %% Top row: bird-eye-view bathymetry plot
            [ll, tt, zz, llons, llats] = ...
                bathymetryprofile2d([20000 20000], [81 81], ...
                [obs_struct.metadata.STLO(ic(ii_cases(ii))) ...
                obs_struct.metadata.STLA(ic(ii_cases(ii)))], ...
                obs_struct.metadata.BAZ(ic(ii_cases(ii)))+180);

            % This is the mean elevation used for SPECFEM2D runs in the
            % past
            stats.mean(ii_stats, 1) = mean(zz(41, :));
            stats.std(ii_stats, 1) = std(zz(41, :));

            if ~strcmpi(bath, 'GEBCO')
                ddir = fullfile(getenv('REMOTE3D'), '20250717_MERMAID_INSTASEIS', ...
                    sprintf('LAYERED_OC_MERMAID_%s_%02d', keyfolders{kk}, ...
                    ii_cases(ii)));
                itfs = loadinterfacefiles3d(fullfile(ddir, 'DATA', ...
                    'meshfem3D_files', 'interfaces.dat'));
                zz = itfs{1}.Z';

                % "mean" bathymetry for plotting is the edge elevation
                mean_zz = zz(1,1);
                std_zz = std(zz(:));
            else
                mean_zz = mean(zz(:));
                std_zz = std(zz(:));
            end
        
            subplot('Position', [0.030+xposshift(ii) 0.65 0.16 0.2562])
            imagesc([-10 10], [-10 10], rot90(zz'));
            colormap(kelicol);
            clim(mean_zz + std_zz * [-3 3])
            % cb = colorbar('Location', 'southoutside');
            % cb.Label.String = 'elevation (m)';
            % cb.TickDirection = 'out';
            grid on;
            box on;
            hold on;
            plot([-10 10], [0 0], 'LineWidth', 1, 'Color', 'k')
            axis image;
            axis xy;
            xticks(-10:5:10)
            yticks(-10:5:10)
            if ii > 1
                nolabels(gca, 3);  % only keep y-axis label for leftmost plot
            else
                nolabels(gca, 1);
                yticklabels(-10:5:10)
                ylabel('transverse position (km)')
            end
            set(gca, 'TickDir', 'out', 'FontSize', 9, ...
                'Position', [0.033+xposshift(ii) 0.68 0.16 0.2562])
            title(sprintf('%d-%s', obs_struct.metadata.USER7(ic(ii_cases(ii))), ...
                obs_struct.metadata.KSTNM{ic(ii_cases(ii))}), ...
                'FontWeight', 'normal')
            subtitle(' ', 'FontSize', 2)
        
            % add mean and std info the corners
            ax1_mean = axes('Position', [0.038+xposshift(ii) 0.68 0.08 0.025]);
            axeslabel(ax1_mean, 0.5, 0.5, sprintf('med  =  %.3f km', ...
                round(mean_zz)/1000), 'FontSize', 8, ...
                'HorizontalAlignment', 'center');
            nolabels(ax1_mean, 3)
            set(ax1_mean, 'Box', 'on', 'XTick', [], 'YTick', [])
    
            ax1_std = axes('Position', [0.128+xposshift(ii) 0.68 0.06 0.025]);
            axeslabel(ax1_std, 0.5, 0.5, sprintf('std  =  %d m', ...
                round(std_zz)), 'FontSize', 8, ...
                'HorizontalAlignment', 'center');
            nolabels(ax1_std, 3)
            set(ax1_std, 'Box', 'on', 'XTick', [], 'YTick', [])
    
            % boxed label
            ax1b = axes('Position', [0.005+xposshift(ii) 0.94 0.014*1.2 0.024*1.2]);
            axeslabel(ax1b, 0.5, 0.5, letterlabel{ii}, 'FontSize', 10, ...
                'HorizontalAlignment', 'center');
            set(ax1b, 'Box', 'on', 'XTick', [], 'YTick', [])
            nolabels(ax1b, 3)
    
            %% Middle row: bathymetry profile along great-circle path
            ah(ii) = subplot('Position', [0.038+xposshift(ii) 0.48 0.15 0.18]);
            hold on
            for jj = 1:size(zz,1)
                hold on
                if jj ~= 41
                    % plot(ll(1,:)/1000, zz(jj,:)/1000, 'LineWidth', 0.25, ...
                    %     'Color', [0.7 0.7 0.7]);
                    plot(ll(1,:)/1000, zz(jj,:)/1000, 'LineWidth', 0.25, ...
                        'Color', [0.95-0.4/80*(jj-1) 0.8 0.55+0.4/80*(jj-1)]);
                else
                    plot(ll(1,:)/1000, zz(jj,:)/1000, 'LineWidth', 1, ...
                        'Color', 'k');
                end
            end
            %plot(ll(1,:)/1000, zz(41,:)/1000, 'LineWidth', 1, 'Color', 'k');
            grid on
            box on
            xticks(-10:5:10)
            xticklabels(-10:5:10)
            xlabel('radial position (km)')
            ylabel('elevation (km)')
            ylim(mean_zz/1000 + [-3 3] * std_zz / 1000)
            ylim(mean_zz/1000 + 1.1*[-1 1] * max(abs(zz(:)-mean_zz))/ 1000);
            getzm(ii)=mean_zz/1000;
            getzs(ii)=max(abs(zz(:)-mean_zz))/1000;
           
            
            %% Bottom row: seismograms
            ax3 = subplot('Position', [0.038+xposshift(ii) 0.07 0.15 0.26]);
        
            % identify observed pressure file
            eventid = obs_struct.metadata.USER7(ic(ii_cases(ii)));
            stationid = indeks(obs_struct.metadata.KSTNM{ic(ii_cases(ii))}, '2:5');
            try
                obsfile = cindeks(ls2cell(sprintf('%s%d/*.%s_*.sac', ...
                    obsmasterdir, eventid, stationid), 1), 1);
            catch ME
                if strcmp(ME.message, 'This directory or file does not exist')
                    obsfile = cindeks(ls2cell(sprintf('%s%d/*.%s_*.sac', ...
                        obsmasterdir, eventid, stationid(end-1:end)), 1), 1);
                end
            end
        
            % read observed pressure file
            [seis_o, hdr_o] = readsac(obsfile);
            [dt_ref_o, ~, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);
            t_relative = seconds(dts_o - dt_ref_o) - hdr_o.T0;
            
            % remove instrument response
            seis_o = counts2pa(seis_o, fs_o, [0.01 0.02 5 10], [], [], false);
    
            % corner frequency
            fc = obs_struct.fcorners(ic(ii_cases(ii)), :);
    
            % badpass filtering
            seis_o2 = bandpass(detrend(seis_o), fs_o, fc(1), fc(2), 4, 2, ...
                'butter', 'linear');
    
            % plot observed seismogram
            wh = and(t_relative >= -10, t_relative <= 30);
            lo = plot(t_relative, seis_o2 / max(abs(seis_o2(wh))), 'LineWidth', 1, 'Color', 'k');
            xlim([-10 30])
            grid on
            hold on
            
            % identify synthetic pressure file
            stationid = indeks(obs_struct.metadata.KSTNM{ic(ii_cases(ii))}, '2:5');
            presfile = cindeks(ls2cell(sprintf('%s/%d/*_%s_0_*.sac', ...
                presmasterdir, eventid, stationid), 1), 1);
            
            % read synthetic pressure file
            [seis_p, hdr_p] = readsac(presfile);
            [dt_ref_p, ~, ~, fs_p, ~, dts_p] = gethdrinfo(hdr_p);
    
            % bandpass filtering
            seis_p2 = bandpass(detrend(seis_p), fs_p, fc(1), fc(2), 4, 2, ...
                'butter', 'linear');
    
            % normalize and plot
            if have_stats_input
                t_shift = stats.t_shifts(ii_stats, 1);
            else
                t_shift = obs_struct.t_shifts(ic(ii_cases(ii)), 2);
            end
            wh_s = and(t_relative + t_shift >= -10, t_relative + t_shift <= 30);
            l2d = plot(t_relative + t_shift, seis_p2 / max(abs(seis_p2(wh_s))) + 2.5, 'LineWidth', 0.75, 'Color', 'r');
            
            % window for cross correlation
            wh3 = and(t_relative >= -5, t_relative <= 5);
            seis_o3 = seis_o2(wh3);
            t_relative_o3 = t_relative(wh3);
            dt_now = datetime('now');
    
            % ask the user if they like the automatic correlation
            % Y - Yes, N - No, Q - Quit
            if handpick
                prompt = "Do you like this? Y/N/Q [Y]: ";
                txt = upper(input(prompt,"s"));
                if isempty(txt)
                    txt = 'N';
                end
                % Q means Quit
                if strcmp(txt, 'Q')
                    return
                elseif strcmp(txt, 'K')
                    keyboard
                end
                
            else
                txt = 'Y';
            end
            if ~strcmp(txt, 'Y')
                [~, ~, lags0, cc0, ~, s0] = ccscale(seis_o3, seis_p2, dt_now + seconds(t_relative_o3(1)), dt_now + seconds(t_relative(1)), fs_o, seconds(20), 'hard', false, false);
                cc1 = cc0;
                cc1(cc1 < 0.25 * max(cc1)) = 0;
                [pks, locs] = findpeaks(cc1);
                ii_locs = ceil(length(locs) / 2);
                while ~strcmp(txt, 'Y') && ii_locs < length(locs)
                    t_shift = lags0(locs(ii_locs));
                    CCmax0 = pks(ii_locs);
                    Smax0 = s0(locs(ii_locs));
        
                    % redraw the figure
                    delete(l2d);
                    wh_s = and(t_relative + t_shift >= -10, t_relative + t_shift <= 30);
                    l2d = plot(t_relative + t_shift, seis_p2 / max(abs(seis_p2(wh_s))) + 2.5, 'LineWidth', 0.75, 'Color', 'r');
    
                    % ask if they like this
                    prompt = sprintf("Do you like this? (%d/%d) Y/N/P/Q [Y]: ", ii_locs, length(locs));
                    txt = upper(input(prompt,"s"));
                    if isempty(txt)
                        txt = 'N';
                    end
    
                    if strcmp(txt, 'Q')
                        return
                    elseif strcmp(txt, 'K')
                        keyboard
                    end
        
                    % move to the next candidate
                    if strcmp(txt, 'P')
                        ii_locs = max(ii_locs - 1, 1);
                    elseif strcmp(txt, 'Y')
                    else
                        ii_locs = ii_locs + 1;
                    end
                end
            else
                if have_stats_input
                    CCmax0 = stats.ccmaxs(ii_stats, 1);
                    t_shift = stats.t_shifts(ii_stats, 1);
                else
                    CCmax0 = obs_struct.CCmaxs(ic(ii_cases(ii)), 2);
                    t_shift = obs_struct.t_shifts(ic(ii_cases(ii)), 2);
                end
            end    
    
            set(gca, 'FontSize', 9, 'TickDir', 'out')
            ax3.Children = ax3.Children([2 1]);
    
            % read Instaseis z-displacement at the ocean bottom
            synfile = cindeks(ls2cell(sprintf('%s%d/*_%s_*.sac', synmasterdir, eventid, stationid), 1), 1);
            [seis_s, hdr_s, ~, ~, tims_s] = readsac(synfile);
            [~, ~, ~, fs_s] = gethdrinfo(hdr_s);
    
             % apply bandpass filter
            seis_sf = bandpass(detrend(seis_s), fs_s, fc(1), fc(2), 4, ...
                2, 'butter', 'linear');
    
            ddir = fullfile(getenv('REMOTE3D'), '20250717_MERMAID_INSTASEIS', ...
                sprintf('LAYERED_OC_MERMAID_%s_%02d', keyfolders{kk}, ...
                ii_cases(ii)));
    
            try
                % read modeled acoustic pressure seismograms from SPECFEM3D
                [t,x] = read_seismogram(fullfile(ddir, 'OUTPUT_FILES/MH.P0009.CXP.semp'));
                % read the locations of the stations
                [~, ~, ~, ~, ~, z] = readstations3d(fullfile(ddir, 'DATA', 'STATIONS'));
                % Calculate the time of the simulation relative to estimated
                % P-wave arrival time from Simon et al. (2022) by
                % first calculate the time in the SPECFEM3D simulation when
                % P-wave arrive the ocean bottom below the float for the first
                % time
                t0 = calculatearrivaltime(ddir);
                % tTauP = indeks(tauptime('mod', 'ak135', ...
                %     'dep', obs_struct.metadata.EVDP(ic(ii_cases(ii))), ...
                %     'ph', 'p,P,PKP,PKIKP', ...
                %     'deg', obs_struct.metadata.GCARC(ic(ii_cases(ii))), ...
                %     'stdp', -z(1)/1000), 1).time;
                tSPECFEM = t - t0 + (hdr_o.USER8 - hdr_o.T0) - (hdr_s.USER8 - hdr_s.T0); % + tTauP;
                fs_SPECFEM = (length(t) - 1) / (t(end) - t(1));
        
                % interpolate the seismogram to align samples with the observed
                % seismogram
                x = lowpass(detrend(x) .* shanning(length(x), 0.2), fs_SPECFEM, 10, 2, 2, 'butter', 'linear');
                x = shannon(tSPECFEM, x, t_relative);
        
                % apply bandpass filter
                xf = bandpass(detrend(x) .* shanning(length(x), 0.2), fs_o, fc(1), fc(2), 4, 2, ...
                    'butter', 'linear');
        
                % read z-displacement at the ocean bottom from SPECFEM3D
                [tb,xb] = read_seismogram(fullfile(ddir, 'OUTPUT_FILES/AA.OBS01.CXZ.semd'));
        
                % downsample down to 20 hz
                xb = lowpass(detrend(xb) .* shanning(length(xb), 0.2), fs_SPECFEM, 10, 2, 2, 'butter', 'linear');
                xb = downsample(xb, floor(fs_SPECFEM/20));
                tb = downsample(tb, floor(fs_SPECFEM/20));
        
                % apply bandpass filter
                xbf = bandpass(detrend(xb) .* shanning(length(xb), 0.2), fs_SPECFEM / floor(fs_SPECFEM/20), ...
                    fc(1), fc(2), 4, 2, 'butter', 'linear');
        
                % amplitude scaling for SPECFEM seismogram
                amp_xbf = max(abs(xbf(and(tb >= t0-2, tb <= t0+2))));
                amp_seis_sf = max(abs(seis_sf(and(tims_s >= hdr_s.T0 - 10, ...
                    tims_s <= hdr_s.T0 + 30))));
        
                % rescale
                xf = xf * amp_seis_sf / amp_xbf;
        
                % compute the cross-correlation of the envelope
                if have_stats_input
                    CCmax2 = stats.ccmaxs(ii_stats, 2);
                    t_shift3D = stats.t_shifts(ii_stats, 2);
                else
                    [t_shift1, CCmax1, lags1, cc1, Smax1, s1] = ccscale(seis_o3, xf, dt_now + seconds(t_relative_o3(1)), dt_now + seconds(t_relative(1)), fs_o, seconds(15), 'soft', true, false);
                    [t_shift2, CCmax2, lags2, cc2, Smax2, s2] = ccscale(seis_o3, xf, dt_now + seconds(t_relative_o3(1)), dt_now + seconds(t_relative(1) +t_shift1), fs_o, seconds(2), 'hard', false, false);
                    t_shift3D = t_shift1 + t_shift2;
                end
        
                wh_s3d = and(t_relative + t_shift3D >= -10, t_relative + t_shift3D <= 30);
                l3d = plot(t_relative + t_shift3D, xf / max(abs(xf(wh_s3d))) - 2.5, 'LineWidth', 0.75, 'Color', [0 0.2 0.9]);
        
                xlabel('time since picked P-wave (s)')
                legend(sprintf('2D (cc: %.2f, \\Delta\\tau: %.2f s)', ...
                    CCmax0, t_shift), 'observed', ...
                    sprintf('3D (cc: %.2f, \\Delta\\tau: %.2f s)', ...
                    CCmax2, t_shift3D), 'Location', 'southoutside')
                yticks([-2.5 0 2.5]);
                nolabels(gca, 2);
                set(gca, 'Position', [0.038+xposshift(ii)  0.14 0.15 0.26], ...
                    'YLim', [-4.5 4])
        
                % ask the user if they like the automatic correlation
                if handpick
                    prompt = "Do you like this? Y/N/Q [Y]: ";
                    txt = upper(input(prompt,"s"));
                    if isempty(txt)
                        txt = 'N';
                    end
            
                    if strcmp(txt, 'Q')
                        return
                    elseif strcmp(txt, 'K')
                        keyboard
                    end
            
                    % handpick system
                    [t_shift3, CCmax3, lags3, cc3, Smax3, s3] = ccscale(seis_o3, xf, dt_now + seconds(t_relative_o3(1)), dt_now + seconds(t_relative(1)), fs_o, seconds(20), 'hard', false, false);
                    cc4 = cc3;
                    cc4(cc4 < 0.25 * max(cc4)) = 0;
                    [pks, locs] = findpeaks(cc4);
                    ii_locs = ceil(length(locs) / 2);
                    while ~strcmp(txt, 'Y') && ii_locs < length(locs)
                        t_shift3D = lags3(locs(ii_locs));
                        CCmax2 = pks(ii_locs);
                        Smax2 = s3(locs(ii_locs));
            
                        % redraw the figure
                        delete(l3d);
                        wh_s3d = and(t_relative + t_shift3D >= -10, t_relative + t_shift3D <= 30);
                        l3d = plot(t_relative + t_shift3D, xf / max(abs(xf(wh_s3d))) - 2.5, 'LineWidth', 0.75, 'Color', [0 0.2 0.9]);
                        
                        xlabel('time since picked P-wave (s)')
                        legend(sprintf('2D (cc: %.2f, \\Delta\\tau: %.2f s)', ...
                            CCmax0, t_shift), 'observed', ...
                            sprintf('3D (cc: %.2f, \\Delta\\tau: %.2f s)', ...
                            CCmax2, t_shift3D), 'Location', 'southoutside')
                        yticks([-2.5 0 2.5]);
                        nolabels(gca, 2);
            
                        % ask if they like this
                        prompt = sprintf("Do you like this? (%d/%d) Y/N/P/Q [Y]: ", ii_locs, length(locs));
                        txt = upper(input(prompt,"s"));
                        if isempty(txt)
                            txt = 'N';
                        end
            
                        if strcmp(txt, 'Q')
                            return
                        elseif strcmp(txt, 'K')
                            keyboard
                        end
            
                        % move to the next candidate
                        if strcmp(txt, 'P')
                            ii_locs = max(ii_locs - 1, 1);
                        elseif strcmp(txt, 'Y')
                        else
                            ii_locs = ii_locs + 1;
                        end
                    end
                end
                % calculate the correlation coefficient between 2D and 3D runs
                wh2d = and(t_relative + t_shift >= -10, t_relative + t_shift <= 30);
                wh3d = and(t_relative + t_shift3D >= -10, t_relative + t_shift3D <= 30);
            
                seis_p2 = seis_p2(wh2d);
                seis_p3 = xf(wh3d);
            
                length_p = min(length(seis_p2), length(seis_p3));
            
                [R, L] = xcorr(seis_p2(1:length_p), seis_p3(1:length_p), 'coeff');
                CCmax_2D_3D = max(R);

                title(sprintf('2D-3D CC: %.2f', CCmax_2D_3D), 'FontWeight', 'normal')
                subtitle(' ', 'FontSize', 2)

                % calculate the correlation coefficient to the observed
                % seismogram over the entire window (-10 to 30 seconds).
                seis_o4 = seis_o2(wh);
                M = corrcoef(seis_o4, seis_p2);
                stats.ccmaxs(ii_stats, 4) = M(1, 2);
                text(-8, 3, sprintf('%.2f', M(1, 2)), 'FontSize', 9, 'Color', 'r', 'VerticalAlignment', 'middle')
                M = corrcoef(seis_o4, seis_p3);
                stats.ccmaxs(ii_stats, 5) = M(1, 2);
                text(-8, -2, sprintf('%.2f', M(1, 2)), 'FontSize', 9, 'Color', [0 0.2 0.9], 'VerticalAlignment', 'middle')
            catch ME
                t_shift3D = nan;
                CCmax2 = nan;
                CCmax_2D_3D = nan;

                plot([-10 30], [-2.5 -2.5], 'LineWidth', 0.75, 'Color', [0 0.2 0.9], 'LineStyle', '--')
                xlabel('time since picked P-wave (s)')
                legend(sprintf('2D (cc: %.2f, \\Delta\\tau: %.2f s)', ...
                    CCmax0, t_shift), 'observed', ...
                    sprintf('3D (cc: %s, \\Delta\\tau: %s s)', ...
                    'n/a', 'n/a'), 'Location', 'southoutside')
                yticks([-2.5 0 2.5]);
                nolabels(gca, 2);
                set(gca, 'Position', [0.038+xposshift(ii)  0.14 0.15 0.26], ...
                    'YLim', [-4.5 4])
                title(sprintf('2D-3D CC: %s', 'n/a'), 'FontWeight', 'normal')
            end

            % Highlight [-5 5] window
            [xbox, ybox] = boxcorner([-5 5], ax3.YLim);
            pgon = polyshape(xbox, ybox);
            bx = plot(ax3, pgon, 'FaceColor', [1 0.9 0.4], ...
                'FaceAlpha', 0.4, 'EdgeColor', 'none', 'EdgeAlpha', 1, ...
                'HandleVisibility', 'off');
            uistack(bx, 'bottom');
    
            % Bandpass corner frequency information
            ax3_fc = axes('Position', [0.128+xposshift(ii) 0.14 0.06 0.025]);
            axeslabel(ax3_fc, 0.5, 0.5, ...
                sprintf('$$%.2f-%.2f \\textnormal{ Hz}$$', fc(1), fc(2)), ...
                'FontSize', 8, 'HorizontalAlignment', 'center', ...
                'Interpreter', 'latex');
            nolabels(ax3_fc, 3)
            set(ax3_fc, 'Box', 'on', 'XTick', [], 'YTick', [])

            % gather pick stats
            stats.is_low(ii_stats) = (kk == 1);
            stats.id(ii_stats) = ic(ii_cases(ii));
            stats.mean(ii_stats, 2) = mean_zz;
            stats.std(ii_stats, 2) = std_zz;
            stats.t_shifts(ii_stats, 1) = t_shift;
            stats.t_shifts(ii_stats, 2) = t_shift3D;
            stats.ccmaxs(ii_stats, 1) = CCmax0;
            stats.ccmaxs(ii_stats, 2) = CCmax2;
            stats.ccmaxs(ii_stats, 3) = CCmax_2D_3D;
        end
    
        % Run through and set the same axes
    
        for ii=1:5
            axes(ah(ii))
            ylim(getzm(ii)+1.1*[-1 1]*max(getzs))
        end
    
        set(gcf, 'Renderer', 'painters')
        figdisp(sprintf('%s_%s_%02d-%02d', mfilename, keyfolders{kk}, ii_cases(1), ii_cases(end)), [], [], 2, [], 'epstopdf');
    end
end
end
