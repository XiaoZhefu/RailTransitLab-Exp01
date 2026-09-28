clear
clc
close all

postDir = fileparts(mfilename('fullpath'));
sourceDir = fullfile(postDir, '..', 'source');
fileNames = {'01.txt', '02.txt', '03.txt'};

%% 报告采用的频率与衰减指数
centerFrequencyHz = [8.9093; 6.4688; 5.3623];
decayRatePerS = [0.462125284875350; 0.291824351293369; 0.225818382564419];

%% 读取波形
[timeByCase, accelerationByCase, deltaTByCase] = readWaveforms(sourceDir, fileNames);

%% 绘制时域波形
redrawTimeWaveforms(timeByCase, accelerationByCase, postDir);

%% 绘制频谱并提取频率
redrawFrequencySpectra(timeByCase, accelerationByCase, deltaTByCase, postDir);

%% 绘制峰值并拟合衰减指数
redrawDecayPeaks(timeByCase, accelerationByCase, deltaTByCase, ...
    centerFrequencyHz, postDir);

%% 由两组频率比分别求等效质量
solveDeltaM(centerFrequencyHz);

%% 联合拟合质量、刚度和阻尼
format long g
fitMassStiffnessDamping(centerFrequencyHz, decayRatePerS);


function [timeByCase, accelerationByCase, deltaTByCase] = readWaveforms(sourceDir, fileNames)
    timeByCase = cell(numel(fileNames), 1);
    accelerationByCase = cell(numel(fileNames), 1);
    deltaTByCase = zeros(numel(fileNames), 1);

    for i = 1:numel(fileNames)
        %% 导入波形文件并按时间戳去重
        opts = delimitedTextImportOptions("NumVariables", 2);
        opts.Delimiter = "\t";
        opts.VariableNames = ["time", "Y"];
        opts.VariableTypes = ["string", "double"];
        opts.EmptyLineRule = "skip";
        data = readtable(fullfile(sourceDir, fileNames{i}), opts);

        timeText = strtrim(data.time);
        yAll = data.Y;
        deltat = yAll(find(timeText == "delta t", 1));
        sampleRows = contains(timeText, "/") & isfinite(yAll);
        sampleTime = timeText(sampleRows);
        sampleY = yAll(sampleRows);
        [~, firstIndex] = unique(sampleTime, 'stable');
        x = sampleY(firstIndex);

        deltaTByCase(i) = deltat;
        accelerationByCase{i} = x;
        timeByCase{i} = (0:numel(x)-1).' * deltat;
    end
end


function redrawTimeWaveforms(timeByCase, accelerationByCase, postDir)
    %% 图像设置
    fig = figure(1);
    set(fig, 'Units', 'centimeters', 'Position', [3, 3, 25, 15], ...
        'PaperPositionMode', 'auto', 'Color', 'w')
    layout = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    matlabBlue = [0, 0.4470, 0.7410];
    fontSize = 12;

    for i = 1:numel(accelerationByCase)
        t = timeByCase{i};
        x = accelerationByCase{i};

        %% 绘制时域波形
        ax = nexttile(layout, i);
        plot(ax, t, x, 'Color', matlabBlue, 'LineWidth', 0.75)
        grid(ax, 'on')
        set(ax, 'FontName', 'Times New Roman', 'FontSize', fontSize)
        title(ax, ['\fontname{SimSun}加速度波形图 \fontname{Times New Roman}- ' ...
            num2str(i) ' kg'], 'Interpreter', 'tex', ...
            'FontName', 'SimSun', 'FontSize', fontSize)
        xlabel(ax, '时间/秒', 'FontName', 'SimSun', 'FontSize', fontSize)
        ylabel(ax, '\fontname{SimSun}加速度/\fontname{Times New Roman}\it g', ...
            'Interpreter', 'tex', 'FontName', 'SimSun', 'FontSize', fontSize)
        xlim(ax, [t(1), t(end)])
        if i == 2
            ylim(ax, [-0.5, 0.5])
            yticks(ax, [-0.5, -0.25, 0, 0.25, 0.5])
        end
    end

    %% 保存图像
    print(fig, fullfile(postDir, '04_time_waveforms.png'), '-dpng', '-r300');
end


function redrawFrequencySpectra( ...
    timeByCase, accelerationByCase, deltaTByCase, postDir)
    %% 频谱分析设置
    analysisWindowsSec = [0, 12; 0, 12; 0, 12];
    displayBandHz = [0, 10];
    peakSearchBandHz = [1, 50];
    centerHalfWidthHz = 5;

    %% 图像设置
    fig = figure(2);
    set(fig, 'Units', 'centimeters', 'Position', [3, 3, 25, 8], ...
        'PaperPositionMode', 'auto', 'Color', 'w')
    layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    matlabBlue = [0, 0.4470, 0.7410];
    fontSize = 12;

    peakFrequencyHz = zeros(3, 1);
    weightedCenterHz = zeros(3, 1);
    modeBandLowHz = zeros(3, 1);
    modeBandHighHz = zeros(3, 1);

    for i = 1:numel(accelerationByCase)
        t = timeByCase{i};
        x = accelerationByCase{i};
        deltat = deltaTByCase(i);

        %% 计算单边幅值谱
        timeWindow = analysisWindowsSec(i, :);
        selected = t >= timeWindow(1) & t <= timeWindow(2);
        xSegment = x(selected);
        N = numel(xSegment);
        xSegment = xSegment - mean(xSegment);
        hannWindow = 0.5 - 0.5*cos(2*pi*(0:N-1).' / N);
        xWindowed = xSegment .* hannWindow;

        fs = 1 / deltat;
        spectrum = fft(xWindowed);
        positiveBins = floor(N / 2) + 1;
        frequencyHz = (0:positiveBins-1).' * (fs / N);
        amplitudeG = abs(spectrum(1:positiveBins)) / sum(hannWindow);
        if mod(N, 2) == 0
            amplitudeG(2:end-1) = 2 * amplitudeG(2:end-1);
        else
            amplitudeG(2:end) = 2 * amplitudeG(2:end);
        end

        %% 提取主峰与加权中心频率
        peakBins = frequencyHz >= peakSearchBandHz(1) & ...
            frequencyHz <= peakSearchBandHz(2);
        peakFrequencies = frequencyHz(peakBins);
        [~, peakIndex] = max(amplitudeG(peakBins));
        peakFrequencyHz(i) = peakFrequencies(peakIndex);

        modeBandLowHz(i) = max(peakSearchBandHz(1), ...
            peakFrequencyHz(i) - centerHalfWidthHz);
        modeBandHighHz(i) = min(peakSearchBandHz(2), ...
            peakFrequencyHz(i) + centerHalfWidthHz);
        modeBins = frequencyHz >= modeBandLowHz(i) & ...
            frequencyHz <= modeBandHighHz(i);
        spectralIntensity = amplitudeG(modeBins).^2;
        weightedCenterHz(i) = sum(frequencyHz(modeBins) .* spectralIntensity) ...
            / sum(spectralIntensity);

        %% 绘制频谱
        visible = frequencyHz >= displayBandHz(1) & frequencyHz <= displayBandHz(2);
        ax = nexttile(layout, i);
        stem(ax, frequencyHz(visible), amplitudeG(visible), ...
            'Marker', 'none', 'Color', matlabBlue, 'LineWidth', 0.75)
        grid(ax, 'on')
        set(ax, 'FontName', 'Times New Roman', 'FontSize', fontSize)
        title(ax, ['\fontname{SimSun}加速度频谱图 \fontname{Times New Roman}- ' ...
            num2str(i) ' kg'], 'Interpreter', 'tex', ...
            'FontName', 'SimSun', 'FontSize', fontSize)
        xlabel(ax, '\fontname{SimSun}频率/\fontname{Times New Roman}Hz', ...
            'Interpreter', 'tex', 'FontName', 'SimSun', 'FontSize', fontSize)
        ylabel(ax, '\fontname{SimSun}幅值/\fontname{Times New Roman}\it g', ...
            'Interpreter', 'tex', 'FontName', 'SimSun', 'FontSize', fontSize)
        xlim(ax, displayBandHz)
    end

    %% 输出频率
    results = table((1:3).', analysisWindowsSec(:, 1), analysisWindowsSec(:, 2), ...
        peakFrequencyHz, weightedCenterHz, modeBandLowHz, modeBandHighHz, ...
        'VariableNames', {'Mass_kg', 'Start_s', 'End_s', 'Peak_Hz', ...
        'WeightedCenter_Hz', 'BandLow_Hz', 'BandHigh_Hz'});
    disp(results)

    %% 保存图像
    print(fig, fullfile(postDir, '05_frequency_spectra_0_10Hz.png'), '-dpng', '-r300');
end


function redrawDecayPeaks( ...
    timeByCase, accelerationByCase, deltaTByCase, centerFrequencyHz, postDir)
    %% 峰值选取设置
    fitWindowsSec = [0.1, 5.0; 0.1, 8.0; 0.1, 8.0];

    %% 图像设置
    fig = figure(3);
    set(fig, 'Units', 'centimeters', 'Position', [3, 3, 25, 8], ...
        'PaperPositionMode', 'auto', 'Color', 'w')
    layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    matlabBlue = [0, 0.4470, 0.7410];
    markerOrange = [0.8500, 0.3250, 0.0980];
    fitGray = [0.2, 0.2, 0.2];
    fontSize = 12;
    peakCounts = zeros(3, 1);
    decayRates = zeros(3, 1);
    peakTimesByCase = cell(3, 1);
    logPeaksByCase = cell(3, 1);
    fitParameters = zeros(3, 2);

    for i = 1:numel(accelerationByCase)
        t = timeByCase{i};
        x = accelerationByCase{i};
        deltat = deltaTByCase(i);

        %% 提取峰值并拟合衰减指数
        baseline = median(x(t >= t(end)-2));
        x = x - baseline;
        noiseStd = std(x(t >= t(end)-2), 1);
        minDistanceSamples = floor(0.75 / (centerFrequencyHz(i) * deltat));
        [peakValues, peakLocations] = findpeaks(x, ...
            'MinPeakDistance', minDistanceSamples);
        peakTimes = t(peakLocations);
        selectedPeaks = peakTimes >= fitWindowsSec(i, 1) & ...
            peakTimes <= fitWindowsSec(i, 2) & peakValues > 5*noiseStd;
        peakTimes = peakTimes(selectedPeaks);
        peakValues = peakValues(selectedPeaks);
        peakCounts(i) = numel(peakTimes);
        fitCoefficients = polyfit(peakTimes, log(peakValues), 1);
        decayRates(i) = -fitCoefficients(1);
        peakTimesByCase{i} = peakTimes;
        logPeaksByCase{i} = log(peakValues);
        fitParameters(i, :) = fitCoefficients;

        %% 绘制波形、峰值与包络
        visible = t >= fitWindowsSec(i, 1) & t <= fitWindowsSec(i, 2);
        ax = nexttile(layout, i);
        plot(ax, t(visible), x(visible), 'Color', matlabBlue, 'LineWidth', 0.75)
        hold(ax, 'on')
        plot(ax, t(visible), exp(polyval(fitCoefficients, t(visible))), '--', ...
            'Color', fitGray, 'LineWidth', 1.2)
        plot(ax, peakTimes, peakValues, 'o', 'LineStyle', 'none', ...
            'MarkerEdgeColor', markerOrange, 'MarkerFaceColor', 'w', ...
            'MarkerSize', 4, 'LineWidth', 1.1)
        grid(ax, 'on')
        set(ax, 'FontName', 'Times New Roman', 'FontSize', fontSize)
        title(ax, ['\fontname{SimSun}加速度衰减波形图 \fontname{Times New Roman}- ' ...
            num2str(i) ' kg'], 'Interpreter', 'tex', ...
            'FontName', 'SimSun', 'FontSize', fontSize)
        xlabel(ax, '\fontname{SimSun}时间/\fontname{Times New Roman}s', ...
            'Interpreter', 'tex', 'FontName', 'SimSun', 'FontSize', fontSize)
        ylabel(ax, '\fontname{SimSun}加速度/\fontname{Times New Roman}\it g', ...
            'Interpreter', 'tex', 'FontName', 'SimSun', 'FontSize', fontSize)
        xlim(ax, fitWindowsSec(i, :))
    end

    %% 输出衰减指数
    disp(table((1:3).', peakCounts, decayRates, ...
        'VariableNames', {'Mass_kg', 'PeakCount', 'DecayRate_per_s'}))

    %% 保存图像
    print(fig, fullfile(postDir, '06_decay_peaks.png'), '-dpng', '-r300');

    %% 绘制峰值对数拟合
    fitFig = figure(4);
    set(fitFig, 'Units', 'centimeters', 'Position', [3, 3, 25, 8], ...
        'PaperPositionMode', 'auto', 'Color', 'w')
    fitLayout = tiledlayout(fitFig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:numel(accelerationByCase)
        ax = nexttile(fitLayout, i);
        plot(ax, peakTimesByCase{i}, logPeaksByCase{i}, 'o', ...
            'LineStyle', 'none', 'MarkerEdgeColor', markerOrange, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 1.1)
        hold(ax, 'on')
        fitTime = linspace(fitWindowsSec(i, 1), fitWindowsSec(i, 2), 200);
        plot(ax, fitTime, polyval(fitParameters(i, :), fitTime), '--', ...
            'Color', fitGray, 'LineWidth', 1.2)
        grid(ax, 'on')
        set(ax, 'FontName', 'Times New Roman', 'FontSize', fontSize)
        title(ax, ['\fontname{SimSun}峰值对数线性拟合 \fontname{Times New Roman}- ' ...
            num2str(i) ' kg'], 'Interpreter', 'tex', ...
            'FontName', 'SimSun', 'FontSize', fontSize)
        xlabel(ax, '\fontname{SimSun}时间/\fontname{Times New Roman}s', ...
            'Interpreter', 'tex', 'FontName', 'SimSun', 'FontSize', fontSize)
        ylabel(ax, '\fontname{Times New Roman}ln(\it A\rm /g)', ...
            'Interpreter', 'tex', 'FontName', 'Times New Roman', 'FontSize', fontSize)
        xlim(ax, fitWindowsSec(i, :))
        slope = fitParameters(i, 1);
        intercept = fitParameters(i, 2);
        if intercept < 0
            equationText = sprintf('ln(\\it A\\rm /g) = %.4f\\it t\\rm - %.4f', ...
                slope, abs(intercept));
        else
            equationText = sprintf('ln(\\it A\\rm /g) = %.4f\\it t\\rm + %.4f', ...
                slope, intercept);
        end
        text(ax, 0.04, 0.06, equationText, 'Units', 'normalized', ...
            'Interpreter', 'tex', 'FontName', 'Times New Roman', ...
            'FontSize', 10, 'BackgroundColor', 'w', 'Margin', 2)
    end

    %% 保存图像
    print(fitFig, fullfile(postDir, '07_log_peak_linear_fit.png'), '-dpng', '-r300');
end


function solveDeltaM(centerFrequencyHz)
    m = 1.0; % kg
    a = [2, 3];
    ratio = centerFrequencyHz(1) ./ centerFrequencyHz(2:3).';

    deltaM = (a - ratio.^2) ./ (ratio.^2 - 1) * m;
    fprintf('第一式：Delta m = %.4f kg\n', deltaM(1));
    fprintf('第二式：Delta m = %.4f kg\n', deltaM(2));
end


function fitMassStiffnessDamping(centerFrequencyHz, decayRatePerS)
    %% 修正无阻尼固有频率
    massKg = [1; 2; 3];
    naturalFrequencyHz = sqrt(centerFrequencyHz.^2 + (decayRatePerS/(2*pi)).^2);

    fprintf('Inputs and frequency correction:\n')
    disp(table(massKg, centerFrequencyHz, decayRatePerS, ...
        naturalFrequencyHz-centerFrequencyHz, naturalFrequencyHz, ...
        'VariableNames', {'Mass_kg', 'CenterFrequency_Hz', 'DecayRate_per_s', ...
        'FrequencyCorrection_Hz', 'NaturalFrequency_Hz'}))

    %% 固定 Delta m 时线性拟合 beta
    initialDeltaMass = 0.1; % kg，仅作为一维优化初值
    basisAt = @(deltaM) 1./sqrt(massKg+deltaM);
    scaleAt = @(deltaM) (basisAt(deltaM).'*naturalFrequencyHz) ...
        /(basisAt(deltaM).'*basisAt(deltaM));
    frequencyModel = @(deltaM) scaleAt(deltaM)*basisAt(deltaM);
    sumSquaredError = @(logDeltaM) sum( ...
        (naturalFrequencyHz-frequencyModel(exp(logDeltaM))).^2);
    options = optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-14);
    [logDeltaMass, SSE, exitFlag, optimizationInfo] = fminsearch( ...
        sumSquaredError, log(initialDeltaMass), options);

    deltaMassKg = exp(logDeltaMass);
    modelBasis = basisAt(deltaMassKg);
    scaleHzSqrtKg = scaleAt(deltaMassKg);
    stiffnessNPerM = (2*pi*scaleHzSqrtKg)^2;
    fittedFrequencyHz = scaleHzSqrtKg*modelBasis;
    frequencyResidualHz = naturalFrequencyHz-fittedFrequencyHz;
    equivalentMassKg = massKg+deltaMassKg;
    dampingRatio = decayRatePerS./(2*pi*naturalFrequencyHz);
    dampingNsPerM = 2*equivalentMassKg.*decayRatePerS;

    %% 输出拟合结果
    fprintf('Nonlinear least-squares fit:\n')
    fprintf('  Initial Delta m = %.4f kg\n', initialDeltaMass)
    fprintf('  Profiled beta = %.10f Hz*sqrt(kg)\n', scaleHzSqrtKg)
    fprintf('  Fitted k = %.6f N/m, Delta m = %.6f kg\n', stiffnessNPerM, deltaMassKg)
    fprintf('  SSE = %.10g Hz^2, RMSE = %.10g Hz\n', ...
        SSE, sqrt(mean(frequencyResidualHz.^2)))
    fprintf('  exitFlag = %d, iterations = %d, function evaluations = %d\n', ...
        exitFlag, optimizationInfo.iterations, optimizationInfo.funcCount)

    fprintf('Per-condition fitted values, residuals, and damping:\n')
    disp(table(massKg, modelBasis, naturalFrequencyHz, fittedFrequencyHz, ...
        frequencyResidualHz, equivalentMassKg, dampingRatio, dampingNsPerM, ...
        'VariableNames', {'Mass_kg', 'Basis_per_sqrt_kg', 'MeasuredNatural_Hz', 'FittedNatural_Hz', ...
        'Residual_Hz', 'EquivalentMass_kg', 'DampingRatio', 'Damping_N_s_per_m'}))
    fprintf('Mean damping coefficient = %.6f N*s/m\n', mean(dampingNsPerM))
end
