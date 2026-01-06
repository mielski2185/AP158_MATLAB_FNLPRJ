%% Topological Analysis: Multi-Window Visualization
%  Generates separate figures for Physics, Fluence, Cell Visuals, and Networks.

clear; clc; close all;
addpath("BCT\2019_03_03_BCT");
%% --- 1. CONFIGURATION ---

% A. RADIATION & PHYSICS SETUP
radType = 'Photon';         % Try: 'Photon', 'Proton', 'Alpha'
surfaceDoseVals = 0:2:16;   % Discrete energy levels (Gy) applied at surface
targetDepth = 4.0;          % Tumor location (cm)
tissueThickness = 10;       % Max depth (cm)

% B. CELL FACTORS
cellProps.name = 'Deep Tumor (Glioblastoma)';
cellProps.alpha = 0.04;      
cellProps.beta  = 0.025;     
cellProps.attenuation = 0.2; 
baseColonySize = 150;
numSimulations = 50;                                % Reduced slightly for speed
netThreshold = 20;                                  % Pixel distance for network connection

fprintf('Simulating %s transport to depth %.1f cm...\n', radType, targetDepth);

%% --- 2. PHYSICS ENGINE ---

% Depth Axis for Physics Plots
depthAxis = linspace(0, tissueThickness, 200);

% 1. Dose Profile (Energy Deposition)
doseProfile = getDepthDoseProfile(radType, depthAxis, targetDepth);
doseModAtTarget = interp1(depthAxis, doseProfile, targetDepth);

% 2. Fluence Profile (Particle Count)
fluenceProfile = getFluenceProfile(radType, depthAxis, targetDepth);

% RBE
rbe = getRBE(radType);

%% --- 3. MAIN SIMULATION LOOP ---

results = struct();
avg_b0 = zeros(length(surfaceDoseVals), 1);
avg_surv = zeros(length(surfaceDoseVals), 1);

avg_ecc  = zeros(length(surfaceDoseVals),1);
avg_area = zeros(length(surfaceDoseVals),1);
avg_per  = zeros(length(surfaceDoseVals),1);
avg_rnd  = zeros(length(surfaceDoseVals),1);

fprintf('Running Simulations');
% Store all simulation metrics for statistical analysis
all_C = zeros(numSimulations, length(surfaceDoseVals));
all_L = zeros(numSimulations, length(surfaceDoseVals));
all_sigma = zeros(numSimulations, length(surfaceDoseVals));

for d = 1:length(surfaceDoseVals)
    fprintf('.');
    surfDose = surfaceDoseVals(d);
    
    % Effective Bio Dose = Surface * DepthFactor * RBE
    effBioDose = surfDose * doseModAtTarget * rbe;
    
    temp_C = zeros(numSimulations,1);
    temp_L = zeros(numSimulations,1);
    temp_sigma = zeros(numSimulations,1);

    temp_ecc  = zeros(numSimulations,1);
    temp_area = zeros(numSimulations,1);
    temp_per  = zeros(numSimulations,1);
    temp_rnd  = zeros(numSimulations,1);

    temp_b0 = zeros(numSimulations, 1);
    temp_surv = zeros(numSimulations, 1);
    
    
    
    for sim = 1:numSimulations
        % Generate & Irradiate
        [mask, dens] = generateColony(baseColonySize);
        survMask = applyBioResponse(mask, dens, effBioDose, cellProps);

        % Analyze Topology
        cc = bwconncomp(survMask, 8);
        temp_b0(sim) = cc.NumObjects;
        temp_surv(sim) = sum(survMask(:)) / max(sum(mask(:)), 1);

        nodeMetrics = computeNodeAndShapeMetrics(survMask);

        temp_ecc(sim)  = nodeMetrics.meanEccentricity;
        temp_area(sim) = nodeMetrics.meanArea;
        temp_per(sim)  = nodeMetrics.meanPerimeter;
        temp_rnd(sim)  = nodeMetrics.meanRoundness;


        A = buildAdjacency(survMask, netThreshold);
        netMetrics = computeNetworkMetrics(A);
        
        temp_C(sim) = netMetrics.C;
        temp_L(sim) = netMetrics.L;
        temp_sigma(sim) = netMetrics.sigma;

        
        % Save first sim image for visualization
        if sim == 1
            results(d).Image = survMask;
            results(d).Dose = surfDose;
        end
    end

    all_C(:, d) = temp_C;
    all_L(:, d) = temp_L;
    all_sigma(:, d) = temp_sigma;


    avg_b0(d) = mean(temp_b0);
    avg_surv(d) = mean(temp_surv);

    avg_C(d) = nanmean(temp_C);
    avg_L(d) = nanmean(temp_L);
    avg_sigma(d) = nanmean(temp_sigma);

    avg_ecc(d)  = nanmean(temp_ecc);
    avg_area(d) = nanmean(temp_area);
    avg_per(d)  = nanmean(temp_per);
    avg_rnd(d)  = nanmean(temp_rnd);

end
fprintf(' Done!\n');



%% --- 4. VISUALIZATION (SEPARATE FIGURES) ---

% === FIGURE 1: Target Location & Depth Dose ===
figure('Name', 'Fig 1: Dose & Target Location', 'Color', 'w', 'Position', [50, 500, 600, 400]);
plot(depthAxis, doseProfile, 'LineWidth', 3, 'Color', [0.8 0.2 0.2]); % Red Line
hold on;
xline(targetDepth, '--k', 'LineWidth', 2); % Target Line
scatter(targetDepth, doseModAtTarget, 100, 'k', 'filled'); % Target Dot

title(['Radiation Dose vs Depth (' radType ')']);
xlabel('Depth in Tissue (cm)');
ylabel('Relative Energy Deposition (Normalized)');
legend('Dose Profile', 'Tumor Location', 'Dose at Tumor');
grid on;
subtitle(['Target Depth: ' num2str(targetDepth) ' cm']);


% === FIGURE 2: Fluence Profile ===
figure('Name', 'Fig 2: Particle Fluence', 'Color', 'w', 'Position', [700, 500, 600, 400]);
area(depthAxis, fluenceProfile, 'FaceColor', [0.2 0.6 0.2], 'FaceAlpha', 0.3, 'EdgeColor', [0.2 0.6 0.2], 'LineWidth', 2);
title(['Particle Fluence vs Depth (' radType ')']);
xlabel('Depth in Tissue (cm)');
ylabel('Particle Count / Flux (Normalized)');
xline(targetDepth, '--k', 'Target Depth');
grid on;
subtitle('Shows how many particles survive to reach the depth');


% === FIGURE 3: Cancer Cells at Different Energy Levels ===
figure('Name', 'Fig 3: Cell Morphology at Different Levels', 'Color', 'w', 'Position', [50, 50, 1200, 350]);
t = tiledlayout(1, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
indices = floor(linspace(1, length(surfaceDoseVals), 6)); % Pick 6 evenly spaced samples

for i = 1:6
    idx = indices(i);
    nexttile;
    % Show inverted image (Black cells on white background)
    imshow(~results(idx).Image); 
    title([num2str(results(idx).Dose) ' Gy']);
    if i == 1; ylabel('Cell Structure'); end
end
title(t, 'Evolution of Cancer Colony Structure under Radiation', 'FontSize', 14);


% === FIGURE 4: Network Topology Analysis ===
figure('Name', 'Fig 4: Topological Networks', 'Color', 'w', 'Position', [800, 50, 800, 500]);

% Generate a high-dose example for the network demo
[netMask, netDens] = generateColony(baseColonySize);
doseForNet = surfaceDoseVals(end) * doseModAtTarget * rbe;
netSurvivor = applyBioResponse(netMask, netDens, doseForNet, cellProps);

subplot(1, 2, 1);
plotCellNetwork(netMask, netThreshold, 'Pre-Radiation Network');

subplot(1, 2, 2);
plotCellNetwork(netSurvivor, netThreshold, 'Post-Radiation Network');


% === FIGURE 5: Quantitative Metrics (Beta-0 & Survival) ===
figure('Name', 'Fig 5: Quantitative Metrics', 'Color', 'w', 'Position', [1000, 500, 500, 400]);
yyaxis left
plot(surfaceDoseVals, avg_b0, '-o', 'LineWidth', 2);
ylabel('Fragmentation (\beta_0 Count)');

yyaxis right
plot(surfaceDoseVals, avg_surv * 100, '--', 'LineWidth', 2);
ylabel('Survival Fraction (%)');

xlabel('Surface Dose (Gy)');
title('Topological & Biological Statistics');
legend('Fragmentation', 'Survival');
grid on;

% === FIGURE 6: Topological Data ===
figure('Name','Fig 6: Network Metrics','Color','w','Position', [1100, 400, 1000, 400]);
tiledlayout(1,3);

nexttile
plot(surfaceDoseVals, avg_C, '-o','LineWidth',2);
ylabel('Clustering Coefficient');
xlabel('Surface Dose (Gy)');
grid on
axis square

nexttile
plot(surfaceDoseVals, avg_L, '-s','LineWidth',2);
ylabel('Characteristic Path Length');
xlabel('Surface Dose (Gy)');
grid on
axis square

nexttile
plot(surfaceDoseVals, avg_sigma, '-^','LineWidth',2);
ylabel('Small-World Coefficient');
xlabel('Surface Dose (Gy)');
grid on
axis square

title(t,'Radiation-Induced Network Topology Changes');

% === FIGURE 7: Morphological Data ===
figure('Name','Fig 7: Morphology','Color','w','Position', [1200, 300, 1000, 400]);
tiledlayout(1,4);

nexttile
plot(surfaceDoseVals, avg_ecc, '-o','LineWidth',2);
ylabel('Mean Node Eccentricity');
xlabel('Dose (Gy)');
grid on
axis square

nexttile
plot(surfaceDoseVals, avg_area, '-s','LineWidth',2);
ylabel('Mean Area (px)');
xlabel('Dose (Gy)');
grid on
axis square

nexttile
plot(surfaceDoseVals, avg_per, '-^','LineWidth',2);
ylabel('Mean Perimeter (px)');
xlabel('Dose (Gy)');
grid on
axis square

nexttile
plot(surfaceDoseVals, avg_rnd, '-d','LineWidth',2);
ylabel('Mean Roundness');
xlabel('Dose (Gy)');
grid on
axis square

figure('Name','Fig 6: Network Metrics with Significance','Color','w','Position', [1100, 400, 1000, 400]);
tiledlayout(1,3);

nexttile
plot(surfaceDoseVals, avg_C, '-o','LineWidth',2); hold on;
sigIdx = find(p_C < 0.05);
scatter(surfaceDoseVals(sigIdx), avg_C(sigIdx), 100, 'r', 'filled');
ylabel('Clustering Coefficient'); xlabel('Surface Dose (Gy)'); grid on; axis square

nexttile
plot(surfaceDoseVals, avg_L, '-s','LineWidth',2); hold on;
sigIdx = find(p_L < 0.05);
scatter(surfaceDoseVals(sigIdx), avg_L(sigIdx), 100, 'r', 'filled');
ylabel('Characteristic Path Length'); xlabel('Surface Dose (Gy)'); grid on; axis square

nexttile
plot(surfaceDoseVals, avg_sigma, '-^','LineWidth',2); hold on;
sigIdx = find(p_sigma < 0.05);
scatter(surfaceDoseVals(sigIdx), avg_sigma(sigIdx), 100, 'r', 'filled');
ylabel('Small-World Coefficient'); xlabel('Surface Dose (Gy)'); grid on; axis square

title(tiledlayout,'Radiation-Induced Network Topology Changes (Significant Doses Marked)');


%% --- 5. HELPER FUNCTIONS ---

function profile = getFluenceProfile(type, depth, range)
    switch type
        case 'Photon', profile = exp(-0.15 .* depth); % Attenuates
        case 'Proton', profile = 1 ./ (1 + exp(20 * (depth - range))); % Stops at range
        case 'Alpha',  profile = double(depth <= 0.5); % Stops immediately
        case 'Electron', profile = double(depth < (range+1)); 
    end
end

function profile = getDepthDoseProfile(type, depth, range)
    switch type
        case 'Photon'
            profile = exp(-0.15 .* depth);
            profile(depth < 0.5) = linspace(0.5, 1, sum(depth < 0.5));
        case 'Proton'
            sigma = 0.5; base = 0.3;
            peak = exp(-((depth - range).^2) / (2*sigma^2));
            profile = base + (1-base)*peak;
            profile(depth > (range + sigma*2)) = 0;
        case 'Alpha'
            profile = double(depth <= 0.5);
    end
    profile = profile / max(profile); % Normalize to 1.0 peak
end

function rbe = getRBE(type)
    switch type
        case 'Photon', rbe = 1.0;
        case 'Proton', rbe = 1.1;
        case 'Alpha', rbe = 20.0;
    end
end

function [mask, density] = generateColony(numCells)
    sz = 200; mask = zeros(sz);
    pts = [sz/2 + randn(numCells,1)*25, sz/2 + randn(numCells,1)*25];
    [x, y] = meshgrid(1:sz, 1:sz);
    for i = 1:numCells
        mask = mask + (sqrt((x-pts(i,1)).^2 + (y-pts(i,2)).^2) < 5);
    end
    density = imgaussfilt(mask, 8); density = density/max(density(:));
    mask = mask > 0;
end

function mask = applyBioResponse(mask, density, effectiveDose, props)
    if effectiveDose <= 0; return; end
    shielding = props.attenuation * density;
    localDose = effectiveDose .* (1 - shielding);
    survivalProb = exp(-props.alpha.*localDose - props.beta.*localDose.^2);
    mask = mask & (rand(size(mask)) < survivalProb);
    mask = bwareaopen(mask, 2); 
end

function plotCellNetwork(mask, threshold, plotTitle)
    imshow(~mask); hold on; title(plotTitle);
    s = regionprops(mask, 'Centroid');
    if isempty(s), return; end
    nodes = cat(1, s.Centroid);
    plot(nodes(:,1), nodes(:,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    dists = pdist2(nodes, nodes);
    numNodes = size(nodes, 1);
    for i = 1:numNodes
        for j = i+1:numNodes
            if dists(i,j) < threshold
                line([nodes(i,1), nodes(j,1)], [nodes(i,2), nodes(j,2)], ...
                     'Color', [0.2 0.2 0.8 0.5], 'LineWidth', 1);
            end
        end
    end
end

function A = buildAdjacency(mask, threshold)
    s = regionprops(mask, 'Centroid');
    if numel(s) < 3
        A = [];
        return;
    end
    
    nodes = cat(1, s.Centroid);
    D = pdist2(nodes, nodes);
    
    A = (D < threshold) & (D > 0);
    A = double(A);
end

function metrics = computeNetworkMetrics(A)

    metrics.C = NaN;
    metrics.L = NaN;
    metrics.sigma = NaN;

    if isempty(A) || size(A,1) < 5
        return;
    end

    % --- Extract giant component ---
    G = graph(A);
    bins = conncomp(G);
    giant = bins == mode(bins);

    if sum(giant) < 5
        return;
    end

    A = A(giant, giant);

    % --- Clustering ---
    C = mean(clustering_coef_bu(A));

    % --- Path length ---
    D = distance_bin(A);
    finiteD = D(~isinf(D) & D > 0);
    L = mean(finiteD);

    % --- Random graph (FAST version) ---
    p = nnz(A) / (numel(A)-size(A,1));
    Arand = rand(size(A)) < p;
    Arand = triu(Arand,1);
    Arand = Arand + Arand';

    Crand = mean(clustering_coef_bu(Arand));
    Drand = distance_bin(Arand);
    finiteDr = Drand(~isinf(Drand) & Drand > 0);
    Lrand = mean(finiteDr);

    % --- Small-world coefficient ---
    if L > 0 && Lrand > 0 && Crand > 0
        metrics.sigma = (C/Crand) / (L/Lrand);
    end

    metrics.C = C;
    metrics.L = L;
end

function metrics = computeNodeAndShapeMetrics(mask)

    metrics.meanEccentricity = NaN;
    metrics.meanArea        = NaN;
    metrics.meanPerimeter   = NaN;
    metrics.meanRoundness   = NaN;

    cc = bwconncomp(mask, 8);

    if cc.NumObjects < 2
        return;
    end

    % --- Shape metrics ---
    props = regionprops(cc, 'Area', 'Perimeter', 'Centroid');

    areas = [props.Area];
    perims = [props.Perimeter];

    roundness = 4*pi*areas ./ (perims.^2 + eps);

    metrics.meanArea      = mean(areas);
    metrics.meanPerimeter = mean(perims);
    metrics.meanRoundness = mean(roundness);

    % --- Build graph for eccentricity ---
    nodes = cat(1, props.Centroid);
    D = pdist2(nodes, nodes);

    threshold = 20;   % same as netThreshold
    A = (D < threshold) & (D > 0);

    G = graph(A);
    bins = conncomp(G);
    giant = bins == mode(bins);

    if sum(giant) < 2
        return;
    end

    Gg = subgraph(G, find(giant));

    distMat = distances(Gg);
    ecc = max(distMat, [], 2, 'omitnan');

    metrics.meanEccentricity = mean(ecc(~isinf(ecc)));
end


