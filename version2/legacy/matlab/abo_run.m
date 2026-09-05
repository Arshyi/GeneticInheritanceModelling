%% ABO-ONLY: Prompt → Simulate → Plot (Genotype + Phenotype)
% Run this file: abo_run
% Accepts parent inputs like: AA, AO, OO, BO, BB, AB (will ignore any trailing Rh like 'AAdd')

function abo_run
    fprintf('=== ABO ONLY (no Rh) ===\n');
    p1 = strtrim(input('Parent 1 ABO genotype (AA/AO/OO/BO/BB/AB or AAdd etc.): ','s'));
    p2 = strtrim(input('Parent 2 ABO genotype (AA/AO/OO/BO/BB/AB or AAdd etc.): ','s'));
    G  = input('Number of generations to simulate [default 12]: ');
    if isempty(G), G = 12; end

    out = abo_pipeline_from_parents(p1, p2, G);
    plot_abo_results(out, p1, p2);
end

%% ===== Pipeline (ABO only) =====
function out = abo_pipeline_from_parents(p1, p2, G)
    abo1 = extract_abo_from_parent(p1);
    abo2 = extract_abo_from_parent(p2);

    genoABO0 = abo_offspring_from_parents_str_aboOnly(abo1, abo2);  % 1x6
    trajABO = abo_sim(genoABO0, G);
    trajABOpheno = abo_geno_to_pheno_series(trajABO);

    out = struct('trajABO',trajABO,'trajABOpheno',trajABOpheno);
end

%% ===== Core ABO simulation (random mating) =====
function traj = abo_sim(geno0, G)
    geno0 = geno0(:).'/sum(geno0);
    traj = zeros(G+1, 6);
    traj(1,:) = geno0;
    for t = 1:G
        traj(t+1,:) = abo_step_random_mating(traj(t,:));
    end
end

function geno_next = abo_step_random_mating(geno)
    geno = geno(:).'; 
    geno = geno / sum(geno);
    P = geno.' * geno;  % 6x6 ordered pair probabilities

    % Gamete allele distributions for each ABO genotype, columns=[A B O]
    GAM = [ 0,   0,   1  ;  % OO
            0.5, 0,   0.5;  % AO
            1,   0,   0  ;  % AA
            0,   0.5, 0.5;  % BO
            0,   1,   0  ;  % BB
            0.5, 0.5, 0 ];  % AB

    geno_next = zeros(1,6);
    for i = 1:6
        gi = GAM(i,:);
        for j = 1:6
            pij = P(i,j);
            if pij==0, continue; end
            gj = GAM(j,:);
            off = abo_offspring_from_gametes(gi, gj); % 1x6
            geno_next = geno_next + pij * off;
        end
    end
    geno_next = geno_next / sum(geno_next);
end

function off = abo_offspring_from_gametes(g1, g2)
    pA1=g1(1); pB1=g1(2); pO1=g1(3);
    pA2=g2(1); pB2=g2(2); pO2=g2(3);

    OO = pO1*pO2;
    AO = pA1*pO2 + pO1*pA2;
    AA = pA1*pA2;
    BO = pB1*pO2 + pO1*pB2;
    BB = pB1*pB2;
    AB = pA1*pB2 + pB1*pA2;

    off = [OO AO AA BO BB AB];
    off = off / sum(off);
end

%% ===== Phenotype mapping =====
function P = abo_geno_to_pheno_series(traj6)
    OO=traj6(:,1); AO=traj6(:,2); AA=traj6(:,3);
    BO=traj6(:,4); BB=traj6(:,5); AB=traj6(:,6);

    O   = OO;
    A   = AO + AA;
    B   = BO + BB;
    ABp = AB;

    P = [O A B ABp];
end

%% ===== Parental parsing (ABO only) =====
function aboStr = extract_abo_from_parent(parentStr)
    s = upper(strtrim(parentStr));
    dpos = regexp(s,'D','once');
    if ~isempty(dpos)
        s = s(1:dpos-1); % strip any Rh suffix like DD/Dd/dd
    end
    if isequal(sort(s), ['A','O']), aboStr='AO';
    elseif isequal(sort(s), ['B','O']), aboStr='BO';
    elseif isequal(sort(s), ['A','B']), aboStr='AB';
    else
        aboStr = s; % 'AA','BB','OO' OK
    end
end

function genoABO = abo_offspring_from_parents_str_aboOnly(abo1, abo2)
    g1 = abo_gamete_from_genostr(abo1);
    g2 = abo_gamete_from_genostr(abo2);
    genoABO = abo_offspring_from_gametes(g1, g2);
end

function g = abo_gamete_from_genostr(gs)
    switch upper(gs)
        case 'AA', g=[1 0 0];
        case 'AO', g=[0.5 0 0.5];
        case 'OO', g=[0 0 1];
        case 'BO', g=[0 0.5 0.5];
        case 'BB', g=[0 1 0];
        case 'AB', g=[0.5 0.5 0];
        otherwise, error('Bad ABO genotype: %s', gs);
    end
end

%% ===== Pretty plotting (Genotypes + Phenotypes) =====
function plot_abo_results(out, p1, p2)
    G = size(out.trajABO,1)-1;
    gens = 0:G;

    % Fixed order + labels
    genoLabels = {'OO','AO','AA','BO','BB','AB'};
    phenoLabels= {'O','A','B','AB'};

    % Color coding (consistent palette)
    genoColors = [ ...
        0.20 0.20 0.70;   % OO
        0.10 0.60 0.80;   % AO
        0.90 0.30 0.25;   % AA
        0.90 0.60 0.10;   % BO
        0.20 0.65 0.20;   % BB
        0.55 0.35 0.70];  % AB

    phenoColors = [ ...
        0.20 0.20 0.70;   % O
        0.90 0.30 0.25;   % A
        0.20 0.65 0.20;   % B
        0.55 0.35 0.70];  % AB

    % Tiled layout: Genotypes | Phenotypes
    t = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % --- Genotypes
    nexttile; hold on;
    set(gca,'ColorOrder',genoColors,'NextPlot','replacechildren');
    plot(gens, out.trajABO, 'LineWidth',1.8);
    xlabel('Generation'); ylabel('Frequency');
    title({'ABO Genotypes',sprintf('Parents: %s × %s | Random mating, no selection', upper(p1), upper(p2))});
    legend(genoLabels,'Location','bestoutside'); grid on; ylim([0 1]);

    % --- Phenotypes
    nexttile; hold on;
    set(gca,'ColorOrder',phenoColors,'NextPlot','replacechildren');
    plot(gens, out.trajABOpheno, 'LineWidth',1.8);
    xlabel('Generation'); ylabel('Frequency');
    title({'ABO Phenotypes', 'Phenotypes: O, A, B, AB'});
    legend(phenoLabels,'Location','bestoutside'); grid on; ylim([0 1]);

    % "Subtitle" across the figure
    sgtitle('ABO Population Evolution (Random Mating, No Selection)');
end
