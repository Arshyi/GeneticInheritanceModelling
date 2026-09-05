%% ABO + Rh: Prompt → Simulate → Plot (Genotype + Phenotype)
% Run this file: aborh_run
% Accepts parents like: AAdd, AOdd, ABdd, BOdd, BBdd, OOdd, AODd, ABDD, etc.

function aborh_run
    fprintf('=== ABO + Rh (unlinked, no selection) ===\n');
    p1 = strtrim(input('Parent 1 (e.g., AAdd, AOdd, ABDD): ','s'));
    p2 = strtrim(input('Parent 2 (e.g., AAdd, AOdd, ABDD): ','s'));
    G  = input('Number of generations to simulate [default 12]: ');
    if isempty(G), G = 12; end

    out = aborh_pipeline_from_parents(p1, p2, G);
    plot_aborh_results(out, p1, p2);
end

%% ===== Pipeline (ABO + Rh) =====
function out = aborh_pipeline_from_parents(p1, p2, G)
    [abo1, rh1] = parse_parent_str_aborh(p1);
    [abo2, rh2] = parse_parent_str_aborh(p2);

    genoABO0 = abo_offspring_from_parents_str(abo1, abo2); % 1x6
    genoRh0  = rh_offspring_from_parents_str(rh1, rh2);    % 1x3

    trajABO = abo_sim(genoABO0, G);
    trajRh  = rh_sim(genoRh0,  G);

    trajABORh = zeros(G+1, 18);
    for t = 1:G+1
        trajABORh(t,:) = kron(trajABO(t,:), trajRh(t,:)); % ABO-major ⊗ Rh-minor
    end

    trajABOpheno    = abo_geno_to_pheno_series(trajABO);
    trajRhpheno     = rh_geno_to_pheno_series(trajRh);
    trajABORhpheno  = aborh_geno_to_pheno_series(trajABORh);

    out = struct('trajABO',trajABO,'trajRh',trajRh,'trajABORh',trajABORh,...
                 'trajABOpheno',trajABOpheno,'trajRhpheno',trajRhpheno,...
                 'trajABORhpheno',trajABORhpheno);
end

%% ===== Parent parsing (ABO + Rh) =====
function [aboStr, rhStr] = parse_parent_str_aborh(parentStr)
    s = upper(strtrim(parentStr));
    dpos = regexp(s,'D');
    if isempty(dpos)
        if length(s) < 2
            error('Parent must include ABO (e.g., AO) and optionally Rh (e.g., dd).');
        end
        aboRaw = s(1:end-2); rhStr = 'dd'; % default dd if Rh missing
    else
        aboRaw = s(1:dpos(1)-1);
        rhRaw  = s(dpos(1):end);
        if strcmpi(rhRaw,'DD'), rhStr='DD';
        elseif strcmpi(rhRaw,'Dd'), rhStr='Dd';
        elseif strcmpi(rhRaw,'dd'), rhStr='dd';
        else, error('Bad Rh string: %s', rhRaw);
        end
    end
    if isequal(sort(aboRaw), ['A','O']), aboStr='AO';
    elseif isequal(sort(aboRaw), ['B','O']), aboStr='BO';
    elseif isequal(sort(aboRaw), ['A','B']), aboStr='AB';
    else
        aboStr = aboRaw; % 'AA','BB','OO'
    end
end

%% ===== ABO core =====
function traj = abo_sim(geno0, G)
    geno0 = geno0(:).'/sum(geno0);
    traj = zeros(G+1, 6); traj(1,:) = geno0;
    for t = 1:G, traj(t+1,:) = abo_step_random_mating(traj(t,:)); end
end

function geno_next = abo_step_random_mating(geno)
    geno = geno(:).'; geno = geno / sum(geno);
    P = geno.' * geno;
    GAM = [ 0,   0,   1  ;
            0.5, 0,   0.5;
            1,   0,   0  ;
            0,   0.5, 0.5;
            0,   1,   0  ;
            0.5, 0.5, 0 ];
    geno_next = zeros(1,6);
    for i = 1:6
        gi = GAM(i,:);
        for j = 1:6
            pij = P(i,j); if pij==0, continue; end
            gj = GAM(j,:);
            off = abo_offspring_from_gametes(gi, gj);
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
    off = [OO AO AA BO BB AB]; off = off / sum(off);
end

function genoABO = abo_offspring_from_parents_str(abo1, abo2)
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

function P = abo_geno_to_pheno_series(traj6)
    OO=traj6(:,1); AO=traj6(:,2); AA=traj6(:,3);
    BO=traj6(:,4); BB=traj6(:,5); AB=traj6(:,6);
    O   = OO; A = AO + AA; B = BO + BB; ABp = AB;
    P = [O A B ABp];
end

%% ===== Rh core =====
function traj = rh_sim(geno0, G)
    geno0 = geno0(:).'/sum(geno0);
    traj = zeros(G+1,3); traj(1,:) = geno0;
    for t = 1:G, traj(t+1,:) = rh_step_random_mating(traj(t,:)); end
end

function geno_next = rh_step_random_mating(geno)
    geno = geno(:).'; geno = geno / sum(geno);
    P = geno.' * geno; % 3x3
    GAM = [ 0,   1 ;   % dd -> all d
            0.5, 0.5; % Dd
            1,   0 ]; % DD
    geno_next = zeros(1,3);
    for i = 1:3
        gi = GAM(i,:);
        for j = 1:3
            pij = P(i,j); if pij==0, continue; end
            gj = GAM(j,:);
            off = rh_offspring_from_gametes(gi, gj); % [dd Dd DD]
            geno_next = geno_next + pij * off;
        end
    end
    geno_next = geno_next / sum(geno_next);
end

function off = rh_offspring_from_gametes(g1, g2)
    pD1=g1(1); pd1=g1(2); pD2=g2(1); pd2=g2(2);
    DD = pD1*pD2;
    dd = pd1*pd2;
    Dd = pD1*pd2 + pd1*pD2;
    off = [dd Dd DD]; off = off / sum(off);
end

function genoRh = rh_offspring_from_parents_str(rh1, rh2)
    g1 = rh_gamete_from_genostr(rh1); 
    g2 = rh_gamete_from_genostr(rh2);
    genoRh = rh_offspring_from_gametes(g1, g2);
end

function g = rh_gamete_from_genostr(gs)
    if     strcmpi(gs,'DD'), g=[1 0];
    elseif strcmpi(gs,'Dd'), g=[0.5 0.5];
    elseif strcmpi(gs,'dd'), g=[0 1];
    else, error('Bad Rh genotype: %s', gs);
    end
end

function P = rh_geno_to_pheno_series(traj3)
    dd=traj3(:,1); Dd=traj3(:,2); DD=traj3(:,3);
    P = [dd, Dd+DD]; % [Rh- Rh+]
end

function P8 = aborh_geno_to_pheno_series(traj18)
    T = size(traj18,1); P8 = zeros(T,8);
    for t=1:T
        row = traj18(t,:);
        row = reshape(row, [3,6]).'; % 6x3 blocks: [dd Dd DD]
        Oblock  = row(1,:);
        Ablock  = row(2,:) + row(3,:);
        Bblock  = row(4,:) + row(5,:);
        ABblock = row(6,:);
        P8(t,:) = [Oblock(1), sum(Oblock(2:3)), ...
                   Ablock(1), sum(Ablock(2:3)), ...
                   Bblock(1), sum(Bblock(2:3)), ...
                   ABblock(1),sum(ABblock(2:3))];
    end
end

%% ===== Pretty plotting (Genotypes + Phenotypes) =====
function plot_aborh_results(out, p1, p2)
    G = size(out.trajABO,1)-1;
    gens = 0:G;

    % Labels / fixed order
    aboGenoLabels = {'OO','AO','AA','BO','BB','AB'};
    aboPhenoLabels= {'O','A','B','AB'};
    rhGenoLabels  = {'dd','Dd','DD'};
    rhPhenoLabels = {'Rh-','Rh+'};
    bothPhenoLabels = {'O-','O+','A-','A+','B-','B+','AB-','AB+'};

    % Color palettes (consistent)
    aboGenoColors = [ ...
        0.20 0.20 0.70; 0.10 0.60 0.80; 0.90 0.30 0.25; ...
        0.90 0.60 0.10; 0.20 0.65 0.20; 0.55 0.35 0.70 ];
    aboPhenoColors = [ ...
        0.20 0.20 0.70; 0.90 0.30 0.25; 0.20 0.65 0.20; 0.55 0.35 0.70 ];
    rhGenoColors = [0.25 0.25 0.25; 0.55 0.55 0.55; 0.85 0.85 0.85];
    rhPhenoColors= [0.25 0.25 0.25; 0.70 0.70 0.70];
    bothPhenoColors = [ ...
        0.25 0.25 0.55; 0.45 0.45 0.85; ... % O- O+
        0.75 0.30 0.30; 0.95 0.55 0.55; ... % A- A+
        0.20 0.60 0.25; 0.50 0.85 0.55; ... % B- B+
        0.50 0.35 0.70; 0.75 0.60 0.90];    % AB- AB+

    % Figure 1: ABO Genotypes | ABO Phenotypes
    figure('Name','ABO Evolution','Color','w');
    t1 = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    set(gca,'ColorOrder',aboGenoColors,'NextPlot','replacechildren');
    plot(gens, out.trajABO, 'LineWidth',1.8);
    xlabel('Generation'); ylabel('Frequency'); ylim([0 1]); grid on;
    title({'ABO Genotypes',sprintf('Parents: %s × %s', upper(p1), upper(p2))});
    legend(aboGenoLabels,'Location','bestoutside');

    nexttile; hold on;
    set(gca,'ColorOrder',aboPhenoColors,'NextPlot','replacechildren');
    plot(gens, out.trajABOpheno, 'LineWidth',1.8);
    xlabel('Generation'); ylabel('Frequency'); ylim([0 1]); grid on;
    title({'ABO Phenotypes','O, A, B, AB'});
    legend(aboPhenoLabels,'Location','bestoutside');

    sgtitle('ABO Population Evolution (Random Mating, No Selection)');

    % Figure 2: Rh Genotypes | Rh Phenotypes
    figure('Name','Rh Evolution','Color','w');
    t2 = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    set(gca,'ColorOrder',rhGenoColors,'NextPlot','replacechildren');
    plot(gens, out.trajRh, 'LineWidth',1.8);
    xlabel('Generation'); ylabel('Frequency'); ylim([0 1]); grid on;
    title('Rh Genotypes'); legend(rhGenoLabels,'Location','bestoutside');

    nexttile; hold on;
    set(gca,'ColorOrder',rhPhenoColors,'NextPlot','replacechildren');
    plot(gens, out.trajRhpheno, 'LineWidth',1.8);
    xlabel('Generation'); ylabel('Frequency'); ylim([0 1]); grid on;
    title('Rh Phenotypes'); legend(rhPhenoLabels,'Location','bestoutside');

    sgtitle('Rh Population Evolution (Random Mating, No Selection)');

    % Figure 3: Combined ABO×Rh Phenotypes (8 series)
    figure('Name','ABO×Rh Phenotypes','Color','w');
    set(gca,'ColorOrder',bothPhenoColors,'NextPlot','replacechildren'); hold on;
    plot(gens, out.trajABORhpheno, 'LineWidth',1.8);
    xlabel('Generation'); ylabel('Frequency'); ylim([0 1]); grid on;
    title({'ABO × Rh Phenotypes','O-, O+, A-, A+, B-, B+, AB-, AB+'});
    legend(bothPhenoLabels,'Location','bestoutside');
end
