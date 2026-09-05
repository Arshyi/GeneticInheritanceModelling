function bloodgroup18_simulation()
% bloodgroup18_simulation.m
% Simple ABO+Rh simulator via self-cross transition matrix

  % 1) Define genotype classes
  ABO = {'AA','AO','BB','BO','AB','OO'};
  Rh  = {'DD','Dd','dd'};
  classes = cell(18,1);
  k = 0;
  for i = 1:6
    for j = 1:3
      k = k + 1;
      classes{k} = [ABO{i} '/' Rh{j}];
    end
  end

  % 2) User input
  dad = input('Father genotype (e.g. AO/Dd): ','s');
  mom = input('Mother genotype (e.g. AO/Dd): ','s');
  n   = input('Number of generations: ');

  % 3) Compute generation 1
  v = zeros(18,n);
  v(:,1) = punnett18(dad,mom,ABO,Rh);

  % 4) Iterate generations via random mating
  for g = 2:n
    Xprev = v(:,g-1);
    Xnew  = zeros(18,1);
    for i = 1:18
      for j = 1:18
        w = Xprev(i) * Xprev(j);
        if w > 0
          Xnew = Xnew + w * punnett18(classes{i}, classes{j}, ABO, Rh);
        end
      end
    end
    v(:,g) = Xnew;
  end

  % 6) Display Gen-n genotype frequencies
  fprintf('Generation %d genotype frequencies:\n', n);
  for c = 1:18
    fprintf('  %5s: %5.2f%%\n', classes{c}, v(c,n)*100);
  end

  % 7) Compute phenotypes (8 categories)
  phen = zeros(8,n);
  labels = {'A+','A-','B+','B-','AB+','AB-','O+','O-'};
  for g = 1:n
    for c = 1:18
      parts = strsplit(classes{c}, '/');
      abo = parts{1}; rhg = parts{2};
      switch abo
        case {'AA','AO'}, base = 1;
        case {'BB','BO'}, base = 3;
        case 'AB',       base = 5;
        case 'OO',       base = 7;
      end
      idx = base + strcmp(rhg,'dd');
      phen(idx,g) = phen(idx,g) + v(c,g);
    end
  end

  % 8) Display Gen-n phenotype frequencies
  fprintf('Generation %d phenotype frequencies:\n', n);
  for p = 1:8
    fprintf('  %3s: %5.2f%%\n', labels{p}, phen(p,n)*100);
  end

  % 9) Plot genotype trajectories
  gens = 1:n;
  figure; hold on;
  for c = 1:18
    plot(gens, v(c,:)*100, 'DisplayName', classes{c});
  end
  legend('Location','eastoutside'); xlabel('Gen'); ylabel('Percent');
  title('Genotype Frequencies over Generations');

  % 10) Plot phenotype trajectories
  figure; hold on;
  for p = 1:8
    plot(gens, phen(p,:)*100, 'DisplayName', labels{p});
  end
  legend('Location','eastoutside'); xlabel('Gen'); ylabel('Percent');
  title('Phenotype Frequencies over Generations');
end

function p = punnett18(g1,g2,ABO,Rh)
  % Split ABO/Rh
  parts1 = strsplit(g1,'/'); a1 = parts1{1}; r1 = parts1{2};
  parts2 = strsplit(g2,'/'); a2 = parts2{1}; r2 = parts2{2};
  % ABO cross
  pA = punnett_simple(a1,a2,ABO);
  % Rh cross
  pR = punnett_rh(r1,r2,Rh);
  % Joint dist
  p = kron(pA,pR);
end

function p = punnett_simple(a1,a2,ABO)
  p = zeros(6,1);
  for i = 1:2
    for j = 1:2
      pair = sort([a1(i) a2(j)]);
      s = char(pair);
      if strcmp(s,'OA'), s='AO'; end
      if strcmp(s,'OB'), s='BO'; end
      idx = find(strcmp(ABO,s));
      p(idx) = p(idx) + 0.25;
    end
  end
end

function p = punnett_rh(r1,r2,Rh)
  p = zeros(3,1);
  for i = 1:2
    for j = 1:2
      pair = [r1(i) r2(j)];
      s = char(pair);
      if strcmp(s,'dD'), s='Dd'; end
      idx = find(strcmp(Rh,s));
      p(idx) = p(idx) + 0.25;
    end
  end
end
