function abo_simulation()
% abo_simulation.m
% Simulate ABO genotype & phenotype frequencies for n generations
% - Gen 1 from user-specified parents via Punnett
% - Gen 2..n via random mating among all six genotypes

  clc; clear;

  % 1) User input
  valid = {'AA','AO','BB','BO','AB','OO'};
  dad = upper(strtrim(input('Enter father genotype (AA, AO, BB, BO, AB, OO): ','s')));
  mom = upper(strtrim(input('Enter mother genotype (AA, AO, BB, BO, AB, OO): ','s')));
  if ~ismember(dad,valid) || ~ismember(mom,valid)
      error('Genotype must be one of: AA, AO, BB, BO, AB, OO');
  end
  n = input('Enter number of generations to simulate (integer ≥1): ');
  if n<1 || floor(n)~=n, error('n must be a positive integer'); end

  % 3) Prepare
  types = {'AA','AO','BB','BO','AB','OO'};
  history = zeros(6,n);

  % Gen 1
  history(:,1) = punnett(dad,mom);

  % 4) Random‐mating update for gens 2..n
  for g=2:n
    Xprev = history(:,g-1);
    Xnew = zeros(6,1);
    for i=1:6
      for j=1:6
        p_pair = Xprev(i)*Xprev(j);
        Xnew = Xnew + p_pair*punnett(types{i}, types{j});
      end
    end
    history(:,g) = Xnew;
  end

  % 5) Display Gen-n results
  fprintf('\nGen %d genotype frequencies:\n',n);
  for t=1:6
    fprintf('  %-2s: %5.2f%%\n', types{t}, history(t,n)*100);
  end

  phen = zeros(4,n);
  % A = AA+AO, B = BB+BO, AB = AB, O = OO
  phen(1,:) = history(1,:)+history(2,:);
  phen(2,:) = history(3,:)+history(4,:);
  phen(3,:) = history(5,:);
  phen(4,:) = history(6,:);

  fprintf('\nGen %d phenotype frequencies:\n',n);
  labels = {'A','B','AB','O'};
  for p=1:4
    fprintf('  %-2s: %5.2f%%\n', labels{p}, phen(p,n)*100);
  end

  % 6) Plot over time
  gens = 1:n;
  figure; hold on; box on;
  colors = ['r','b','g','m','c','k'];
  for t=1:6
    plot(gens, history(t,:)*100, ['-' colors(t)], 'LineWidth',1.5);
  end
  legend(types,'Location','Best');
  xlabel('Generation'); ylabel('Genotype %');
  title('ABO Genotype Trajectories');
  grid on;

  figure; hold on; box on;
  phcols = ['r','b','g','k'];
  for p=1:4
    plot(gens, phen(p,:)*100, ['-' phcols(p)], 'LineWidth',1.5);
  end
  legend(labels,'Location','Best');
  xlabel('Generation'); ylabel('Phenotype %');
  title('ABO Phenotype Trajectories');
  grid on;
end

% -------------------------------------------------------------------------
% Local function must come after the main function in a function file.
function v = punnett(g1,g2)
  % Returns a 6×1 vector of offspring genotype probs for parents g1×g2.
  % g1, g2 are strings: 'AA','AO','BB','BO','AB','OO'
  kids = strings(4,1);
  idx = 0;
  for i=1:2
    for j=1:2
      idx = idx+1;
      a = g1(i);  b = g2(j);
      pair = sort([a b]);  % e.g. 'O'+'A' -> ['A','O']
      kids(idx) = pair;
    end
  end
  types = ["AA","AO","BB","BO","AB","OO"];
  v = zeros(6,1);
  for k=1:4
    str = kids(k);
    if str=="OA", str="AO"; end
    if str=="OB", str="BO"; end
    for t=1:6
      if types(t)==str
        v(t) = v(t) + 1/4;
      end
    end
  end
end
