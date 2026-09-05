% sickle_cell_model.m
%
% Assumptions:
% 1. Alleles: HbA, HbS → genotypes HbA/HbA, HbA/HbS, HbS/HbS.
% 2. Canada (no malaria advantage).
% 3. HbS/HbS individuals survive but do NOT reproduce.
% 4. Parents must be one of:
%       'HbAHbA', 'HbAHbS', 'HbSHbA', 'HbSHbS'
% 5. Gen-1 via Punnett; gens 2…n via a 3×3 transition matrix.

clc; clear;

%% 1. Prompt for inputs
dad_genotype = input(...
  'Enter father genotype (HbAHbA, HbAHbS, HbSHbA, or HbSHbS): ', 's');
mom_genotype = input(...
  'Enter mother genotype (HbAHbA, HbAHbS, HbSHbA, or HbSHbS): ', 's');
n = input('Enter number of generations to simulate (positive integer): ');

%% 2. Generation-1 via Punnett logic
pair = sort({dad_genotype, mom_genotype});
key  = strjoin(pair, '_');
switch key
  case 'HbAHbA_HbAHbA'
    X = [1; 0; 0];
  case {'HbAHbA_HbAHbS','HbAHbA_HbSHbA'}
    X = [0.5; 0.5; 0];
  case 'HbAHbA_HbSHbS'
    X = [0; 1; 0];
  case {'HbAHbS_HbAHbS','HbAHbS_HbSHbA','HbSHbA_HbSHbA'}
    X = [0.25; 0.5; 0.25];
  case {'HbAHbS_HbSHbS','HbSHbA_HbSHbS'}
    X = [0; 0.5; 0.5];
  case 'HbSHbS_HbSHbS'
    X = [0; 0; 1];
  otherwise
    error('Invalid genotype. Use HbAHbA, HbAHbS, HbSHbA, or HbSHbS');
end

%% 3. Transition matrix
M = [1,   1/2, 1/4;
     0,   1/2, 1/2;
     0,   0,   1/4];

%% 4. Simulate up to generation n
history = zeros(3, n);
history(:,1) = X;
for gen = 2:n
  X_next = M * X;
  history(:,gen) = X_next;   % record raw
  X_next(3) = 0;             % HbS/HbS do not reproduce
  X = X_next;                % carry forward
end

%% 5. Display final generation
fprintf('\nResults at generation %d:\n', n);
fprintf('  HbA/HbA: %.3f%%\n', history(1,n)*100);
fprintf('  HbA/HbS: %.3f%%\n', history(2,n)*100);
fprintf('  HbS/HbS: %.3f%%\n\n', history(3,n)*100);

%% 6. Plot and highlight generation n
gens = 1:n;
figure; hold on; box on;
plot(gens, history(1,:)*100, '-g','LineWidth',2);
plot(gens, history(2,:)*100, '-b','LineWidth',2);
plot(gens, history(3,:)*100, '-r','LineWidth',2);
scatter(n, history(1,n)*100, 80,'g','filled');
scatter(n, history(2,n)*100, 80,'b','filled');
scatter(n, history(3,n)*100, 80,'r','filled');
xlabel('Generation');
ylabel('Genotype Probability (%)');
title(sprintf('Genotype Distribution up to Generation %d', n));
legend('HbA/HbA','HbA/HbS','HbS/HbS','Location','Best');
grid on;
