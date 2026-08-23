import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


# ============================================================
# Helpers
# ============================================================

def canonical_pair(a, b):
    """Return genotype pair in a consistent sorted form."""
    return tuple(sorted((a, b)))


def punnett_distribution(parent1, parent2):
    """
    Compute offspring genotype probabilities for one gene.

    parent1, parent2: tuples like ('A','O') or ('+','-')
    returns: dict {child_genotype_pair: probability}
    """
    outcomes = defaultdict(float)
    for a1 in parent1:
        for a2 in parent2:
            child = canonical_pair(a1, a2)
            outcomes[child] += 0.25
    return dict(outcomes)


# ============================================================
# ABO system
# ============================================================

ABO_ALLELES = ['A', 'B', 'O']
ABO_GENOTYPES = [
    ('A', 'A'),
    ('A', 'O'),
    ('A', 'B'),
    ('B', 'B'),
    ('B', 'O'),
    ('O', 'O')
]


def abo_phenotype(genotype):
    g = canonical_pair(*genotype)
    if g in [('A', 'A'), ('A', 'O')]:
        return 'A'
    if g in [('B', 'B'), ('B', 'O')]:
        return 'B'
    if g == ('A', 'B'):
        return 'AB'
    return 'O'


# ============================================================
# Rh system
# ============================================================

RH_ALLELES = ['+', '-']
RH_GENOTYPES = [
    ('+', '+'),
    ('+', '-'),
    ('-', '-')
]


def rh_phenotype(genotype):
    return '+' if '+' in genotype else '-'


# ============================================================
# Simplified eye-color model
# ============================================================
# IMPORTANT:
# This is only a toy Mendelian approximation.
# Real eye color is polygenic and influenced by multiple genes.

EYE_ALLELES = ['B', 'b']   # B = brown-dominant, b = blue-recessive
EYE_GENOTYPES = [
    ('B', 'B'),
    ('B', 'b'),
    ('b', 'b')
]


def eye_phenotype(genotype):
    return 'Brown' if 'B' in genotype else 'Blue'


# ============================================================
# Combine independent loci
# ============================================================

def combine_locus_distributions(*dists):
    """
    Combine multiple independent locus distributions.
    Each dist is a dict {genotype_tuple: prob}.
    Returns dict {combined_genotype_tuple: prob}
    """
    combined = {(): 1.0}

    for dist in dists:
        new_combined = defaultdict(float)
        for prefix, p1 in combined.items():
            for geno, p2 in dist.items():
                new_combined[prefix + (geno,)] += p1 * p2
        combined = dict(new_combined)

    return combined


# ============================================================
# ABO + Rh simulation for a fixed parental pairing pattern
# ============================================================
# This matches the style of your earlier sickle-cell project:
# generation 1 comes from the actual parent pair,
# then each generation is produced by mating within the same
# genotype-distribution pattern from the previous generation.

class MultiTraitInheritanceModel:
    def __init__(self, include_eye=False):
        self.include_eye = include_eye
        self.state_genotypes = self._build_state_space()
        self.index = {g: i for i, g in enumerate(self.state_genotypes)}

    def _build_state_space(self):
        states = []
        if self.include_eye:
            for abo in ABO_GENOTYPES:
                for rh in RH_GENOTYPES:
                    for eye in EYE_GENOTYPES:
                        states.append((abo, rh, eye))
        else:
            for abo in ABO_GENOTYPES:
                for rh in RH_GENOTYPES:
                    states.append((abo, rh))
        return states

    def initial_distribution_from_parents(
        self,
        parent1_abo, parent2_abo,
        parent1_rh, parent2_rh,
        parent1_eye=None, parent2_eye=None
    ):
        """
        Generation 1 offspring distribution from two parents.
        """
        abo_dist = punnett_distribution(parent1_abo, parent2_abo)
        rh_dist = punnett_distribution(parent1_rh, parent2_rh)

        if self.include_eye:
            if parent1_eye is None or parent2_eye is None:
                raise ValueError("Eye genotypes must be provided when include_eye=True.")
            eye_dist = punnett_distribution(parent1_eye, parent2_eye)
            combined = combine_locus_distributions(abo_dist, rh_dist, eye_dist)
        else:
            combined = combine_locus_distributions(abo_dist, rh_dist)

        X0 = np.zeros(len(self.state_genotypes), dtype=float)
        for genotype, prob in combined.items():
            X0[self.index[genotype]] = prob

        return X0

    def next_generation(self, current_distribution):
        """
        Produce the next generation by random mating within the current generation.

        This is more general than using a hand-written transition matrix:
        it builds offspring probabilities from every possible parent-pair,
        weighted by current genotype frequencies.
        """
        next_dist = np.zeros_like(current_distribution)

        for i, g1 in enumerate(self.state_genotypes):
            p1 = current_distribution[i]
            if p1 == 0:
                continue

            for j, g2 in enumerate(self.state_genotypes):
                p2 = current_distribution[j]
                if p2 == 0:
                    continue

                pair_weight = p1 * p2

                if self.include_eye:
                    abo1, rh1, eye1 = g1
                    abo2, rh2, eye2 = g2

                    abo_child = punnett_distribution(abo1, abo2)
                    rh_child = punnett_distribution(rh1, rh2)
                    eye_child = punnett_distribution(eye1, eye2)

                    child_dist = combine_locus_distributions(abo_child, rh_child, eye_child)
                else:
                    abo1, rh1 = g1
                    abo2, rh2 = g2

                    abo_child = punnett_distribution(abo1, abo2)
                    rh_child = punnett_distribution(rh1, rh2)

                    child_dist = combine_locus_distributions(abo_child, rh_child)

                for child_genotype, child_prob in child_dist.items():
                    next_dist[self.index[child_genotype]] += pair_weight * child_prob

        return next_dist

    def simulate(self, initial_distribution, generations):
        """
        Simulate genotype distributions across generations.

        Returns history array of shape (num_states, generations)
        """
        if generations < 1:
            raise ValueError("generations must be at least 1")

        history = np.zeros((len(self.state_genotypes), generations), dtype=float)
        X = initial_distribution.copy()

        for g in range(generations):
            history[:, g] = X
            X = self.next_generation(X)

        return history

    def phenotype_summary(self, history):
        """
        Aggregate genotype history into phenotype history.
        """
        generations = history.shape[1]
        summary = defaultdict(lambda: np.zeros(generations))

        for i, genotype in enumerate(self.state_genotypes):
            if self.include_eye:
                abo, rh, eye = genotype
                label = f"{abo_phenotype(abo)}{rh_phenotype(rh)} | {eye_phenotype(eye)}"
            else:
                abo, rh = genotype
                label = f"{abo_phenotype(abo)}{rh_phenotype(rh)}"

            summary[label] += history[i]

        return dict(summary)

    def plot_phenotypes(self, history, title="Phenotype Distribution Across Generations"):
        summary = self.phenotype_summary(history)

        plt.figure(figsize=(11, 6))
        for label, values in sorted(summary.items()):
            plt.plot(np.arange(1, len(values) + 1), values * 100, label=label)

        plt.xlabel("Generation")
        plt.ylabel("Frequency (%)")
        plt.title(title)
        plt.grid(True, alpha=0.3)
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.show()


# ============================================================
# Example usage
# ============================================================

if __name__ == "__main__":
    # -------- ABO + Rh only --------
    model = MultiTraitInheritanceModel(include_eye=False)

    # Example parents:
    # Parent 1: AO, +-
    # Parent 2: BO, --
    X0 = model.initial_distribution_from_parents(
        parent1_abo=('A', 'O'),
        parent2_abo=('B', 'O'),
        parent1_rh=('+', '-'),
        parent2_rh=('-', '-')
    )

    history = model.simulate(X0, generations=10)
    model.plot_phenotypes(history, title="ABO + Rh Across Generations")

    # -------- ABO + Rh + simplified eye color --------
    model_eye = MultiTraitInheritanceModel(include_eye=True)

    # Example:
    # Parent 1 eye = Bb
    # Parent 2 eye = bb
    X0_eye = model_eye.initial_distribution_from_parents(
        parent1_abo=('A', 'O'),
        parent2_abo=('B', 'O'),
        parent1_rh=('+', '-'),
        parent2_rh=('-', '-'),
        parent1_eye=('B', 'b'),
        parent2_eye=('b', 'b')
    )

    history_eye = model_eye.simulate(X0_eye, generations=10)
    model_eye.plot_phenotypes(
        history_eye,
        title="ABO + Rh + Simplified Eye Color Across Generations"
    )