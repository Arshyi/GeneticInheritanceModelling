import numpy as np
import matplotlib.pyplot as plt


def sickle_cell_model(dad_genotype, mom_genotype, n, plot=True):
    """
    Simulate sickle-cell genotype frequencies across n generations.

    Genotype labels accepted:
        'HbAHbA'  -> normal
        'HbAHbS'  -> carrier
        'HbSHbS'  -> sickle-cell disease

    Assumptions mirrored from the MATLAB writeup:
    - Generation 1 is computed directly from the parents using Punnett logic.
    - From generation 2 onward, a fixed 3x3 transition matrix is used.
    - HbSHbS individuals are counted each generation but do NOT reproduce into the next generation.

    Parameters
    ----------
    dad_genotype : str
    mom_genotype : str
    n : int
        Number of generations to simulate (must be >= 1).
    plot : bool
        Whether to display the plot.

    Returns
    -------
    history : np.ndarray, shape (3, n)
        Rows correspond to:
            0 -> HbAHbA
            1 -> HbAHbS
            2 -> HbSHbS
        Columns correspond to generations 1..n
    """

    if n < 1:
        raise ValueError("n must be at least 1.")

    valid = {"HbAHbA", "HbAHbS", "HbSHbS"}
    if dad_genotype not in valid or mom_genotype not in valid:
        raise ValueError("Genotypes must be one of: 'HbAHbA', 'HbAHbS', 'HbSHbS'.")

    # --- 1) Initial distribution from Punnett-square logic ---
    pair = tuple(sorted([dad_genotype, mom_genotype]))

    # X = [P(HbAHbA), P(HbAHbS), P(HbSHbS)] for generation 1
    punnett_map = {
        ("HbAHbA", "HbAHbA"): np.array([1.0, 0.0, 0.0]),
        ("HbAHbA", "HbAHbS"): np.array([0.5, 0.5, 0.0]),
        ("HbAHbA", "HbSHbS"): np.array([0.0, 1.0, 0.0]),
        ("HbAHbS", "HbAHbS"): np.array([0.25, 0.5, 0.25]),
        ("HbAHbS", "HbSHbS"): np.array([0.0, 0.5, 0.5]),
        ("HbSHbS", "HbSHbS"): np.array([0.0, 0.0, 1.0]),
    }

    X = punnett_map[pair].copy()

    # --- 2) Transition matrix from the essay / MATLAB code ---
    # columns correspond to parent-pair classes:
    #   AAxAA, AAxAS, ASxAS
    # rows correspond to offspring genotype distribution:
    #   AA, AS, SS
    M = np.array([
        [1.0, 0.5, 0.25],
        [0.0, 0.5, 0.5],
        [0.0, 0.0, 0.25]
    ])

    # --- 3) Simulate up to generation n ---
    history = np.zeros((3, n), dtype=float)
    history[:, 0] = X

    for gen in range(1, n):
        X_next = M @ X              # raw next-generation distribution
        history[:, gen] = X_next    # record before removing reproduction
        X_next[2] = 0.0             # HbSHbS do not reproduce
        X = X_next                  # carry forward

    # --- 4) Display results for generation n ---
    print(f"Generation {n}:")
    print(f"  HbA/HbA: {history[0, n - 1] * 100:.3f}%")
    print(f"  HbA/HbS: {history[1, n - 1] * 100:.3f}%")
    print(f"  HbS/HbS: {history[2, n - 1] * 100:.3f}%")

    # --- 5) Plot trajectories and highlight generation n ---
    if plot:
        gens = np.arange(1, n + 1)

        plt.figure(figsize=(9, 5))
        plt.plot(gens, history[0, :] * 100, label="HbA/HbA", linewidth=2)
        plt.plot(gens, history[1, :] * 100, label="HbA/HbS", linewidth=2)
        plt.plot(gens, history[2, :] * 100, label="HbS/HbS", linewidth=2)

        plt.scatter(n, history[0, n - 1] * 100, s=70)
        plt.scatter(n, history[1, n - 1] * 100, s=70)
        plt.scatter(n, history[2, n - 1] * 100, s=70)

        plt.xlabel("Generation")
        plt.ylabel("Genotype Probability (%)")
        plt.title(f"Genotype Distribution up to Generation {n}")
        plt.legend(loc="best")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    return history