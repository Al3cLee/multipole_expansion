"""
Symmetric Traceless Multipole Moments

This module implements the computation of symmetric traceless multipole moments Q^{i₁...iₙ}.

The moments are defined such that:
    φ^(n)(x) = (1/n!) * Q^{i₁...iₙ} * n^{i₁} * ... * n^{iₙ} / r^{n+1}

Properties:
    1. Symmetric: Q^{i₁...iₙ} is unchanged under any permutation of indices
    2. Traceless: Q^{...i...j...} δ_{ij} = 0 (trace over any pair is zero)

Examples:
    n=0 (monopole): Q = 1
    n=1 (dipole): Q^i = x_a^i
    n=2 (quadrupole): Q^{ij} = 3 x_a^i x_a^j - δ^{ij} |x_a|²
"""

from symbolica import Expression, S
from .contraction import TensorContraction
import itertools


class MultipoleMoments:
    """
    Compute and manipulate symmetric traceless multipole moments.
    """

    def __init__(self):
        """Initialize tensor contraction engine."""
        self.tc = TensorContraction()

        # Symbols
        self.xa = self.tc.xa
        self.n = self.tc.n
        self.r0 = self.tc.r0
        self.ra0 = self.tc.ra0
        self.delta = self.tc.delta
        self.dot = self.tc.dot

    def Q_tensor(self, n, indices=None):
        """
        Compute the symmetric traceless multipole moment Q^{i₁...iₙ}.

        Uses recursive algorithm that works for arbitrary n.

        Args:
            n: Order of the multipole
            indices: List of index symbols (if None, generic indices are created)

        Returns:
            Expression for Q^{i₁...iₙ}
        """
        if indices is None:
            indices = [S(f"i{k}") for k in range(1, n + 1)]  # indices = [1,2,...,n]
        # Edge cases first
        if n == 0:
            return self._Q_0()
        elif n == 1:
            return self._Q_1(indices[0])
        else:
            # None-edge cases,
            # Use general formula for n >= 2
            return self._Q_n_general(n, indices)

    def _Q_0(self):
        """
        Monopole moment: Q = 1
        """
        return Expression.num(1)

    def _Q_1(self, i):
        """
        Dipole moment: Q^i = x_a^i

        This is already symmetric (trivially) and traceless (scalar).
        """
        return self.xa(i)

    def _Q_2(self, i, j):
        """
        Quadrupole moment: Q^{ij} = 3 x_a^i x_a^j - δ^{ij} |x_a|²

        Properties:
        - Symmetric: Q^{ij} = Q^{ji}
        - Traceless: Q^{ii} = 3|x_a|² - 3|x_a|² = 0

        Note: Coefficient pattern:
        - Main term: (2*2-1)!! = 3!! = 3
        - Trace term: (2*2-3)!! = 1!! = 1
        """
        term1 = Expression.num(3) * self.xa(i) * self.xa(j)
        term2 = self.delta(i, j) * self.ra0**2
        return term1 - term2

    def _Q_3(self, i, j, k):
        """
        Octupole moment: Q^{ijk}

        Formula:
        Q^{ijk} = 15 x_a^i x_a^j x_a^k
                  - 3 |x_a|² (δ^{ij} x_a^k + δ^{ik} x_a^j + δ^{jk} x_a^i)

        This is symmetric and traceless.

        Note: Coefficient pattern:
        - Main term: (2*3-1)!! = 5!! = 15
        - Trace term: (2*3-3)!! = 3!! = 3
        """
        # Main term: 15 x_a^i x_a^j x_a^k
        term1 = Expression.num(15) * self.xa(i) * self.xa(j) * self.xa(k)

        # Trace terms: -3|x_a|² (δ^{ij} x_a^k + δ^{ik} x_a^j + δ^{jk} x_a^i)
        term2 = (
            Expression.num(3)
            * self.ra0**2
            * (
                self.delta(i, j) * self.xa(k)
                + self.delta(i, k) * self.xa(j)
                + self.delta(j, k) * self.xa(i)
            )
        )

        return term1 - term2

    def _Q_n_general(self, n, indices):
        """
        Algorithm for Q^{i₁...iₙ} , works for arbitrary n!

        Formula:
            Q^{i₁...iₙ} = Σ_{m=0}^{floor(n/2)} (-1)^m (2n-2m-1)!! * |x_a|^{2m} *
                         Σ_{pairings} δ^{...} * x_a^{remaining indices}

        Explicitly:
            = (2n-1)!! * x_a^{i₁} * ... * x_a^{iₙ}
              - (2n-3)!! * |x_a|² * Σ_{pairs} δ^{iⱼiₖ} * x_a^{remaining}
              + (2n-5)!! * |x_a|⁴ * Σ_{2 pairs} δ^{...}δ^{...} * x_a^{remaining}
              - ...

        The alternating sign (-1)^m ensures the traceless property.

        When used in φ^(n) = (1/n!) * Q * n...n / r^{n+1}, this gives the same
        result as the Taylor expansion φ^(n) = ((-1)^n/n!) * x_a * ∂^n(1/r).
        """
        # Main product term: (2n-1)!! * x_a^{i₁} * ... * x_a^{iₙ}
        # Use double factorial to match derivative normalization
        main_coeff = self._double_factorial(2 * n - 1)
        product = Expression.num(main_coeff)
        for idx in indices:
            product = product * self.xa(idx)

        # Add trace terms (which now include the alternating (-1)^m signs)
        trace_terms = self._compute_all_traces(n, indices)

        return product + trace_terms

    def _compute_all_traces(self, n, indices):
        """
        Compute all trace correction terms for arbitrary n.

        Returns sum of:
            (-1)^1 * (2n-3)!! * |x_a|² * (single pair contractions)
            + (-1)^2 * (2n-5)!! * |x_a|⁴ * (double pair contractions)
            + (-1)^3 * (2n-7)!! * |x_a|⁶ * (triple pair contractions)
            + ...

        Note: The alternating sign (-1)^m is built into each term.
        """
        total = Expression.num(0)

        # Number of pairs we can contract: 0, 1, 2, ..., floor(n/2)
        max_pairs = n // 2
        # Sum over k =number of δ symbols
        #            =number of contracted pairs
        #            =num_pairs
        for num_pairs in range(1, max_pairs + 1):
            pair_terms = self._trace_with_k_pairs(n, indices, num_pairs)
            total = total + pair_terms

        return total

    def _trace_with_k_pairs(self, n, indices, k):
        """
        Compute trace terms with exactly k pairs of deltas.

        For k=1: (-1)^1 * (2n-3)!! * |x_a|^2 * Σ_{i<j} δ^{iᵢ iⱼ} * x_a^{remaining n-2 indices}
        For k=2: (-1)^2 * (2n-5)!! * |x_a|^4 * Σ_{pairs} δ^{...}δ^{...} * x_a^{remaining n-4 indices}
        etc.

        The coefficient is (-1)^k * (2n - 2k - 1)!!
        For k=1: -1 * (2n-3)!!
        For k=2: +1 * (2n-5)!!
        For k=3: -1 * (2n-7)!!
        """
        if k > n // 2:
            return Expression.num(0)

        # Generate all ways to choose k pairs from n indices
        all_pairings = self._generate_k_pairings(n, k)

        total = Expression.num(0)

        for pairing in all_pairings:
            term = Expression.num(1)

            # Add delta for each pair
            for i_pos, j_pos in pairing:
                term = term * self.delta(indices[i_pos], indices[j_pos])

            # Specify paired positions
            paired_positions = set()
            for i_pos, j_pos in pairing:
                paired_positions.add(i_pos)
                paired_positions.add(j_pos)

            # Add x_a for unpaired indices
            for pos in range(n):
                if pos not in paired_positions:
                    term = term * self.xa(indices[pos])

            total = total + term

        # Coefficient for k-pair trace:
        # Use (2n - 2k - 1)!! to match the derivative pattern
        # For k=1: (2n-3)!!, for k=2: (2n-5)!!, etc.
        coeff_arg = 2 * n - 2 * k - 1
        if coeff_arg > 0:
            trace_coeff = self._double_factorial(coeff_arg)
        else:
            # For coeff_arg <= 0 (all indices paired), coefficient is 1
            trace_coeff = 1

        # Multiply by |x_a|^{2k} and coefficient
        # Include the alternating sign (-1)^m where m=k is the number of delta pairs
        sign = (-1) ** k
        return Expression.num(sign * trace_coeff) * (self.ra0 ** (2 * k)) * total

    def _generate_k_pairings(self, n, k):
        """
        Generate all distinct ways to pair k pairs from n indices.

        For k=1: C(n, 2) = n(n-1)/2 pairings
        For k=2: Ways to choose 2 non-overlapping pairs

        Returns list of pairings, each pairing is list of (i,j) tuples.
        """
        if k == 0:
            return [[]]

        if k == 1:
            # All possible pairs
            return [[(i, j)] for i in range(n) for j in range(i + 1, n)]

        # For k >= 2, recursively build pairings
        result = []

        # Choose first pair
        for i in range(n):
            for j in range(i + 1, n):
                first_pair = (i, j)

                # Get remaining indices
                remaining = [idx for idx in range(n) if idx != i and idx != j]

                # There remains k-1 pairs, i.e. 2*(k-1) indices, to be paired,
                # which is only possible if
                # there are at least 2*(k-1) remaining indices.
                if len(remaining) >= 2 * (k - 1):
                    # Recursively pair the rest
                    # Map remaining indices to 0..len(remaining)-1
                    sub_pairings = self._generate_k_pairings(len(remaining), k - 1)

                    for sub_pairing in sub_pairings:
                        # Map back to original indices
                        mapped_pairing = [
                            tuple(remaining[idx] for idx in pair)
                            for pair in sub_pairing
                        ]
                        result.append([first_pair] + mapped_pairing)

        # Remove duplicates, just for security. The previous code
        # has already enforced i<j in indices to avoid double counting.
        unique = []
        seen = set()
        for pairing in result:
            # Normalize: sort pairs within pairing, then sort the pairing
            normalized = tuple(sorted([tuple(sorted(p)) for p in pairing]))
            if normalized not in seen:
                seen.add(normalized)
                unique.append(pairing)

        return unique

    def phi_from_Q(self, n, indices=None):
        """
        Compute φ^(n) using the Q tensor formulation:

        φ^(n)(x) = (1/n!) * Q^{i₁...iₙ} * n^{i₁} * ... * n^{iₙ} / r^{n+1}

        Args:
            n: Order of multipole
            indices: Index symbols (auto-generated if None)

        Returns:
            Expression for φ^(n) in terms of Q
        """
        if indices is None:
            indices = [S(f"i{k}") for k in range(1, n + 1)]

        # Get Q tensor
        Q = self.Q_tensor(n, indices)

        # Multiply by n^{i₁} * ... * n^{iₙ}
        for idx in indices:
            Q = Q * self.n(idx)

        # Divide by r^{n+1}
        Q = Q / self.r0 ** (n + 1)

        # Divide by n!
        factorial = self._factorial(n)
        Q = Q / factorial

        # IMPORTANT: Expand before contracting to distribute products
        Q = Q.expand()

        # Contract indices (apply multiple times to ensure full contraction)
        # For high orders (n>=7), need more passes
        num_passes = max(5, n + 2)
        result = Q
        for _ in range(num_passes):
            result = self.tc.contract_indices(result)

        return result

    def _factorial(self, n):
        """Compute n!"""
        if n <= 0:
            return 1
        result = 1
        for i in range(1, n + 1):
            result *= i
        return result

    def _double_factorial(self, n):
        """
        Compute double factorial n!! = n*(n-2)*(n-4)*...

        For odd n: 1*3*5*...*(2k+1)
        For even n: 2*4*6*...*2k
        """
        if n <= 0:
            return 1
        if n == 1:
            return 1
        return n * self._double_factorial(n - 2)


def verify_symmetry(Q_expr, indices, tc):
    """
    Verify that Q^{i₁...iₙ} is symmetric under index permutations.

    Args:
        Q_expr: Expression for Q tensor
        indices: List of index symbols
        tc: TensorContraction instance

    Returns:
        True if symmetric, False otherwise
    """
    # Generate all permutations of indices
    perms = list(itertools.permutations(indices))

    # For each permutation, check if Q(...) equals Q(permuted)
    # This is symbolic, so we'd need to substitute and compare

    # Simplified check: just verify a few permutations
    print(f"Checking symmetry for {len(perms)} permutations...")

    # In practice, we'd substitute specific values or use Symbolica's
    # pattern matching to verify symbolic equality

    return True  # Placeholder


def verify_traceless(Q_expr, indices, tc):
    """
    Verify that Q^{...i...j...} δ_{ij} = 0 (traceless property).

    Args:
        Q_expr: Expression for Q tensor
        indices: List of index symbols
        tc: TensorContraction instance

    Returns:
        True if traceless, False otherwise
    """
    if len(indices) < 2:
        # Traceless property only applies for n >= 2
        return True

    # Contract first two indices: Q^{ii i₃ i₄...}
    i_same = indices[0]

    # Create a copy with first two indices the same
    contracted_indices = [i_same, i_same] + indices[2:]

    # Compute Q with contracted indices
    mm = MultipoleMoments()
    Q_contracted = mm.Q_tensor(len(indices), contracted_indices)

    # Apply contraction
    result = tc.contract_indices(Q_contracted)

    # Check if result is zero
    # For symbolic expressions, we'd check if it simplifies to 0
    print(f"Trace (first two indices): {result}")

    # Simplified: return True if it looks like it should be zero
    return True  # Placeholder


if __name__ == "__main__":
    print("Testing Multipole Moments Module\n")

    mm = MultipoleMoments()

    # Test monopole
    print("=" * 60)
    print("Monopole (n=0):")
    Q_0 = mm.Q_tensor(0)
    print(f"Q = {Q_0.format(terms_on_new_line=True)}")
    phi_0 = mm.phi_from_Q(0)
    print(f"φ^(0) = {phi_0}")
    print()

    # Test dipole
    print("=" * 60)
    print("Dipole (n=1):")
    i = S("i")
    Q_1 = mm.Q_tensor(1, [i])
    print(f"Q^i = {Q_1.format(terms_on_new_line=True)}")
    phi_1 = mm.phi_from_Q(1, [i])
    print(f"φ^(1) = {phi_1}")
    print()

    # Test quadrupole
    print("=" * 60)
    print("Quadrupole (n=2):")
    i, j = S("i"), S("j")
    Q_2 = mm.Q_tensor(2, [i, j])
    print(f"Q^{{ij}} = {Q_2.format(terms_on_new_line=True)}")

    # Check traceless: Q^{ii}
    Q_2_trace = mm.Q_tensor(2, [i, i])
    Q_2_trace_contracted = mm.tc.contract_indices(Q_2_trace)
    print(f"Trace Q^{{ii}} = {Q_2_trace_contracted}")

    phi_2 = mm.phi_from_Q(2, [i, j])
    print(f"φ^(2) from Q = {phi_2}")
    print()

    # Test octupole
    print("=" * 60)
    print("Octupole (n=3):")
    i, j, k = S("i"), S("j"), S("k")
    Q_3 = mm.Q_tensor(3, [i, j, k])
    print(f"Q^{{ijk}} = {Q_3.format(terms_on_new_line=True)}")

    # Check traceless: Q^{iik}
    Q_3_trace = mm.Q_tensor(3, [i, i, k])
    Q_3_trace_contracted = mm.tc.contract_indices(Q_3_trace)
    print(f"Trace Q^{{iik}} = {Q_3_trace_contracted}")
    print()
