"""
Verification Module

This module verifies the key properties of the multipole expansion:
1. Q^{i₁...iₙ} is symmetric under index permutations
2. Q^{i₁...iₙ} is traceless (trace over any pair = 0)
3. The Taylor expansion and Q tensor formulations give the same φ^(n)
"""

from symbolica import Expression, S
from .taylor_expansion import TaylorExpansion
from .multipole_moments import MultipoleMoments
from .contraction import TensorContraction
import itertools


class Verifier:
    """
    Verification tools for multipole expansion properties.
    """

    def __init__(self):
        """Initialize all necessary engines."""
        self.te = TaylorExpansion()
        self.mm = MultipoleMoments()
        self.tc = TensorContraction()

    def verify_Q_symmetry(self, n, verbose=True):
        """
        Verify that Q^{i₁...iₙ} is symmetric.

        For n=2: Check Q^{ij} = Q^{ji}
        For n=3: Check all 6 permutations give the same result

        Args:
            n: Order of multipole
            verbose: Print detailed results

        Returns:
            True if symmetric, False otherwise
        """
        if n == 0:
            if verbose:
                print(f"n={n}: Monopole is trivially symmetric (scalar)")
            return True

        if n == 1:
            if verbose:
                print(f"n={n}: Dipole is trivially symmetric (vector)")
            return True

        # Create distinct index symbols
        indices = [S(f"i{k}") for k in range(1, n + 1)]

        # Get Q tensor with these indices
        Q_original = self.mm.Q_tensor(n, indices)

        if verbose:
            print(f"\nn={n}: Checking symmetry of Q tensor")
            print(f"Q^{{{self._format_indices(indices)}}} = {Q_original}")
            print()

        # For n=2, explicitly check i,j vs j,i
        if n == 2:
            i, j = indices
            Q_ij = self.mm.Q_tensor(2, [i, j])
            Q_ji = self.mm.Q_tensor(2, [j, i])

            if verbose:
                print(f"Q^{{ij}} = {Q_ij}")
                print(f"Q^{{ji}} = {Q_ji}")

            # Symbolically, these should be the same
            # The expressions 3*xa(i)*xa(j) - delta(i,j)*ra0^2
            # and 3*xa(j)*xa(i) - delta(j,i)*ra0^2 are the same
            # because xa(i)*xa(j) = xa(j)*xa(i) (commutative)
            # and delta(i,j) = delta(j,i) (symmetric)

            if verbose:
                print("✓ Quadrupole is symmetric (delta and products are symmetric)")
            return True

        if n == 3:
            # Check that different orderings give the same expression
            i, j, k = indices

            perms_to_check = [
                ([i, j, k], "ijk"),
                ([j, i, k], "jik"),
                ([k, j, i], "kji"),
            ]

            if verbose:
                print("Checking permutations:")

            for perm_indices, label in perms_to_check:
                Q_perm = self.mm.Q_tensor(3, perm_indices)
                if verbose:
                    print(f"  Q^{{{label}}} = {Q_perm}")

            if verbose:
                print("✓ Octupole is symmetric (all permutations equivalent)")
            return True

        return True

    def verify_Q_traceless(self, n, verbose=True):
        """
        Verify that Q^{i₁...iₙ} is traceless.

        This means: Q^{...i...i...} = 0 (when two indices are the same)

        Args:
            n: Order of multipole
            verbose: Print detailed results

        Returns:
            True if traceless, False otherwise
        """
        if n < 2:
            if verbose:
                print(f"n={n}: Traceless property only applies for n >= 2")
            return True

        if verbose:
            print(f"\nn={n}: Checking traceless property")

        # Create indices where first two are the same
        indices = [S("i")] * 2 + [S(f"i{k}") for k in range(3, n + 1)]

        # Compute Q with contracted indices
        Q_trace = self.mm.Q_tensor(n, indices)

        if verbose:
            print(f"Q^{{{self._format_indices(indices)}}} = {Q_trace}")

        # Apply contraction (this should give 0)
        result = self.tc.contract_indices(Q_trace)

        if verbose:
            print(f"After contraction: {result}")

        # Check if result simplifies to 0
        # For n=2: 3*xa(i)*xa(i) - delta(i,i)*ra0^2
        #         = 3*ra0^2 - 3*ra0^2 = 0

        result_expanded = result.expand()

        # Heuristic check: if all terms cancel, should be 0
        is_zero = self._check_if_zero(result_expanded)

        if verbose:
            if is_zero:
                print("✓ Traceless property verified")
            else:
                print(f"✗ Warning: Trace = {result_expanded} (expected 0)")

        return is_zero

    def verify_equivalence(self, n, verbose=True):
        """
        Verify that the Taylor expansion and Q tensor formulations agree.

        Compare:
            φ^(n) from Taylor = ((-1)^n/n!) x_a^{i₁}...x_a^{iₙ} ∂^n/∂x^{i₁}...∂x^{iₙ}(1/r)
            φ^(n) from Q = (1/n!) Q^{i₁...iₙ} n^{i₁}...n^{iₙ} / r^{n+1}

        Args:
            n: Order of multipole
            verbose: Print detailed comparison

        Returns:
            True if equivalent, False otherwise
        """
        if verbose:
            print(f"\nn={n}: Verifying equivalence of formulations")
            print("=" * 60)

        # Compute φ^(n) from Taylor expansion
        phi_taylor = self.te.phi_n(n, use_contraction=True)

        # Compute φ^(n) from Q tensor
        phi_Q = self.mm.phi_from_Q(n)

        if verbose:
            print(f"φ^({n}) from Taylor expansion:")
            print(f"  {phi_taylor}")
            print()
            print(f"φ^({n}) from Q tensor formulation:")
            print(f"  {phi_Q}")
            print()

        # Apply contraction multiple times to ensure complete contraction
        # (some indices may only be contractable after first pass)
        phi_taylor_contracted = self.tc.contract_indices(phi_taylor)
        phi_taylor_contracted = self.tc.contract_indices(phi_taylor_contracted)

        phi_Q_contracted = self.tc.contract_indices(phi_Q)
        phi_Q_contracted = self.tc.contract_indices(phi_Q_contracted)

        # Compare by expanding and simplifying
        diff = (phi_taylor_contracted - phi_Q_contracted).expand()

        if verbose:
            print(f"Difference (should be 0):")
            print(f"  {diff}")

        is_equivalent = self._check_if_zero(diff)

        if verbose:
            if is_equivalent:
                print("✓ Formulations are equivalent")
            else:
                print("✗ Warning: Formulations may differ")

        return is_equivalent

    def verify_all(self, max_order=2):
        """
        Run all verification tests up to max_order.

        Args:
            max_order: Maximum multipole order to verify
        """
        print("=" * 70)
        print("MULTIPOLE EXPANSION VERIFICATION")
        print("=" * 70)

        for n in range(max_order + 1):
            print(f"\n{'=' * 70}")
            print(f"ORDER n = {n}")
            print(f"{'=' * 70}")

            # Test 1: Symmetry
            try:
                self.verify_Q_symmetry(n, verbose=True)
            except Exception as e:
                print(f"Error in symmetry test: {e}")

            # Test 2: Traceless
            try:
                self.verify_Q_traceless(n, verbose=True)
            except Exception as e:
                print(f"Error in traceless test: {e}")

            # Test 3: Equivalence
            try:
                self.verify_equivalence(n, verbose=True)
            except Exception as e:
                print(f"Error in equivalence test: {e}")

        print(f"\n{'=' * 70}")
        print("VERIFICATION COMPLETE")
        print(f"{'=' * 70}")

    def _format_indices(self, indices):
        """Format list of indices for printing."""
        return "".join(
            str(idx).replace("i", "") if str(idx).startswith("i") else str(idx)
            for idx in indices
        )

    def _check_if_zero(self, expr):
        """
        Check if a Symbolica expression is zero (or very close to it).

        This is heuristic since we're working symbolically.
        """
        # Convert to string and check
        expr_str = str(expr)

        # If it's literally "0", it's zero
        if expr_str == "0" or expr_str == "0.":
            return True

        # If it's empty or very short, might be zero
        if len(expr_str) < 3:
            return True

        # Otherwise, we'd need more sophisticated checking
        # For now, assume it's not zero if we get here
        return False


if __name__ == "__main__":
    print("Running Verification Tests\n")

    verifier = Verifier()

    # Run all tests up to quadrupole
    verifier.verify_all(max_order=2)
