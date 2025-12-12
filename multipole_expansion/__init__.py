"""
Multipole Expansion Package

A pure Symbolica implementation for computing multipole moments in electrodynamics.

Main Components:
    - TensorContraction: Handle Einstein summation and index contraction
    - DerivativeEngine: Compute derivatives of 1/r^n
    - TaylorExpansion: Multipole Taylor series
    - MultipoleMoments: Symmetric traceless Q tensors
    - Verifier: Verification of symmetry, traceless, and equivalence
    - SphericalExpansion: Convert to spherical harmonics representation

Quick Start:
    >>> from multipole_expansion import MultipoleExpansion
    >>> mp = MultipoleExpansion()
    >>> mp.compute_all(max_order=2)
"""

from .contraction import TensorContraction, create_indexed_product
from .derivatives import DerivativeEngine, generate_derivative_table
from .taylor_expansion import TaylorExpansion
from .multipole_moments import MultipoleMoments
from .verification import Verifier
from .spherical_expansion import SphericalExpansion

__version__ = "1.0.0"
__author__ = "Multipole Expansion Project"

__all__ = [
    "TensorContraction",
    "DerivativeEngine",
    "TaylorExpansion",
    "MultipoleMoments",
    "Verifier",
    "SphericalExpansion",
    "MultipoleExpansion",
    "create_indexed_product",
    "generate_derivative_table",
]


class MultipoleExpansion:
    """
    Main interface for multipole expansion computations.

    This class provides a unified interface to all the multipole expansion
    functionality.
    """

    def __init__(self, max_order=3):
        """
        Initialize the multipole expansion framework.

        Args:
            max_order: Maximum order to compute (default: 3)
        """
        self.max_order = max_order

        # Initialize all engines
        self.tc = TensorContraction()
        self.de = DerivativeEngine()
        self.te = TaylorExpansion()
        self.mm = MultipoleMoments()
        self.verifier = Verifier()
        self.se = SphericalExpansion()

    def taylor_term(self, n):
        """
        Compute φ^(n) using Taylor expansion.

        Args:
            n: Multipole order

        Returns:
            Expression for φ^(n)
        """
        return self.te.phi_n(n, use_contraction=True)

    def compute_Q_tensor(self, n, indices=None):
        """
        Compute the symmetric traceless Q^{i₁...iₙ} tensor.

        Args:
            n: Multipole order
            indices: Index symbols (auto-generated if None)

        Returns:
            Expression for Q tensor
        """
        return self.mm.Q_tensor(n, indices)

    def Q_formulation(self, n):
        """
        Compute φ^(n) using Q tensor formulation.

        Args:
            n: Multipole order

        Returns:
            Expression for φ^(n)
        """
        return self.mm.phi_from_Q(n)

    def verify_equivalence(self, n):
        """
        Verify that Taylor and Q formulations agree for order n.

        Args:
            n: Multipole order
        """
        return self.verifier.verify_equivalence(n, verbose=True)

    def verify_Q_properties(self, n):
        """
        Verify that Q is symmetric and traceless.

        Args:
            n: Multipole order
        """
        print(f"\nVerifying Q^{{i₁...i_{n}}} properties:")
        print("=" * 60)

        sym = self.verifier.verify_Q_symmetry(n, verbose=True)
        trace = self.verifier.verify_Q_traceless(n, verbose=True)

        return sym and trace

    def spherical_moment(self, l, m):
        """
        Compute spherical multipole moment q_{lm}.

        Args:
            l: Angular momentum
            m: Magnetic quantum number

        Returns:
            Expression for q_{lm}
        """
        return self.se.q_lm(l, m)

    def compute_all(self, max_order=None, run_verification=True):
        """
        Compute all multipole moments up to max_order and optionally verify.

        Args:
            max_order: Maximum order (uses self.max_order if None)
            run_verification: Run verification tests
        """
        if max_order is None:
            max_order = self.max_order

        print("=" * 70)
        print("MULTIPOLE EXPANSION COMPUTATION")
        print("=" * 70)

        # Compute Taylor expansion terms
        print("\n1. TAYLOR EXPANSION FORMULATION")
        print("-" * 70)
        for n in range(max_order + 1):
            try:
                phi_n = self.taylor_term(n)
                order_name = (
                    ["Monopole", "Dipole", "Quadrupole", "Octupole"][n]
                    if n < 4
                    else f"n={n}"
                )
                print(f"{order_name}: φ^({n}) = {phi_n}")
            except Exception as e:
                print(f"φ^({n}): Error - {e}")

        # Compute Q tensors
        print("\n2. Q TENSOR FORMULATION")
        print("-" * 70)
        for n in range(max_order + 1):
            try:
                Q_n = self.compute_Q_tensor(n)
                phi_n = self.Q_formulation(n)
                order_name = (
                    ["Monopole", "Dipole", "Quadrupole", "Octupole"][n]
                    if n < 4
                    else f"n={n}"
                )
                print(f"{order_name}:")
                print(f"  Q^{{...}} = {Q_n}")
                print(f"  φ^({n}) = {phi_n}")
            except Exception as e:
                print(f"Q^({n}): Error - {e}")

        # Spherical moments
        print("\n3. SPHERICAL HARMONICS FORMULATION")
        print("-" * 70)
        self.se.print_spherical_moments(max_l=min(max_order, 2))

        # Run verification
        if run_verification:
            print("\n4. VERIFICATION")
            print("-" * 70)
            self.verifier.verify_all(max_order=min(max_order, 2))

    def summary(self):
        """Print a summary of available methods."""
        print("=" * 70)
        print("MULTIPOLE EXPANSION - AVAILABLE METHODS")
        print("=" * 70)
        print("""
1. Compute multipole terms:
   - taylor_term(n)      : φ^(n) from Taylor expansion
   - compute_Q_tensor(n) : Symmetric traceless Q tensor
   - Q_formulation(n)    : φ^(n) from Q tensor
   - spherical_moment(l,m): q_{lm} in spherical basis

2. Verification:
   - verify_equivalence(n)  : Check Taylor = Q formulation
   - verify_Q_properties(n) : Check symmetry & traceless

3. Complete computation:
   - compute_all(max_order) : Everything at once
        """)


if __name__ == "__main__":
    # Example usage
    print("Multipole Expansion Package\n")

    mp = MultipoleExpansion(max_order=2)

    # Show summary
    mp.summary()

    # Run complete computation
    mp.compute_all(max_order=2, run_verification=True)
