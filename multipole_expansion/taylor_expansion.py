"""
Taylor Expansion for Multipole Moments

This module implements the Taylor series expansion:
    φ(x) = 1/|x - x_a| = Σ φ^(n)(x)

where:
    φ^(n)(x) = ((-1)^n / n!) * x_a^{i₁} ... x_a^{iₙ} * ∂^n/∂x^{i₁}...∂x^{iₙ} (1/r)

The key challenge is handling the sum over all index combinations with
proper contraction of repeated indices.
"""

from symbolica import Expression, S
from .contraction import TensorContraction, create_indexed_product
from .derivatives import DerivativeEngine
import itertools


class TaylorExpansion:
    """
    Implements the multipole Taylor expansion.
    """

    def __init__(self):
        """Initialize with contraction and derivative engines."""
        self.tc = TensorContraction()
        self.de = DerivativeEngine()

        # Symbols
        self.xa = self.tc.xa
        self.n = self.tc.n
        self.r0 = self.tc.r0
        self.ra0 = self.tc.ra0
        self.dot = self.tc.dot
        self.delta = self.tc.delta

    def phi_n(self, n, use_contraction=True):
        """
        Compute the n-th term in the multipole expansion.

        φ^(n)(x) = ((-1)^n / n!) * x_a^{i₁} ... x_a^{iₙ} * ∂^n/∂x^{i₁}...∂x^{iₙ} (1/r)

        Args:
            n: Order of the multipole (0=monopole, 1=dipole, 2=quadrupole, ...)
            use_contraction: If True, apply index contraction rules

        Returns:
            Expression for φ^(n)
        """
        if n == 0:
            return self._phi_0()
        elif n == 1:
            return self._phi_1(use_contraction)
        elif n == 2:
            return self._phi_2(use_contraction)
        else:
            return self._phi_n_general(n, use_contraction)

    def _phi_0(self):
        """
        Monopole term: φ^(0) = 1/r
        """
        return Expression.num(1) / self.r0

    def _phi_1(self, use_contraction=True):
        """
        Dipole term: φ^(1) = -x_a^i * ∂/∂x^i (1/r)

        Result (after contraction): dot(xa, n) / r^2 = (x_a · x) / r^3
        """
        # Create index
        i = S("i")

        # x_a^i
        xa_i = self.xa(i)

        # ∂/∂x^i (1/r) = -x^i / r^3
        deriv = self.de.derivative_1_over_r(i)

        # -x_a^i * (-x^i / r^3) = x_a^i * x^i / r^3
        result = -xa_i * deriv

        if use_contraction:
            # Contract x_a^i * x^i -> x_a^i * n^i * r = dot(xa, n) * r
            # But we have x^i not n^i, so: x_a^i * x^i = dot(xa, n) * r
            # Actually, we need to be more careful here
            # Let's expand x^i = n^i * r
            result = result.replace(self.tc.x(i), self.n(i) * self.r0)
            result = result.expand()
            result = self.tc.contract_indices(result)
            result = self.tc.contract_indices(result)

        return result

    def _phi_2(self, use_contraction=True):
        """
        Quadrupole term: φ^(2) = (1/2) * x_a^i * x_a^j * ∂²/∂x^i∂x^j (1/r)

        Result: (3(x_a·n)² - x_a²) / (2r³)
        """
        # Create indices
        i, j = S("i"), S("j")

        # x_a^i * x_a^j
        xa_product = self.xa(i) * self.xa(j)

        # ∂²/∂x^i∂x^j (1/r) = (3 x^i x^j - δ^{ij} r²) / r^5
        deriv = self.de.second_derivative_1_over_r(i, j)

        # (1/2) * x_a^i * x_a^j * deriv
        result = Expression.num(1) / 2 * xa_product * deriv

        if use_contraction:
            # Replace x^i with n^i * r
            result = result.replace(self.tc.x(i), self.n(i) * self.r0)
            result = result.replace(self.tc.x(j), self.n(j) * self.r0)

            # Expand to distribute products
            result = result.expand()

            # Contract indices (multiple passes for complete contraction)
            result = self.tc.contract_indices(result)
            result = self.tc.contract_indices(result)
            result = self.tc.contract_indices(result)

        return result

    def _phi_n_general(self, n, use_contraction=True):
        """
        General n-th order term (placeholder for n >= 3).

        For full implementation, we would need to:
        1. Generate all n indices
        2. Compute the n-th derivative
        3. Contract all repeated indices
        """
        # Create n different indices
        indices = [S(f"i{k}") for k in range(1, n + 1)]

        # Create x_a^{i₁} * ... * x_a^{iₙ}
        xa_product = Expression.num(1)
        for idx in indices:
            xa_product = xa_product * self.xa(idx)

        # Compute n-th derivative (if implemented)
        try:
            deriv = self.de.nth_derivative_1_over_r(n, indices)
        except NotImplementedError:
            # Return symbolic form
            deriv_sym = S("deriv_n")
            return Expression.parse(f"((-1)^{n} / {n}!) * xa_product * deriv_{n}")

        # Sign and factorial
        sign = (-1) ** n
        factorial = self._factorial(n)

        result = Expression.num(sign) / factorial * xa_product * deriv

        if use_contraction:
            # Replace x^i with n^i * r for all indices
            for idx in indices:
                result = result.replace(self.tc.x(idx), self.n(idx) * self.r0)

            # Expand and contract multiple times
            result = result.expand()
            for _ in range(5):  # More passes for higher orders
                result = self.tc.contract_indices(result)

        return result

    def _factorial(self, n):
        """Compute factorial."""
        if n <= 0:
            return 1
        result = 1
        for i in range(1, n + 1):
            result *= i
        return result

    def multipole_series(self, max_order=3):
        """
        Compute the multipole series up to a given order.

        φ(x) ≈ φ^(0) + φ^(1) + φ^(2) + ... + φ^(max_order)

        Args:
            max_order: Maximum order to compute

        Returns:
            Sum of all terms up to max_order
        """
        result = Expression.num(0)
        for n in range(max_order + 1):
            try:
                phi_n = self.phi_n(n, use_contraction=True)
                result = result + phi_n
                print(f"φ^({n}) = {phi_n}")
            except Exception as e:
                print(f"Warning: Could not compute order {n}: {e}")
                break
        return result


def compare_to_exact(taylor_approx, num_terms):
    """
    Compare the Taylor approximation to the exact result.

    The exact result is: 1/|x - x_a|

    For small x_a/x ratio, the Taylor series should converge to this.
    """
    # This would require numerical evaluation
    # Placeholder for now
    pass


if __name__ == "__main__":
    print("Testing Taylor Expansion Module\n")

    te = TaylorExpansion()

    # Test monopole
    print("=" * 60)
    print("Monopole (n=0):")
    phi_0 = te.phi_n(0)
    print(f"φ^(0) = {phi_0}")
    print()

    # Test dipole
    print("=" * 60)
    print("Dipole (n=1):")
    phi_1 = te.phi_n(1)
    print(f"φ^(1) = {phi_1}")
    print()

    # Test quadrupole
    print("=" * 60)
    print("Quadrupole (n=2):")
    phi_2 = te.phi_n(2)
    print(f"φ^(2) = {phi_2}")
    print()

    # Full series
    print("=" * 60)
    print("Full multipole series (up to n=2):")
    series = te.multipole_series(max_order=2)
    print(f"\nTotal: {series}")
    print()
