"""
Spherical Harmonics Expansion

This module converts Cartesian multipole moments Q^{i₁...iₙ} to
spherical multipole moments q_{lm}.

The spherical expansion is:
    φ(x) = Σ_{l,m} (4π/(2l+1)) * q_{lm}/r^{l+1} * Y_{lm}(θ,φ)

where q_{lm} = r_a^l * Y_{lm}^*(θ_a, φ_a)

For the connection to Cartesian moments, we need to relate:
    - Cartesian components x_a^i to spherical components (r_a, θ_a, φ_a)
    - Cartesian Q^{i₁...iₗ} to spherical q_{lm}

This is done through the spherical harmonics addition theorem and
the traceless symmetric tensor decomposition.
"""

from symbolica import Expression, S
from .contraction import TensorContraction
from .multipole_moments import MultipoleMoments
import math


class SphericalExpansion:
    """
    Convert Cartesian multipole moments to spherical form.
    """

    def __init__(self):
        """Initialize engines."""
        self.tc = TensorContraction()
        self.mm = MultipoleMoments()

        # Cartesian components (source position)
        self.xa = self.tc.xa
        self.ra0 = self.tc.ra0

        # Spherical symbols
        self.theta_a = S("theta_a")  # Polar angle of source
        self.phi_a = S("phi_a")  # Azimuthal angle of source

        # Observation point spherical coordinates
        self.theta = S("theta")
        self.phi = S("phi")
        self.r = self.tc.r0

    def q_lm(self, l, m):
        """
        Compute spherical multipole moment q_{lm}.

        For a point source at x_a:
            q_{lm} = r_a^l * Y_{lm}^*(θ_a, φ_a)

        where Y_{lm} are spherical harmonics.

        Args:
            l: Angular momentum quantum number (l >= 0)
            m: Magnetic quantum number (-l <= m <= l)

        Returns:
            Expression for q_{lm} in terms of Cartesian components
        """
        if abs(m) > l:
            raise ValueError(f"m must satisfy -l <= m <= l, got l={l}, m={m}")

        # For now, we provide the connection formulas for low l values
        if l == 0:
            return self._q_00()
        elif l == 1:
            return self._q_1m(m)
        elif l == 2:
            return self._q_2m(m)
        else:
            raise NotImplementedError(f"q_{{lm}} for l={l} not yet implemented")

    def _q_00(self):
        """
        Monopole: q_{0,0} = √(1/4π)

        In normalized form: q_{0,0} = 1 (charge)
        """
        # Y_{0,0} = √(1/4π), so q_{0,0} = √(1/4π) for unit charge
        # We'll use unnormalized form: q_{0,0} = 1
        return Expression.num(1)

    def _q_1m(self, m):
        """
        Dipole moments: q_{1,m} for m = -1, 0, 1

        Connection to Cartesian:
            q_{1,0} ∝ x_a^3 = z_a
            q_{1,1} ∝ -(x_a^1 + i*x_a^2) = -(x_a + i*y_a)
            q_{1,-1} ∝ (x_a^1 - i*x_a^2) = (x_a - i*y_a)

        Here we use convention: x_a^1=x, x_a^2=y, x_a^3=z
        """
        # Define Cartesian components
        # xa(1) = x, xa(2) = y, xa(3) = z

        if m == 0:
            # q_{1,0} ∝ z_a
            return self.xa(S("z"))
        elif m == 1:
            # q_{1,1} ∝ -(x + iy)
            # We'll represent complex numbers symbolically
            i_unit = S("I")  # Imaginary unit
            return -(self.xa(S("x")) + i_unit * self.xa(S("y")))
        elif m == -1:
            # q_{1,-1} ∝ (x - iy)
            i_unit = S("I")
            return self.xa(S("x")) - i_unit * self.xa(S("y"))
        else:
            raise ValueError(f"Invalid m={m} for l=1")

    def _q_2m(self, m):
        """
        Quadrupole moments: q_{2,m} for m = -2, -1, 0, 1, 2

        Connection to Cartesian Q^{ij}:
            q_{2,0} ∝ Q^{33} - (Q^{11} + Q^{22})/2 ∝ (3z² - r²)
            q_{2,±1} ∝ Q^{13} ± iQ^{23} ∝ xz ± iyz
            q_{2,±2} ∝ (Q^{11} - Q^{22}) ± 2iQ^{12} ∝ (x² - y²) ± 2ixy

        These are related to the traceless Q tensor.
        """
        i_unit = S("I")

        # Get Cartesian components
        x = self.xa(S("x"))
        y = self.xa(S("y"))
        z = self.xa(S("z"))

        if m == 0:
            # q_{2,0} ∝ 3z² - r²
            # From Q^{33}: 3*z*z - ra0^2 (but we need to be careful with indices)
            # Actually: 2z² - x² - y²
            return 2 * z**2 - x**2 - y**2

        elif m == 1:
            # q_{2,1} ∝ -xz - iyz = -(x + iy)z
            return -(x + i_unit * y) * z

        elif m == -1:
            # q_{2,-1} ∝ xz - iyz = (x - iy)z
            return (x - i_unit * y) * z

        elif m == 2:
            # q_{2,2} ∝ (x² - y²) + 2ixy = (x + iy)²
            return (x + i_unit * y) ** 2

        elif m == -2:
            # q_{2,-2} ∝ (x² - y²) - 2ixy = (x - iy)²
            return (x - i_unit * y) ** 2

        else:
            raise ValueError(f"Invalid m={m} for l=2")

    def Q_to_q_conversion(self, l, verbose=True):
        """
        Show the conversion from Cartesian Q^{i₁...iₗ} to spherical q_{lm}.

        This prints the formulas for all m values at a given l.

        Args:
            l: Angular momentum
            verbose: Print the conversion formulas
        """
        if verbose:
            print(f"Conversion from Q^{{i₁...i_{l}}} to q_{{{l},m}}:")
            print("=" * 60)

        for m in range(-l, l + 1):
            try:
                q = self.q_lm(l, m)
                if verbose:
                    print(f"q_{{{l},{m:+d}}} = {q}")
            except Exception as e:
                if verbose:
                    print(f"q_{{{l},{m:+d}}} : {e}")

        if verbose:
            print()

    def spherical_series(self, max_l=2):
        """
        Generate the full spherical multipole series up to l=max_l.

        φ(x) = Σ_{l=0}^{max_l} Σ_{m=-l}^{l} (4π/(2l+1)) q_{lm}/r^{l+1} Y_{lm}(θ,φ)

        Returns a symbolic representation.
        """
        # This is complex as it requires implementing spherical harmonics
        # For now, we return the moments themselves
        result = {}

        for l in range(max_l + 1):
            result[l] = {}
            for m in range(-l, l + 1):
                try:
                    result[l][m] = self.q_lm(l, m)
                except Exception as e:
                    result[l][m] = f"Error: {e}"

        return result

    def print_spherical_moments(self, max_l=2):
        """Print all spherical moments up to max_l in a nice format."""
        print("=" * 70)
        print("SPHERICAL MULTIPOLE MOMENTS")
        print("=" * 70)

        for l in range(max_l + 1):
            print(
                f"\nl = {l} ({'monopole' if l == 0 else 'dipole' if l == 1 else 'quadrupole' if l == 2 else 'multipole'}):"
            )
            print("-" * 70)

            for m in range(-l, l + 1):
                try:
                    q = self.q_lm(l, m)
                    print(f"  q_{{{l},{m:+2d}}} = {q}")
                except Exception as e:
                    print(f"  q_{{{l},{m:+2d}}} : Not implemented")
            print()


def compare_cartesian_spherical(n=2):
    """
    Compare the Cartesian and spherical representations for order n.

    For n=2 (quadrupole), show how Q^{ij} components relate to q_{2,m}.
    """
    print("=" * 70)
    print(f"CARTESIAN vs SPHERICAL COMPARISON (n={n})")
    print("=" * 70)

    mm = MultipoleMoments()
    se = SphericalExpansion()

    # Cartesian Q tensor
    print("\nCartesian Q tensor:")
    print("-" * 70)
    i, j = S("i"), S("j")
    Q = mm.Q_tensor(n, [i, j])
    print(f"Q^{{ij}} = {Q}")
    print()

    # Specific components
    print("Specific components:")
    for idx1 in ["x", "y", "z"]:
        for idx2 in ["x", "y", "z"]:
            Q_comp = mm.Q_tensor(n, [S(idx1), S(idx2)])
            print(f"  Q^{{{idx1}{idx2}}} = {Q_comp}")
    print()

    # Spherical moments
    print("Spherical moments:")
    print("-" * 70)
    se.Q_to_q_conversion(n, verbose=True)


if __name__ == "__main__":
    print("Testing Spherical Expansion Module\n")

    se = SphericalExpansion()

    # Print all spherical moments
    se.print_spherical_moments(max_l=2)

    # Compare Cartesian and spherical
    print()
    compare_cartesian_spherical(n=2)
