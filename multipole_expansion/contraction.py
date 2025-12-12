"""
Tensor Index Contraction Module

This module implements Einstein summation convention for tensor expressions.
When an index appears twice in a product, it represents a sum over that index.

Examples:
    x_a^i * n^i = x_a · n (dot product)
    x_a^i * x_a^i = |x_a|^2 (magnitude squared)
    δ^{ij} * A^j = A^i (Kronecker delta contraction)
"""

from symbolica import Expression, S


class TensorContraction:
    """
    Handles tensor index contraction using pattern matching.

    Convention:
    - xa(i) represents x_a^i (source position component i)
    - x(i) represents x^i (observation position component i)
    - n(i) represents n^i = x^i/r (unit vector component i)
    - delta(i,j) represents δ^{ij} (Kronecker delta)
    """

    def __init__(self):
        """Initialize symbols for tensor operations."""
        # Vector components as functions
        # The S function maps a string to a symbolica symbol.
        # E.g., self.delta=S("delta") means
        # when symbolica sees a string containing "delta",
        # it knows that this string means some symbol self.delta,
        # whose behavior will be defined later in this file:
        # see line 125, the function _contract_delta.
        self.xa = S("xa")  # Source position
        self.x = S("x")  # Observation position
        self.n = S("n")  # Unit vector n = x/r
        self.delta = S("delta")  # Kronecker delta

        # Scalar quantities
        self.r0 = S("r0")  # |x| magnitude
        self.ra0 = S("ra0")  # |x_a| magnitude
        self.dot = S("dot")  # Dot product function

        # Index wildcards (for pattern matching)
        # This `x_` syntax is basically exactly the same as
        # the un-named function or pattern matching
        # of Wolfram Mathematica.
        self.i_ = S("i_")
        self.j_ = S("j_")
        self.k_ = S("k_")
        self.l_ = S("l_")

    def contract_indices(self, expr):
        """
        Apply Einstein summation convention to contract repeated indices.

        This searches for patterns where the same index appears twice
        and replaces them with the appropriate contracted form.

        Args:
            expr: Symbolica Expression containing indexed tensors

        Returns:
            Expression with contracted indices replaced
        """
        result = expr

        # Contract xa(i) * xa(i) -> ra0^2
        result = self._contract_xa_xa(result)

        # Contract xa(i) * n(i) -> dot(xa, n)
        result = self._contract_xa_n(result)

        # Contract x(i) * x(i) -> r0^2
        result = self._contract_x_x(result)

        # Contract x(i) * n(i) -> r0
        result = self._contract_x_n(result)

        # Contract n(i) * n(i) -> 1
        result = self._contract_n_n(result)

        # Contract delta(i,j) with vectors
        result = self._contract_delta(result)

        return result

    def _contract_xa_xa(self, expr):
        """Contract xa(i) * xa(i) = ra0^2"""
        # Pattern: xa(i_) * xa(i_) where same index appears twice
        pattern = self.xa(self.i_) * self.xa(self.i_)
        replacement = self.ra0**2
        return expr.replace(pattern, replacement)

    def _contract_xa_n(self, expr):
        """Contract xa(i) * n(i) = dot(xa, n)"""
        pattern = self.xa(self.i_) * self.n(self.i_)
        replacement = self.dot(self.xa, self.n)
        # Also handle reverse order
        result = expr.replace(pattern, replacement)
        pattern_rev = self.n(self.i_) * self.xa(self.i_)
        result = result.replace(pattern_rev, replacement)
        return result

    def _contract_x_x(self, expr):
        """Contract x(i) * x(i) = r0^2"""
        pattern = self.x(self.i_) * self.x(self.i_)
        replacement = self.r0**2
        return expr.replace(pattern, replacement)

    def _contract_x_n(self, expr):
        """Contract x(i) * n(i) = r0 (since n = x/r)"""
        pattern = self.x(self.i_) * self.n(self.i_)
        replacement = self.r0
        result = expr.replace(pattern, replacement)
        pattern_rev = self.n(self.i_) * self.x(self.i_)
        result = result.replace(pattern_rev, replacement)
        return result

    def _contract_n_n(self, expr):
        """Contract n(i) * n(i) = 1"""
        pattern = self.n(self.i_) * self.n(self.i_)
        replacement = Expression.num(1)
        return expr.replace(pattern, replacement)

    def _contract_delta(self, expr):
        """
        Contract Kronecker delta with vectors.
        delta(i,j) * A(j) = A(i)
        """
        # delta(i_, j_) * xa(j_) -> xa(i_)
        pattern = self.delta(self.i_, self.j_) * self.xa(self.j_)
        replacement = self.xa(self.i_)
        result = expr.replace(pattern, replacement)

        # Handle reverse order
        pattern_rev = self.xa(self.j_) * self.delta(self.i_, self.j_)
        result = result.replace(pattern_rev, replacement)

        # Same for n(i): contract n_i δ_ij and δ_ij n_j
        pattern = self.delta(self.i_, self.j_) * self.n(self.j_)
        replacement = self.n(self.i_)
        result = result.replace(pattern, replacement)

        pattern_rev = self.n(self.j_) * self.delta(self.i_, self.j_)
        result = result.replace(pattern_rev, replacement)

        # Because we only have position vector and normal vectors
        # throughout the calculation, the kronecker delta only needs to
        # take care of these replacement rules.

        # delta(i_, i_) -> 3 (trace in 3D)
        pattern = self.delta(self.i_, self.i_)
        replacement = Expression.num(3)
        result = result.replace(pattern, replacement)

        # IMPORTANT: Handle xa(i_)*delta(i_, k_) where i_ is repeated and k_ is free
        # When i is a dummy index (appears twice), xa(i)*delta(i,k) = xa(k)
        # This is the key for traceless verification of higher order tensors!
        # Note: This only applies when i appears elsewhere (making it a dummy index)

        # Try the contraction (though without knowing dummy vs free, we leave it)
        # The mathematical simplification xa(i)*delta(i,k) = xa(k) is valid
        # when i is summed over, which happens in traceless verification

        return result

    def expand_dot_products(self, expr):
        """
        Expand dot products back to component form if needed.
        dot(xa, n) -> xa(i) * n(i)

        This is useful when we need to work with explicit indices again.
        """
        pattern = self.dot(self.xa, self.n)
        # Use a fresh index 'i1' for the expansion
        i1 = S("i1")
        replacement = self.xa(i1) * self.n(i1)
        return expr.replace(pattern, replacement)


def create_indexed_product(tensor_fn, indices):
    """
    Create a product of indexed tensors.

    Args:
        tensor_fn: Function symbol (e.g., xa, n)
        indices: List of index symbols

    Returns:
        Product tensor_fn(i1) * tensor_fn(i2) * ...

    Example:
        >>> create_indexed_product(xa, ['i1', 'i2'])
        xa(i1) * xa(i2)
    """
    result = Expression.num(1)
    for idx in indices:
        idx_sym = S(idx) if isinstance(idx, str) else idx
        result = result * tensor_fn(idx_sym)
    return result


if __name__ == "__main__":
    # Test the contraction module
    print("Testing Tensor Contraction Module\n")

    tc = TensorContraction()

    # Test 1: xa(i) * xa(i) -> ra0^2
    test1 = tc.xa(S("i")) * tc.xa(S("i"))
    print("Test 1: xa(i) * xa(i)")
    print("Before:", test1)
    result1 = tc.contract_indices(test1)
    print("After:", result1)
    print()

    # Test 2: xa(i) * n(i) -> dot(xa, n)
    test2 = tc.xa(S("i")) * tc.n(S("i"))
    print("Test 2: xa(i) * n(i)")
    print("Before:", test2)
    result2 = tc.contract_indices(test2)
    print("After:", result2)
    print()

    # Test 3: More complex expression
    i = S("i")
    test3 = tc.xa(i) * tc.n(i) * tc.r0 ** (-3) + tc.xa(i) * tc.xa(i) * tc.r0 ** (-5)
    print("Test 3: xa(i)*n(i)/r0^3 + xa(i)*xa(i)/r0^5")
    print("Before:", test3)
    result3 = tc.contract_indices(test3)
    print("After:", result3)
    print()
