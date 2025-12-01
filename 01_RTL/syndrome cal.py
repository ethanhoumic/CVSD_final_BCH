#!/usr/bin/env python3
"""
BCH Decoder Verification Tool

This script provides complete BCH decoding:
- Syndrome calculation
- Berlekamp-Massey algorithm for error locator polynomial
- Chien search for finding error locations
"""

import numpy as np

# Primitive polynomials for each field (represented as integers, LSB = constant term)
PRIMITIVE_POLYS = {
    6: 0b1000011,      # 1 + X + X^6
    8: 0b100011101,    # 1 + X^2 + X^3 + X^4 + X^8
    10: 0b10000001001  # 1 + X^3 + X^10
}

BCH_PARAMS = {
    1: {'n': 63, 'k': 51, 'm': 6, 't': 2},
    2: {'n': 255, 'k': 239, 'm': 8, 't': 2},
    3: {'n': 1023, 'k': 983, 'm': 10, 't': 4}
}


class GaloisField:
    """Galois Field GF(2^m) arithmetic"""
    
    def __init__(self, m, primitive_poly):
        self.m = m
        self.primitive_poly = primitive_poly
        self.n = (1 << m) - 1  # 2^m - 1
        
        # Build lookup tables
        self.exp_table = [0] * (2 * self.n + 1)
        self.log_table = [0] * (self.n + 1)
        
        self._build_tables()
    
    def _build_tables(self):
        """Build exponent and logarithm tables"""
        x = 1
        for i in range(self.n):
            self.exp_table[i] = x
            self.log_table[x] = i
            x <<= 1
            if x & (1 << self.m):
                x ^= self.primitive_poly
        
        for i in range(self.n, 2 * self.n + 1):
            self.exp_table[i] = self.exp_table[i - self.n]
    
    def multiply(self, a, b):
        """Multiply two elements in GF(2^m)"""
        if a == 0 or b == 0:
            return 0
        log_sum = self.log_table[a] + self.log_table[b]
        if log_sum >= self.n:
            log_sum -= self.n
        return self.exp_table[log_sum]
    
    def add(self, a, b):
        """Add two elements in GF(2^m) - XOR operation"""
        return a ^ b
    
    def power(self, a, p):
        """Compute a^p in GF(2^m)"""
        if a == 0:
            return 0
        log_a = self.log_table[a]
        log_result = (log_a * p) % self.n
        return self.exp_table[log_result]
    
    def alpha_power(self, i):
        """Return alpha^i"""
        i = i % self.n
        return self.exp_table[i]
    
    def inverse(self, a):
        """Compute multiplicative inverse of a"""
        if a == 0:
            raise ValueError("Cannot invert zero")
        return self.exp_table[self.n - self.log_table[a]]


def parse_llr_input(lines, code_type):
    """Parse LLR input data and extract received polynomial r(X)"""
    params = BCH_PARAMS[code_type]
    n = params['n']
    
    llr_values = []
    for line in lines:
        line = line.strip()
        if len(line) != 64:
            continue
        for i in range(8):
            byte_str = line[i*8:(i+1)*8]
            val = int(byte_str, 2)
            if val >= 128:
                val -= 256
            llr_values.append(val)
    
    r = [0] * n
    
    for i in range(1, n + 1):
        if i < len(llr_values):
            llr = llr_values[i]
            bit = 0 if llr >= 0 else 1
            r[n - i] = bit
    
    return r, llr_values


def compute_syndromes(r, gf, t):
    """
    Compute syndromes S_1, S_2, ..., S_2t
    
    S_i = r(alpha^i) = sum_{j=0}^{n-1} r_j * alpha^(i*j)
    """
    n = len(r)
    syndromes = []
    
    for i in range(1, 2*t + 1):
        S_i = 0
        for j in range(n):
            if r[j] == 1:
                alpha_power = gf.alpha_power(i * j)
                S_i = gf.add(S_i, alpha_power)
        syndromes.append(S_i)
    
    return syndromes


def berlekamp_massey(syndromes, gf, t):
    """
    Berlekamp-Massey algorithm to find error locator polynomial
    
    Returns:
        sigma: Error locator polynomial coefficients [σ_0, σ_1, ..., σ_l]
               where σ(X) = σ_0 + σ_1*X + σ_2*X^2 + ... + σ_l*X^l
        l: Degree of error locator polynomial
    
    Algorithm:
    - Initialize: σ^(-1)(X) = 1, l_{-1} = 0, d_{-1} = 1
    - For μ = 0 to 2t-1:
        1. Compute discrepancy: d_μ = S_{μ+1} + Σ(i=1 to l_μ) σ_i^(μ) * S_{μ+1-i}
        2. If d_μ = 0: σ^(μ+1) = σ^(μ), l_{μ+1} = l_μ
        3. If d_μ ≠ 0: 
           - Find ρ: the most recent iteration where d_ρ ≠ 0 and (ρ - l_ρ) is maximum
           - σ^(μ+1) = σ^(μ) + (d_μ / d_ρ) * X^(μ-ρ) * σ^(ρ)
           - l_{μ+1} = max(l_μ, l_ρ + μ - ρ)
    """
    # Initialize
    sigma = [[0] * (t + 1) for _ in range(2 * t + 1)]
    l = [0] * (2 * t + 1)
    d = [0] * (2 * t + 1)
    
    # σ^(-1)(X) = 1
    sigma[-1 + 1][0] = 1  # Shift index by 1
    l[-1 + 1] = 0
    d[-1 + 1] = 1
    
    # σ^(0)(X) = 1
    sigma[0][0] = 1
    l[0] = 0
    # Compute d_0 = S_1
    d[0] = syndromes[0]
    
    # Store candidates for correction
    cand_idx = -1
    cand_max = -1
    
    for mu in range(2 * t):
        if d[mu] == 0:
            # No correction needed
            sigma[mu + 1] = sigma[mu].copy()
            l[mu + 1] = l[mu]
        else:
            # Find the best candidate ρ
            if mu == 0:
                # First iteration with d_0 ≠ 0
                # Use σ^(-1) as reference
                rho = -1 + 1  # Index shift
                sigma[mu + 1] = sigma[mu].copy()
                sigma[mu + 1][1] = gf.multiply(d[mu], gf.inverse(d[rho]))
                l[mu + 1] = 1
            else:
                # Find ρ where d_ρ ≠ 0 and (ρ - l_ρ) is maximum
                rho = -1
                max_val = -1
                for i in range(mu):
                    if d[i] != 0 and (i - l[i]) > max_val:
                        rho = i
                        max_val = i - l[i]
                
                if rho == -1:
                    # No previous non-zero discrepancy, use initial
                    rho = -1 + 1
                
                # Compute correction factor
                d_mu_inv_d_rho = gf.multiply(d[mu], gf.inverse(d[rho]))
                
                # Update σ^(μ+1) = σ^(μ) + (d_μ / d_ρ) * X^(μ-ρ) * σ^(ρ)
                sigma[mu + 1] = sigma[mu].copy()
                shift = mu - rho
                for i in range(t + 1):
                    if i + shift <= t and sigma[rho][i] != 0:
                        correction = gf.multiply(d_mu_inv_d_rho, sigma[rho][i])
                        sigma[mu + 1][i + shift] = gf.add(sigma[mu + 1][i + shift], correction)
                
                # Update degree
                l[mu + 1] = max(l[mu], l[rho] + shift)
        
        # Compute next discrepancy d_{μ+1}
        if mu + 1 < 2 * t:
            d[mu + 1] = syndromes[mu + 1]
            for i in range(1, l[mu + 1] + 1):
                if sigma[mu + 1][i] != 0:
                    d[mu + 1] = gf.add(d[mu + 1], 
                                       gf.multiply(sigma[mu + 1][i], syndromes[mu + 1 - i]))
    
    return sigma[2 * t], l[2 * t]


def chien_search(sigma, l, gf, n):
    """
    Chien search to find error locations
    
    Tests if σ(α^(-j)) = 0 for j = 0, 1, ..., n-1
    If σ(α^(-j)) = 0, then α^j is a root, meaning error at position j
    
    Returns:
        error_locations: List of error positions (sorted in ascending order)
    """
    error_locations = []
    
    # For each possible error location j
    for j in range(n):
        # Evaluate σ(α^(-j))
        alpha_inv_j = gf.alpha_power((gf.n - j) % gf.n)  # α^(-j) = α^(2^m - 1 - j)
        
        result = sigma[0]  # σ_0
        alpha_inv_j_power = alpha_inv_j
        
        for i in range(1, l + 1):
            if sigma[i] != 0:
                term = gf.multiply(sigma[i], alpha_inv_j_power)
                result = gf.add(result, term)
            alpha_inv_j_power = gf.multiply(alpha_inv_j_power, alpha_inv_j)
        
        if result == 0:
            error_locations.append(j)
    
    return sorted(error_locations)


def correct_errors(r, error_locations):
    """
    Correct errors in received polynomial
    
    Returns:
        corrected: The corrected codeword
    """
    corrected = r.copy()
    for loc in error_locations:
        corrected[loc] ^= 1  # Flip the bit
    return corrected


def decode_bch(r, gf, t, n, verbose=True):
    """
    Complete BCH decoding process
    """
    if verbose:
        print("\n" + "="*70)
        print("BCH DECODING PROCESS")
        print("="*70)
    
    # Step 1: Compute syndromes
    syndromes = compute_syndromes(r, gf, t)
    
    if verbose:
        print("\nStep 1: Syndrome Calculation")
        print("-" * 70)
        for i, S in enumerate(syndromes, 1):
            if S == 0:
                print(f"  S_{i} = 0")
            else:
                log_S = gf.log_table[S]
                print(f"  S_{i} = {S} (0x{S:0{(gf.m+3)//4}x}, 0b{S:0{gf.m}b}) = α^{log_S}")
    
    # Check if all syndromes are zero
    if all(s == 0 for s in syndromes):
        if verbose:
            print("\n  => All syndromes are 0: No errors detected!")
        return True, [], r
    
    if verbose:
        print("\n  => Non-zero syndromes detected: Errors present")
    
    # Step 2: Berlekamp-Massey algorithm
    sigma, l = berlekamp_massey(syndromes, gf, t)
    
    if verbose:
        print("\nStep 2: Berlekamp-Massey Algorithm")
        print("-" * 70)
        print(f"  Error locator polynomial degree: l = {l}")
        print_polynomial(sigma, l, gf, "σ(X)")
        
        # Also print in binary coefficient form for hardware verification
        print(f"\n  Binary coefficient representation:")
        print(f"    σ(X) coefficients = [", end="")
        coef_strs = []
        for i in range(l + 1):
            coef_strs.append(f"0x{sigma[i]:0{(gf.m+3)//4}x}")
        print(", ".join(coef_strs) + "]")
    
    # Step 3: Chien search
    error_locations = chien_search(sigma, l, gf, n)
    
    if verbose:
        print("\nStep 3: Chien Search")
        print("-" * 70)
        print(f"  Number of roots found: {len(error_locations)}")
        print(f"  Error locations: {error_locations}")
    
    # Check if decoding is successful
    if len(error_locations) != l:
        if verbose:
            print(f"\n  => DECODING FAILED: Found {len(error_locations)} roots but degree is {l}")
        return False, error_locations, None
    
    if len(error_locations) > t:
        if verbose:
            print(f"\n  => DECODING FAILED: Too many errors ({len(error_locations)} > t={t})")
        return False, error_locations, None
    
    # Step 4: Correct errors
    corrected = correct_errors(r, error_locations)
    
    if verbose:
        print("\nStep 4: Error Correction")
        print("-" * 70)
        print(f"  Corrected {len(error_locations)} error(s)")
        print(f"\n  => DECODING SUCCESSFUL!")
    
    return True, error_locations, corrected


def syndrome_to_str(syndrome, m):
    """Convert syndrome value to various representations"""
    if syndrome == 0:
        return "0"
    return f"{syndrome} (0x{syndrome:0{(m+3)//4}x}, 0b{syndrome:0{m}b})"

def print_polynomial(sigma, l, gf, poly_name="σ(X)"):
    """
    Print error locator polynomial in a readable format
    Shows both polynomial form and power form
    """
    print(f"  {poly_name} coefficients:")
    for i in range(l + 1):
        if sigma[i] != 0:
            if sigma[i] == 1:
                print(f"    σ_{i} = 1")
            else:
                log_val = gf.log_table[sigma[i]]
                print(f"    σ_{i} = {sigma[i]} = α^{log_val}")
    
    # Print as polynomial expression
    print(f"  {poly_name} = ", end="")
    terms = []
    for i in range(l + 1):
        if sigma[i] != 0:
            # Show field element in power notation
            if sigma[i] == 1:
                coef_str = "1"
            else:
                log_val = gf.log_table[sigma[i]]
                coef_str = f"α^{log_val}"
            
            if i == 0:
                terms.append(coef_str)
            elif i == 1:
                if coef_str == "1":
                    terms.append("X")
                else:
                    terms.append(f"{coef_str}·X")
            else:
                if coef_str == "1":
                    terms.append(f"X^{i}")
                else:
                    terms.append(f"{coef_str}·X^{i}")
    print(" + ".join(terms) if terms else "0")

def main():
    print("=" * 70)
    print("BCH Decoder Verification Tool")
    print("=" * 70)
    
    # Example input
    example_input = """0000000000000000100000000000000000000000100000000000000010000000
0000000000000000000000000000000010000000100000001000000010000000
1000000010000000000000001000000010000000000000001000000000000000
0000000010000000000000001000000000000000100000000000000000000000
0000000010000000100000001000000000000000000000001000000010000000
0000000010000000000000001000000010000000100000000000000010000000
0000000010000000000000001000000010000000000000001000000010000000
0000000000000000000000001000000000000000100000000000000000000000
1000000010000000000000000000000000000000000000001000000000000000
0000000000000000100000001000000000000000000000001000000000000000
1000000010000000000000000000000000000000100000000000000010000000
0000000000000000000000000000000010000000000000001000000000000000
1000000010000000000000001000000010000000000000001000000000000000
0000000000000000000000001000000010000000100000001000000000000000
0000000010000000100000000000000000000000000000001000000010000000
0000000000000000100000001000000010000000100000001000000010000000
0000000000000000000000001000000010000000000000000000000010000000
1000000010000000000000000000000000000000000000000000000000000000
1000000010000000000000001000000000000000100000001000000000000000
1000000010000000000000000000000010000000100000000000000000000000
1000000010000000100000000000000010000000000000000000000010000000
0000000010000000000000001000000010000000000000000000000000000000
1000000010000000000000001000000010000000000000001000000010000000
0000000000000000100000001000000000000000000000001000000010000000
1000000010000000100000001000000000000000100000000000000000000000
1000000010000000000000001000000010000000100000000000000000000000
1000000010000000100000000000000000000000000000000000000000000000
1000000010000000100000000000000010000000100000000000000010000000
0000000000000000000000001000000010000000000000000000000010000000
0000000000000000000000000000000010000000100000001000000010000000
0000000000000000000000000000000000000000000000000000000000000000
1000000000000000000000001000000010000000000000000000000010000000
1000000010000000000000000000000010000000100000000000000000000000
1000000010000000100000001000000000000000100000000000000000000000
1000000000000000100000000000000010000000000000000000000000000000
0000000010000000000000000000000000000000000000001000000010000000
1000000010000000000000000000000000000000000000000000000010000000
1000000000000000100000000000000000000000100000001000000000000000
1000000010000000000000001000000010000000100000000000000000000000
1000000000000000100000000000000010000000000000001000000010000000
1000000000000000100000001000000000000000100000001000000000000000
0000000010000000100000000000000010000000100000000000000010000000
1000000000000000100000000000000010000000000000001000000000000000
0000000010000000100000001000000010000000100000001000000010000000
1000000000000000100000000000000010000000000000001000000000000000
0000000000000000100000000000000000000000100000001000000000000000
0000000010000000100000000000000000000000100000000000000000000000
1000000000000000000000000000000010000000100000001000000010000000
1000000000000000100000001000000000000000000000001000000010000000
0000000010000000100000001000000000000000100000000000000010000000
0000000010000000100000000000000000000000100000001000000000000000
1000000010000000000000001000000000000000000000000000000010000000
1000000010000000000000000000000010000000000000001000000000000000
0000000010000000100000000000000000000000100000001000000000000000
1000000010000000000000000000000000000000000000000000000010000000
0000000000000000000000001000000000000000000000000000000010000000
1000000010000000100000001000000010000000000000001000000000000000
0000000000000000100000001000000010000000000000001000000000000000
1000000010000000000000000000000010000000000000000000000010000000
0000000000000000100000001000000000000000000000000000000000000000
1000000000000000000000000000000010000000000000001000000000000000
0000000010000000100000001000000010000000000000001000000000000000
0000000010000000000000001000000000000000100000000000000000000000
0000000010000000100000001000000000000000100000001000000010000000
1000000000000000000000001000000000000000000000000000000000000000
0000000000000000000000001000000000000000100000000000000000000000
0000000000000000000000001000000010000000000000000000000000000000
0000000000000000000000001000000010000000100000000000000000000000
0000000000000000000000000000000010000000000000000000000010000000
0000000000000000000000001000000000000000100000000000000010000000
0000000000000000000000001000000000000000000000000000000010000000
0000000000000000100000000000000000000000000000000000000000000000
0000000000000000000000000000000010000000000000000000000000000000
0000000000000000100000000000000000000000000000000000000010000000
1000000000000000100000001000000010000000000000001000000000000000
0000000000000000100000000000000010000000000000000000000000000000
0000000010000000100000000000000000000000000000000000000000000000
1000000000000000100000000000000000000000000000000000000000000000
0000000000000000100000000000000000000000000000000000000000000000
1000000000000000000000000000000010000000100000001000000010000000
0000000000000000000000001000000000000000000000001000000010000000
1000000000000000100000001000000010000000000000000000000000000000
1000000000000000000000001000000010000000100000000000000000000000
1000000010000000000000001000000010000000100000000000000000000000
0000000010000000000000001000000010000000000000000000000010000000
0000000010000000100000000000000010000000000000001000000010000000
0000000000000000000000001000000000000000100000001000000010000000
1000000010000000000000001000000000000000000000001000000000000000
1000000010000000000000000000000000000000100000001000000000000000
1000000010000000000000000000000010000000100000001000000000000000
1000000010000000000000000000000010000000000000001000000010000000
0000000000000000000000000000000000000000000000000000000000000000
1000000010000000100000000000000010000000100000001000000010000000
1000000000000000100000001000000000000000100000000000000010000000
1000000010000000000000000000000010000000000000000000000010000000
1000000010000000000000001000000000000000100000000000000000000000
0000000010000000100000001000000000000000100000000000000010000000
1000000000000000000000000000000010000000000000000000000010000000
0000000010000000000000000000000000000000100000000000000010000000
1000000000000000000000001000000000000000000000001000000010000000
1000000010000000100000000000000010000000000000001000000010000000
1000000000000000000000000000000010000000000000000000000000000000
0000000010000000100000001000000010000000000000001000000010000000
0000000010000000000000000000000000000000000000001000000010000000
0000000000000000100000000000000010000000000000000000000000000000
0000000010000000100000001000000010000000000000000000000000000000
1000000010000000000000001000000010000000100000001000000000000000
1000000010000000100000000000000000000000100000000000000000000000
0000000010000000100000001000000000000000000000000000000010000000
1000000000000000100000000000000000000000000000000000000000000000
1000000000000000100000000000000010000000100000001000000000000000
1000000010000000000000001000000010000000100000000000000000000000
0000000000000000100000001000000000000000100000001000000000000000
1000000010000000100000001000000010000000000000000000000010000000
1000000000000000000000001000000010000000000000001000000000000000
0000000000000000100000001000000000000000100000000000000010000000
0000000010000000100000000000000000000000100000001000000000000000
1000000010000000100000000000000010000000000000001000000000000000
1000000010000000100000000000000010000000000000000000000010000000
1000000010000000100000001000000000000000000000000000000000000000
1000000000000000000000000000000000000000000000001000000010000000
0000000000000000000000001000000000000000100000001000000000000000
1000000000000000000000001000000010000000000000001000000000000000
0000000010000000000000000000000010000000100000001000000010000000
0000000000000000000000000000000000000000000000001000000010000000
1000000000000000100000001000000010000000000000001000000000000000
1000000010000000000000001000000010000000100000001000000010000000
1000000010000000000000001000000000000000100000001000000010000000"""
    
    lines = example_input.strip().split('\n')
    
    # Test for each code type
    for code_type in [1, 2, 3]:
        params = BCH_PARAMS[code_type]
        n, k, m, t = params['n'], params['k'], params['m'], params['t']
        
        print(f"\n{'='*70}")
        print(f"BCH Code ({n}, {k}): m={m}, t={t}")
        print(f"Primitive polynomial: {bin(PRIMITIVE_POLYS[m])}")
        print(f"{'='*70}")
        
        # Create Galois Field
        gf = GaloisField(m, PRIMITIVE_POLYS[m])
        
        # Parse input
        r, llr_values = parse_llr_input(lines, code_type)
        
        # Print received polynomial info
        error_positions = [i for i in range(n) if r[i] == 1]
        print(f"\nReceived polynomial r(X):")
        print(f"  Number of 1s: {sum(r)}")
        if len(error_positions) <= 20:
            print(f"  Positions with 1: {error_positions}")
        
        # Perform complete decoding
        success, error_locs, corrected = decode_bch(r, gf, t, n, verbose=True)
        
        if success and corrected is not None:
            # Verify correction by computing syndromes of corrected word
            corrected_syndromes = compute_syndromes(corrected, gf, t)
            all_zero = all(s == 0 for s in corrected_syndromes)
            
            print("\nVerification:")
            print("-" * 70)
            if all_zero:
                print("  ✓ Corrected codeword has all-zero syndromes")
            else:
                print("  ✗ Warning: Corrected codeword still has non-zero syndromes")


def test_with_known_errors():
    """Test with known error patterns"""
    print("\n" + "=" * 70)
    print("Testing BCH Decoder with Known Error Patterns")
    print("=" * 70)
    
    for code_type in [1, 2, 3]:
        params = BCH_PARAMS[code_type]
        n, k, m, t = params['n'], params['k'], params['m'], params['t']
        
        print(f"\n--- ({n}, {k}) BCH Code ---")
        
        gf = GaloisField(m, PRIMITIVE_POLYS[m])
        
        # Test 1: No errors
        print("\nTest 1: No errors")
        r = [0] * n
        success, locs, corrected = decode_bch(r, gf, t, n, verbose=False)
        assert success and len(locs) == 0, "Should detect no errors"
        print("  PASS")
        
        # Test 2: Single error
        print(f"\nTest 2: Single error at position 5")
        r = [0] * n
        r[5] = 1
        success, locs, corrected = decode_bch(r, gf, t, n, verbose=False)
        assert success and locs == [5], f"Should find error at position 5, got {locs}"
        assert corrected == [0] * n, "Should correct to all zeros"
        print("  PASS")
        
        # Test 3: Two errors
        print(f"\nTest 3: Two errors at positions 3 and 10")
        r = [0] * n
        r[3] = 1
        r[10] = 1
        success, locs, corrected = decode_bch(r, gf, t, n, verbose=False)
        assert success and set(locs) == {3, 10}, f"Should find errors at 3 and 10, got {locs}"
        assert corrected == [0] * n, "Should correct to all zeros"
        print("  PASS")
        
        # Test 4: Maximum correctable errors
        print(f"\nTest 4: {t} errors (maximum correctable)")
        r = [0] * n
        error_positions = list(range(t))
        for pos in error_positions:
            r[pos] = 1
        success, locs, corrected = decode_bch(r, gf, t, n, verbose=False)
        assert success and set(locs) == set(error_positions), f"Should find all {t} errors"
        assert corrected == [0] * n, "Should correct to all zeros"
        print("  PASS")


if __name__ == "__main__":
    # Run tests first
    test_with_known_errors()
    
    # Process example input
    main()