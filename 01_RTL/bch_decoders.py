#!/usr/bin/env python3
"""
BCH Decoder Verification Tool

This script provides complete BCH decoding:
- Syndrome calculation
- Berlekamp-Massey algorithm for error locator polynomial
- Chien search for finding error locations
- Chase algorithm for soft-decision decoding
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
    
    for mu in range(2 * t):
        if d[mu] == 0:
            # No correction needed
            sigma[mu + 1] = sigma[mu].copy()
            l[mu + 1] = l[mu]
        else:
            # Find the best candidate ρ
            if mu == 0:
                # First iteration with d_0 ≠ 0
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
    
    for j in range(n):
        alpha_inv_j = gf.alpha_power((gf.n - j) % gf.n)
        
        result = sigma[0]
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
    Complete BCH hard-decision decoding process
    """
    if verbose:
        print("\n" + "="*70)
        print("BCH HARD-DECISION DECODING PROCESS")
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


# ============================================================================
# SOFT-DECISION DECODING (Chase Algorithm)
# ============================================================================

def find_least_reliable_bits(llr_values, n, p):
    """
    Find the p least reliable bit positions
    LLR with smaller absolute value = less reliable
    
    Args:
        llr_values: List of LLR values (index 0 is don't care, index 1 to n are valid)
                    LLR[i] corresponds to r[n-i], i.e., LLR[1]->r[n-1], LLR[n]->r[0]
        n: Codeword length
        p: Number of least reliable bits to find
    
    Returns:
        List of p least reliable bit positions (in polynomial coefficient order)
    """
    reliability = []
    
    # 遍歷所有多項式係數位置 r[0] 到 r[n-1]
    for poly_pos in range(n):
        llr_idx = n - poly_pos  # r[poly_pos] 對應 LLR[n - poly_pos]
        if llr_idx < len(llr_values) and llr_idx >= 1:
            reliability.append((poly_pos, abs(llr_values[llr_idx])))
    
    # Sort by reliability (ascending = least reliable first)
    reliability.sort(key=lambda x: x[1])
    
    # Return the p least reliable positions
    least_reliable = [pos for pos, _ in reliability[:p]]
    
    return least_reliable


def generate_test_patterns(r, least_reliable_positions, p):
    """
    Generate 2^p test patterns by flipping combinations of least reliable bits
    
    Args:
        r: Original received polynomial (hard decision)
        least_reliable_positions: List of p least reliable bit positions
        p: Number of least reliable bits
    
    Returns:
        List of 2^p test patterns
    """
    test_patterns = []
    
    for i in range(2 ** p):
        pattern = r.copy()
        for j in range(p):
            if (i >> j) & 1:  # Check if j-th bit should be flipped
                pos = least_reliable_positions[j]
                pattern[pos] ^= 1
        test_patterns.append(pattern)
    
    return test_patterns


def decode_test_patterns(test_patterns, gf, t, n):
    """
    Decode all test patterns using hard-decision decoding
    
    Args:
        test_patterns: List of test patterns to decode
        gf: Galois Field object
        t: Error correction capability
        n: Codeword length
    
    Returns:
        List of successfully decoded results
    """
    results = []
    
    for idx, pattern in enumerate(test_patterns):
        success, error_locs, corrected = decode_bch(pattern, gf, t, n, verbose=False)
        if success and corrected is not None:
            results.append({
                'pattern_idx': idx,
                'corrected': corrected,
                'error_locations': error_locs
            })
    
    return results


def calculate_correlation(corrected, llr_values, n):
    """
    Calculate correlation value for a decoded result
    
    Formula: sum_{i=0}^{n-1} l_i * (1 - 2 * c_i)
    where l_i is the LLR of r_i, and c_i is the decoded bit
    
    c_i = 0 => (1 - 2*c_i) = 1
    c_i = 1 => (1 - 2*c_i) = -1
    
    Args:
        corrected: Corrected codeword
        llr_values: Original LLR values
        n: Codeword length
    
    Returns:
        Correlation value (higher is better)
    """
    correlation = 0
    
    for i in range(n):
        llr_idx = n - i  # Convert polynomial position to LLR index
        if llr_idx < len(llr_values):
            l_i = llr_values[llr_idx]
            c_i = corrected[i]
            correlation += l_i * (1 - 2 * c_i)
    
    return correlation


def select_best_result(results, llr_values, n):
    """
    Select the decoded result with maximum correlation value
    
    Args:
        results: List of successfully decoded results
        llr_values: Original LLR values
        n: Codeword length
    
    Returns:
        Best result (with highest correlation)
    """
    best_result = None
    max_correlation = float('-inf')
    
    for result in results:
        corr = calculate_correlation(result['corrected'], llr_values, n)
        result['correlation'] = corr
        
        if corr > max_correlation:
            max_correlation = corr
            best_result = result
    
    return best_result


def decode_bch_soft_decision(r, llr_values, gf, t, n, p=2, verbose=True):
    """
    BCH Soft-Decision Decoding using Chase Algorithm
    
    Args:
        r: Received polynomial (hard decision from LLR)
        llr_values: Original LLR values
        gf: Galois Field object
        t: Error correction capability
        n: Codeword length
        p: Number of least reliable bits to test (default: 2)
        verbose: Print detailed output
    
    Returns:
        success: Boolean indicating if decoding succeeded
        error_locations: List of error positions (relative to original r)
        corrected: Corrected codeword (or None if failed)
    """
    if verbose:
        print("\n" + "="*70)
        print("BCH SOFT-DECISION DECODING (Chase Algorithm)")
        print(f"Testing p = {p} least reliable bits -> {2**p} test patterns")
        print("="*70)
    
    # Step 1: Find p least reliable bits
    least_reliable = find_least_reliable_bits(llr_values, n, p)
    
    if verbose:
        print(f"\nStep 1: Find Least Reliable Bits")
        print("-" * 70)
        print(f"  Least reliable bit positions: {least_reliable}")
        print(f"  Their LLR values: ", end="")
        llr_vals = []
        for pos in least_reliable:
            llr_idx = n - pos
            if llr_idx < len(llr_values):
                llr_vals.append(f"r{pos}={llr_values[llr_idx]}")
        print(", ".join(llr_vals))
    
    # Step 2: Generate 2^p test patterns
    test_patterns = generate_test_patterns(r, least_reliable, p)
    
    # 記錄每個 pattern 翻轉了哪些位置
    pattern_flipped_bits = []
    for idx in range(2 ** p):
        flipped = []
        for j in range(p):
            if (idx >> j) & 1:
                flipped.append(least_reliable[j])
        pattern_flipped_bits.append(flipped)
    
    if verbose:
        print(f"\nStep 2: Generate Test Patterns")
        print("-" * 70)
        print(f"  Generated {len(test_patterns)} test patterns")
        for idx, flipped in enumerate(pattern_flipped_bits):
            flip_str = ", ".join([f"r{pos}" for pos in flipped]) if flipped else "none"
            print(f"    Pattern {idx}: flip [{flip_str}]")
    
    # Step 3: Decode all test patterns
    if verbose:
        print(f"\nStep 3: Decode Test Patterns")
        print("-" * 70)
    
    results = []
    for idx, pattern in enumerate(test_patterns):
        if verbose:
            print(f"\n  --- Pattern {idx} ---")
            flip_str = ", ".join([f"r{pos}" for pos in pattern_flipped_bits[idx]]) if pattern_flipped_bits[idx] else "none"
            print(f"  Flipped bits: [{flip_str}]")
        
        # Compute syndromes
        syndromes = compute_syndromes(pattern, gf, t)
        
        if verbose:
            print(f"  Syndromes:")
            for i, S in enumerate(syndromes, 1):
                if S == 0:
                    print(f"    S_{i} = 0")
                else:
                    log_S = gf.log_table[S]
                    print(f"    S_{i} = {S} (0x{S:0{(gf.m+3)//4}x}, 0b{S:0{gf.m}b}) = α^{log_S}")
        
        # Check if all syndromes are zero
        if all(s == 0 for s in syndromes):
            if verbose:
                print(f"  => No errors in this pattern")
            # No errors, the pattern itself is a valid codeword
            flipped_set = set(pattern_flipped_bits[idx])
            error_from_original_r = sorted(list(flipped_set))
            results.append({
                'pattern_idx': idx,
                'corrected': pattern.copy(),
                'flipped_bits': pattern_flipped_bits[idx],
                'errors_in_pattern': [],
                'errors_from_r': error_from_original_r
            })
            continue
        
        # Berlekamp-Massey
        sigma, l = berlekamp_massey(syndromes, gf, t)
        
        if verbose:
            print(f"  Error locator polynomial:")
            print(f"    Degree: l = {l}")
            
            # Print coefficients with tuple representation
            print(f"    Coefficients:")
            for i in range(l + 1):
                if sigma[i] == 0:
                    print(f"      σ_{i} = 0")
                elif sigma[i] == 1:
                    print(f"      σ_{i} = 1")
                else:
                    log_val = gf.log_table[sigma[i]]
                    print(f"      σ_{i} = {sigma[i]} (0x{sigma[i]:0{(gf.m+3)//4}x}, 0b{sigma[i]:0{gf.m}b}) = α^{log_val}")
            
            # Print polynomial expression
            print(f"    σ(X) = ", end="")
            terms = []
            for i in range(l + 1):
                if sigma[i] != 0:
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
        
        # Chien search
        error_locs_in_pattern = chien_search(sigma, l, gf, n)
        
        if verbose:
            print(f"  Chien search roots: {error_locs_in_pattern}")
            if error_locs_in_pattern:
                print(f"  Roots LLR values:")
                for pos in error_locs_in_pattern:
                    llr_idx = n - pos
                    if llr_idx < len(llr_values) and llr_idx >= 1:
                        print(f"    r{pos}: LLR[{llr_idx}] = {llr_values[llr_idx]}")
                    else:
                        print(f"    r{pos}: LLR out of range")
        
        # Check if decoding successful
        if len(error_locs_in_pattern) != l or len(error_locs_in_pattern) > t:
            if verbose:
                print(f"  => DECODING FAILED for this pattern")
            continue
        
        # Correct errors
        corrected = correct_errors(pattern, error_locs_in_pattern)
        
        # 計算相對於原始 r 的錯誤位置
        flipped_set = set(pattern_flipped_bits[idx])
        error_in_pattern_set = set(error_locs_in_pattern)
        
        # XOR: 對稱差集
        error_from_original_r = flipped_set.symmetric_difference(error_in_pattern_set)
        error_from_original_r = sorted(list(error_from_original_r))
        
        if verbose:
            print(f"  => DECODING SUCCESS")
            print(f"     Errors in pattern: {error_locs_in_pattern}")
            print(f"     Errors from original r: {error_from_original_r}")
        
        results.append({
            'pattern_idx': idx,
            'corrected': corrected,
            'flipped_bits': pattern_flipped_bits[idx],
            'errors_in_pattern': error_locs_in_pattern,
            'errors_from_r': error_from_original_r
        })
    
    if verbose:
        print(f"\n" + "-" * 70)
        print(f"Summary: Successfully decoded {len(results)}/{len(test_patterns)} patterns")
    
    if not results:
        if verbose:
            print(f"\n  => DECODING FAILED: No test pattern decoded successfully")
        return False, [], None
    
    # Step 4: Calculate correlation and select best result
    best = select_best_result(results, llr_values, n)
    
    if verbose:
        print(f"\nStep 4: Evaluate Correlation Values")
        print("-" * 70)
        for result in results:
            marker = " <-- BEST" if result == best else ""
            print(f"  Pattern {result['pattern_idx']}: correlation = {result['correlation']}{marker}")
    
    final_error_locations = best['errors_from_r']
    
    if verbose:
        print(f"\nFinal Result:")
        print("-" * 70)
        print(f"  Selected pattern: {best['pattern_idx']}")
        print(f"  Error locations (from original r): {final_error_locations}")
        print(f"\n  => DECODING SUCCESSFUL!")
    
    return True, final_error_locations, best['corrected']


# ============================================================================
# UNIFIED DECODER INTERFACE
# ============================================================================

def decode_bch_unified(r, gf, t, n, mode='hard', llr_values=None, p=2, verbose=True):
    """
    Unified BCH decoder interface supporting both hard and soft decision
    
    Args:
        r: Received polynomial
        gf: Galois Field object
        t: Error correction capability
        n: Codeword length
        mode: 'hard' for hard-decision, 'soft' for soft-decision
        llr_values: LLR values (required for soft-decision)
        p: Number of least reliable bits for Chase algorithm (default: 2)
        verbose: Print detailed output
    
    Returns:
        success: Boolean indicating if decoding succeeded
        error_locations: List of error positions
        corrected: Corrected codeword (or None if failed)
    """
    if mode == 'hard':
        return decode_bch(r, gf, t, n, verbose=verbose)
    elif mode == 'soft':
        if llr_values is None:
            raise ValueError("LLR values are required for soft-decision decoding")
        return decode_bch_soft_decision(r, llr_values, gf, t, n, p=p, verbose=verbose)
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'hard' or 'soft'")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def syndrome_to_str(syndrome, m):
    """Convert syndrome value to various representations"""
    if syndrome == 0:
        return "0"
    return f"{syndrome} (0x{syndrome:0{(m+3)//4}x}, 0b{syndrome:0{m}b})"


def print_polynomial(sigma, l, gf, poly_name="σ(X)"):
    """Print error locator polynomial in a readable format"""
    print(f"  {poly_name} coefficients:")
    for i in range(l + 1):
        if sigma[i] != 0:
            if sigma[i] == 1:
                print(f"    σ_{i} = 1")
            else:
                log_val = gf.log_table[sigma[i]]
                print(f"    σ_{i} = {sigma[i]} = α^{log_val}")
    
    print(f"  {poly_name} = ", end="")
    terms = []
    for i in range(l + 1):
        if sigma[i] != 0:
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


# ============================================================================
# MAIN FUNCTIONS
# ============================================================================

def main(mode='hard', code_type=None):
    """
    Main function to demonstrate BCH decoding
    
    Args:
        mode: 'hard' for hard-decision, 'soft' for soft-decision
        code_type: 1 for (63,51) m=6, 2 for (255,239) m=8, 3 for (1023,983) m=10
                   None for all codes
    """
    print("=" * 70)
    print(f"BCH Decoder Verification Tool - {mode.upper()}-DECISION MODE")
    if code_type is not None:
        params = BCH_PARAMS[code_type]
        print(f"Code: ({params['n']}, {params['k']}), m={params['m']}, t={params['t']}")
    print("=" * 70)
    
    # Example input
    example_input = """0000000001111101011110101111000110000101100011011011011110110111
1001010110100010101110000101011010111010010001110100100101000111
0111101101100000101111110101001101001111100001110110000010100100
1010010101111101101100111001100110011001010111100001110110010110
0111110101001001011001111010011101011000000000100110110101010101
1011101101011110011001010101001010001100011101101001100001011110
0110111001100011010101000100111111000000010110000100000001101010
0010000010110111101111110110111101000111010001100110110001111100"""
    
    lines = example_input.strip().split('\n')
    
    # 決定要處理哪些 code types
    if code_type is None:
        code_types = [1, 2, 3]  # 全部執行
    else:
        code_types = [code_type]  # 只執行指定的
    
    # Test for each code type
    for ct in code_types:
        params = BCH_PARAMS[ct]
        n, k, m, t = params['n'], params['k'], params['m'], params['t']
        
        print(f"\n{'='*70}")
        print(f"BCH Code ({n}, {k}): m={m}, t={t}")
        print(f"Primitive polynomial: {bin(PRIMITIVE_POLYS[m])}")
        print(f"{'='*70}")
        
        # Create Galois Field
        gf = GaloisField(m, PRIMITIVE_POLYS[m])
        
        # Parse input
        r, llr_values = parse_llr_input(lines, ct)
        
        # Print received polynomial info
        error_positions = [i for i in range(n) if r[i] == 1]
        print(f"\nReceived polynomial r(X):")
        print(f"  Number of 1s: {sum(r)}")
        if len(error_positions) <= 20:
            print(f"  Positions with 1: {error_positions}")
        
        # Perform decoding based on mode
        if mode == 'hard':
            success, error_locs, corrected = decode_bch(r, gf, t, n, verbose=True)
        else:  # soft
            success, error_locs, corrected = decode_bch_soft_decision(
                r, llr_values, gf, t, n, p=2, verbose=True
            )
        
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


def test_with_known_errors(mode='hard'):
    """Test with known error patterns"""
    print("\n" + "=" * 70)
    print(f"Testing BCH Decoder with Known Error Patterns - {mode.upper()}-DECISION")
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


def test_soft_decision():
    """Test soft-decision decoding with cases where hard-decision fails"""
    print("\n" + "=" * 70)
    print("Testing Soft-Decision Decoding (Chase Algorithm)")
    print("=" * 70)
    
    for code_type in [1, 2, 3]:
        params = BCH_PARAMS[code_type]
        n, k, m, t = params['n'], params['k'], params['m'], params['t']
        
        print(f"\n--- ({n}, {k}) BCH Code, t={t} ---")
        
        gf = GaloisField(m, PRIMITIVE_POLYS[m])
        
        # Create a scenario with t+1 errors where soft-decision can help
        # The idea: t+1 errors, but one of them is at a "least reliable" position
        print(f"\nTest: {t+1} errors (exceeds hard-decision capability)")
        
        # Create LLR values - most are reliable (high absolute value)
        llr_values = [0]  # Index 0 is don't care
        for i in range(1, n + 1):
            llr_values.append(100)  # High reliability, bit = 0
        
        # Create error pattern with t+1 errors
        error_positions = list(range(t + 1))
        r = [0] * n
        for pos in error_positions:
            r[pos] = 1
            llr_idx = n - pos
            llr_values[llr_idx] = -100  # Flip LLR sign for error bits
        
        # Make one error position "unreliable" (low absolute LLR)
        # This allows Chase algorithm to potentially flip it back
        unreliable_pos = error_positions[0]
        llr_idx = n - unreliable_pos
        llr_values[llr_idx] = -5  # Low absolute value = unreliable
        
        # Hard decision should fail
        print("  Hard-decision decoding:")
        success_hard, _, _ = decode_bch(r, gf, t, n, verbose=False)
        print(f"    Result: {'SUCCESS' if success_hard else 'FAILED (expected)'}")
        
        # Soft decision might succeed
        print("  Soft-decision decoding:")
        success_soft, locs, corrected = decode_bch_soft_decision(
            r, llr_values, gf, t, n, p=2, verbose=False
        )
        print(f"    Result: {'SUCCESS' if success_soft else 'FAILED'}")
        if success_soft:
            print(f"    Error locations found: {locs}")


if __name__ == "__main__":
    import sys
    
    # Default settings
    mode = 'hard'
    code_type = None  # None = 全部執行
    
    # Parse command line arguments
    # Usage: python bch_decoder.py [hard|soft] [1|2|3]
    #        python bch_decoder.py test
    
    args = sys.argv[1:]
    
    for arg in args:
        if arg in ['hard', 'soft']:
            mode = arg
        elif arg in ['1', '2', '3']:
            code_type = int(arg)
        elif arg == 'test':
            # Run tests
            test_with_known_errors('hard')
            test_soft_decision()
            sys.exit(0)
        elif arg in ['-h', '--help', 'help']:
            print(f"Usage: {sys.argv[0]} [mode] [code_type]")
            print()
            print("Mode:")
            print("  hard    Hard-decision decoding (default)")
            print("  soft    Soft-decision decoding (Chase algorithm)")
            print("  test    Run test cases")
            print()
            print("Code Type:")
            print("  1       (63, 51) BCH code, m=6, t=2")
            print("  2       (255, 239) BCH code, m=8, t=2")
            print("  3       (1023, 983) BCH code, m=10, t=4")
            print("  (none)  Run all three codes")
            print()
            print("Examples:")
            print(f"  {sys.argv[0]} hard 1      # Hard-decision, (63,51) only")
            print(f"  {sys.argv[0]} soft 2      # Soft-decision, (255,239) only")
            print(f"  {sys.argv[0]} soft        # Soft-decision, all codes")
            print(f"  {sys.argv[0]} 3           # Hard-decision (default), (1023,983) only")
            sys.exit(0)
        else:
            print(f"Unknown argument: {arg}")
            print(f"Use '{sys.argv[0]} --help' for usage information.")
            sys.exit(1)
    
    # Run main with selected mode and code type
    main(mode=mode, code_type=code_type)