#!/usr/bin/env python3
"""
GF(2^m) Polynomial to Power Representation Converter

Input: A polynomial representation (binary string) and field parameter m
Output: The power representation (alpha^k)

Polynomial representation: bits represent coefficients of 1, α, α², α³, ...
  e.g., "000011" = 1 + α (LSB is coefficient of 1, next bit is coefficient of α)
  
Power representation: α^k for some k
  e.g., α^6

Note: In GF(2^m), every non-zero element can be expressed as α^k for some k.
"""

# Primitive polynomials for each field
PRIMITIVE_POLYS = {
    6: 0b1000011,      # 1 + X + X^6
    8: 0b100011101,    # 1 + X^2 + X^3 + X^4 + X^8
    10: 0b10000001001  # 1 + X^3 + X^10
}


class GaloisField:
    """Galois Field GF(2^m) arithmetic"""
    
    def __init__(self, m, primitive_poly=None):
        self.m = m
        if primitive_poly is None:
            if m in PRIMITIVE_POLYS:
                primitive_poly = PRIMITIVE_POLYS[m]
            else:
                raise ValueError(f"No default primitive polynomial for m={m}")
        self.primitive_poly = primitive_poly
        self.n = (1 << m) - 1  # 2^m - 1
        
        # Build lookup tables
        self.exp_table = [0] * (2 * self.n + 1)  # alpha^i -> polynomial form
        self.log_table = [0] * (self.n + 1)      # polynomial form -> exponent
        
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
    
    def poly_to_power(self, poly_value):
        """
        Convert polynomial representation to power representation
        
        Args:
            poly_value: integer representing the polynomial
                        (bit i = coefficient of α^i)
        
        Returns:
            Power k such that poly_value = α^k, or None if poly_value is 0
        """
        if poly_value == 0:
            return None  # 0 has no power representation
        
        if poly_value > self.n:
            raise ValueError(f"Value {poly_value} exceeds field size {self.n}")
        
        return self.log_table[poly_value]
    
    def power_to_poly(self, power):
        """
        Convert power representation to polynomial representation
        
        Args:
            power: the exponent k in α^k
        
        Returns:
            Integer representing the polynomial form
        """
        power = power % self.n
        return self.exp_table[power]
    
    def poly_str(self, poly_value):
        """Convert polynomial value to readable string like '1 + α + α^3'"""
        if poly_value == 0:
            return "0"
        
        terms = []
        for i in range(self.m):
            if poly_value & (1 << i):
                if i == 0:
                    terms.append("1")
                elif i == 1:
                    terms.append("α")
                else:
                    terms.append(f"α^{i}")
        
        return " + ".join(terms)


def parse_poly_input(poly_str):
    """
    Parse polynomial input string to integer value
    
    Supports formats:
    - Binary string: "000011", "0b000011"
    - Hex string: "0x03"
    - Decimal: "3" or just the number
    """
    poly_str = poly_str.strip()
    
    if poly_str.startswith("0x") or poly_str.startswith("0X"):
        return int(poly_str, 16)
    elif poly_str.startswith("0b") or poly_str.startswith("0B"):
        return int(poly_str, 2)
    elif all(c in '01' for c in poly_str) and len(poly_str) > 1:
        # Assume binary if only 0s and 1s and length > 1
        return int(poly_str, 2)
    else:
        return int(poly_str)


def poly_to_power(poly_input, m):
    """
    Main function: Convert polynomial representation to power representation
    
    Args:
        poly_input: polynomial as binary string, hex string, or integer
        m: field parameter (GF(2^m))
    
    Returns:
        String describing the power representation
    """
    if m not in PRIMITIVE_POLYS:
        raise ValueError(f"Unsupported field parameter m={m}. Supported: {list(PRIMITIVE_POLYS.keys())}")
    
    gf = GaloisField(m)
    
    # Parse input
    if isinstance(poly_input, str):
        poly_value = parse_poly_input(poly_input)
    else:
        poly_value = int(poly_input)
    
    # Convert
    power = gf.poly_to_power(poly_value)
    
    if power is None:
        return "0 (no power representation)"
    else:
        return f"α^{power}"


def main():
    """Interactive mode and examples"""
    print("=" * 60)
    print("GF(2^m) Polynomial to Power Representation Converter")
    print("=" * 60)
    
    # Examples
    examples = [
        # (poly_binary, m, description)
        ("000011", 6, "1 + α in GF(2^6)"),
        ("000010", 6, "α in GF(2^6)"),
        ("000001", 6, "1 in GF(2^6)"),
        ("100000", 6, "α^5 in GF(2^6)"),
        ("110000", 6, "α^5 + α^4 in GF(2^6)"),
        ("00000010", 8, "α in GF(2^8)"),
        ("00011101", 8, "1 + α^2 + α^3 + α^4 in GF(2^8)"),
        ("0000000010", 10, "α in GF(2^10)"),
    ]
    
    print("\nExamples:")
    print("-" * 60)
    
    for poly_str, m, desc in examples:
        gf = GaloisField(m)
        poly_value = int(poly_str, 2)
        power = gf.poly_to_power(poly_value)
        poly_readable = gf.poly_str(poly_value)
        
        if power is None:
            power_str = "0"
        else:
            power_str = f"α^{power}"
        
        print(f"Input: {poly_str} (m={m})")
        print(f"  Polynomial: {poly_readable}")
        print(f"  Power:      {power_str}")
        print()
    
    # Print element tables for reference
    print("\n" + "=" * 60)
    print("GF(2^6) Element Table (first 20 elements)")
    print("=" * 60)
    gf6 = GaloisField(6)
    print(f"{'Power':<12} {'Decimal':<10} {'Binary':<10} {'Polynomial':<20}")
    print("-" * 52)
    print(f"{'0':<12} {'0':<10} {'000000':<10} {'0':<20}")
    for i in range(20):
        poly = gf6.exp_table[i]
        poly_str = gf6.poly_str(poly)
        print(f"α^{i:<10} {poly:<10} {poly:06b}     {poly_str:<20}")
    
    print("\n" + "=" * 60)
    print("GF(2^8) Element Table (first 20 elements)")
    print("=" * 60)
    gf8 = GaloisField(8)
    print(f"{'Power':<12} {'Decimal':<10} {'Binary':<12} {'Polynomial':<25}")
    print("-" * 60)
    print(f"{'0':<12} {'0':<10} {'00000000':<12} {'0':<25}")
    for i in range(20):
        poly = gf8.exp_table[i]
        poly_str = gf8.poly_str(poly)
        print(f"α^{i:<10} {poly:<10} {poly:08b}     {poly_str:<25}")
    
    print("\n" + "=" * 60)
    print("GF(2^10) Element Table (first 20 elements)")
    print("=" * 60)
    gf10 = GaloisField(10)
    print(f"{'Power':<12} {'Decimal':<10} {'Binary':<14} {'Polynomial':<30}")
    print("-" * 66)
    print(f"{'0':<12} {'0':<10} {'0000000000':<14} {'0':<30}")
    for i in range(20):
        poly = gf10.exp_table[i]
        poly_str = gf10.poly_str(poly)
        print(f"α^{i:<10} {poly:<10} {poly:010b}     {poly_str:<30}")


def interactive():
    """Run interactive conversion"""
    print("\n" + "=" * 60)
    print("Interactive Mode")
    print("=" * 60)
    print("Enter 'q' to quit")
    
    while True:
        try:
            m_str = input("\nEnter m (6, 8, or 10): ").strip()
            if m_str.lower() == 'q':
                break
            m = int(m_str)
            
            poly_str = input("Enter polynomial (binary, hex, or decimal): ").strip()
            if poly_str.lower() == 'q':
                break
            
            result = poly_to_power(poly_str, m)
            
            gf = GaloisField(m)
            poly_value = parse_poly_input(poly_str)
            poly_readable = gf.poly_str(poly_value)
            
            print(f"\nInput:      {poly_str}")
            print(f"Decimal:    {poly_value}")
            print(f"Binary:     {poly_value:0{m}b}")
            print(f"Polynomial: {poly_readable}")
            print(f"Power:      {result}")
            
        except ValueError as e:
            print(f"Error: {e}")
        except EOFError:
            break


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) == 3:
        # Command line mode: python script.py <poly> <m>
        poly_input = sys.argv[1]
        m = int(sys.argv[2])
        
        gf = GaloisField(m)
        poly_value = parse_poly_input(poly_input)
        power = gf.poly_to_power(poly_value)
        poly_readable = gf.poly_str(poly_value)
        
        print(f"Input:      {poly_input}")
        print(f"Decimal:    {poly_value}")
        print(f"Binary:     {poly_value:0{m}b}")
        print(f"Polynomial: {poly_readable}")
        if power is None:
            print(f"Power:      0 (no power representation)")
        else:
            print(f"Power:      α^{power}")
    
    elif len(sys.argv) == 2 and sys.argv[1] == "-i":
        # Interactive mode
        interactive()
    
    else:
        # Default: show examples
        main()