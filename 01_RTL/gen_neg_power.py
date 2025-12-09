def get_alpha_inverse_poly(m, i):
    """
    計算 alpha^(-i) 在 GF(2^m) 中的多項式表示法。
    支援 m = 6, 8, 10 (對應您的 BCH Project 規格)。
    """
    
    # 定義本原多項式 (Primitive Polynomials)
    # 這裡只儲存最高次項以外的部分 (Mask)
    # 例如 m=6, p(x) = 1 + x + x^6 => x^6 = 1 + x (binary: 000011 -> 0x03)
    polys = {
        6: 0x03,   # 1 + x + x^6
        8: 0x1D,   # 1 + x^2 + x^3 + x^4 + x^8
        10: 0x09   # 1 + x^3 + x^10
    }
    
    if m not in polys:
        return f"錯誤: 目前不支援 m={m}，請手動新增該 Field 的本原多項式。"
    
    prim_poly = polys[m]
    # 有限域的週期 (Order) 是 2^m - 1
    period = (1 << m) - 1
    
    # 計算 -i 在模 (2^m - 1) 下的等效正指數
    # alpha^(-i) = alpha^(period - i)
    target_exp = (-i) % period
    
    # 開始計算 alpha^target_exp
    # reg 代表目前的數值，初始化為 1 (即 alpha^0)
    # 我們使用 LFSR 的方式模擬乘法
    reg = 1
    
    for _ in range(target_exp):
        # 檢查最高位 (MSB) 是否為 1
        is_msb_set = (reg >> (m - 1)) & 1
        
        # 左移一位 (相當於乘 alpha)
        reg = (reg << 1) & period # 確保不超過 m bits
        
        # 如果原本最高位是 1，則需要 XOR 本原多項式 (模運算)
        if is_msb_set:
            reg ^= prim_poly
            
    # 將結果轉換為多項式字串格式回傳
    return format_poly(reg)

def format_poly(val):
    """將整數轉換為多項式字串表示法 (例如: alpha^5 + alpha^2 + 1)"""
    if val == 0:
        return "0"
    
    terms = []
    # 檢查每一個 bit
    # 假設我們不超過 20 bits，這裡寫個迴圈檢查
    for k in range(20, -1, -1):
        if (val >> k) & 1:
            if k == 0:
                terms.append("1")
            elif k == 1:
                terms.append("alpha")
            else:
                terms.append(f"alpha^{k}")
                
    return " + ".join(terms)

# --- 使用範例 ---

m = 6
i = 16
print(f"GF(2^{m}), alpha^(-{i}) = {get_alpha_inverse_poly(m, i)}")

m = 6
i = 24
print(f"GF(2^{m}), alpha^(-{i}) = {get_alpha_inverse_poly(m, i)}")

m = 6
i = 32
print(f"GF(2^{m}), alpha^(-{i}) = {get_alpha_inverse_poly(m, i)}")
