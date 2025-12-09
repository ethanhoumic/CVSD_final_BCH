import os
import random
from typing import List, Optional, Tuple

# ------------------------------------------------------------
# BCH 規格表：m -> (m, n, k, t, g(x) 的非零次方)
# ------------------------------------------------------------
BCH_SPECS = {
    6: {
        "m": 6,
        "n": 63,
        "k": 51,
        "t": 2,
        # g(X) = 1 + X^3 + X^4 + X^5 + X^8 + X^10 + X^12
        "g_exponents": [0, 3, 4, 5, 8, 10, 12],
    },
    8: {
        "m": 8,
        "n": 255,
        "k": 239,
        "t": 2,
        # g(X) = 1 + X + X^5 + X^6 + X^8 + X^9 + X^10
        #        + X^11 + X^13 + X^14 + X^16
        "g_exponents": [0, 1, 5, 6, 8, 9, 10, 11, 13, 14, 16],
    },
    10: {
        "m": 10,
        "n": 1023,
        "k": 983,
        "t": 4,
        # g(X) = 1 + X + X^3 + X^4 + X^7 + X^9 + X^10 + X^11 + X^12
        #      + X^16 + X^19 + X^21 + X^22 + X^23 + X^24 + X^25
        #      + X^27 + X^29 + X^30 + X^31 + X^33 + X^39 + X^40
        "g_exponents": [
            0, 1, 3, 4, 7, 9, 10, 11, 12,
            16, 19, 21, 22, 23, 24, 25,
            27, 29, 30, 31, 33, 39, 40
        ],
    },
}

# 對應到 pcode.txt 的兩位碼
PCODE_MAP = {
    6:  "01",
    8:  "10",
    10: "11",
}


# ------------------------------------------------------------
# 多項式工具：用 int 表示 GF(2) 多項式
# ------------------------------------------------------------
def build_gen_bits(g_exponents: List[int]) -> List[int]:
    """由非零次方列表建立 generator polynomial 的係數 bits[degree]."""
    deg = max(g_exponents)
    bits = [0] * (deg + 1)
    for e in g_exponents:
        bits[e] = 1
    return bits


def poly_mod_int(dividend: int, divisor: int) -> int:
    """GF(2) 多項式整除：回傳 dividend / divisor 的餘數（int 型態）。"""
    if divisor == 0:
        raise ValueError("divisor polynomial is zero")
    deg_div = divisor.bit_length() - 1
    r = dividend
    while r and (r.bit_length() - 1) >= deg_div:
        shift = (r.bit_length() - 1) - deg_div
        r ^= divisor << shift
    return r


# ------------------------------------------------------------
# BCH 編碼 (systematic)：m(x) * x^(n-k) mod g(x)
# ------------------------------------------------------------
def encode_bch(msg_bits: List[int], spec: dict) -> List[int]:
    """
    msg_bits: 長度 k，index i 對應 x^i 的係數
    回傳 codeword bits: 長度 n，index j 對應 x^j 的係數
    """
    n = spec["n"]
    k = spec["k"]
    g_bits = build_gen_bits(spec["g_exponents"])

    # message m(x)
    m_int = sum((bit & 1) << i for i, bit in enumerate(msg_bits))
    # generator g(x)
    g_int = sum((bit & 1) << i for i, bit in enumerate(g_bits))

    # systematic encoding: c(x) = m(x) * x^(n-k) + r(x)
    shifted = m_int << (n - k)
    rem = poly_mod_int(shifted, g_int)
    cw_int = shifted ^ rem

    cw_bits = [(cw_int >> i) & 1 for i in range(n)]
    return cw_bits


# ------------------------------------------------------------
# LLR byte 產生
# ------------------------------------------------------------
def llr_byte_from_bit_hard(bit: int) -> str:
    """
    硬判決版本（帶隨機 magnitude）：
      bit = 0 -> 正值或 0（MSB=0）
      bit = 1 -> 負值      （MSB=1）
    """
    if bit == 0:
        # MSB=0 → signed: 0 ~ +127
        mag = random.randint(0, 127)
        val = mag
    else:
        # MSB=1 → signed: -1 ~ -128  (完整 signed 8-bit)
        mag = random.randint(1, 128)     # absolute value
        val = (-mag) & 0xFF              # 二補數
    return f"{val:08b}"



def llr_byte_from_bit_soft(bit: int, mag: int) -> str:
    """
    soft 版本：
      bit = 0 -> +mag   (0xxxxxxx)，mag 會被限制在 0~127
      bit = 1 -> -mag   (1xxxxxxx，二補數)，mag 可到 128 → -128
    """
    assert 0 <= mag <= 127

    if bit == 0:
        # 正值：避免 mag=128 變成 1000_0000 (看起來像負數)
        if mag > 127:
            mag = 127
        val = mag                  # 0..127, MSB=0
    else:
        # 負值：允許 mag=1..128 → -1..-127
        if mag == 0:
            mag = 1
        elif mag > 127:
            mag = 127
        val = (-mag) & 0xFF        # 8-bit 二補數，MSB=1
    return f"{val:08b}"



def llr_rows_from_llr_array(llr_array: List[str], m: int) -> List[str]:
    """
    llr_array[e]: X^e 的 LLR (字串 8 bits)，e = 0..2^m-1
    輸出：每列 64 bits，順序從 X^(2^m-1) -> X^0，每 8 個 LLR 一列。
    """
    M = 1 << m
    assert len(llr_array) == M

    rows: List[str] = []
    cur: List[str] = []

    # 從最高次 X^(2^m-1) 到最低次 X^0
    for e in range(M - 1, -1, -1):
        cur.append(llr_array[e])
        if len(cur) == 8:
            rows.append("".join(cur))
            cur = []

    assert not cur  # 剛好切整除
    return rows


# ------------------------------------------------------------
# hard / soft 版本的 codeword -> LLR rows
# ------------------------------------------------------------
def cw_bits_to_llr_rows_hard(cw_bits: List[int], m: int) -> List[str]:
    """
    硬判決：完全照你原本的 p100.txt 格式，只是自動產生。
    """
    n = (1 << m) - 1
    assert len(cw_bits) == n
    M = 1 << m

    llr_array = ["00000000"] * M

    # X^0..X^(n-1) 的係數
    for e in range(n):
        llr_array[e] = llr_byte_from_bit_hard(cw_bits[e])

    # dummy X^n：係數固定 0，LLR 也要固定 00000000
    llr_array[n] = "00000000"

    return llr_rows_from_llr_array(llr_array, m)


def cw_bits_to_llr_rows_soft(
    cw_bits: List[int],
    mags: List[int],
    m: int,
) -> List[str]:
    """
    soft 判決：LLR 由 cw_bits 的 sign + mags 的絕對值組合而成。
    mags[e] 是 X^e 的 |LLR| (0..127)，包含 dummy 那一個。
    """
    n = (1 << m) - 1
    assert len(cw_bits) == n
    M = 1 << m
    assert len(mags) == M  # 包含 dummy

    llr_array = ["00000000"] * M

    for e in range(n):
        bit = cw_bits[e]
        mag = mags[e]
        llr_array[e] = llr_byte_from_bit_soft(bit, mag)

    # dummy bit X^n：係數無意義且 LLR 必須是 00000000
    llr_array[n] = "00000000"

    return llr_rows_from_llr_array(llr_array, m)


# ------------------------------------------------------------
# 產生「單一 codeword」的 p_rows / pa_lines
# ------------------------------------------------------------
def generate_single_bch_codeword(
    m: int,
    soft: bool,
) -> Tuple[List[str], List[str]]:
    """
    產生一個 codeword 的測資。

    回傳:
      p_rows_one : 該 codeword 對應的 p.txt 多行（每行 64 bits）
      pa_one     : 該 codeword 對應的 pa.txt 多行（每行 10 bits）
    """
    if m not in BCH_SPECS:
        raise ValueError(f"Unsupported m={m}")

    spec = BCH_SPECS[m]
    n, k, t = spec["n"], spec["k"], spec["t"]

    # 1. 隨機產生 message bits
    msg_bits = [random.randint(0, 1) for _ in range(k)]

    # 2. BCH systematic encoding -> 真實 codeword
    cw_clean = encode_bch(msg_bits, spec)

    pa_lines_one: List[str] = []

    if not soft:
        # -------- Hard decision 測資：最多 t 個錯誤 --------
        num_err = random.randint(0, t)
        cw_bits = cw_clean[:]
        error_positions: List[int] = []

        if num_err > 0:
            error_positions = random.sample(range(n), num_err)
            for pos in error_positions:
                cw_bits[pos] ^= 1

            for pos in sorted(error_positions):
                pa_lines_one.append(f"{pos:010b}")
        else:
            pa_lines_one.append(f"{1023:010b}")

        p_rows_one = cw_bits_to_llr_rows_hard(cw_bits, m)
        return p_rows_one, pa_lines_one

    # -------- Soft decision 測資：最多 t+2 錯誤 --------
    cw_bits = cw_clean[:]
    M = 1 << m

    # 先挑兩個「最不可靠」的位置：硬體會先挑這兩個去 flip
    s1, s2 = random.sample(range(n), 2)
    soft_candidates = [s1, s2]

    # 初始化每個位置的 |LLR|：
    # 先全部給「大範圍」[64, 127]，再覆蓋特殊的
    mags = [0] * M
    for e in range(n):
        mags[e] = random.randint(64, 127)
    # dummy bit X^n：係數無意義，mags[n] 也不用理，最後 LLR 會是 00000000
    mags[n] = 0

    # 這兩個 soft 候選：給最小 |LLR| 範圍 [1, 16]
    for s in soft_candidates:
        mags[s] = random.randint(1, 16)

    # 先對 soft_candidates 決定要不要翻（最多 2 個錯）
    soft_err = set()
    for s in soft_candidates:
        if random.choice([True, False]):
            cw_bits[s] ^= 1
            soft_err.add(s)

    # 剩下 bits 再選 0..t 個錯誤
    remaining_indices = [i for i in range(n) if i not in soft_candidates]
    num_err_hard = random.randint(0, t)
    hard_err = set()
    if num_err_hard > 0 and remaining_indices:
        chosen = random.sample(
            remaining_indices,
            min(num_err_hard, len(remaining_indices)),
        )
        for pos in chosen:
            cw_bits[pos] ^= 1
            hard_err.add(pos)

        # 這些「剩下翻錯的 bit」：給第二小 |LLR| 範圍 [17, 32]
        for pos in chosen:
            mags[pos] = random.randint(17, 32)

    # 實際錯誤位置 = soft_err ∪ hard_err  (最多 t+2 個)
    error_positions = sorted(soft_err.union(hard_err))

    if error_positions:
        for pos in error_positions:
            pa_lines_one.append(f"{pos:010b}")
    else:
        pa_lines_one.append(f"{1023:010b}")

    # 用 cw_bits + mags 產生 LLR rows
    p_rows_one = cw_bits_to_llr_rows_soft(cw_bits, mags, m)
    return p_rows_one, pa_lines_one


# ------------------------------------------------------------
# 固定一種 (m, soft/hard) 的版本（舊功能，保留）
# ------------------------------------------------------------
def write_bch_test_files(
    m: int,
    num_codewords: int,
    p_filename: str = "p.txt",
    pa_filename: str = "pa.txt",
    seed: Optional[int] = None,
    soft: bool = False,
) -> None:
    """
    產生只含單一 (m, mode) 的測資，寫到 pattern/ 裡：
      p.txt  / pa.txt

    soft = False -> hard decision 測資 (最多 t 個錯)
    soft = True  -> soft decision 測資 (最多 t+2 個錯)
    """
    os.makedirs("pattern", exist_ok=True)
    if seed is not None:
        random.seed(seed)

    all_p_rows: List[str] = []
    all_pa_lines: List[str] = []

    for _ in range(num_codewords):
        p_rows_one, pa_one = generate_single_bch_codeword(m, soft)
        all_p_rows.extend(p_rows_one)
        all_pa_lines.extend(pa_one)

    with open(os.path.join("pattern", p_filename), "w") as f:
        f.write("\n".join(all_p_rows))
        f.write("\n")

    with open(os.path.join("pattern", pa_filename), "w") as f:
        f.write("\n".join(all_pa_lines))
        f.write("\n")


# ------------------------------------------------------------
# ✅ 混和測資版本：每個 codeword 隨機選 (m, mode)
# ------------------------------------------------------------
def write_bch_mixed_test_files(
    num_codewords: int,
    p_filename: str = "p.txt",
    pa_filename: str = "pa.txt",
    pmode_filename: str = "pmode.txt",
    pcode_filename: str = "pcode.txt",
    seed: Optional[int] = None,
) -> None:
    """
    產生「混和規格」的測資：
      - 每個 codeword 有自己的 (m, mode)
      - mode: 0 = hard, 1 = soft      -> 寫到 pmode.txt
      - m   : 01=6, 10=8, 11=10       -> 寫到 pcode.txt

    p.txt  : 所有 codeword 的 LLR rows 串在一起
    pa.txt : 所有 codeword 的錯誤位置 (10 bits) 串在一起
    """
    os.makedirs("pattern", exist_ok=True)
    if seed is not None:
        random.seed(seed)

    all_p_rows: List[str] = []
    all_pa_lines: List[str] = []
    pmode_lines: List[str] = []
    pcode_lines: List[str] = []

    for _ in range(num_codewords):
        # 隨機選 mode: 0=hard, 1=soft
        soft = bool(random.getrandbits(1))
        pmode_lines.append("1" if soft else "0")

        # 隨機選 m ∈ {6,8,10}
        m_choice = random.choice([6, 8, 10])
        pcode_lines.append(PCODE_MAP[m_choice])

        # 產生這個 codeword 的測資
        p_rows_one, pa_one = generate_single_bch_codeword(m_choice, soft)

        all_p_rows.extend(p_rows_one)
        all_pa_lines.extend(pa_one)

    # 寫檔
    with open(os.path.join("pattern", p_filename), "w") as f:
        f.write("\n".join(all_p_rows))
        f.write("\n")

    with open(os.path.join("pattern", pa_filename), "w") as f:
        f.write("\n".join(all_pa_lines))
        f.write("\n")

    with open(os.path.join("pattern", pmode_filename), "w") as f:
        f.write("\n".join(pmode_lines))
        f.write("\n")

    with open(os.path.join("pattern", pcode_filename), "w") as f:
        f.write("\n".join(pcode_lines))
        f.write("\n")


# ------------------------------------------------------------
# 範例：直接執行時
# ------------------------------------------------------------
if __name__ == "__main__":
    # 例 1：單一規格 
    write_bch_test_files(
        m=6,
        num_codewords=50,
        p_filename="p100.txt",
        pa_filename="p100a.txt",
        seed=77,
        soft=False,
    )
    
    write_bch_test_files(
        m=8,
        num_codewords=50,
        p_filename="p200.txt",
        pa_filename="p200a.txt",
        seed=9,
        soft=False,
    )
    
    write_bch_test_files(
        m=10,
        num_codewords=50,
        p_filename="p300.txt",
        pa_filename="p300a.txt",
        seed=7,
        soft=False,
    )
    
    write_bch_test_files(
        m=6,
        num_codewords=50,
        p_filename="p400.txt",
        pa_filename="p400a.txt",
        seed=777,
        soft=True,
    )
    
    write_bch_test_files(
        m=8,
        num_codewords=50,
        p_filename="p500.txt",
        pa_filename="p500a.txt",
        seed=79,
        soft=True,
    )
    
    write_bch_test_files(
        m=10,
        num_codewords=50,
        p_filename="p600.txt",
        pa_filename="p600a.txt",
        seed=779,
        soft=True,
    )

    # write_bch_test_files(
    #     m=10,
    #     num_codewords=50,
    #     p_filename="p600.txt",
    #     pa_filename="p600a.txt",
    #     seed=124,
    #     soft=True,
    # )

    # 例 2：混和規格 (hard/soft x m=6/8/10)
    # write_bch_mixed_test_files(
    #     num_codewords=50,
    #     p_filename="p_mix.txt",
    #     pa_filename="pa_mix.txt",
    #     pmode_filename="pmode_mix.txt",
    #     pcode_filename="pcode_mix.txt",
    #     seed=33,
    # )
