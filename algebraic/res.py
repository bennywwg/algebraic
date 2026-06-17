import sys
from typing import List
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

CHECK_CORRECTNESS = True

def find_index(arr, pred):
    return next((i for i, x in enumerate(arr) if pred(x)), -1)

def inv3_pow2(k):
    x = 1
    for _ in range(k.bit_length()):
        x = x * (2 - 3*x) & ((1 << k) - 1)
    return x

# lowest high bit index
def lhb(x):
    return (x & -x).bit_length() - 1

# run collatz on x, and return the number of shifts done on this step
# a 3x+1 step will return 0 shifts
def col(x):
    if x <= 0:
        return 0, -1
    s = lhb(x)
    if s == 0:
        x = x * 3 + 1
        return x, 0
    else:
        return x >> s, s

# split an integer into two sections, one with the n lowest bits, and one with the other bits, shifted by n
# returns (low_bits, high_bits shifted)
def trim(n, x):
    trimmed = x & ((1 << n) - 1)
    return trimmed, ((x - trimmed) >> n)

class Res:
    val = 0 # full value after expansion
    n = 0 # number of steps
    def __init__(self, n, val):
        self.val = val
        self.n = n
        if CHECK_CORRECTNESS:
            assert val < 3**n
            assert val >= 0
    def __hash__(self):
        return hash((self.n, self.val))
    def __eq__(self, other):
        return self.n == other.n and self.val == other.val
    def __repr__(self):
        return f"{self.val}|{self.n}"
    def __lt__(self, other):
        if self.n == other.n:
            return self.val < other.val
        return self.n < other.n

class Block:
    n = 0
    val = 0
    def __init__(self, n, val):
        self.n = n
        self.val = val
        if CHECK_CORRECTNESS:
            assert val < 2**n
            assert val >= 0
    def __hash__(self):
        return hash((self.n, self.val))
    def __eq__(self, other):
        return self.n == other.n and self.val == other.val
    def __repr__(self):
        as_bin = format(self.val, f"0{self.n}b")
        return f"{as_bin}|{self.n}"
    def __lt__(self, other):
        if self.n == other.n:
            return self.val < other.val
        return self.n < other.n
    def res_cz(self) -> Res:
        res = Res(0, 0)
        curr_val = self.val
        n_remaining = self.n
        while n_remaining > 0:
            curr_val, s = col(curr_val)
            if s == -1:
                break
            curr_val, carried = trim(n_remaining, curr_val)
            if s == 0:
                res.val *= 3
                res.val += carried
                res.n += 1
            n_remaining -= s
        return res
    def res_m(self, in_res: Res):
        self_val = self.val * pow(3, in_res.n)
        block_val, res_val = trim(self.n, self_val + in_res.val)
        return Res(in_res.n, res_val), Block(self.n, block_val)



# enumerate all residual sets up to n values long
def enumerate_all_res(n):
    for length in range(n + 1):
        for val in range(3**length):
            yield Res(length, val)
            
# enumerate all residual lists exactly of length n
def enumerate_all_res_exactly(n):
    for val in range(3**n):
        yield Res(n, val)
 
# enumerate 
def enumerate_all_blocks(n):
    for val in range(1 << n):
        yield Block(n, val)
    


def compute_initial_block(res: Res, b_final: Block):
    m = 2**b_final.n
    inv3 = inv3_pow2(b_final.n)
    base = b_final.val * pow(inv3, res.n, m)
    return Block(b_final.n, (base - res.val) % m)

# Compute all the residual lists up to length max_len that convert the value of b_initial to b_final
def find_to_shortest(b_initial: Block, b_final: Block, max_len):
    final_res = []
    for res_len in range(max_len):
        res_list = make_all_res_exactly(res_len)
        for res in res_list:
            _, out_val = b_initial.res_m(res)
            if out_val == b_final:
                final_res.append(res)
    return final_res

# Compute all the shortest residual lists needed to convert the value of b_initial to b_final
def find_to_shortest_algo(b_initial: Block, b_final: Block):
    m = 2**b_initial.n

    seq_len = 0
    while True:
        pow3 = 3**seq_len
        base_val = b_initial.val * pow3
        seq_extent = pow3 - 1
        diff = (b_final.val - base_val) % m
        if diff <= seq_extent:
            return Res(seq_len, diff)
        seq_len += 1

def expand_solution(bits, expr):
    seen_vars = set()

    def find_next_var(expr):
        pattern = r'(x|xi|m)(\d)(\d)?'
        for match in re.finditer(pattern, expr):
            g3 = match.group(3)
            res = (
                match.start(),
                len(match.group(0)),
                match.group(1),
                int(match.group(2)),
                int(g3) if g3 is not None else None,
            )
            if res[2] == "x" and res[4] == 1:
                seen_vars.add(match.group(0))
                continue
            return res
        return None

    original_expr = expr
    while True:
        next_var = find_next_var(expr)
        if next_var is None:
            break

        replacement = None
        if next_var[2] == 'm':
            myt = next_var[3]
            replacement = f"3**CzM(x{myt}{myt})"

        if next_var[2] == 'x':
            myindex = next_var[3]
            myt = next_var[4]
            if myindex < myt:
                raise ValueError(f"Invalid variable x{myindex}{myt}")
            replacement = f"(xi{next_var[3]}{next_var[4]-1} % {2**bits})"

        if next_var[2] == 'xi':
            myindex = next_var[3]
            myt = next_var[4]
            if myindex == myt:
                replacement = f"(Cz(x{myindex}{myt}) << {bits})"
            elif myindex > myt:
                replacement = f"(x{myindex}{myt} * m{myt} + (xi{myindex-1}{myt} >> {bits}))"
            else:
                raise ValueError(f"Invalid variable xi{myindex}{myt}")

        if replacement is None:
            break

        if replacement is not None:
            start = next_var[0]
            length = next_var[1]
            expr = expr[:start] + replacement + expr[start + length:]
    
    def Cz(x):
        return Block(bits, x).res_cz().val
    
    def CzM(x):
        return Block(bits, x).res_cz().n

    seen_vars_list = sorted(seen_vars)

    lambda_str = "lambda " + ", ".join(seen_vars_list) + ": " + expr

    print(f"{original_expr} = {expr}")
    print("lambda form:", lambda_str)

    return eval(f"{lambda_str}", {"Cz": Cz, "CzM": CzM})


def solution(bits, x1):
    x = []
    xi = []

    def row(arr, time):
        while time - 1 >= len(arr):
            arr.append([])
        return arr[time - 1]

    def get(row_val, index):
        while index - 1 >= len(row_val):
            row_val.append(0)
        return row_val[index - 1]
    
    def set(row_val, index, val):
        get(row_val, index)
        row_val[index - 1] = val
    
    def randval():
        return 0#random.getrandbits(bits)

    set(row(x, 1), 1, x1)

    def compute_intermediates(t):
        xt = row(x, t) # this row
        Mt = 3**Block(bits, get(xt, t)).res_cz().n
        # first compute the intermediate xi values
        index  = t
        while True:
            val = None
            if index == t:
                val = Block(bits, get(xt, index)).res_cz().val * 2**bits
            else:
                val = get(row(x, t), index) * Mt + (get(row(xi, t), index - 1) >> bits)
            if index >= len(xt) and val == 0:
                break
            set(row(xi, t), index, val)
            index += 1

            
    def compute_actuals(t):
        m = 2**bits
        xti = row(xi, t - 1)
        xt = row(x, t)
        index = t
        while True:
            low, high = trim(bits, get(xti, index))
            set(xt, index, low)
            index += 1
            if index >= len(xt) and high == 0:
                break

    for t in range(1, 10):
        print(f"--- TIME {t} ---")
        xrow = row(x, t)
        xirow = row(xi, t)

        compute_intermediates(t)
        compute_actuals(t + 1)

        for i in range(t + 5, t - 1, -1):
            print(f"X{i}{t} ={get(xrow, i)}".ljust(14), end="")
        print()
        for i in range(t + 5, t - 1, -1):
            print(f"Xi{i}{t}={get(xirow, i)}".ljust(14), end="")
        print("\n")


    #while True:
    #    print(" Ÿ ")
    


def main():
    if sys.argv[1] == "r":
        n = int(sys.argv[2])
        points = []
        for block in enumerate_all_blocks(n):
            res = block.res_cz()
            points.append((block.val, res.val, res.val / max(block.val, 1)))

        # sort the points by y, leave only the top t, then make a bar plot labeling the x vals
        if False:
            p = 20
            #points = sorted(points, key=lambda x: x[1], reverse=True)[:p]
            for point in points:
                # print point[1] in binary, left padded to 64 bits, and using □ as 0 and ■ as 1
                res_label = ''.join(['■' if c == '1' else '□' for c in format(point[0], f"0{n}b")])
                res_label = res_label.rjust(24, '□')
                print(f"block {point[0]} -> {res_label} -> {point[1]}")

        xs = [x[0] for x in points]
        ys = [y[2] for y in points]
        plt.bar(range(len(xs)), ys)
        plt.xlabel("block value")
        plt.ylabel("residual value / (block value or 1)")
        plt.show()
        return
    
    if sys.argv[1] == "e":
        bits = int(sys.argv[2])
        target_val = int(sys.argv[3])
        m = 2**bits
        r_fn = expand_solution(bits, "x44")
        g_fn = expand_solution(bits, "x55")
        b_fn = expand_solution(bits, "x66")

        fig, ax = plt.subplots()

        img = np.zeros((m, m, 3), dtype=np.float32)
        im = ax.imshow(img)
        ax.axis("off")

        def update(frame):
            target_val = frame
            r = np.fromfunction(np.vectorize(lambda i, j: 1 if (r_fn(target_val, int(i), int(j), 0) == target_val) else 0), (m, m), dtype=int)
            g = np.fromfunction(np.vectorize(lambda i, j: 1 if (g_fn(target_val, int(i), int(j), 0, 0) == int(i)) else 0), (m, m), dtype=int)
            b = np.fromfunction(np.vectorize(lambda i, j: 1 if (b_fn(target_val, int(i), int(j), 0, 0, 0) == int(j)) else 0), (m, m), dtype=int)
            im.set_data(np.stack((r, g, b), axis=-1).astype(np.float32))
            return (im,)
        
        FuncAnimation(fig, update, frames=range(m), interval=1, blit=True)
        plt.show()
        return
    
    if sys.argv[1] == "s":
        bits = int(sys.argv[2])
        x1 = int(sys.argv[3])
        solution(bits, x1)
        return
    
    if sys.argv[1] == "ft":
        n = int(sys.argv[2])
        from_block = Block(n, int(sys.argv[3]))
        to_block = Block(n, int(sys.argv[4]))
        print(f"computing residuals needed to convert block {from_block} to {to_block}")
        final_res = find_to_shortest_algo(from_block, to_block)[0]
        
        print(f"{final_res} meets condition")
        tmp_block = from_block
        for res in (final_res):
            new_val = tmp_block.res_m([res])[1]
            print(f"{tmp_block} <- {res} = {new_val}")
            tmp_block = new_val
        return
    
    if sys.argv[1] == "aft":
        n = int(sys.argv[2])
        m = 2**n

        def from_to(f, t):
            #shortest = find_to_shortest(Block(n, f), Block(n, t), 7)
            shortest = find_to_shortest_algo(Block(n, f), Block(n, t))
            return shortest.val

        a = np.fromfunction(np.vectorize(lambda i, j: from_to(i, j)), (m, m), dtype=int)

        if True:
            plt.imshow(a)
            plt.xlabel("to")
            plt.ylabel("from")
            plt.colorbar()
            plt.show()

        # plot the sum of each column and row
        if False:
            row_sums = np.sum(a, axis=1)
            col_sums = np.sum(a, axis=0)
            plt.plot(range(m), row_sums / m, label="from")
            plt.plot(range(m), col_sums / m, label="to")
            plt.xlabel("block value")
            plt.ylabel("avg of residual lengths")
            plt.legend()
            plt.show()

        return
    
    # linear from to, shows how many residuals it takes to get to a single value from all possible values
    # 2nd param is the target values
    if sys.argv[1] == "lft":
        n = int(sys.argv[2])
        target = int(sys.argv[3])
        xs = []
        ys = []
        for f in range(2**n):
            res = find_to_shortest_algo(Block(n, f), Block(n, target))
            xs.append(f)
            ys.append(res.val)
            print(f"{Block(n, f)} -> {Block(n, target)} via {res}")
        # curve
        plt.plot(xs, ys)
        plt.show()
        return
    
    n = int(sys.argv[1])
    counts = {}

    if len(sys.argv) <= 2:
        for val in range(2**n):
            b = Block(n, val)
            res = b.res_cz()
            counts[res_key(res)] = counts.get(res_key(res), 0) + 1
            print(f"{b} -> {res}")
        xs = sorted(counts.keys())
        ys = [counts[x] for x in xs]
        #plt.bar(xs, ys)
        #plt.show()
        print(f"Coverage = {len(counts)} / {2**n}, {round(100 * len(counts) / (2**n), 2)}% inputs, {round(100 * len(counts) / (3**n), 2)}% residuals")
    else:
        val = int(sys.argv[2])
        b = Block(n, val)
        res = b.res_cz()
        print(res)


if __name__ == "__main__":
    main()
