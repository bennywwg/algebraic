import sys
from typing import List
import matplotlib.pyplot as plt

def find_index(arr, pred):
    return next((i for i, x in enumerate(arr) if pred(x)), -1)

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
def trim(n, x):
    trimmed = x & ((1 << n) - 1)
    return trimmed, (x - trimmed) >> n

# convert a residual array into a unique number
def res_key(res):
    val = 0
    for r in reversed(res):
        val <<= 3
        val += (r + 4)
    return val

# make a residual from a number, the inverse of the res_key function
def make_res(val):
    res = []
    while val != 0:
        res.append(val & 3)
        val >>= 3
    return res

class Res:
    val = 0
    owner = None
    def __init__(self, val_key, owner):
        self.val = val_key
        self.owner = owner
        if not owner is None:
            found_index = find_index(self.root_residuals, lambda res: res_key(res) == val_key)
            self.val = -found_index if (found_index == -1) else self.val
    def __hash__(self):
        return hash(self.get_proper_int())
    def __eq__(self, other):
        return self.get_proper_int() == other.get_proper_int()
    def __repr__(self):
        return f"{'R' if self.is_root() else ''}{make_res(self.get_proper_int())}"
    def __lt__(self, other):
        if self.val >= 0 and other.val >= 0:
            return self.val < other.val
        else:
            return self.get_proper_int() < other.get_proper_int()
    
    def is_root(self):
        return self.val < 0

    def get_proper_int(self):
        if self.is_root():
            return self.owner.root_residuals[-self.val].val
        else:
            return self.val

class Block:
    n = 0
    val = 0
    def __init__(self, n, val):
        self.n = n
        self.val = val

    def __eq__(self, other):
        return self.n == other.n and self.val == other.val
    
    def __hash__(self):
        return hash((self.n, self.val))
    
    def __repr__(self):
        as_bin = format(self.val, f"0{self.n}b")
        return f"{as_bin}|{self.n}"
    
    def res_cz(self):
        res = []
        curr_val = self.val
        n_remaining = self.n
        while n_remaining > 0:
            #print(f"curr_val = {curr_val}, n_remaining = {n_remaining}")
            curr_val, s = col(curr_val)
            if s < 0:
                break
            #print(f"new_val = {curr_val}, s = {s}")
            curr_val, carried = trim(n_remaining, curr_val)
            if s == 0:
                res.append(carried)
            n_remaining -= s
        return res

    def res_m(self, in_res):
        result = []
        val = self.val
        for res in in_res:
            val = val * 3 + res
            val, carried = trim(self.n, val)
            result.append(carried)
        return result, val
    
# make all residual sets up to n values long
def make_all_res(n):
    res = [[]]
    for l in range(1, n + 1):
        for key in range(3**l):
            iter_res = []
            for _ in range(l):
                iter_res.append(key % 3)
                key //= 3
            res.append(iter_res)
    return res

def make_all_res_exactly(n):
    res = []
    for key in range(3**n):
        iter_res = []
        for _ in range(n):
            iter_res.append(key % 3)
            key //= 3
        res.append(iter_res)
    return res

def remove_front_zeroes(res):
    i = 0
    n = len(res)
    while i < n and res[i] == 0:
        i += 1
    return res[i:]

def remove_front_zeroes_from_all(all_res):
    out = []
    seen = set()
    for lst in all_res:
        i = 0
        n = len(lst)
        while i < n and lst[i] == 0:
            i += 1
        trimmed = lst[i:]
        k = res_key(trimmed)
        if k not in seen:
            seen.add(k)
            out.append(list(trimmed))
    return out
    

def make_all_blocks(n) -> List[Block]:
    res = []
    for val in range(2**n):
        res.append(Block(n, val))
    return res
    
def make_cz_res(n):
    seen = set()
    out = []
    for val in range(1 << n):
        b = Block(n, val)
        r = b.res_cz()
        k = res_key(r)
        if k not in seen:
            seen.add(k)
            out.append(r)
    return out

def make_interior_res(n, initial_res):
    seen = set()
    out = []
    for val in range(1 << n):
        b = Block(n, val)
        for r in initial_res:
            res, block_val = b.res_m(r)
            final_key = res_key(res)
            final_key = final_key << n | block_val
            if final_key not in seen:
                seen.add(final_key)
                out.append((res, block_val))
    return out

# possible residual and block values for a block with initial value 0 only
def make_last_res(n, initial_res):
    seen = set()
    out = []
    b = Block(n, 0)
    for r in initial_res:
        res, block_val = b.res_m(r)
        res = remove_front_zeroes(res)
        final_key = res_key(res)
        final_key = final_key << n | block_val
        if final_key not in seen:
            seen.add(final_key)
            out.append((res, block_val))
    return out

def analyze_residuals(n):
    counts = {}
    total_count = (2**n)
    initial_res = make_all_res(n)
    all_res = make_last_res(n, initial_res)
    for res, block_val in all_res:
        final_key = res_key(res)
        print(f"{Block(n, block_val)} -> {res}")
        #final_key = final_key << n | block_val
        #final_key = block_val
        counts[final_key] = counts.get(final_key, 0) + 1
    print(f"Coverage = {len(counts)} / {total_count}, {round(100 * len(counts) / total_count, 2)}% inputs")
    print(list(map(lambda k: make_res(k), counts.keys())))
    if True:
        xs = sorted(counts.keys())
        ys = [counts[x] for x in xs]
        plt.bar(xs, ys)
        plt.show()

# a0




class Context:
    n = 0
    root_residuals = []
    block_to_index = {}
            

    def __init__(self, n):
        self.n = n
        self.root_residuals = sorted(list(map(lambda r: Res(res_key(r), None), make_cz_res(n))))
        print(self.root_residuals)
        #for block in make_all_blocks(n):
            #self.block_to_index[block] = 
        #print(f"res {len(self.root_residuals)} / all {2**n} -> {round(100 * len(self.root_residuals) / 2**n, 2)}%")
        print(self.root_residuals)
        print(self.block_to_index)

    def find_root_res(self, in_res):
        test_key = res_key(in_res)
        return find_index(self.root_residuals, lambda res: res_key(res) == test_key)
    
    def next(self):
        for block in make_all_blocks(self.n):
            print(block)
            for res in self.root_residuals:
                print(f" {self.find_root_res(block.res_m(res)[0])}")


def analyze_mapping(n):
    Context(n)

    return
    all_res = make_cz_res(n)
    all_res = sorted(all_res, key=res_key)
    all_blocks = make_all_blocks(n)
    #all_res = make_all_res(n)
    for block in all_blocks:
        counts = {}
        for res in all_res:
            _, val = block.res_m(res)
            final_key = val
            counts[final_key] = counts.get(final_key, 0) + 1
        print(f"Block {block.val} -> {len(counts)} ({round(100 * len(counts) / 2**n, 2)}%)")



def main():
    if sys.argv[1] == "r":
        n = int(sys.argv[2])
        analyze_residuals(n)
        return
    
    if sys.argv[1] == "m":
        n = int(sys.argv[2])
        analyze_mapping(n)
        return

    n = int(sys.argv[1])
    counts = {}

    if len(sys.argv) <= 2:
        for val in range(2**n):
            b = Block(n, val)
            res = b.res_cz()
            counts[res_key(res)] = counts.get(res_key(res), 0) + 1
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
