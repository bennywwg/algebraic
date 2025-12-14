import z3
from z3 import Bool, BoolVal, Xor, Or, And, Not, Solver

# a4_1 a3_1 a2_1 a1_1
#      a3_2 a2_2 a1_2
#           a2_3 a1_3
#                a1_4

# c3_2

n = 12

bits = []
carries = []
residuals = []

def ternary(condition, yes, no):
    pass

def initialize_variables():
    total = 0
    for it in range(1, n + 1):
        bits.append([])
        for bit in range(1, n + 2 - it):
            bits[-1].append(Bool(f"a{bit}_{it}"))
            total += 1

    for it in range(2, n + 1):
        tmp = []
        for carry in range(3, n + 1):
            if carry <= (n - it + 1):
                tmp.append(Bool(f"c{carry}_{it}"))
                total += 1
        if len(tmp) > 0:
            carries.append(tmp)
    
    print(f"Total variables defined: {total}, exp = {2**total}")

def get_bit(bit, it):
    return bits[it - 1][bit - 1]

# When computing the result of iteration i on a multiplication op,
# Represents the carry bit used in the addition column of this iteration
# No carry bit exists for columns 1 and 2, because they are always 1
# Also no carry bit exists for row 1 because that is the initial value
def get_carry(carry, it):
    if it < 2:
        raise ValueError("Iteration must be at least 1")
    if carry > (n - it + 1):
        #print(f"Getting carry c{carry}_{it}, None")
        return None
    if carry < 3:
        #print(f"Getting carry c{carry}_{it}, True")
        return BoolVal(True)
    #print(f"Getting carry c{carry}_{it}, {carries[it - 2][carry - 3]}")
    return carries[it - 2][carry - 3]

# Creates the expression assigning the output bit and carry bit
def half_adder(s: Solver, out_bit, out_iteration):
    if out_iteration <= 1:
        raise ValueError("Iteration must be greater than 1")
    if out_bit < 1:
        raise ValueError("Bit index starts from 1")
    
    op_bit_var = get_bit(1, out_iteration - 1)
    out_bit_var = get_bit(out_bit, out_iteration)
    if out_bit == 1:
        s.add(out_bit_var == get_bit(2, out_iteration - 1))
    elif out_bit >= 2:
        in_bit_l      = get_bit  (out_bit + 1, out_iteration - 1)
        in_bit_r      = get_bit  (out_bit    , out_iteration - 1)
        inout_bit_c   = get_carry(out_bit,     out_iteration)
        bit_m_op_exp = Xor(Xor(in_bit_r, in_bit_l), inout_bit_c)
        bit_s_op_exp = in_bit_l
        bit_exp = Or(And(op_bit_var, bit_m_op_exp), And(Not(op_bit_var), bit_s_op_exp))
        s.add(out_bit_var == bit_exp)
        if out_bit >= 3:
            prior_carry_var = get_carry(out_bit - 1, out_iteration)
            prior_bit = get_bit(out_bit - 1, out_iteration - 1)
            carry_m_op_exp = Or(And(in_bit_r, prior_bit), And(prior_bit, prior_carry_var), And(in_bit_r, prior_carry_var))
            s.add(inout_bit_c == And(op_bit_var, carry_m_op_exp));

def construct_system(s):
    for it in range(n, 1, -1):
        for bit in range(1, n + 2 - it):
            #print(f"Constructing half adder for bit {bit}, iteration {it}")
            half_adder(s, bit, it)

    if False:
        for bit in range(1, n + 1):
            bit_var = get_bit(bit, 1)
            exp = bit_var == (BoolVal(True) if bit == 1 else BoolVal(False))
            s.add(exp)

    if True:
        for it in range(1, n + 1):
            if it != 1:
                continue
            bit_var = get_bit(1, it)
            exp = bit_var == BoolVal(True)
            s.add(exp)

initialize_variables()

print(bits)
print(carries)

s = Solver()
construct_system(s)
#print(s)

print("\nSolving...")
print(s.check())
m = s.model()


solutions = []
while False and s.check() == z3.sat:
    m = s.model()
    #print("Found solution m {m}")
    solutions.append(m)

    val = 0
    for bit in range(1, n + 1):
        bit_val = 1 if m.evaluate(get_bit(bit, 1)) else 0
        val += 2**bit * bit_val
        print(bit_val, "", end="")
    
    print(" <-", val)

    block = []
    for d in m:
        v = d()                 # the Z3 variable
        val = m[v]              # True or False
        if z3.is_true(val):
            block.append(v)
        else:
            block.append(z3.Not(v))

    # Prevent this exact combination from occurring again
    s.add(z3.Or([z3.Not(b) for b in block]))



exit()
for it in range(1, n + 1):
    if it != 1:
        for bit in range(n, 0, -1):
            carry_val = get_carry(bit, it)
            if carry_val is None:
                print(f"".ljust(14), end="")
            else:
                print(f"{1 if m.evaluate(carry_val) else 0}({carry_val})".ljust(14), end="")
        print()
    
    for bit in range(n, 0, -1):
        bit_val = get_bit(bit, it) if bit <= (n - it + 1) else None
        val = "____" if (bit_val is None) else f"{1 if m.evaluate(bit_val) else 0}({bit_val})"
        print(val.ljust(14), sep="", end="")
    
    print("\n")