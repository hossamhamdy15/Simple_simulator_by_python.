"""
This is a program written as a single function. The function should accept the netlist
as an input (a path to a text file). The function should return the symbolic and
numerical results
"""

import re
from sympy import symbols, sympify, Matrix
import numpy as np

#-------------------------------------------------------------
# Part 1: Reading the netlist
#-------------------------------------------------------------
netpath = r"Put here netlist name"

with open(netpath, 'r') as f:
    raw_netlist = f.read()

# Clean spaces
netlist = re.sub(r' +', ' ', raw_netlist)
netlist = re.sub(r' I', 'I', netlist)
netlist = re.sub(r' R', 'R', netlist)
netlist = re.sub(r' V', 'V', netlist)
netlist = re.findall(r'[^\n]+', netlist)

#-------------------------------------------------------------
# Part 2: Parsing (you must define parse_netlist)
#-------------------------------------------------------------
def parse_netlist(netlist, element_type):
    """Example parser: expects lines like R1 1 2 1000"""
    names, n1, n2, values = [], [], [], []
    for line in netlist:
        if line.strip().startswith(element_type):
            parts = line.strip().split()
            names.append(parts[0])
            n1.append(parts[1])
            n2.append(parts[2])
            values.append(float(parts[3]))
    return n1, n2, values, names

R_Node_1, R_Node_2, R_Values, R_Names = parse_netlist(netlist, 'R')
V_Node_1, V_Node_2, V_Values, V_Names = parse_netlist(netlist, 'V')
I_Node_1, I_Node_2, I_Values, I_Names = parse_netlist(netlist, 'I')

#-------------------------------------------------------------
# Part 3: Create matrices
#-------------------------------------------------------------
nodes_list = R_Node_1 + R_Node_2 + V_Node_1 + V_Node_2 + I_Node_1 + I_Node_2
nodes_number = max(map(int, nodes_list))
matrices_size = nodes_number + len(V_Names)

#-------------------------------------------------------------
# Z matrix
#-------------------------------------------------------------
z = ['0'] * matrices_size

# stamp current sources
for i, name in enumerate(I_Names):
    n1 = int(I_Node_1[i])
    n2 = int(I_Node_2[i])
    if n1 != 0:
        z[n1 - 1] += f"-{name}"
    if n2 != 0:
        z[n2 - 1] += f"+{name}"

# voltage sources
for i, name in enumerate(V_Names):
    z[nodes_number + i] = name

Z = Matrix([sympify(expr) for expr in z])

#-------------------------------------------------------------
# X matrix
#-------------------------------------------------------------
x = [f"V_{i+1}" for i in range(nodes_number)]
for name in V_Names:
    x.append(f"I_{name}")

X = Matrix(symbols(' '.join(x)))

#-------------------------------------------------------------
# G matrix
#-------------------------------------------------------------
G = [['0' for _ in range(nodes_number)] for _ in range(nodes_number)]
for i, name in enumerate(R_Names):
    n1 = int(R_Node_1[i])
    n2 = int(R_Node_2[i])
    if n1 != 0:
        G[n1-1][n1-1] += f"+1/{name}"
    if n2 != 0:
        G[n2-1][n2-1] += f"+1/{name}"
    if n1 != 0 and n2 != 0:
        G[n1-1][n2-1] += f"-1/{name}"
        G[n2-1][n1-1] += f"-1/{name}"

#-------------------------------------------------------------
# B matrix
#-------------------------------------------------------------
B = [['0'] * len(V_Names) for _ in range(nodes_number)]
for i, name in enumerate(V_Names):
    n1 = int(V_Node_1[i])
    n2 = int(V_Node_2[i])
    if n1 != 0:
        B[n1-1][i] = B[n1-1][i] + '+1'
    if n2 != 0:
        B[n2-1][i] = B[n2-1][i] + '-1'

#-------------------------------------------------------------
# Combine A matrix
#-------------------------------------------------------------
C = list(map(list, zip(*B)))  # transpose
A_top = [G_row + B_row for G_row, B_row in zip(G, B)]
A_bottom = [C_row + ['0'] * len(V_Names) for C_row in C]
A = Matrix([[sympify(e) for e in row] for row in (A_top + A_bottom)])

#-------------------------------------------------------------
# Part 4: Solve
#-------------------------------------------------------------
symbolic_ans = A.LUsolve(Z)

# Assign numeric values
locals_dict = {}

for i, name in enumerate(R_Names):
    locals_dict[name] = R_Values[i]
for i, name in enumerate(V_Names):
    locals_dict[name] = V_Values[i]
for i, name in enumerate(I_Names):
    locals_dict[name] = I_Values[i]

numeric_ans = [ans.evalf(subs=locals_dict) for ans in symbolic_ans]

#-------------------------------------------------------------
# Print results
#-------------------------------------------------------------
for var, val in zip(x, numeric_ans):
    print(f"{var} = {float(val):.6f}")
