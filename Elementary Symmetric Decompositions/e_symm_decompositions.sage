### Code to compute and format TriCat_n(q,t) in terms of elementary symmetric functions

R.<q,t> = PolynomialRing(ZZ)
E.<e1,e2> = PolynomialRing(QQ)


# Compute longest chain from v to all other vertices in graph reachable from v;
# exploits fact that poset is a DAG
def longest_chain_lengths(G, v):
    top_order = G.topological_sort()
    distances = {v: 0}
    for v in top_order:
        if v not in distances:
            continue
        for w in G.neighbor_out_iterator(v):
            if w in distances:
                distances[w] = max(distances.get(w), distances[v] + 1)
            else:
                distances[w] = distances[v] + 1
    return distances


def i(b,n):
    P = posets.TamariLattice(n)
    H = P.hasse_diagram().reverse()
    distances = longest_chain_lengths(H,b)
    f = 0
    for key in distances:
        f = f+q^(distances[key])
    return f


# Converts a Dyck word to a Dyck path
def path(b):
    temp = DyckWord(list(b[0:-1])).to_non_decreasing_parking_function()
    path = []
    for i in temp:
        path.append(i-1)
    return path


# Returns a polynomial in t
def dinv(b,n):
    c = []
    path_b = path(b)
    for i in range(len(path_b)):
        c.append(i-path_b[i])
    count = 0
    for j in range(0,n):
        for i in range(0,j):
            if c[i] - c[j] == 0 or c[i] - c[j] == 1:
                count = count+1
    return t^count


# Computes TriCat_n^{BPR}(q,t)
def TriCat(n):
    P = posets.TamariLattice(n)
    f = 0

    for b in P:
        f = f + i(b,n)*dinv(b,n)
    return f


# Take a symmetric polynomial f in variables q,t and return its expansion in terms of elementary symmetric functions
# (Courtesy of Claude; I've checked this)
def elementary_expansion(f):
    f = R(f)

    if f != f.subs({q: t, t: q}):
        raise ValueError("Input polynomial is not symmetric in q, t")

    result = E(0)

    while f != 0:
        # Leading coefficient with respect to lex ordering (q > t)
        c = f.leading_coefficient()

        # Exponents of the leading monomial with respect to lex ordering (q > t)
        a,b = f.lm().exponents()[0]  # q^a * t^b

        # Since c(q + t)^{a-b} q^b t^b = c\sum_{i=0}^{a-b} \binom{a-b}{i} q^{a-i}t^{b+i} 
        # = c q^a t^b + c\sum_{i=1}^{a-b} \binom{a-b}{i} q^{a-i}t^{b+i}, the difference
        # f(q,t) - c(q + t)^{a-b} (qt)^b strictly smaller leading coefficient that f with
        # respect to the lex-ordering with q > t.
        term = c * e1^(a - b) * e2^b
        result += term

        # Subtract off its expansion back in q,t
        g = c * (q + t)^(a - b) * (q * t)^b
        f = f - g

    return result


# Compile table of monomial coefficients for TriCat_n^{BPR}(q,t) in terms of e1 & e2
def TriCat_elementary_expansion(a,b):
    polynomials = []
    degrees_in_e2 = []
    for n in [a..b]:
        f = elementary_expansion(TriCat(n))
        polynomials.append(f)
        degrees_in_e2.append(f.degree(e2))
    
    max_deg_2 = max(degrees_in_e2)
    
    Rows = []
    for n in [a..b]:
        f = polynomials.pop(0)
        deg_1 = f.degree(e1)

        header = ["n = " + str(n)]
        for j in [0..max_deg_2]:
            header.append("$e_2^{" + str(j) + "}$")

        Rows.append(header)

        for i in [0..deg_1]:
            row = ["$e_1^{" + str(i) + "}$"]
            for j in [0..max_deg_2]:
                # Add coefficient of e1^i * e2^j to table
                row.append(f.monomial_coefficient(e1^i * e2^j))
            Rows.append(row)

    return table(Rows,frame=True,align='center',header_row=True,header_column=True)

# print(latex(TriCat_elementary_expansion(4,6)))