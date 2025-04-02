from sage.combinat.words import *
from sage.graphs.graph_plot import *
from sage.graphs.graph_plot_js import gen_html_code
from sage.misc.viewer import viewer, browser

def Reset_balanced_vector(s):
    """
    Calculate the balanced vector of a word.

    Args:
        self (Word): The input word.

    Returns:
        dict: The balanced vector of the word.
    """
    balanced_vector = {key: 0 for key in s.parent().alphabet()}
    return balanced_vector

def Balanced_pair_list(s, t, n=1000):
    """
    Check if the given string has balanced pairs of infinite words using Word from Sage.

    Args:
        s, t (word): The input infinite words.

    Returns:
        Alphabet: The Balanced Pair of 2 infinite words.
    """
    if s.parent().alphabet() != t.parent().alphabet():
        raise ValueError("The two words must have the same alphabet.")

    min_length = min(s.length(), t.length(), n)
    head = 0
    startpos = 0
    A = []
    s.balanced_vector = Reset_balanced_vector(s)
    t.balanced_vector = Reset_balanced_vector(t)

    for head in range(min_length):
        s.balanced_vector[s[head]] += 1
        t.balanced_vector[t[head]] += 1

        if s.balanced_vector == t.balanced_vector:
            A.append((str(s[startpos: head + 1]), str(t[startpos: head + 1])))
            startpos = head + 1
            s.balanced_vector = Reset_balanced_vector(s)
            t.balanced_vector = Reset_balanced_vector(t)
    return A

def Balanced_pair_alphabet(s, t, n=1000):
    return set(Balanced_pair_list(s, t, n))

def Lifted_Morphism(s, A):
    """
    Lift a morphism to a new morphism.

    Args:
        s (WordMorphism): The input morphism.
        A (alphabet): The set of balanced pairs.

    Returns:
        WordMorphism: The lifted morphism.
    """
    lifted_morphism = {}
    lifted_alphabet = set(list(A))

    for balanced_pair in A:
        lifted_morphism[balanced_pair] = Balanced_pair_list(
            Word(balanced_pair[0]).apply_morphism(s),
            Word(balanced_pair[1]).apply_morphism(s),
        )
        lifted_alphabet = lifted_alphabet.union(set(lifted_morphism[balanced_pair]))

    return (WordMorphism(lifted_morphism), lifted_alphabet)

def Whole_irreducible_balanced_pair(s, alpha, iteration=100):
    """
    Check if the given string has whole irreducible balanced pairs.

    Args:
        s (morphism): The initial morphism.
        iteration (int): The number of iterations to check.

    Returns:
        bool: True if the word has whole irreducible balanced pairs, False otherwise.
    """
    alpha_0 = alpha

    for i in range(iteration):
        (mor, alpha_1) = Lifted_Morphism(s, alpha_0)

        if len(alpha_1.difference(alpha_0)) == 0:
            print("irreducible_balanced_pair iteration:", i)
            return mor

        alpha_0 = alpha_1

    return False

def lifted_morphism_matrix(m):
    """
    Create a matrix of lifted morphisms.

    Args:
        m (WordMorphism): The input morphism.

    Returns:
        Matrix: The matrix of lifted morphisms.
    """
    return m.incidence_matrix()

def lifted_morphism_graph(m):
    """
    Create a graph of lifted morphisms.

    Args:
        m (matrix): The input morphism incidence matrix.

    Returns:
        Graph: The graph of lifted morphisms.
    """
    G = DiGraph(m, loops=True)
    return G

# Example usage:
s = WordMorphism({"0": "01", "1": "02", "2": "0"})
F = s.fixed_point("0")

print("F:", F)
print("s:", s)
print("s is_primitive:", s.is_primitive())

result = Balanced_pair_alphabet(F, F[1:])
print("balanced pair alphabet by shift:", result)

mor = Whole_irreducible_balanced_pair(s, result, 30)
mor._codomain = mor.domain()

print(latex(mor))

incidence = lifted_morphism_matrix(mor)
if not mor.is_endomorphism():
    print("lifted_morphism is not endomorphism")
    diff = set(mor.domain().alphabet()).difference(set(mor.codomain().alphabet()))
    print("lifted_morphism domain difference:", diff)
else:
    print("lifted_morphism is primitive:", mor.is_primitive())

print("incidence:", incidence)

G = lifted_morphism_graph(mor._morph)
filename = gen_html_code(
    G, vertex_labels=True, gravity=0.05, force_spring_layout=True, charge=-500
)

with open(filename, "r") as f:
    data = f.read()

with open("/tmp/dom.html", "w") as f:
    f.write(data)

os.system("sensible-browser file://///wsl.localhost/Ubuntu-seo/tmp/dom.html")

P = GraphPlot(G, {"vertex_size": 1000, "layout": "spring"})
save(P, "/tmp/dom.png")
os.system("display /tmp/dom.png &")