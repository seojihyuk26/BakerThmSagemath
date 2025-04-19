from sage.combinat.words import *
from sage.graphs.graph_plot import *
from sage.graphs.graph_plot_js import gen_html_code
from sage.misc.viewer import viewer, browser
from sage.combinat.words.morphism import WordMorphism


class BalancedPair:
    """
    A class to represent a balanced pair and provide custom LaTeX representation.

    Attributes:
        pair (tuple): A tuple representing the balanced pair.
    """

    def __init__(self, pair):
        """
        Initialize the BalancedPair object.

        Args:
            pair (tuple): A tuple representing the balanced pair.
        """
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise ValueError(
                "BalancedPair must be initialized with a tuple of length 2."
            )
        self.pair = pair

    def latex_layout(self, layout=None):
        r"""
        Get or set the actual latex layout (oneliner vs array).

        INPUT:

        - ``layout`` - string (default: ``None``), can take one of the
          following values:

          - ``None`` - Returns the actual latex layout. By default, the
            layout is ``'array'``
          - ``'oneliner'`` - Set the layout to ``'oneliner'``
          - ``'array'`` - Set the layout to ``'array'``

        EXAMPLES::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: s.latex_layout()
            'array'
            sage: s.latex_layout('oneliner')
            sage: s.latex_layout()
            'oneliner'
        """
        if layout is None:
            # return the layout
            if not hasattr(self, "_latex_layout"):
                self._latex_layout = "array"
            return self._latex_layout
        else:
            # change the layout
            self._latex_layout = layout

    def _latex_(self):
        """
        Custom LaTeX representation for the BalancedPair.

        Returns:
            str: The LaTeX representation of the balanced pair.
        """
        from sage.misc.latex import LatexExpr
        latex_layout = self.latex_layout()
        if latex_layout == "oneliner":
            # One-liner format
            return LatexExpr(r"(%s , %s)" % (self.pair[0], self.pair[1]))
        elif latex_layout == "array":
            # Array format
            s = r"\begin{bmatrix}" +self.pair[0] + r"\\ " + self.pair[1] + r"\end{bmatrix}"
            return LatexExpr(s)
        else:
            raise ValueError(f"Unknown latex_layout(={self._latex_layout})")

    def __repr__(self):
        """
        String representation of the BalancedPair.

        Returns:
            str: The string representation of the balanced pair.
        """
        return f"BalancedPair({self.pair})"

    def __getitem__(self, index):
        """
        Allow access to the elements of the pair using indexing.

        Args:
            index (int): The index of the element to access (0 or 1).

        Returns:
            str: The element at the specified index.
        """
        return self.pair[index]

    def __eq__(self, other):
        """
        Check equality between two BalancedPair objects.

        Args:
            other (BalancedPair): Another BalancedPair object.

        Returns:
            bool: True if the pairs are equal, False otherwise.
        """
        if not isinstance(other, BalancedPair):
            return False
        return self.pair == other.pair

    def __hash__(self):
        """
        Compute the hash value of the BalancedPair object.

        Returns:
            int: The hash value of the pair.
        """
        return hash(self.pair)


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


def Balanced_pair_polynomial(l):
    """
    Calculate the balanced polynomial of a balanced pair.

    Args:
        l (tuple): The input letter of balanced pairs.

    Returns:
        Polynomial: The balanced polynomial of the letter.
    """
    x = PolynomialRing(QQ, "x").gen()
    balanced_polynomial = 0
    # print("balanced pair polynomial:", l[0], l[1])
    # print("balanced pair polynomial length:", len(l[0]))
    for i in range(len(l[0])):
        balanced_polynomial += (int(l[0][i]) - int(l[1][i])) * x ** Integer(-i)
    return balanced_polynomial * x ** Integer(len(l[0])-1)


def Balanced_pair_polynomial_dictionary(s):
    """
    Calculate the balanced polynomial of a list of balanced pairs.

    Args:
        s (set): The input list of balanced pairs.

    Returns:
        dictionary : The dictionary of balanced polynomial of the list.
    """
    balanced_polynomial_dictionary = {}
    for i in s:
        _balanced_pair_polynomial = Balanced_pair_polynomial(i)
        if _balanced_pair_polynomial == 0:
            balanced_polynomial_dictionary[i] = [_balanced_pair_polynomial,0]
        else:
            balanced_polynomial_dictionary[i] = [_balanced_pair_polynomial,_balanced_pair_polynomial.factor()]
    return balanced_polynomial_dictionary


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
            A.append(
                BalancedPair((str(s[startpos : head + 1]), str(t[startpos : head + 1])))
            )
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
            print("irreducible balanced pair iteration:", i, r" \\ ")
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


def format_latex_tuple_to_bmatrix(t):
    return r"\begin{bmatrix}" + r"\\ ".join(t) + r"\end{bmatrix}"

def custom_latex_dictionary(d):
    """
    Custom LaTeX representation for a dictionary.

    Args:
        d (dict): The input dictionary.

    Returns:
        str: The LaTeX representation of the dictionary with multiple lines.
    """
    lines = []
    for key, value in d.items():
        if isinstance(key, BalancedPair):
            key_latex = latex(key)
        else:
            key_latex = str(key)

        if isinstance(value, (int, str)):
            value_latex = str(value)
        else:
            value_latex = latex(value)

        lines.append(f"{key_latex} : {value_latex}")

    # Join lines with LaTeX line breaks
    return (
        r"  \begin{array}{l} " + r"\\ ".join(lines) + r" \end{array}"
    )

def custom_latex_mismatch_set(s):
    """
    Custom LaTeX representation for a set of mismatched elements.

    Args:
        s (set): The input set.

    Returns:
        str: The LaTeX representation of the set.
    """
    return (
        r"$\{" + r", ".join([format_latex_tuple_to_bmatrix(i) for i in s]) + r"\}$ \\"
    )


def custom_latex_morphism(morphism):
    """
    Custom LaTeX representation for WordMorphism.

    Args:
        morphism (WordMorphism): The input morphism.

    Returns:
        str: The LaTeX representation of the morphism.
    """

    latex_representation = []
    for key in morphism.domain().alphabet():
        values = morphism.image(key)
        if isinstance(key, BalancedPair):
            joint_string = r"\\ "
            begin_string = r"\begin{array}{l}"
            end_string = r"\end{array}"
            key_latex = latex(key)
        else:
            joint_string = r", "
            begin_string = ""
            end_string = ""
            key_latex = key
        value_latex = ""
        for value in values:
            if isinstance(value, BalancedPair):
                value_latex += latex(value)
            else:
                value_latex += value

        latex_representation.append(f"{key_latex} \\mapsto {value_latex}" + r"\\")

    return (
        r"$$"
        + begin_string
        + joint_string.join(latex_representation)
        + end_string
        + r"$$ \\"
    )


# Override the _latex_ method for WordMorphism
def custom_latex_method(self):
    return custom_latex_morphism(self)


WordMorphism._latex_ = custom_latex_method

def BalancedPair_morphism_latex(morphism_dict,shift=1):
    """
    Custom LaTeX representation for BalancedPair.

    Args:
        dict (dictionary): The input dictionary for morphism.

    Returns:
        str: The LaTeX representation of the balanced pair.
    """
    s = WordMorphism(morphism_dict)
    F = s.fixed_point("0")

    print("")
    print("s:", latex(s))
    print("s is primitive:", s.is_primitive(), r" \\")
    print("F[0:]: $", latex(F), r" $ \\")
    print("F["+str(shift)+":]: $", latex(F[shift:]), r" $ \\")

    result = Balanced_pair_alphabet(F, F[shift:])
    print("balanced pair alphabet by shift: $$", latex(list(result)), r"$$ \\")

    mor = Whole_irreducible_balanced_pair(s, result, 30)
    mor._codomain = mor.domain()
    print("Balanced Pair to Polynomials: $$",custom_latex_dictionary(Balanced_pair_polynomial_dictionary(mor.codomain().alphabet())), r"$$ \\")
    print("Lifted Morphism: ",latex(mor))

    incidence = lifted_morphism_matrix(mor)
    if not mor.is_endomorphism():
        print("lifted morphism is not endomorphism")
        diff = set(mor.domain().alphabet()).difference(set(mor.codomain().alphabet()))
        print("lifted morphism domain difference:", diff)
    else:
        print("lifted morphism is primitive:", mor.is_primitive(), r" \\")

    # print("incidence:", incidence) mor._morph

    G = lifted_morphism_graph(incidence)
    # filename = gen_html_code(
    #     G, vertex_labels=True, gravity=0.05, force_spring_layout=True, charge=-500
    # )

    # with open(filename, "r") as f:
    #     data = f.read()

    # with open("/tmp/dom.html", "w") as f:
    #     f.write(data)

    # os.system("sensible-browser file://///wsl.localhost/Ubuntu-seo/tmp/dom.html")

    P = GraphPlot(G, {"vertex_size": 1000, "layout": "spring"})
    save(P, "/tmp/dom.png")
    os.system("display /tmp/dom.png &")
    print(r"\noindent\makebox[\linewidth]{\rule{\paperwidth}{0.4pt}}")

def BalancedPair_shift_latex(num,shift):
    print("F[0:]: $", latex(num), r" $ \\")
    print("F["+str(shift)+":]: $", latex(num[shift:]), r" $ \\")
    result = Balanced_pair_alphabet(Word(num,alphabet=[i for i in set(str(num))]), Word(num[shift:],alphabet=[i for i in set(str(num))]))
    print("balanced pair alphabet by shift: $$", latex(list(result)), r"$$ \\")
    print("Balanced Pair to Polynomials: $$",custom_latex_dictionary(Balanced_pair_polynomial_dictionary(result)), r"$$ \\")
# Example usage:
BalancedPair_morphism_latex({"0": "020", "1": "01","2": "01"},1)
BalancedPair_morphism_latex({"0": "020", "1": "01","2": "01"},3)
# x = RealField(200)(3).nth_root(2)
# F = x.str(base=2)[2:]
# BalancedPair_shift_latex(F,1)
# BalancedPair_shift_latex(F,2)
# BalancedPair_shift_latex(F,3)
# BalancedPair_shift_latex(F,4)
# BalancedPair_shift_latex(F,5)
# BalancedPair_shift_latex(F,6)

# print("F[0:]: $", latex(F), r" $ \\")
# print("F[1:]: $", latex(F[1:]), r" $ \\")
# result = Balanced_pair_alphabet(Word(F,alphabet=["0","1"]), Word(F[1:],alphabet=["0","1"]))
# print("balanced pair alphabet by shift: $$", latex(list(result)), r"$$ \\")
# print("$$",custom_latex_dictionary(Balanced_pair_polynomial_dictionary(result)), r"$$ \\")