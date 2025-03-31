from sage.combinat.words import *

def Reset_balanced_vector(s):
    """
    Calculate the balanced vector of a word.
    
    Args:
        self (Word): The input word.
    
    Returns:
        dict: The balanced vector of the word.
    """  
    # print("self:", self)
    balanced_vector = {key: 0 for key in s.parent().alphabet()}
    return balanced_vector

def Balanced_pair_list(s,t, n=1000):
    """
    Check if the given string has balanced pairs of infinite words using Word from Sage.
    
    Args:
        s,t (word): The input infinite words.
    
    Returns:
        Alphabet: The Balanced Pair of 2 infinite words.
    """  
    if s.parent().alphabet() != t.parent().alphabet() :
        raise ValueError("The two words must have the same alphabet.")
    min_length = min(s.length(), t.length(), n)
    # print("min_length:", min_length)
    head = 0
    startpos = 0
    A = []
    s.balanced_vector = Reset_balanced_vector(s) 
    t.balanced_vector = Reset_balanced_vector(t) 
    for head in range(min_length):
            # print("subword_s:", subword_s)
            s.balanced_vector[s[head]] +=1
            t.balanced_vector[t[head]] +=1
            # print("s.balanced_vector:", s.balanced_vector)
            # print("t.balanced_vector:", t.balanced_vector)
            if s.balanced_vector==t.balanced_vector:
                A.append((str(s[startpos:head+1]), str(t[startpos:head+1])))
                # print("A+:", [s[startpos:head+1], t[startpos:head+1]])
                # print("head:" , head)
                startpos = head + 1
                s.balanced_vector = Reset_balanced_vector(s) 
                t.balanced_vector = Reset_balanced_vector(t) 
    return A

def Balanced_pair_alphabet(s,t, n=1000):
    return set(Balanced_pair_list(s,t, n))

def Lifted_Morphism(s,A):
    """
    Lift a morphism to a new morphism.
    
    Args:
        s (WordMorphism): The input morphism.
        A (alphabet): The set of balanced pairs.
    
    Returns:
        WordMorphism: The lifted morphism.
    """  
    # print("s:", s)
    # print("s._letter:", s._letter)
    # print("s._letter:", s._letter)
    # print("s._letter:", s._letter)
    lifted_morphism = {}
    lifted_alphabet = set(list(A))
    for balanced_pair in A:
        # print("balanced_pair:", balanced_pair)
        # print("lifted_alphabet:", lifted_alphabet)
        lifted_morphism[balanced_pair] = Balanced_pair_list(Word(balanced_pair[0]).apply_morphism(s),Word(balanced_pair[1]).apply_morphism(s))
        # print("lifted_morphism[balanced_pair]:", lifted_morphism[balanced_pair])
        # print("set lifted_morphism[balanced_pair]:", set(lifted_morphism[balanced_pair]))
        lifted_alphabet = lifted_alphabet.union(set(lifted_morphism[balanced_pair]))
    return (WordMorphism(lifted_morphism),lifted_alphabet)


# Example usage:
s = WordMorphism({'0': '01', '1': '02', '2': '0'})
F = s.fixed_point("0")
print("F.parent:", F.parent().alphabet())
result = Balanced_pair_alphabet(F, F[1:])
print("balanced pair alphabet:", result)
print("Lifted_Morphism:", Lifted_Morphism(s, result))
mor = s
alpha_0 = result
for i in range(30):
    (mor,alpha_1) = Lifted_Morphism(s, alpha_0)
    # print("alpha_1:", alpha_1)
    # print("alpha_0:", alpha_0)
    if len(alpha_1.difference(alpha_0))==0:
        break
    alpha_0 = alpha_1
    # print("Lifted_Morphism alphabet:", mor.domain())
# print("Lifted_Morphism:", mor)
print(latex(mor))
