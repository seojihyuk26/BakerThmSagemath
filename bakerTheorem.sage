from sage.modules.free_module_integer import IntegerLattice

def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]


# SageMath에서 변수를 정의합니다.
var("x, y, z, b, L, n, h")

L = 3#6
n = 2
h = 1#3
D = 1#3
# over derivate f
d = 0

def bakerAuxiliaryFunction(L, n, h, D):
    c = matrix(
        L + 1,
        L + 1,
        [
            var("c_{}{}".format(i, j), latex_name="c_{{{}{}}}".format(i, j))
            for i in [0..L] for j in [0..L]
        ],
    )
    print("(h+1)*(D+1)^(n-1) < (L+1)^n : ", (h + 1) * (D + 1) ^ (n - 1) < (L + 1) ^ n)
    print(
        "------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    # 주어진 함수를 정의합니다.
    f = sum(
        sum(c[s, t] * 3 ^ (s * x) * 2 ^ (t * x) for s in range(L + 1)) for t in range(L + 1)
    )
    # 함수 f에서 가장 윗 계수 c[L,L]을 -1로 대체합니다.
    f_AntiUnitaire = f.subs(c[L, L] == -1)

    # 함수를 x에 대해 D번까지 미분합니다.
    df = [f_AntiUnitaire.diff(x, i) for i in range(D + 1)]

    # 미분한 결과를 출력합니다.
    print(latex(df))
    print(
        "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )


    # log(3)을 b*log(2)로 대체합니다. (log(3)^3을 b^3**log(2)^3 로 문제없이 잘 대체합니다.)
    dfb = [d.subs(log(3) == b * log(2)) for d in df]

    # 결과를 출력합니다.
    print(latex(dfb))
    print(
        "------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )


    # b에 대한 항과 그렇지 않은 항을 분리합니다.
    dfb_b_terms = [d.coefficients(b) for d in dfb]

    # 결과를 출력합니다.
    # print(latex(dfb_b_terms))
    # print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    # 각 계수를 별도의 함수로 정의합니다.
    cf = [
        [(term[0] / log(2) ^ m).simplify_full().expand() for term in terms]
        for m, terms in enumerate(dfb_b_terms)
    ]

    print("number of equations: ", len(flatten_list(cf)) * (h + 1))
    print("number of variables: ", (L + 1) ** n)
    if len(flatten_list(cf)) * (h + 1) > (L + 1) ** n:
        print("The number of equations > the number of variables.")
        return None
    print(
        "------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    # 결과를 출력합니다.

    # print(latex(cf))
    # print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    #
    # print(latex(c))
    # print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

    # f_b=[]
    # 각 함수를 출력합니다.
    # for i, fi in enumerate(cf):
    # f_b.append((fi/log(2)^m).simplify_full().expand())#[m.expand() for m in .operands()]
    # print(f"f_b^{i}(x) = {[m.expand() for m in (fi/log(2)).simplify_full().expand().operands()]}")
    # for j in range(h+1):
    # print(f_b[i](x=j))#[g(x=j) for g in ])


    # 방정식 시스템을 정의합니다.
    eqs = [
        [[cf[i][j](x=l) == 0 for l in range(h + 1)] for j in range(i + 1)]
        for i in range(D + 1)
    ]
    eqs_flat = flatten_list(flatten_list(eqs))

    # 결과를 출력합니다.
    print(eqs)
    print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    # print(latex(eqs_flat))
    # print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------")


    def linear_equations_to_matrix(eqns, c_list):
        c_list= c_list[:-1]
        # 각 방정식을 재배열합니다.
        print([eq.lhs() for eq in eqns])
        v=[]
        for eq in eqns:
            substituting = dict(zip(c_list, [0] * len(c_list)))
            v.append(-eq.lhs().subs(substituting))
        return (matrix(ZZ, [[eq.lhs().coefficient(v) for v in c_list] for eq in eqns]),vector(ZZ,v), c_list)

    (A, v,u) = linear_equations_to_matrix(eqs_flat, flatten_list(c))

    # Now, 'A' is the coefficient matrix and 'b' is the constant vector

    print("Coefficient matrix, A:")
    print(A)#latex()
    print("Constant vector, u:")
    print(u)#latex()
    print("Right side vector, v:")
    print(v)#latex()
    print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    B, U = A.transpose().LLL(transformation=True)
    nz = sum(1 for r in B.rows() if r==0)   # number of zero rows in B
    B = B.delete_rows(range(nz))            # trimming first nz rows of B
    U = U.delete_rows(range(nz))            # trimming first nz rows of U
    assert U*A.transpose() == B      # the key property of U
    Lattice = IntegerLattice(B, lll_reduce=False)  # our basis is already reduced and should not be altered
    # print(Lattice)
    # print(U)
    assert v in Lattice                        # just in case checking that r belongs to Lattice
    solution = Lattice.coordinate_vector(v) * U
    assert solution*A.transpose() == v                  # have we got correct result?
    # print(solution)
    # 방정식 시스템을 풉니다.
    # solution = solve(eqs_flat, flatten_list(c), solution_dict=True)
    # variableLiberty =  (L + 1) ** n - len(flatten_list(cf)) * (h + 1)
    # r1, r2, ..., rn까지의 변수 리스트 생성
    # r_vars = [var(f'r{i}') for i in range(1, variableLiberty+1)]  

    OneDimensionalSolution = [{u[i]:solution[i]} for i in range(len(u))]
    # 해를 출력합니다.
    print(solution)
    # print(latex(solution))
    print(
        "------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    logsbyderivation = [(log(i))^y for i in A[1]]
    TwoDimensionalSolution = [{u[i]:solution[i]*logsbyderivation[i]} for i in range(len(u))]
    TwoDimensionalSolution.append({c[L,L]:(-1)*(L*log(6))^y})
    # 해를 함수로 변환합니다.
    f_sol = f_AntiUnitaire.subs(OneDimensionalSolution)
    df_sol = f.subs(TwoDimensionalSolution)
    return (f_sol,df_sol)


(f_sol,fxy_sol) = bakerAuxiliaryFunction(L, n, h, D)
if f_sol is not None:
    
    # 해를 함수로 변환합니다.
    df_sol = [f_sol.diff(x, i) for i in range(D + d + 1)]
    # print(latex(df_sol))
    print(
        "------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    print(fxy_sol)
    # 확인합니다.
    for i in range(D + d + 1):
        for l in range(-d-1, h + d+1):
            print(f"f^({i})({l}) =~ ",df_sol[i](x=l).n() , " = ", df_sol[i](x=l))

    # 그래픽을 출력합니다.
    # g = Graphics()
    # for i in range(D + d + 1):
        # g += plot(df_sol[i](x=x), (x, -(0.01), h + 0.01), rgbcolor=hue(i * pi))
    # g.show()

    # 3D plot
    P = implicit_plot3d(z==fxy_sol(x=x,y=y), (x,-1 , 2*h + 0.01), (y, 0, D + d + 1), (z,-1,1))
    P.show()
