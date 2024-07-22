import inflect, re, symfem, sympy, sys

def get_basis(geom, order, dx, dy, dz):

    elem = symfem.create_element(geom, "Nedelec", order)
    basis = elem.get_basis_functions()

    if geom == "triangle" or geom == "quadrilateral":
        for e in range(len(elem.reference.edges)):
            edge_dofs = elem.dof_plot_positions()[e*order : (e+1)*order]
            edge_funs = elem.get_basis_functions()[e*order : (e+1)*order]
            basis[e*order : (e+1)*order] = [f for _, f in sorted(zip(edge_dofs, edge_funs), key = lambda x: x[0])]

    if geom == "triangle":
        basis[: order] = basis[order - 1 : : -1]
        basis[order : 2*order] = basis[2*order - 1 : order - 1 : -1]
        basis[order : 2*order] = [-f for f in basis[order : 2*order]]
        basis = basis[2*order : 3*order] + basis[: 2*order] + basis[3*order : ]
    elif geom == "quadrilateral":
        basis[order : 2*order] = basis[2*order - 1 : order - 1 : -1]
        basis[order : 2*order] = [-f for f in basis[order : 2*order]]
        basis[3*order : 4*order] = basis[4*order - 1 : 3*order - 1 : -1]
        basis[3*order : 4*order] = [-f for f in basis[3*order : 4*order]]
        basis = basis[ : order] + basis[2*order : 4*order] + basis[order : 2*order] + basis[4*order : ]

    if geom == "tetrahedron":
        for e in range(len(elem.reference.edges)):
            edge = elem.reference.edges[e]
            edge_dofs = elem.dof_plot_positions()[e*order : (e+1)*order]
            edge_funs = elem.get_basis_functions()[e*order : (e+1)*order]
            basis[e*order : (e+1)*order] = [f for _, f in sorted(zip(edge_dofs, edge_funs), key = lambda x: x[0][edge[1]-1])]
        
        temp = basis

        basis[        :   order] = temp[5*order : 6*order]
        basis[  order : 2*order] = temp[2*order : 3*order]
        basis[2*order : 3*order] = temp[4*order : 5*order]
        basis[3*order : 4*order] = temp[3*order : 4*order]
        basis[4*order : 5*order] = temp[  order : 2*order]
        basis[5*order : 6*order] = temp[        :   order]

        basis = basis[: 6*order] + basis[6*order : ]

    if geom == "hexahedron":

        edge_oreintation_dict = dict([((0, 1), 0), 
                                    ((0, 2), 1), 
                                    ((0, 4), 2), 
                                    ((1, 3), 1), 
                                    ((1, 5), 2), 
                                    ((2, 3), 0), 
                                    ((2, 6), 2), 
                                    ((3, 7), 2), 
                                    ((4, 5), 0), 
                                    ((4, 6), 1), 
                                    ((5, 7), 1), 
                                    ((6, 7), 0)])

        for e in range(len(elem.reference.edges)):
            edge = elem.reference.edges[e]
            edge_dofs = elem.dof_plot_positions()[e*order : (e+1)*order]
            edge_funs = elem.get_basis_functions()[e*order : (e+1)*order]
            basis[e*order : (e+1)*order] = [f for _, f in sorted(zip(edge_dofs, edge_funs), key = lambda x: x[0][edge_oreintation_dict[edge]])]

            temp = basis

            basis[         :    order] = temp[         :   order]
            basis[  order  :  2*order] = temp[3*order  : 4*order]
            basis[2*order  :  3*order] = temp[5*order  : 6*order]
            basis[2*order  :  3*order] = basis[3*order - 1 : 2*order - 1 : -1]
            basis[2*order  :  3*order] = [-f for f in basis[2*order  :  3*order]]
            basis[3*order  :  4*order] = temp[  order  : 2*order]
            basis[4*order  :  5*order] = temp[2*order  : 3*order]
            basis[5*order  :  6*order] = temp[4*order  : 5*order]
            basis[6*order  :  7*order] = temp[7*order  : 8*order]
            basis[7*order  :  8*order] = temp[6*order  : 7*order]
            basis[8*order  :  9*order] = temp[8*order  : 9*order]
            basis[9*order  : 10*order] = temp[10*order : 11*order]
            basis[10*order : 11*order] = temp[11*order : 12*order]
            basis[10*order : 11*order] = basis[11*order - 1 : 10*order - 1 : -1]
            basis[10*order : 11*order] = [-f for f in basis[10*order : 11*order]]
            basis[11*order : 12*order] = temp[ 9*order : 10*order]

            basis = basis[: 12*order] + basis[12*order : ]

    x, y, z = sympy.symbols('x y z')
    if geom == "quadrilateral":
        basis = [f.subs((x, y, z), ((1 + x) / 2, (1 + y) / 2, (1 + z) / 2)) for f in basis]
        basis = [f / sympy.sympify(2) for f in basis]

    basis = [f.diff((x, dx)) for f in basis]
    basis = [f.diff((y, dy)) for f in basis]
    basis = [f.diff((z, dz)) for f in basis]

    xi, eta, zeta = sympy.symbols('xi eta zeta')
    basis = [f.subs((x, y, z), (xi, eta, zeta)) for f in basis]

    for o in range(order, 1, -1):
        basis = [f.subs(xi**o, sympy.UnevaluatedExpr(sympy.sympify(('xi*'*o)[:-1], locals={"xi": xi}, evaluate = False))) for f in basis]
        basis = [f.subs(eta**o, sympy.UnevaluatedExpr(sympy.sympify(('eta*'*o)[:-1], locals={"eta": eta}, evaluate = False))) for f in basis]
        basis = [f.subs(zeta**o, sympy.UnevaluatedExpr(sympy.sympify(('zeta*'*o)[:-1], locals={"zeta": zeta}, evaluate = False))) for f in basis]
        basis = [f.subs((xi + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(xi + 1)*'*o)[:-1], locals={"xi": xi}, evaluate = False))) for f in basis]
        basis = [f.subs((eta + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(eta + 1)*'*o)[:-1], locals={"eta": eta}, evaluate = False))) for f in basis]
        basis = [f.subs((zeta + 1)**o, sympy.UnevaluatedExpr(sympy.sympify(('(zeta + 1)*'*o)[:-1], locals={"zeta": zeta}, evaluate = False))) for f in basis]
        basis = [f.subs(((xi + 1)/2)**o, sympy.UnevaluatedExpr(sympy.sympify(('(xi + 1)/2*'*o)[:-1], locals={"xi": xi}, evaluate = False))) for f in basis]
        basis = [f.subs(((eta + 1)/2)**o, sympy.UnevaluatedExpr(sympy.sympify(('(eta + 1)/2*'*o)[:-1], locals={"eta": eta}, evaluate = False))) for f in basis]
        basis = [f.subs(((zeta + 1)/2)**o, sympy.UnevaluatedExpr(sympy.sympify(('(zeta + 1)/2*'*o)[:-1], locals={"zeta": zeta}, evaluate = False))) for f in basis]

    p = re.compile(r'(\d+)')
    basis = [p.sub(r'\1.', str(f)) for f in basis]

    return basis

dim = int(sys.argv[1])
order = int(sys.argv[2])
derivatives = int(sys.argv[3])

p = inflect.engine()

print("case " + p.number_to_words(p.ordinal(order)).upper() + ":\n"
      "  {\n"
      "    switch (elem->type())\n"
      "      {")

for geom in ["quadrilateral", "triangle"]:

    if geom == "triangle":
        print("      case TRI6:\n"
              "      case TRI7:\n"
              "        {")
    elif geom == "quadrilateral":
        print("      case QUAD8:\n"
              "      case QUAD9:\n"
              "        {")

    elem = symfem.create_reference(geom)

    if derivatives:
        print("          switch (j)\n"
              "            {")

    for d in range(derivatives + 1):
        dx = derivatives - d
        dy = d

        basis = get_basis(geom, order, dx, dy, 0)
        spaces = 6 * " " if derivatives else ""

        if derivatives:
            print(f"              // d" + f"^{derivatives}" * (derivatives > 1) + "()/" +
                                    "dxi" * (dx > 0) + f"^{dx}" * (dx > 1) +
                                    "deta" * (dy > 0) + f"^{dy}" * (dy > 1) + "\n"
                  f"            case {d}:\n"
                   "              {")

        print(spaces + "          switch(ii)\n" +
              spaces + "            {")

        for f in range(len(elem.edges) * order):
            print(spaces + f"            case {f}:\n" +
                  spaces + f"              return sign * RealGradient{basis[f]};")

        for f in range(len(elem.edges) * order, len(basis)):
            print(spaces + f"            case {f}:\n" +
                  spaces + f"              return RealGradient{basis[f]};")

        print(spaces + "            default:\n" +
              spaces + "              libmesh_error_msg(\"Invalid i = \" << i);\n" +
              spaces + "            }")

        if derivatives:
            print(f"              }} // j = {d}\n")

    if derivatives:
        print("            default:\n"
              "              libmesh_error_msg(\"Invalid j = \" << j);\n"
              "            }")

    print("        }\n")

print( "      default:\n"
      f"        libmesh_error_msg(\"ERROR: Unsupported {dim}D element type!: \" << Utility::enum_to_string(elem->type()));\n"
       "      } // end switch (type)\n"
       "  } // end case " + p.number_to_words(p.ordinal(order)).upper())