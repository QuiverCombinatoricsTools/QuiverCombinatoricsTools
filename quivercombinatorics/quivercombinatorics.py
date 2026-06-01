import random
import dot2tex

from sage.all import Integer, Poset
from sage.combinat.partition import Partitions
from sage.combinat.posets.hasse_diagram import HasseDiagram
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.matrix.special import zero_matrix
from sage.misc.latex import LatexExpr
from sage.misc.latex_standalone import TikzPicture
from sage.modules.free_module_element import vector
from sage.structure.element import Element

from quiver import *
from quiver import Quiver as BaseQuiver

def quiver_from_cartan_matrix(C):
    r"""Returns the quiver :math:`Q` given by a Cartan matrix :math:`C`
        
        INPUT:

        - ``C`` -- a square matrix with entries in :math:`\mathbb{Z}`

        OUTPUT: The quiver :math:`Q` given by a Cartan matrix :math:`C`

        EXAMPLE::
        
            sage: from quivercombinatorics import *
            sage: C = [[2, -1, 0], [-1, 2, -2], [0, -2, 2]]
            sage: quiver_from_cartan_matrix(C)
            a quiver with 3 vertices and 3 arrows
        
        """
    for i in range(len(C)):
        C[i][i] /= 2                        
        for j in range(i):
            C[i][j] = 0                     
    E = matrix.identity(len(C)) - matrix(C)
    return Quiver(E)

def random_quiver(vertices, max_arrows_per_edge):
    r"""Returns a randomly generated quiver
        
        INPUT:

        - ``vertices`` -- number of vertices you want in the quiver, in :math:`\mathbb{N}_0`
        - ``max_arrows_per_edge`` -- maximum of number of arrows you want in the quiver, in :math:`\mathbb{N}_0`

        OUTPUT: A randomly generated quiver with the prescribed vertices and maximum number of edges

        EXAMPLES::

            sage: from quivercombinatorics import *
            sage: random_quiver(3, 2)
            a quiver with 3 vertices and 6 arrows

            sage: from quivercombinatorics import *
            sage: random_quiver(46, 25)
            a quiver with 46 vertices and 26388 arrows
            
        
        """
    A = [[random.randint(0, max_arrows_per_edge) for _ in range(vertices)] for _ in range(vertices)]
    G = DiGraph(matrix(A))
    return Quiver(A)

def N_set(S, v):
    r"""For a list of vectors :math:`S` in :math:`\mathbb{Z}Q_0`, returns all possible sums of vectors in :math:`S` as a list, up to the upper bound :math:`v`. Though not directly about quivers, this is a helper function to define :math:`\Sigma_\lambda`
        
        INPUT:

        - ``S`` -- a list of vectors :math:`S` in :math:`\mathbb{Z}Q_0`
        - ``v`` -- vector in :math:`\mathbb{Z}Q_0`

        OUTPUT: A list of all possible sums of vectors in :math:`S` as a list, up to the upper bound :math:`v`

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: N_set([[0, 1]], [0, 3])
            [(0, 1), (0, 2), (0, 3)]
        
        """
    if isinstance(v, (int, Integer)):
            v = (v,)
    if isinstance(v, int):
        v = vector([v])
    else:
        v = vector(v)
    reachable = {tuple(a) for a in S}
    changed = True
    while changed:
        new_reachable = set(reachable)
        for a in reachable:
            a_vec = vector(a)
            for b in S:
                b_vec = vector(b)
                s = a_vec + b_vec
                if all(s[i] <= v[i] for i in range(len(v))):
                    new_reachable.add(tuple(s))
        changed = (new_reachable != reachable)
        reachable = new_reachable
    return [vector(a) for a in reachable]

def vector_decomposition(x, S):
    r"""For a list of vectors :math:`S` in :math:`\mathbb{Z}Q_0`, returns all possible sums of vectors in :math:`S` as a list that sum up to :math:`x`. This is a helper function in order to define CB-decompositions
        
        INPUT:

        - ``x`` -- a vector in :math:`\mathbb{Z}Q_0`
        - ``S`` -- a list of vectors :math:`S` in :math:`\mathbb{Z}Q_0`

        OUTPUT: All possible sums of vectors in :math:`S` as a list that sum up to :math:`x`

        EXAMPLE::
        
            sage: from quivercombinatorics import *
            sage: vector_decomposition((4, 5), [(0, 1), (1, 0), (1, 1)])
            [[[(0, 1), 1], [(1, 1), 4]],
            [[(0, 1), 2], [(1, 0), 1], [(1, 1), 3]],
            [[(0, 1), 3], [(1, 0), 2], [(1, 1), 2]],
            [[(0, 1), 4], [(1, 0), 3], [(1, 1), 1]],
            [[(0, 1), 5], [(1, 0), 4]]]
        
        """
    x = vector(x)
    if all(i == 0 for i in x):
        return [[]]
    if not S:
        return []
    s = vector(S[0])
    decomps = vector_decomposition(x, S[1:]) 
    n = 1
    while all(i >= 0 for i in x - n*s):
        current = vector_decomposition(x - n*s, S[1:])
        current = [[[s,n]] + item for item in current]
        decomps = decomps + current
        n += 1
    return sorted(decomps)

def small_decomposition(v, n):
    r"""For a vector :math:`v` in :math:`\mathbb{Z}Q_0`, it returns all possible partitions of the representation type :math:`[v,n]`
        
        INPUT:

        - ``v`` -- a vector in :math:`\mathbb{Z}Q_0`
        - ``n`` -- a natural number in :math:`\mathbb{N}`

        OUTPUT: All possible partitions of the representation type :math:`[v,n]`

        EXAMPLE::
        
            sage: from quivercombinatorics import *
            sage: small_decomposition((1, 1), 3)
            [[[(1, 1), 1], [(1, 1), 1], [(1, 1), 1]],
            [[(1, 1), 2], [(1, 1), 1]],
            [[(1, 1), 3]]]
        
        """
    output = []
    for partition in Partitions(n).list():
        temp = []
        for i in partition:
            temp = temp + [[v,i]]
        output = output + [temp]
    return sorted(output)

def ext_dimension_vector(tau):
    """For a representation type :math:`\\tau`, returns the associated dimension vector

        INPUT:

        - ``tau`` -- a representation type, i.e., a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}Q_0`, and the second is in :math:`\mathbb{N}`

        OUTPUT: a vector in :math:`\mathbb{N}Q_0`

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: tau = [[(1, 1), 2], [(1, 1), 1], [(1, 1), 1], [(1, 1), 1]]
            sage: ext_dimension_vector(tau)
            (2, 1, 1, 1)
    """
    return vector([pair[1] for pair in tau])

def D_map(m, tau):
    r"""Evaluates a map :math:`D:\mathbb{Z}^k\to\mathbb{Z}^n` from dimension vectors of the :math:`\mathrm{ext}`-quiver to dimension vectors of the original quiver by :math:`D(m_1,\dots,m_k):=\sum_{i=1}^k m_i\beta^{(i)}`

        INPUT:
        
        - ``m`` -- a vector in :math:`\mathbb{Z}\operatorname{ext}(Q)_0`
        - ``tau`` -- a representation type of the original quiver, i.e., a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}Q_0`, and the second is in :math:`\mathbb{N}`

        OUTPUT: a vector in :math:`\mathbb{Z}Q_0`

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: tau = [[(1, 1), 2], [(1, 1), 1], [(1, 1), 1], [(1, 1), 1]]
            sage: D_map((1, 2, 8, -3), tau)
            (8, 8)
    
    """
    return sum(m[i]*vector(tau[i][0]) for i in range(len(m)))

def D_lifting(tau, L):
    r"""Applies ``D_map`` to a representation type of the :math:`\mathrm{ext}`-quiver

        INPUT:

        - ``L`` -- a representation type of the :math:`\mathrm{ext}`-quiver, i.e., a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}\mathrm{ext}(Q)_0`, and the second is in :math:`\mathbb{N}`
        - ``tau`` -- a representation type, i.e., a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}Q_0`, and the second is in :math:`\mathbb{N}`

        OUTPUT: a representation type of the :math:`\mathrm{ext}`-quiver

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: tau = [[(1, 1), 2], [(1, 1), 1], [(1, 1), 1], [(1, 1), 1]]
            sage: D_lifting([[(1, 2, 8, -3), 3]], tau)
            [[(8, 8), 3]]

    """
    if len(tau) > 1 and isinstance(tau[1], str):
        return (sorted([[D_map(pair[0], L), pair[1]] for pair in tau[0]]), tau[1])
    else:
        return sorted([[D_map(pair[0], L), pair[1]] for pair in tau])
    
def is_direct_successor(epsilon, rho):
    r"""Verifies whether ``epsilon`` is the direct successor of ``rho``
        
        INPUT:

        - ``epsilon`` -- an element of :math:`\mathbb{N}^n`
        - ``rho`` -- an element of :math:`\mathbb{N}^n`

        OUTPUT: ``True`` if ``epsilon`` is the direct successor of ``rho``, and ``False`` otherwise

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: rho = [((1, 0), 2), ((1, 0), 3)]
            sage: epsilon = [((1, 0), 5)]
            sage: is_direct_successor(epsilon, rho)
            True
        
        """
    epsilon_copy = list(epsilon)
    rho_copy = list(rho)
    for pair_e in list(epsilon_copy):
        for pair_r in list(rho_copy):
            if vector(pair_e[0]) == vector(pair_r[0]) and pair_e[1] == pair_r[1]:
                epsilon_copy.remove(pair_e)
                rho_copy.remove(pair_r)
                break 
    if len(epsilon_copy) == 2 and len(rho_copy) == 1:
        alpha_r, m_r = rho_copy[0]
        gamma_r, p_r = epsilon_copy[0]
        gamma_r1, p_r1 = epsilon_copy[1]
        return (m_r == p_r == p_r1) and (gamma_r + gamma_r1 == alpha_r)
    elif len(epsilon_copy) == 1 and len(rho_copy) == 2:
        alpha_r, m_r = rho_copy[0]
        alpha_r1, m_r1 = rho_copy[1]
        gamma_r1, p_r1 = epsilon_copy[0]
        return (m_r + m_r1 == p_r1) and (alpha_r == alpha_r1 == gamma_r1)
    return False

class Quiver(BaseQuiver):
    def p_function(self, x):
        r"""Outputs the function :math:`p(x) = 1 - \frac{1}{2}(x, x)`, where :math:`(x, x)` is the symmetrized Euler form. Note that :math:`p(x)\geq 0` if :math:`x` is a root, and :math:`p(x) = 0` if and only if :math:`x` is a real root
        
        INPUT:

        - ``x`` -- an element of :math:`\mathbb{Z}Q_0`

        OUTPUT: :math:`1 - \frac{1}{2}(x, x)`

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: Q.p_function((2, 3))
            0.0
        
        """
        return 1 - 0.5 * self.symmetrized_euler_form(x, x)

    def R_lambda_plus(self, l, v):
        r"""Returns a list of elements of :math:`R_\lambda^+`, the set of positive roots :math:`\alpha` with :math:`\alpha\cdot\lambda=0`, up to the upper bound :math:`v`


        INPUT:

        - ``l`` -- an element of :math:`\mathbb{Z}Q_0`
        - ``v`` -- an element of :math:`\mathbb{N}Q_0`

        OUTPUT: A list of elements of :math:`R_\lambda^+`, where ``l`` is :math:`\lambda`, the set of positive roots :math:`\alpha` with :math:`\alpha\cdot\lambda=0`, up to the upper bound :math:`v`

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: Q.R_lambda_plus((1, -1), (5, 5))
            [(1, 1), (2, 2), (3, 3), (4, 4)]
        
        """
        return list(
            filter(
                lambda a: self.is_root(a) and a.dot_product(vector(l)) == 0,
                self.all_subdimension_vectors(v, proper=True, nonzero=True)
            )
        )

    def sigma_lambda(self, l, v):
        r"""Returns a list of elements of :math:`\Sigma_\lambda`, up to the upper bound :math:`v`


        INPUT:

        - ``l`` -- an element of :math:`\mathbb{Z}Q_0`
        - ``v`` -- an element of :math:`\mathbb{N}Q_0`

        OUTPUT: A list of elements of :math:`\Sigma_\lambda`, where ``l`` is :math:`\lambda`, up to the upper bound :math:`v`

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: Q.sigma_lambda((1, -1), (5, 5))
            [(1, 1)]
        
        """
        v = self._coerce_dimension_vector(v)
        l = self._coerce_vector(l)
        NR_lambda_plus = N_set(self.R_lambda_plus(l, v), v)
        return list(
            filter(
                lambda a: not any(self.symmetrized_euler_form(b, a - b) > -2 and all(x >= 0 for x in a - b) and any(x > 0 for x in a - b) for b in NR_lambda_plus),
                NR_lambda_plus
            )
        )

    def all_representation_types(self, l, x):
        r"""Returns a list of representation types of a quiver, with respect to :math:`x`


        INPUT:

        - ``l`` -- an element of :math:`\mathbb{Z}Q_0`
        - ``x`` -- an element of :math:`\mathbb{N}Q_0`

        OUTPUT: A list of representation types of a quiver, with respect to :math:`x`, and where ``l`` is :math:`\lambda`. Each representation type is stored as a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}Q_0`, and the second is in :math:`\mathbb{N}`

        EXAMPLES::

            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: Q.all_representation_types((1, -1), (5, 5))
            [[[(1, 1), 1], [(1, 1), 1], [(1, 1), 1], [(1, 1), 1], [(1, 1), 1]],
            [[(1, 1), 1], [(1, 1), 1], [(1, 1), 1], [(1, 1), 2]],
            [[(1, 1), 1], [(1, 1), 1], [(1, 1), 3]],
            [[(1, 1), 1], [(1, 1), 2], [(1, 1), 2]],
            [[(1, 1), 1], [(1, 1), 4]],
            [[(1, 1), 2], [(1, 1), 3]],
            [[(1, 1), 5]]]

            sage: from quivercombinatorics import *
            sage: C = [[2, -1, 0], [-1, 2, -2], [0, -2, 2]]
            sage: Q = quiver_from_cartan_matrix(C)
            sage: Q.all_representation_types((0, 0, 0), (2, 4, 3))
            [[[(0, 0, 1), 1],
            [(0, 1, 0), 2],
            [(0, 1, 1), 1],
            [(0, 1, 1), 1],
            [(1, 0, 0), 2]],
            [[(0, 0, 1), 1], [(0, 1, 0), 2], [(0, 1, 1), 2], [(1, 0, 0), 2]],
            [[(0, 0, 1), 2], [(0, 1, 0), 3], [(0, 1, 1), 1], [(1, 0, 0), 2]],
            [[(0, 0, 1), 3], [(0, 1, 0), 4], [(1, 0, 0), 2]],
            [[(0, 1, 0), 1],
            [(0, 1, 1), 1],
            [(0, 1, 1), 1],
            [(0, 1, 1), 1],
            [(1, 0, 0), 2]],
            [[(0, 1, 0), 1], [(0, 1, 1), 1], [(0, 1, 1), 2], [(1, 0, 0), 2]],
            [[(0, 1, 0), 1], [(0, 1, 1), 3], [(1, 0, 0), 2]],
            [[(2, 4, 3), 1]]]
                
        """
        x = self._coerce_dimension_vector(x)
        all_decomps = vector_decomposition(x, self.sigma_lambda(l, x))
        all_reps = []
        for decomp in all_decomps:
            current = [[]]
            for pair in decomp:
                next_current = []
                if self.is_imaginary_root(pair[0]):
                    expansions = small_decomposition(pair[0], pair[1])
                    for item in current:
                        for temp in expansions:
                            next_current.append(sorted(item + temp))
                else:
                    for item in current:
                        next_current.append(sorted(item + [pair]))
                current = next_current
            all_reps.extend(current)
        return sorted(all_reps)
    
    def symplectic_leaf_dimension(self, tau):
        r"""Returns the dimension of the symplectic leaf corresponding to the representation type :math:`\tau`


        INPUT:

        - ``tau`` -- a representation type, i.e., a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}Q_0`, and the second is in :math:`\mathbb{N}`

        OUTPUT: The corresponding symplectic leaf dimension, i.e., :math:`2p` applied to every element of :math:`\tau`

        EXAMPLE::
            
            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: Q.symplectic_leaf_dimension([[(1, 1), 2], [(1, 1), 3]])
            4
        
        """
        d = 0
        for pair in tau:
            d += 2 - self.symmetrized_euler_form(pair[0], pair[0])
        return d

    def CB_decomposition(self, l, x):
        r"""Returns the CB decomposition, which is the representation type maximizing :math:`p`


        INPUT:

        - ``l`` -- an element of :math:`\mathbb{Z}Q_0`
        - ``x`` -- an element of :math:`\mathbb{N}Q_0`

        OUTPUT: The CB decomposition with respect to :math:`x`, where ``l`` is :math:`\lambda`

        EXAMPLE::
            
            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: Q.CB_decomposition((1, -1), (5, 5))
            [[(1, 1), 1], [(1, 1), 1], [(1, 1), 1], [(1, 1), 1], [(1, 1), 1]]
        
        """
        all_reps = self.all_representation_types(l, x)
        d_max = 0
        max_rep = []
        for tau in all_reps:
            d = self.symplectic_leaf_dimension(tau)
            if d >= d_max:
                d_max = d
                max_rep = tau
        return max_rep
    
    def quiver_variety_dimension(self, l, x):
        r"""Returns the dimension of the quiver variety


        INPUT:

        - ``l`` -- an element of :math:`\mathbb{Z}Q_0`
        - ``x`` -- an element of :math:`\mathbb{N}Q_0`

        OUTPUT: The dimension of the corresponding quiver variety

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: Q.quiver_variety_dimension((1, -1), (5, 5))
            10
        
        """
        return self.symplectic_leaf_dimension(self.CB_decomposition(l, x))
    
    def codimension_two_leaves(self, l, x):
        r"""Returns the codimension 2 leaves of the quiver variety

        INPUT:

        - ``l`` -- an element of :math:`\mathbb{Z}Q_0`
        - ``x`` -- an element of :math:`\mathbb{N}Q_0`

        OUTPUT: The representation types of the codimension 2 leaves of the quiver variety

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: Q.codimension_two_leaves((0, 0), (5, 5))
            [[[(0, 1), 1],
            [(1, 0), 1],
            [(1, 1), 1],
            [(1, 1), 1],
            [(1, 1), 1],
            [(1, 1), 1]],
            [[(1, 1), 1], [(1, 1), 1], [(1, 1), 1], [(1, 1), 2]]]
        
        """
        d = self.quiver_variety_dimension(l, x)
        all_reps = self.all_representation_types(l, x)
        leaves = []
        for tau in all_reps:
            if self.symplectic_leaf_dimension(tau) + 2 == d:
                leaves = leaves + [tau]
        return leaves
    
    def ext_quiver(self, tau):
        r"""Given a representation type, returns the ext-quiver

        INPUT:

        - ``tau`` -- a representation type, i.e., a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}Q_0`, and the second is in :math:`\mathbb{N}`

        OUTPUT: The :math:`\mathrm{ext}`-quiver

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0, 1], [1, 0]])
            sage: tau = [[(1, 1), 2], [(1, 1), 1], [(1, 1), 1], [(1, 1), 1]]
            sage: Q.ext_quiver(tau)
            a quiver with 4 vertices and 4 arrows

        """
        k = len(tau)
        A = [[0 for i in range(k)] for j in range(k)]
        for i in range(k):
            A[i][i] = self.p_function(tau[i][0])
            for j in range(i):
                A[i][j] = -self.symmetrized_euler_form(tau[i][0], tau[j][0])
        return Quiver(A)

    def is_minimal_imaginary_root(self, x):
        r"""Tests whether a vector is a vector is a minimal imaginary root

        INPUT:

        - ``x`` -- a vector in :math:`\mathbb{Z}Q_0`

        OUTPUT: ``True`` if ``x`` is a minimal imaginary root, and ``False`` otherwise

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(2)
            sage: Q.is_minimal_imaginary_root((1, 1))
            True

        """
        return self.is_imaginary_root(x) and not any([self.is_root(b) and self.is_imaginary_root(b) for b in self.all_subdimension_vectors(x, proper=True, nonzero=True)])

    def all_minimal_imaginary_positive_roots(self, v):
        r"""Returns all minimal imaginary positive roots up to a bound :math:`v`

        INPUT:

        - ``v`` -- a vector in :math:`\mathbb{N}Q_0`

        OUTPUT: A list of all minimal imaginary positive roots up to a bound :math:`v`

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(2)
            sage: Q.all_minimal_imaginary_positive_roots((3, 3))
            [(1, 1)]

        """
        return list(
            filter(
                lambda a: self.is_minimal_imaginary_root(a),
                self.all_subdimension_vectors(v, proper=True, nonzero=True)       
            )
        )
    
    def all_subminimal_representation_types(self, v):
        r"""Returns all subminimal representation types up to a bound :math:`v`, and classifies them by affine Dynkin type, or :math:`a_{r-1}`, :math:`c_{g}`, :math:`m_{g}`

        INPUT:

        - ``v`` -- a vector in :math:`\mathbb{N}Q_0`

        OUTPUT: A list of 2-tuples, the first element is a subminimal representation type up to a bound :math:`v`, the second element is its corresponding affine Dynkin type, or :math:`a_{r-1}`, :math:`c_{g}`, :math:`m_{g}`

        EXAMPLES::

            sage: from quivercombinatorics import *
            sage: Q = CyclicQuiver(3)
            sage: Q.all_subminimal_representation_types((3, 3, 3))
            [([[(0, 0, 1), 2], [(0, 1, 0), 2], [(1, 0, 0), 2], [(1, 1, 1), 1]], 'A_{2}'),
            ([[(0, 0, 1), 1], [(0, 1, 0), 1], [(1, 0, 0), 1], [(1, 1, 1), 2]], 'A_{2}'),
            ([[(1, 1, 1), 3]], 'A_{2}')]

            sage: from quivercombinatorics import *
            sage: Q = LoopQuiver(3)
            sage: Q.all_subminimal_representation_types((5))
            [([[(1), 1], [(1), 4]], 'm_{3}'), ([[(1), 2], [(1), 3]], 'm_{3}')]

            sage: from quivercombinatorics import *
            sage: A = [[0, 1, 0, 3], [1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 0]]
            sage: Q = Quiver(A)
            sage: Q.all_subminimal_representation_types((1, 1, 4, 1))
            [([[(0, 0, 0, 1), 1], [(0, 0, 1, 0), 4], [(1, 1, 0, 0), 1]], 'A_{1}'),
            ([[(0, 0, 1, 0), 4], [(0, 1, 0, 0), 1], [(1, 0, 0, 1), 1]], 'a_{2}'),
            ([[(0, 0, 0, 1), 1],
            [(0, 0, 1, 0), 1],
            [(0, 0, 1, 0), 3],
            [(0, 1, 0, 0), 1],
            [(1, 0, 0, 0), 1]],
            'm_{2}'),
            ([[(0, 0, 0, 1), 1],
            [(0, 0, 1, 0), 2],
            [(0, 0, 1, 0), 2],
            [(0, 1, 0, 0), 1],
            [(1, 0, 0, 0), 1]],
            'c_{2}')]

        """
        rep_types_with_labels = []
        if isinstance(v, (int, Integer)):
            v = (v,)
        if isinstance(v, int):
            v = vector([v])
        else:
            v = vector(v)
        n = len(v)
        A = self.adjacency_matrix()
        G = DiGraph(matrix(A))
        G.remove_loops()
        Q_without_loops = Quiver.from_digraph(G)
        minimal_imaginary_positive_roots_up_to_v = list(                                  
            filter(
                lambda a: Q_without_loops.is_minimal_imaginary_root(a),
                Q_without_loops.all_subdimension_vectors(v, proper=False, nonzero=True)       
            )
        )
        for delta in minimal_imaginary_positive_roots_up_to_v:
            k = 1
            while all(v[i] - k * delta[i] >= 0 for i in range(n)):
                rep_type = [[delta, k]]
                for i in range(n):
                    if v[i] - k * delta[i]:
                        rep_type = rep_type + [[vector([1 if j == i else 0 for j in range(n)]), v[i] - k * delta[i]]]
                supp = Q_without_loops.support(delta)
                subquiver = Q_without_loops.full_subquiver(supp)
                label = ""
                if len(supp) == 2:
                    r = subquiver.adjacency_matrix()[0][1] + subquiver.adjacency_matrix()[1][0]
                    if r == 2:
                        label = f"A_{{1}}"
                elif len(supp) >= 3 and CartanMatrix(subquiver.cartan_matrix()).is_affine():
                    subquiver_type = CartanMatrix(subquiver.cartan_matrix()).cartan_type()
                    label = f"{subquiver_type.type()}_{{{subquiver_type.rank() - 1}}}"
                if label: rep_types_with_labels.append((sorted(rep_type),label))
                k = k + 1
        for i in range(n):
            for j in range(i):
                r = Q_without_loops.adjacency_matrix()[i][j] + Q_without_loops.adjacency_matrix()[j][i]
                if r >= 3 and v[i] > 0 and v[j] > 0:
                    delta = vector([1 if l == i or l == j else 0 for l in range(n)])
                    k = 1
                    while v[i] >= k and v[j] >= k:
                        rep_type = [[delta, k]]
                        for l in range(n):
                            if v[l] - k * delta[l]:
                                rep_type = rep_type + [[vector([1 if s == l else 0 for s in range(n)]), v[l] - k * delta[l]]]
                        rep_types_with_labels.append((sorted(rep_type),f"a_{{{r-1}}}"))
                        k = k + 1    
        for i in range(n):
            g = self.adjacency_matrix()[i][i]
            if g >= 1:
                for a in range(1, v[i] // 2 + 1):
                    b = v[i] - a
                    rep_type = [[vector([1 if j == i else 0 for j in range(n)]), a], [vector([1 if j == i else 0 for j in range(n)]), b]]
                    for k in range(n):
                        if k != i and v[i]:
                            rep_type = rep_type + [[vector([1 if j == k else 0 for j in range(n)]), v[k]]]
                        
                    label = ""
                    if a == b:
                        label = f"c_{{{g}}}"
                    else:
                        label = f"m_{{{g}}}"
                    rep_types_with_labels.append((sorted(rep_type),label))
        return rep_types_with_labels
    
    def minimal_degenerations(self, L):
        r"""Returns all minimal degenerations of a a symplectic leaf, corresponding to the representation type :math:`L`

        INPUT:

        - ``L`` -- a representation type, i.e., a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}Q_0`, and the second is in :math:`\mathbb{N}`

        OUTPUT: A list of representation types

        EXAMPLES::

            sage: from quivercombinatorics import *
            sage: A = [[0, 1, 0, 3], [1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 0]]
            sage: Q = Quiver(A)
            sage: Q.minimal_degenerations([[(0, 0, 0, 1), 1], [(0, 0, 1, 0), 4], [(1, 1, 0, 0), 1]])
            [([[(0, 0, 1, 0), 4], [(1, 1, 0, 1), 1]], '$a_{2}(1)$'),
            ([[(0, 0, 0, 1), 1], [(0, 0, 1, 0), 1], [(0, 0, 1, 0), 3], [(1, 1, 0, 0), 1]],
            '$m_{2}(1)$'),
            ([[(0, 0, 0, 1), 1], [(0, 0, 1, 0), 2], [(0, 0, 1, 0), 2], [(1, 1, 0, 0), 1]],
            '$c_{2}(1)$')]
        
        """
        Qtilde = self.ext_quiver(L)
        n = ext_dimension_vector(L)
        rep_types = Qtilde.all_subminimal_representation_types(n)
        lifted_rep_types_with_counts = []
        for tau in rep_types:
            rep_type = D_lifting(tau, L)
            if lifted_rep_types_with_counts and lifted_rep_types_with_counts[-1][0] == rep_type:
                lifted_rep_types_with_counts[-1][1] += 1
            else:
                lifted_rep_types_with_counts.append([rep_type, 1])
        return [(sorted(tau[0][0]),fr'{tau[0][1]}({tau[1]})') for tau in lifted_rep_types_with_counts]
    
    def all_decompositions(self, v):
        r"""Constructs all decompositions of a given dimension vector :math:`v`

        INPUT:

        - ``v`` -- an element of :math:`\mathbb{N}Q_0`

        OUTPUT: A list of decompositions of `v`

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: A = [[0, 1, 0, 3], [1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 0]]
            sage: Q = Quiver(A)
            sage: Q.all_decompositions((1, 0, 2, 0))
            [[[(0, 0, 1, 0), 1], [(0, 0, 1, 0), 1], [(1, 0, 0, 0), 1]],
            [[(0, 0, 1, 0), 1], [(1, 0, 1, 0), 1]],
            [[(0, 0, 1, 0), 2], [(1, 0, 0, 0), 1]],
            [[(0, 0, 2, 0), 1], [(1, 0, 0, 0), 1]],
            [[(1, 0, 2, 0), 1]]]

        """
        all_decomps = vector_decomposition(v, self.all_subdimension_vectors(v, proper=False, nonzero=True))
        all_reps = []
        for decomp in all_decomps:
            current = [[]]
            for pair in decomp:
                current = [sorted(temp + item) for item in current for temp in small_decomposition(pair[0], pair[1])]
            all_reps = all_reps + current
        return sorted(all_reps)

    def get_Hasse_diagram(self, l, v, method=1):
        r"""Applies Method 1 or Method 2 to obtain data for the Hasse diagram of minimal degenerations

        INPUT:

        - ``l`` -- an element of :math:`\mathbb{Z}Q_0`
        - ``v`` -- an element of :math:`\mathbb{N}Q_0`
        - ``method`` -- When set to ``1``, Method 1 will be used to obtain the Hasse diagram of minimal degenerations. When set to ``2``, Method 1 will be used to obtain the Hasse diagram of minimal degenerations.

        OUTPUT: 
        
        - ``leaves_poset`` -- the underlying poset of the Hasse diagram of minimal degenerations
        - ``all_leaves`` -- a list of all symplectic leaves, with respect to their representation types, so a list of representation types, i.e., a list whose elements are 2-tuples, the first element is in :math:`\mathbb{N}Q_0`, and the second is in :math:`\mathbb{N}`
        - ``dimensions`` -- a list of dimensions of all the symplectic leaves, using ``symplectic_leaf_dimension``
        - ``edge_labels`` -- a list of 3-tuples, the first two entries of each tuple define the edge, and the last entry is the edge label of the corresponding minimal degeneration, classified by ``all_subminimal_representation_types``

        EXAMPLE::

            sage: from quivercombinatorics import *
            sage: C = [[2, -1, 0], [-1, 2, -2], [0, -2, 2]]
            sage: Q = quiver_from_cartan_matrix(C)
            sage: Q.get_Hasse_diagram((0, 0, 0), (2, 4, 3))
            (Finite poset containing 8 elements,
            [[[(0, 0, 1), 1],
            [(0, 1, 0), 2],
            [(0, 1, 1), 1],
            [(0, 1, 1), 1],
            [(1, 0, 0), 2]],
            [[(0, 0, 1), 1], [(0, 1, 0), 2], [(0, 1, 1), 2], [(1, 0, 0), 2]],
            [[(0, 0, 1), 2], [(0, 1, 0), 3], [(0, 1, 1), 1], [(1, 0, 0), 2]],
            [[(0, 0, 1), 3], [(0, 1, 0), 4], [(1, 0, 0), 2]],
            [[(0, 1, 0), 1],
            [(0, 1, 1), 1],
            [(0, 1, 1), 1],
            [(0, 1, 1), 1],
            [(1, 0, 0), 2]],
            [[(0, 1, 0), 1], [(0, 1, 1), 1], [(0, 1, 1), 2], [(1, 0, 0), 2]],
            [[(0, 1, 0), 1], [(0, 1, 1), 3], [(1, 0, 0), 2]],
            [[(2, 4, 3), 1]]],
            [4, 2, 2, 0, 6, 4, 2, 8],
            [(0, 4, 'A_{1}(1)'),
            (1, 5, 'A_{1}(1)'),
            (1, 0, 'c_{1}(1)'),
            (2, 0, 'A_{1}(1)'),
            (2, 5, 'A_{1}(1)'),
            (3, 2, 'A_{1}(1)'),
            (3, 1, 'A_{1}(1)'),
            (3, 6, 'A_{1}(1)'),
            (4, 7, 'D_{4}(1)'),
            (5, 4, 'c_{1}(1)'),
            (6, 5, 'm_{1}(1)')])
        
        """
        if method != 1 and method != 2:
            raise ValueError("Method can only take values 1 or 2.")
        if isinstance(l, (int, Integer)):
            l = (l,)
        if isinstance(v, (int, Integer)):
            v = (v,)
        all_leaves = self.all_representation_types(l, v)
        dimensions = [self.symplectic_leaf_dimension(tau) for tau in all_leaves]
        num_of_leaves = len(all_leaves)
        edge_labels = []
        relations = []
        if method == 1:
            elements = list(range(num_of_leaves))
            for i in range(num_of_leaves):
                for L in self.minimal_degenerations(all_leaves[i]):
                    j = all_leaves.index(L[0])
                    relations.append((elements[i], elements[j]))
                    edge_labels.append((elements[i], elements[j], L[1]))
            leaves_poset = Poset((elements, relations), cover_relations = False)
        elif method == 2:
            decomps = self.all_decompositions(v)
            indices = [decomps.index(x) for x in all_leaves]
            num_of_decomps = len(decomps)
            elements = range(num_of_leaves)
            for i in range(num_of_decomps):
                for j in range(i):
                    if is_direct_successor(decomps[i], decomps[j]):
                        relations.append((i, j))
                    elif is_direct_successor(decomps[j], decomps[i]):
                        relations.append((j, i))
            P = Poset((elements, relations), cover_relations=False)
            relabelling = {indices[i]: i for i in range(num_of_leaves)}
            leaves_poset = P.subposet(indices).relabel(relabelling)    
            for i in range(num_of_leaves):
                for tau, label in self.minimal_degenerations(all_leaves[i]):
                    if tau in all_leaves:
                        j = all_leaves.index(tau)
                        if leaves_poset.covers(i, j):
                            edge_labels.append((i, j, label))  
        return leaves_poset, all_leaves, dimensions, edge_labels

    def plot_Hasse_diagram(self, l, v, method=1, format="tikz", filename="output"):
        r"""Applies Method 1 or Method 2 to plot the Hasse diagram for minimal degenerations

        INPUT:

        - ``l`` -- an element of :math:`\mathbb{Z}Q_0`
        - ``v`` -- an element of :math:`\mathbb{N}Q_0`
        - ``method`` -- When set to ``1``, Method 1 will be used to obtain the Hasse diagram of minimal degenerations. When set to ``2``, Method 1 will be used to obtain the Hasse diagram of minimal degenerations.
        - ``format`` -- Takes values "tikz", "dot", and "sage", depending on whether you want a ``.tex`` file with *tikzpicture*, a ``.dot`` file, or a sage Hasse diagram output
        - ``filename`` -- filename of output

        OUTPUT: 
        
        - the Hasse diagram for minimal degenerations

        EXAMPLE::

            sage: C = [[2, -1, 0], [-1, 2, -2], [0, -2, 2]]
            sage: Q = quiver_from_cartan_matrix(C)
            sage: Q.plot_Hasse_diagram((0, 0, 0), (2, 4, 3))
            \documentclass[tikz]{standalone}
            \begin{document}

            \begin{tikzpicture}[>=latex,line join=bevel,]
            %%
            \begin{scope}
                \pgfsetstrokecolor{black}
            ---
            77 lines not printed (5057 characters in total).
            Use print to see the full content.
            ---
                \draw (221.5bp,185.5bp) node {$c_{1}(1)$};
                \draw [->] (6) ..controls (293.59bp,110.99bp) and (289.95bp,119.79bp)  .. (284.0bp,126.0bp) .. controls (277.86bp,132.41bp) and (269.81bp,137.25bp)  .. (5);
                \draw (320.5bp,118.5bp) node {$m_{1}(1)$};
            %
            \end{tikzpicture}
            \end{document}
        
        """
        if format not in ["dot", "tikz", "sage"]:
            raise ValueError("Possible formats are sage, dot and tikz.")
        P, _, dimensions, edge_labels = self.get_Hasse_diagram(l, v, method)
        H = P.hasse_diagram()
        possible_dims = sorted(list(set(dimensions)))
        if format == "sage":
            heights = {
                j: [
                    i 
                    for i in range(len(dimensions)) 
                    if dimensions[i] == possible_dims[j]
                ] 
                for j in range(len(possible_dims))
            }
            pos = H.layout_acyclic_dummy(heights=heights)
            H.set_pos(pos)
            for i, j, label in edge_labels:
                H.set_edge_label(i, j, label)
            return H
        elif format == "dot" or format == "tikz":
            lines = []
            lines.append('digraph G {')
            lines.append('    graph [splines=true, overlap=false, rankdir="BT", d2toptions="-e utf8", labeldistance=0.5, nodesep = 0.1, ranksep = 0.2]')
            lines.append('    node [shape=none]')
            string = '    '
            for dim in possible_dims:
                lines.append(f'    dim{dim} [texlbl="$\dim={dim}$"]')
                string = string + f'dim{dim} -> '
            string = string[:-3] + '[style=invis]'
            lines.append(string)
            for i in range(len(dimensions)):
                lines.append(f'    {i} [texlbl="$L_{{{i}}}$"]')
            for dim in possible_dims:
                string = '    {rank=same;'+f' dim{dim}'
                for i in range(len(dimensions)):
                    if dimensions[i] == dim:
                        string += f'; {i}'
                string += '}'
                lines.append(string)
            for i, j, label in edge_labels:
                if H.has_edge(i, j):
                    lines.append(f'    {i} -> {j} [texlbl="${label}$", label="{label}", lp="0,0"]')
            lines.append('}')
            s = '\n'.join(lines)
            if format == "tikz":
                t = TikzPicture(dot2tex.dot2tex(s, format='tikz', figonly='True', prog='dot', rankdir='down'))
                _ = t.tex(filename+".tex")
                return t
            else:
                with open(filename+".dot", "w") as f:
                    f.write(s)
                return s
    
    BaseQuiver.p_function = p_function
    BaseQuiver.R_lambda_plus = R_lambda_plus
    BaseQuiver.N_set = N_set
    BaseQuiver.sigma_lambda = sigma_lambda
    BaseQuiver.vector_decomposition = vector_decomposition
    BaseQuiver.small_decomposition = small_decomposition
    BaseQuiver.all_representation_types = all_representation_types
    BaseQuiver.symplectic_leaf_dimension = symplectic_leaf_dimension
    BaseQuiver.CB_decomposition = CB_decomposition
    BaseQuiver.quiver_variety_dimension = quiver_variety_dimension
    BaseQuiver.codimension_two_leaves = codimension_two_leaves
    BaseQuiver.ext_quiver = ext_quiver
    BaseQuiver.ext_dimension_vector = ext_dimension_vector
    BaseQuiver.is_minimal_imaginary_root = is_minimal_imaginary_root
    BaseQuiver.all_minimal_imaginary_positive_roots = all_minimal_imaginary_positive_roots
    BaseQuiver.all_subminimal_representation_types = all_subminimal_representation_types
    BaseQuiver.D_map = D_map
    BaseQuiver.D_lifting = D_lifting
    BaseQuiver.minimal_degenerations = minimal_degenerations
    BaseQuiver.all_decompositions = all_decompositions
    BaseQuiver.is_direct_successor = is_direct_successor
    BaseQuiver.get_Hasse_diagram = get_Hasse_diagram
    BaseQuiver.plot_Hasse_diagram = plot_Hasse_diagram