#Packages required by QuiverTools
from sage.arith.misc import gcd
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.matrix.special import zero_matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.structure.element import Element
from quiver import Quiver as BaseQuiver

#Extra packages needed for QuiverToolsCombinatorics
from sage.combinat.posets.hasse_diagram import HasseDiagram
from sage.misc.latex import LatexExpr
from collections import defaultdict

def quiver_from_cartan_matrix(C):
    r"""Returns the quiver `Q` given by a Cartan matrix `C`.
        
        INPUT:

        - ``C`` -- a square matrix with entries in :math:`\mathbb{Z}`.

        OUTPUT: The quiver `Q` given by a Cartan matrix `C`.

        EXAMPLE::
        
            sage: from quivercombinatorics import *
            sage: C = [[2,-1,0], [-1,2,-2], [0,-2,2]]
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
    r"""Returns a randomly generated quiver.
        
        INPUT:

        - ``vertices`` -- number of vertices you want in the quiver, in :math:`\mathbb{N}_0`.
        - ``max_arrows_per_edge`` -- maximum of number of arrows you want in the quiver, in :math:`\mathbb{N}_0`.

        OUTPUT: A randomly generated quiver with the prescribed vertices and maximum number of edges.

        EXAMPLES::
        
            sage: random_quiver(3, 2)
            a quiver with 3 vertices and 6 arrows

            sage: random_quiver(46, 25)
            a quiver with 46 vertices and 26388 arrows
            
        
        """
    A = [[random.randint(0, max_arrows_per_edge) for _ in range(vertices)] for _ in range(vertices)]
    G = DiGraph(matrix(A))
    return Quiver(A)

def N_set(S, v):
    r"""For a list of vectors `S` in :math:`\mathbb{Z}^n`, returns all possible sums of vectors in `S` as a list, up to the upper bound `v`.
        
        INPUT:

        - ``S`` -- a list of vectors `S` in :math:`\mathbb{Z}^n`.
        - ``v`` -- vector in :math:`\mathbb{Z}^n`.

        OUTPUT: All possible sums of vectors in `S` as a list, up to the upper bound `v`.

        EXAMPLE::
        
            sage: from quivercombinatorics import *
            sage: N_set([[0,1]],[0,3])
            [(0, 1), (0, 2), (0, 3)]
        
        """
    if isinstance(v, int):                                                                  #make sure v is of the right type
        v = vector([v])
    else:
        v = vector(v)
    reachable = {tuple(a) for a in S}
    changed = True                      #flag to make sure duplicates don't end up in the generated set
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

class Quiver(BaseQuiver):

    def p_function(self, x):
        r"""Outputs the function p(x)=1-0.5(x,x), where (x,x) is the symmetrized Euler form.
        
        INPUT:

        - ``x`` -- an element of :math:`\mathbb{Z}Q_0`.

        OUTPUT: The `p` function `p(x)=1-0.5(x,x)`

        EXAMPLE:

        The `p` function evaluated at `(2,3)` of the doubled A2 quiver::

            sage: from quivercombinatorics import *
            sage: Q = Quiver([[0,1], [1,0]])
            sage: Q.p_function((2,3))
            0.0
        
        """
        return 1 - 0.5 * self.symmetrized_euler_form(x, x)

    def full_subquiver(self, vertices):
        r"""Returns the full subquiver supported on the given set of vertices

        INPUT:

        - ``vertices``: list of vertices for the subquiver

        OUTPUT: the full subquiver on the specified vertices

        EXAMPLES:

        Some basic examples::

            sage: from quivercombinatorics import *
            sage: Q = ThreeVertexQuiver(2, 3, 4)
            sage: print(Q.full_subquiver([0, 1]))
            full subquiver of an acyclic 3-vertex quiver of type (2, 3, 4)
            adjacency matrix:
            [0 2]
            [0 0]
            sage: print(Q.full_subquiver([0, 2]))
            full subquiver of an acyclic 3-vertex quiver of type (2, 3, 4)
            adjacency matrix:
            [0 3]
            [0 0]

        If we specified a non-standard labeling on the vertices we must use it::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q == ThreeVertexQuiver(2, 3, 4)
            True
            sage: print(Q.full_subquiver(["a", "b"]))
            a quiver with 2 vertices and 2 arrows
            adjacency matrix:
            [0 2]
            [0 0]
            sage: print(Q.full_subquiver(["a", "c"]))
            a quiver with 2 vertices and 3 arrows
            adjacency matrix:
            [0 3]
            [0 0]

        """
        name = None
        if self.get_custom_name():
            name = "full subquiver of " + self.get_custom_name()

        return Quiver.from_digraph(self.graph().subgraph(vertices=vertices), name)

    """
    Dimension vectors and roots
    """

    def zero_vector(self):
        r"""
        Returns the zero dimension vector.

        The output is adapted to the vertices.

        OUTPUT: the zero dimension vector

        EXAMPLES:

        Usually it is an actual vector::

            sage: from quivercombinatorics import *
            sage: KroneckerQuiver(3).zero_vector()
            (0, 0)
            sage: type(KroneckerQuiver(3).zero_vector())
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

        But if the quiver has custom vertex labels it is a dict::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q.zero_vector()
            {'a': 0, 'b': 0, 'c': 0}
        """
        return self.__zero_vector()

    @cached_method
    def __zero_vector(self):
        r"""The cacheable implementation of :meth:`Quiver.zero_vector`

        EXAMPLES:

        The zero vector of the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: KroneckerQuiver(3).zero_vector()
            (0, 0)

        If we specified a non-standard labeling on the vertices it is used::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q._Quiver__zero_vector()
            {'a': 0, 'b': 0, 'c': 0}
        """
        if self.__has_vertex_labels():
            return {i: 0 for i in self.vertices()}

        return vector([0] * self.number_of_vertices(), immutable=True)

    def thin_dimension_vector(self):
        r"""
        Returns the thin dimension vector, i.e., all ones

        The output is adapted to the vertices.

        OUTPUT: the thin dimension vector

        EXAMPLES:

        Usually it is an actual vector::

            sage: from quivercombinatorics import *
            sage: KroneckerQuiver(3).thin_dimension_vector()
            (1, 1)
            sage: type(KroneckerQuiver(3).thin_dimension_vector())
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

        But if the quiver has custom vertex labels it is a dict::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q.thin_dimension_vector()
            {'a': 1, 'b': 1, 'c': 1}

        """
        return self.__thin_dimension_vector()

    @cached_method
    def __thin_dimension_vector(self):
        r"""The cacheable implementation of :meth:`Quiver.thin_dimension_vector`

        EXAMPLES:

        The thin dimension vector of the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: KroneckerQuiver(3).thin_dimension_vector()
            (1, 1)

        If we specified a non-standard labeling on the vertices it is used::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q._Quiver__thin_dimension_vector()
            {'a': 1, 'b': 1, 'c': 1}
        """
        if self.__has_vertex_labels():
            return {i: 1 for i in self.vertices()}

        return vector([1] * self.number_of_vertices(), immutable=True)

    def simple_root(self, i):
        r"""
        Returns the simple root at the vertex ``i``

        The output is adapted to the vertices.

        OUTPUT: the simple root at the vertex ``i``

        EXAMPLES:

        Usually it is an actual vector::

            sage: from quivercombinatorics import *
            sage: KroneckerQuiver(3).simple_root(1)
            (0, 1)
            sage: type(KroneckerQuiver(3).simple_root(1))
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

        But if the quiver has custom vertex labels it is a dict::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q.simple_root("b")
            {'a': 0, 'b': 1, 'c': 0}
        """
        return self.__simple_root(i)

    @cached_method
    def __simple_root(self, i):
        r"""The cacheable implementation of :meth:`Quiver.simple_root`

        EXAMPLES:

        The simple root at the source of the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: KroneckerQuiver(3).simple_root(0)
            (1, 0)

        If we specified a non-standard labeling on the vertices it is used::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q._Quiver__simple_root("a")
            {'a': 1, 'b': 0, 'c': 0}
        """
        if self.__has_vertex_labels():
            root = {i: 0 for i in self.vertices()}
            root[i] = 1

            return root

        root = vector([0] * self.number_of_vertices())
        root[i] = 1
        root.set_immutable()

        return root

    def is_root(self, x) -> bool:
        r"""Checks whether ``x`` is a root of the underlying diagram of the quiver.

        A root is a non-zero vector `x` in :math:`\mathbb{Z}Q_0` such that
        the Tits form of `x` is at most 1.

        INPUT:

        - ``x``: integer vector

        OUTPUT: whether ``x`` is a root

        EXAMPLES:

        Some roots and non-roots for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_root((2, 3))
            True
            sage: Q.is_root(Q.zero_vector())
            False
            sage: Q.is_root((4, 1))
            False

        """
        x = self._coerce_vector(x)

        return any(x) and self.tits_form(x) <= 1

    def is_real_root(self, x) -> bool:
        r"""Checks whether ``x`` is a real root of the underlying diagram of the quiver.

        A root is called real if its Tits form equals 1.

        INPUT:

        - ``x``: integer vector

        OUTPUT: whether ``x`` is a real root

        EXAMPLES:

        Some real and non-real for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_real_root((2, 3))
            False
            sage: Q.is_real_root(Q.zero_vector())
            False
            sage: Q.is_real_root((3, 1))
            True

        """
        x = self._coerce_vector(x)

        return self.tits_form(x) == 1

    def is_imaginary_root(self, x) -> bool:
        r"""Checks whether ``x`` is a imaginary root of the quiver.

        A root is called imaginary if its Tits form is non-positive.

        INPUT:

        - ``x``: integer vector

        OUTPUT: whether ``x`` is an imaginary root

        EXAMPLES:

        Some imaginary roots and non imaginary roots for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_imaginary_root((2, 3))
            True
            sage: Q.is_imaginary_root(Q.zero_vector())
            False
            sage: Q.is_imaginary_root((4, 1))
            False

        """
        x = self._coerce_vector(x)

        return any(x) and self.tits_form(x) <= 0

    def is_schur_root(self, d) -> bool:
        r"""Checks if ``d`` is a Schur root.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: whether ``d`` is an imaginary root

        A Schur root is a dimension vector which admits a Schurian representation,
        i.e., a representation whose endomorphism ring is the field itself.
        It is necessarily indecomposable.

        By MR1162487_ :math:`{\bf d}` is a Schur root if and only if it admits a stable
        representation for the canonical stability parameter.

        .. _MR1162487: https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487

        EXAMPLES:

        The dimension vector `(2, 3)` is Schurian for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.is_schur_root([2, 3])
            True

        Examples from Derksen--Weyman's book (Example 11.1.4)::

            sage: from quivercombinatorics import *
            sage: Q = ThreeVertexQuiver(1, 1, 1)
            sage: Q.is_schur_root((1, 1, 2))
            True
            sage: Q.is_schur_root((1, 2, 1))
            False
            sage: Q.is_schur_root((1, 1, 1))
            True
            sage: Q.is_schur_root((2, 2, 2))
            False

        """
        d = self._coerce_dimension_vector(d)
        theta = self.canonical_stability_parameter(d)

        return self.has_stable_representation(d, theta)

    def slope(self, d, theta=None, denom=sum):
        r"""
        Returns the slope of ``d`` with respect to ``theta``

        The slope is defined as the value of ``theta(d)`` divided by the total dimension
        of `d` ``sum(d)``. It is possible to vary the denominator, to use a function
        more general than the sum.

        INPUT:

        - ``d`` -- dimension vector

        - ``theta`` -- (default: canonical stability parameter) stability parameter

        - ``denom`` -- (default: sum) the denominator function

        OUTPUT: the slope of ``d`` with respect to ``theta`` and optional ``denom``

        EXAMPLES:

        Some slopes for the Kronecker quiver, first for the canonical stability
        parameter, then for some other::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: d = (2, 3)
            sage: Q.slope(d, (9, -6))
            0
            sage: Q.slope(d)
            0
            sage: Q.slope(d, (2, -2))
            -2/5

        We can use for instance a constant denominator::

            sage: constant = lambda di: 1
            sage: Q.slope(d, Q.canonical_stability_parameter(d), denom=constant)
            0

        The only dependence on the quiver is the set of vertices, so if we don't
        use vertex labels, the choice of quiver doesn't matter::

            sage: d, theta = (2, 3), (9, -6)
            sage: KroneckerQuiver(3).slope(d, theta)
            0

        """
        d = self._coerce_dimension_vector(d)
        assert denom(d) > 0, "denominator needs to be strictly positive on ``d``"

        if theta is None:
            theta = self.canonical_stability_parameter(d)
        theta = self._coerce_vector(theta)

        return (theta * d) / denom(d)

    def is_subdimension_vector(self, e, d):
        r"""
        Determine whether ``e`` is a subdimension vector of ``d``

        INPUT:

        -- ``e`` -- dimension vector

        -- ``d`` -- dimension vector

        OUTPUT: whether ``e`` is a subdimension vector of ``d``

        EXAMPLES:

        Some basic examples::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_subdimension_vector((1, 2), (2, 3))
            True
            sage: Q.is_subdimension_vector((2, 3), (2, 3))
            True
            sage: Q.is_subdimension_vector((6, 6), (2, 3))
            False

        We can also work with vertex labels::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: d = {"a" : 3, "b" : 3, "c" : 3}
            sage: e = {"a" : 1, "b" : 2, "c" : 3}
            sage: Q.is_subdimension_vector(e, d)
            True
            sage: Q.is_subdimension_vector(d, e)
            False

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        return all(ei <= di for (ei, di) in zip(e, d))

    def _deglex_key(self, e, b=None) -> int:
        r"""
        An integer representation of a dimension vector

        This is the base-b expansion of a dimension vector.

        This is a function which satisfies

            e <_{deglex} d iff deglex_key(e) < deglex_key(d),

        provided that b >> 0.

        For b >> 0 the deglex order is a _total_ order which extends the usual
        entry-wise partial order on dimension vectors.

        INPUT:

        - ``e`` -- dimension vector

        - ``b`` -- the "base" of the key (default: `max(e)+1`)

        OUTPUT: the base-`b` expansion of the dimension vector

        EXAMPLES:

        If we let `b` be the largest entry plus one we get a good key, at least for
        subdimension vectors of the original one::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: d = (2, 3)
            sage: Q._deglex_key(d, max(d) + 1)
            91
            sage: d = (3, 3)
            sage: Q._deglex_key(d)
            111

        """
        e = self._coerce_dimension_vector(e)
        if b is None:
            b = max(e) + 1

        n = self.number_of_vertices()

        return (
            sum(ei * b ** (n - i - 1) for (i, ei) in enumerate(e))
            + sum(self._coerce_dimension_vector(e)) * b**n
        )

    def all_subdimension_vectors(
        self, d, proper=False, nonzero=False, forget_labels=False
    ):
        r"""
        Returns the list of all subdimension vectors of ``d``.

        INPUT:

        - ``d`` -- dimension vector

        - ``proper`` (default: False) -- whether to exclude ``d``

        - ``nonzero`` (default: False) -- whether to exclude the zero vector

        - ``forget_labels`` (default: False) -- whether to forget the vertex labels

        OUTPUT: all subdimension vectors of ``d`` (maybe excluding zero and/or ``d``)

        EXAMPLES:

        The usual use cases::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.all_subdimension_vectors((2, 3))
                [(0, 0),
                 (0, 1),
                 (0, 2),
                 (0, 3),
                 (1, 0),
                 (1, 1),
                 (1, 2),
                 (1, 3),
                 (2, 0),
                 (2, 1),
                 (2, 2),
                 (2, 3)]
            sage: Q.all_subdimension_vectors((2, 3), proper=True)
                [(0, 0),
                 (0, 1),
                 (0, 2),
                 (0, 3),
                 (1, 0),
                 (1, 1),
                 (1, 2),
                 (1, 3),
                 (2, 0),
                 (2, 1),
                 (2, 2)]
            sage: Q.all_subdimension_vectors((2, 3), nonzero=True)
                [(0, 1),
                 (0, 2),
                 (0, 3),
                 (1, 0),
                 (1, 1),
                 (1, 2),
                 (1, 3),
                 (2, 0),
                 (2, 1),
                 (2, 2),
                 (2, 3)]
            sage: Q.all_subdimension_vectors((2, 3), proper=True, nonzero=True)
                [(0, 1),
                 (0, 2),
                 (0, 3),
                 (1, 0),
                 (1, 1),
                 (1, 2),
                 (1, 3),
                 (2, 0),
                 (2, 1),
                 (2, 2)]

        Some exceptional cases::

            sage: Q.all_subdimension_vectors(Q.zero_vector())
            [(0, 0)]
            sage: Q.all_subdimension_vectors(Q.zero_vector(), proper=True)
            []

        If we work with labeled vertices, then we get a list of dicts::

            sage: Q = Quiver.from_string("a---b", forget_labels=False)
            sage: Q.all_subdimension_vectors((1, 2))
            [{'a': 0, 'b': 0},
             {'a': 0, 'b': 1},
             {'a': 0, 'b': 2},
             {'a': 1, 'b': 0},
             {'a': 1, 'b': 1},
             {'a': 1, 'b': 2}]
            sage: Q.all_subdimension_vectors((1, 2), forget_labels=True)
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

        """
        assert self._is_dimension_vector(d), "``d`` needs to be a dimension vector"

        # if zero dimension vector we deal with it separately
        if sum(self._coerce_dimension_vector(d)) == 0:
            if proper or nonzero:
                return []
            return [d]

        vectors = list(cartesian_product([range(di + 1) for di in d]))

        if proper:
            vectors = vectors[:-1]
        if nonzero:
            vectors = vectors[1:]

        if self.__has_vertex_labels() and not forget_labels:
            return list(map(lambda e: dict(zip(self.vertices(), e)), vectors))

        return list(map(vector, vectors))

    def is_theta_coprime(self, d, theta=None) -> bool:
        r"""Checks if ``d`` is ``theta``-coprime.

        A dimension vector `d` is :math:`\theta`-coprime if
        :math:`\mu_{\theta}(e)\neq \mu_{\theta}(e)`
        for all proper non-zero subdimension vectors e of d.

        The default value for ``theta`` is the canonical stability parameter,
        see :meth:`canonical_stability_parameter`.

        INPUT:

        - ``d`` -- dimension vector

        - ``theta`` -- (default: canonical stability paramter) stability parameter

        EXAMPLES:

        Examples of coprimality::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: d = (2, 3)
            sage: Q.is_theta_coprime(d, Q.canonical_stability_parameter(d))
            True
            sage: Q.is_theta_coprime(d)
            True
            sage: Q.is_theta_coprime((3, 3), (1, -1))
            False

        """
        if theta is None:
            theta = self.canonical_stability_parameter(d)

        assert self._is_dimension_vector(d), "``d`` needs to be a dimension vector"
        assert self._is_vector(theta), "`theta` needs to be a stability parameter"

        vectors = self.all_subdimension_vectors(d, proper=True, nonzero=True)

        return all(self.slope(d, theta) != self.slope(e, theta) for e in vectors)

    def is_indivisible(self, d) -> bool:
        """
        Checks if the gcd of all the entries of ``d`` is 1

        INPUT:

        -- ``d`` -- dimension vector

        OUTPUT: whether the dimension vector is indivisible

        EXAMPLES:

        Two examples with the Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_indivisible((2, 3))
            True
            sage: Q.is_indivisible((2, 2))
            False

        """
        return gcd(self._coerce_dimension_vector(d)) == 1

    def support(self, d):
        r"""Returns the support of the dimension vector.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: subset of vertices in the underlying graph in the support

        The support is the set :math:`\{ i \in Q_0 \mid d_i > 0 \}`.

        EXAMPLES:

        The support is the set of vertices for which the value of the dimension
        vector is nonzero::

            sage: from quivercombinatorics import *
            sage: Q = ThreeVertexQuiver(2, 0, 4)
            sage: d = (1, 1, 1)
            sage: Q.support(d)
            [0, 1, 2]
            sage: d = (1, 0, 1)
            sage: Q.support(d)
            [0, 2]

        It takes into account vertex labels::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: d = {"a": 2, "b": 3, "c": 0}
            sage: Q.support(d)
            ['a', 'b']

        """
        assert self._is_dimension_vector(d), "``d`` needs to be a dimension vector"

        return [i for i in self.vertices() if d[i] > 0]

    def in_fundamental_domain(self, d, depth=0):
        r"""Checks if a dimension vector is in the fundamental domain.

        The fundamental domain of :math:`Q` is the set of dimension vectors :math:`d`
        such that

        - :math:`\operatorname{supp}(\mathbf{d})` is connected
        - :math:`\langle d,e_i\rangle + \langle e_i,d\rangle\leq 0` for every simple
          root

        Every :math:`d` in the fundamental domain is an imaginary root and the set of
        imaginary roots is the Weyl group saturation of the fundamental domain.
        If :math:`d` is in the fundamental domain then it is Schurian and a general
        representation of dimension vector :math:`d` is stable for the canonical
        stability parameter.

        The optional parameter ``depth`` allows to make the inequality stricter.

        INPUT:

        - ``d``: dimension vector

        - ``depth`` (default: 0) -- how deep the vector should be in the domain

        OUTPUT: whether ``d`` is in the (interior of) the fundamental domain

        EXAMPLES:

        The fundamental domain of the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.in_fundamental_domain((1, 1))
            True
            sage: Q.in_fundamental_domain((1, 2))
            False
            sage: Q.in_fundamental_domain((2, 3))
            True

        The same calculation now with vertex labels::

            sage: Q = Quiver.from_string("a---b", forget_labels=False)
            sage: Q.in_fundamental_domain({"a" : 1, "b" : 1})
            True
            sage: Q.in_fundamental_domain({"a" : 1, "b" : 2})
            False
            sage: Q.in_fundamental_domain({"a" : 2, "b" : 3})
            True

        We test for dimension vectors in the strict interior, where the depth is
        equal to 1::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.in_fundamental_domain((1, 1), depth=1)
            True
            sage: Q.in_fundamental_domain((2, 3), depth=1)
            False

        """
        assert self._is_dimension_vector(d), "``d`` needs to be a dimension vector"

        # check if `\langle d,e_i\rangle + \langle e_i,d\rangle \leq 0`
        # for all vertices `i\in Q_0`
        inequality = all(
            self.symmetrized_euler_form(d, self.simple_root(i)) <= -depth
            for i in self.vertices()
        )

        # check if the support is connected
        connected = self.full_subquiver(self.support(d)).is_connected()

        return inequality and connected

    def division_order(self, d, e):
        r"""
        Checks if :math:`d\ll e`

        This means that

        - :math:`d_i \leq e_i` for every source `i`,
        - :math:`d_j \geq e_j` for every sink `j`, and
        - :math:`d_k = e_k` for every vertex `k` which is neither a source nor a sink.

        This is used when dealing with Chow rings of quiver moduli, see also
        :meth:`QuiverModuli.chow_ring` and
        :meth:`QuiverModuli._all_minimal_forbidden_subdimension_vectors`.

        EXAMPLES:

        The division order on some dimension vectors for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = (1, 1)
            sage: e = (2, 1)
            sage: f = (2, 2)
            sage: Q.division_order(d, e)
            True
            sage: Q.division_order(e, d)
            False
            sage: Q.division_order(d, f)
            False
            sage: Q.division_order(f, d)
            False
            sage: Q.division_order(e, f)
            False
            sage: Q.division_order(f, e)
            True

        The division order on some dimension vectors for a 3-vertex quiver::

            sage: Q = ThreeVertexQuiver(2, 2, 2)
            sage: d = (1, 1, 1)
            sage: e = (1, 2, 1)
            sage: Q.division_order(d, e)
            False
            sage: Q.division_order(e, d)
            False

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        return (
            all(d[i] <= e[i] for i in self.sources())
            and all(d[i] >= e[i] for i in self.sinks())
            and all(
                d[i] == e[i]
                for i in self.vertices()
                if i not in self.sources() and i not in self.sinks()
            )
        )

    """
    General subdimension vectors and general Hom and Ext
    """

    def is_general_subdimension_vector(self, e, d) -> bool:
        r"""Checks if e is a general subdimension vector of d.

        INPUT:

        - ``e``: dimension vector for the subrepresentation

        - ``d``: dimension vector for the ambient representation

        OUTPUT: whether e is a general subdimension vector of d

        A dimension vector `e` is a general subdimension vector of `d`
        if a general representation of dimension vector `d` possesses
        a subrepresentation of dimension vector `e`.
        By MR1162487_ `e` is a general subdimension vector of `d` if and only if `e` is
        a subdimension vector of `d` and :math:`\langle f,d-e\rangle` is non-negative
        for all general subdimension vectors `f` of `e`.

        .. _MR1162487: https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487

        EXAMPLES:

        Some examples on loop quivers::

            sage: from quivercombinatorics import *
            sage: Q = LoopQuiver(1)
            sage: ds = [vector([i]) for i in range(3)]
            sage: for (e, d) in cartesian_product([ds, ds]):
            ....:     if not Q.is_subdimension_vector(e, d): continue
            ....:     print("{} is general subdimension vector of {}: {}".format(
            ....:         e, d, Q.is_general_subdimension_vector(e,d))
            ....:     )
            (0) is general subdimension vector of (0): True
            (0) is general subdimension vector of (1): True
            (0) is general subdimension vector of (2): True
            (1) is general subdimension vector of (1): True
            (1) is general subdimension vector of (2): True
            (2) is general subdimension vector of (2): True
            sage: Q = LoopQuiver(2)
            sage: for (e, d) in cartesian_product([ds]*2):
            ....:     if not Q.is_subdimension_vector(e, d): continue
            ....:     print("{} is general subdimension vector of {}: {}".format(
            ....:         e, d, Q.is_general_subdimension_vector(e,d))
            ....:     )
            (0) is general subdimension vector of (0): True
            (0) is general subdimension vector of (1): True
            (0) is general subdimension vector of (2): True
            (1) is general subdimension vector of (1): True
            (1) is general subdimension vector of (2): False
            (2) is general subdimension vector of (2): True

        Some examples on generalized Kronecker quivers::

            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: ds = Tuples(range(3), 2)
            sage: for (e, d) in cartesian_product([ds]*2):
            ....:     if not Q.is_subdimension_vector(e, d): continue
            ....:     print("{} is general subdimension vector of {}: {}".format(
            ....:         e, d, Q.is_general_subdimension_vector(e,d))
            ....:     )
            (0, 0) is general subdimension vector of (0, 0): True
            (0, 0) is general subdimension vector of (1, 0): True
            (0, 0) is general subdimension vector of (2, 0): True
            (0, 0) is general subdimension vector of (0, 1): True
            (0, 0) is general subdimension vector of (1, 1): True
            (0, 0) is general subdimension vector of (2, 1): True
            (0, 0) is general subdimension vector of (0, 2): True
            (0, 0) is general subdimension vector of (1, 2): True
            (0, 0) is general subdimension vector of (2, 2): True
            (1, 0) is general subdimension vector of (1, 0): True
            (1, 0) is general subdimension vector of (2, 0): True
            (1, 0) is general subdimension vector of (1, 1): False
            (1, 0) is general subdimension vector of (2, 1): True
            (1, 0) is general subdimension vector of (1, 2): False
            (1, 0) is general subdimension vector of (2, 2): False
            (2, 0) is general subdimension vector of (2, 0): True
            (2, 0) is general subdimension vector of (2, 1): False
            (2, 0) is general subdimension vector of (2, 2): False
            (0, 1) is general subdimension vector of (0, 1): True
            (0, 1) is general subdimension vector of (1, 1): True
            (0, 1) is general subdimension vector of (2, 1): True
            (0, 1) is general subdimension vector of (0, 2): True
            (0, 1) is general subdimension vector of (1, 2): True
            (0, 1) is general subdimension vector of (2, 2): True
            (1, 1) is general subdimension vector of (1, 1): True
            (1, 1) is general subdimension vector of (2, 1): True
            (1, 1) is general subdimension vector of (1, 2): True
            (1, 1) is general subdimension vector of (2, 2): True
            (2, 1) is general subdimension vector of (2, 1): True
            (2, 1) is general subdimension vector of (2, 2): False
            (0, 2) is general subdimension vector of (0, 2): True
            (0, 2) is general subdimension vector of (1, 2): True
            (0, 2) is general subdimension vector of (2, 2): True
            (1, 2) is general subdimension vector of (1, 2): True
            (1, 2) is general subdimension vector of (2, 2): True
            (2, 2) is general subdimension vector of (2, 2): True
            sage: Q = GeneralizedKroneckerQuiver(2)
            sage: for (e, d) in cartesian_product([ds]*2):
            ....:     if not Q.is_subdimension_vector(e, d): continue
            ....:     print("{} is general subdimension vector of {}: {}".format(
            ....:         e, d, Q.is_general_subdimension_vector(e,d))
            ....:     )
            (0, 0) is general subdimension vector of (0, 0): True
            (0, 0) is general subdimension vector of (1, 0): True
            (0, 0) is general subdimension vector of (2, 0): True
            (0, 0) is general subdimension vector of (0, 1): True
            (0, 0) is general subdimension vector of (1, 1): True
            (0, 0) is general subdimension vector of (2, 1): True
            (0, 0) is general subdimension vector of (0, 2): True
            (0, 0) is general subdimension vector of (1, 2): True
            (0, 0) is general subdimension vector of (2, 2): True
            (1, 0) is general subdimension vector of (1, 0): True
            (1, 0) is general subdimension vector of (2, 0): True
            (1, 0) is general subdimension vector of (1, 1): False
            (1, 0) is general subdimension vector of (2, 1): False
            (1, 0) is general subdimension vector of (1, 2): False
            (1, 0) is general subdimension vector of (2, 2): False
            (2, 0) is general subdimension vector of (2, 0): True
            (2, 0) is general subdimension vector of (2, 1): False
            (2, 0) is general subdimension vector of (2, 2): False
            (0, 1) is general subdimension vector of (0, 1): True
            (0, 1) is general subdimension vector of (1, 1): True
            (0, 1) is general subdimension vector of (2, 1): True
            (0, 1) is general subdimension vector of (0, 2): True
            (0, 1) is general subdimension vector of (1, 2): True
            (0, 1) is general subdimension vector of (2, 2): True
            (1, 1) is general subdimension vector of (1, 1): True
            (1, 1) is general subdimension vector of (2, 1): True
            (1, 1) is general subdimension vector of (1, 2): False
            (1, 1) is general subdimension vector of (2, 2): True
            (2, 1) is general subdimension vector of (2, 1): True
            (2, 1) is general subdimension vector of (2, 2): False
            (0, 2) is general subdimension vector of (0, 2): True
            (0, 2) is general subdimension vector of (1, 2): True
            (0, 2) is general subdimension vector of (2, 2): True
            (1, 2) is general subdimension vector of (1, 2): True
            (1, 2) is general subdimension vector of (2, 2): True
            (2, 2) is general subdimension vector of (2, 2): True

        """
        return self.__is_general_subdimension_vector(e, d)

    @cached_method(
        key=lambda self, e, d: (self._coerce_vector(e), self._coerce_vector(d))
    )
    def __is_general_subdimension_vector(self, e, d) -> bool:
        r"""
        The cacheable implementation of :meth:`Quiver.is_general_subdimension_vector`.

        EXAMPLES:

        General subdimension vectors for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q, d, theta = GeneralizedKroneckerQuiver(3), (2, 3), (3, -2)
            sage: for e in Q.all_subdimension_vectors(d):
            ....:     print("{} is general subdimension vector of {}: {}".format(
            ....:         e, d, Q._Quiver__is_general_subdimension_vector(e, d))
            ....:     )
            (0, 0) is general subdimension vector of (2, 3): True
            (0, 1) is general subdimension vector of (2, 3): True
            (0, 2) is general subdimension vector of (2, 3): True
            (0, 3) is general subdimension vector of (2, 3): True
            (1, 0) is general subdimension vector of (2, 3): False
            (1, 1) is general subdimension vector of (2, 3): False
            (1, 2) is general subdimension vector of (2, 3): True
            (1, 3) is general subdimension vector of (2, 3): True
            (2, 0) is general subdimension vector of (2, 3): False
            (2, 1) is general subdimension vector of (2, 3): False
            (2, 2) is general subdimension vector of (2, 3): False
            (2, 3) is general subdimension vector of (2, 3): True
        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        if e == d or all(ei == 0 for ei in e):
            return True

        if not self.is_subdimension_vector(e, d):
            return False

        ds = filter(
            lambda eprime: self.euler_form(eprime, d - e) < 0,
            self.all_subdimension_vectors(e),
        )

        return not any(self.is_general_subdimension_vector(eprime, e) for eprime in ds)

    def all_general_subdimension_vectors(self, d, proper=False, nonzero=False):
        r"""Returns the list of all general subdimension vectors of ``d``.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: list of vectors

        EXAMPLES:

        Some n-Kronecker quivers::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: d = (3, 3)
            sage: Q.all_general_subdimension_vectors(d)
            [(0, 0),
             (0, 1),
             (0, 2),
             (0, 3),
             (1, 1),
             (1, 2),
             (1, 3),
             (2, 2),
             (2, 3),
             (3, 3)]
            sage: Q = GeneralizedKroneckerQuiver(2)
            sage: Q.all_general_subdimension_vectors(d)
            [(0, 0),
             (0, 1),
             (0, 2),
             (0, 3),
             (1, 1),
             (1, 2),
             (1, 3),
             (2, 2),
             (2, 3),
             (3, 3)]
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.all_general_subdimension_vectors(d)
            [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (3, 3)]
            sage: Q.all_general_subdimension_vectors(d, nonzero=True)
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (3, 3)]
            sage: Q.all_general_subdimension_vectors(d, proper=True)
            [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

        """
        d = self._coerce_dimension_vector(d)

        return list(
            filter(
                lambda e: self.is_general_subdimension_vector(e, d),
                self.all_subdimension_vectors(
                    d, proper=proper, nonzero=nonzero, forget_labels=True
                ),
            )
        )

    def general_ext(self, d, e):
        r"""
        Computes :math:`\operatorname{ext}(d, e)`.

        INPUT:

        - ``d``: dimension vector

        - ``e``: dimension vector

        OUTPUT: dimension of the general ext

        According to Theorem 5.4 in Schofield's 'General representations of quivers',
        we have

        .. MATH::

            \operatorname{ext}(a,b) =
            \operatorname{max}\{-\langle c,b\rangle\},

        where :math:`c` runs over the general subdimension vectors of :math:`a`.

        EXAMPLES:

        General ext on the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: ds = [Q.simple_root(0), Q.simple_root(1), Q.thin_dimension_vector()]
            sage: for (d, e) in cartesian_product([ds]*2):
            ....:     print("ext({}, {}) = {}".format(d, e, Q.general_ext(d, e)))
            ext((1, 0), (1, 0)) = 0
            ext((1, 0), (0, 1)) = 3
            ext((1, 0), (1, 1)) = 2
            ext((0, 1), (1, 0)) = 0
            ext((0, 1), (0, 1)) = 0
            ext((0, 1), (1, 1)) = 0
            ext((1, 1), (1, 0)) = 0
            ext((1, 1), (0, 1)) = 2
            ext((1, 1), (1, 1)) = 1

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        return max(
            -self.euler_form(f, e) for f in self.all_general_subdimension_vectors(d)
        )

    def general_hom(self, d, e):
        r"""
        Computes :math:`\operatorname{hom}(d, e)`.

        INPUT:

        - ``d``: dimension vector

        - ``e``: dimension vector

        OUTPUT: dimension of the general hom

        There is a non-empty open subset `U` of :math:`R(Q,d) \times R(Q,e)` such that

        .. MATH::

            \operatorname{dim}(\operatorname{Ext}(M,N)) = \operatorname{ext}(d,e),

        i.e., :math:`\operatorname{dim}(\operatorname{Ext}(M,N))` is minimal for all
        `(M,N)` in `U`.

        Therefore, :math:`\operatorname{dim}(\operatorname{Hom}(M,N)) =
        \langle a,b\rangle + \operatorname{dim}(\operatorname{Ext}(M,N))`
        is minimal, and
        :math:`\operatorname{hom}(a,b) = \langle a,b\rangle + \operatorname{ext}(a,b)`.

        EXAMPLES:

        General hom on the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: ds = [Q.simple_root(0), Q.simple_root(1), Q.thin_dimension_vector()]
            sage: for (d, e) in cartesian_product([ds]*2):
            ....:     print("hom({}, {}) = {}".format(d, e, Q.general_hom(d, e)))
            hom((1, 0), (1, 0)) = 1
            hom((1, 0), (0, 1)) = 0
            hom((1, 0), (1, 1)) = 0
            hom((0, 1), (1, 0)) = 0
            hom((0, 1), (0, 1)) = 1
            hom((0, 1), (1, 1)) = 1
            hom((1, 1), (1, 0)) = 1
            hom((1, 1), (0, 1)) = 0
            hom((1, 1), (1, 1)) = 0

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        return self.euler_form(d, e) + self.general_ext(d, e)

    """
    Harder--Narasimhan types
    """

    @cached_method
    def _all_harder_narasimhan_types(self, d, theta, denom=sum, sorted=True):
        r"""Returns the list of all Harder--Narasimhan types of d.

        INPUT:

        - ``d`` -- dimension vector

        - ``theta` -- stability parameter

        - ``denom`` -- the denominator function (default: sum)

        OUTPUT: tuple of Harder--Narasimhan types

        EXAMPLES:

        The Harder--Narasimhan types for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = (2, 3)
            sage: theta = (3, -2)
            sage: Q._all_harder_narasimhan_types(d, theta)
            (((2, 3),),
             ((1, 1), (1, 2)),
             ((2, 2), (0, 1)),
             ((2, 1), (0, 2)),
             ((1, 0), (1, 3)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 1), (0, 2)),
             ((2, 0), (0, 3)))

        .. NOTE ::

        This is a method of Quiver so that its results can be more efficiently cached.
        See :meth:`QuiverModuli.all_harder_narasimhan_types()` for the better location
        to call it from.

        """
        d = self._coerce_dimension_vector(d)
        theta = self._coerce_vector(theta)

        ds = self.all_subdimension_vectors(d, proper=True, nonzero=True)
        ds = filter(
            lambda e: self.slope(e, theta, denom=denom)
            > self.slope(d, theta, denom=denom),
            ds,
        )
        ds = filter(
            lambda e: self.has_semistable_representation(e, theta, denom=denom),
            ds,
        )
        ds = list(ds)

        if sorted:
            ds.sort(key=(lambda e: self.slope(e, theta, denom=denom)))

        all_types = []
        for e in ds:
            for estar in filter(
                lambda fstar: self.slope(e, theta, denom=denom)
                > self.slope(fstar[0], theta, denom=denom),
                self._all_harder_narasimhan_types(
                    d - e, theta, denom=denom, sorted=sorted
                ),
            ):
                all_types.append((e,) + estar)

        if self.has_semistable_representation(d, theta, denom=denom):
            all_types.insert(0, (d,))

        # because it is a cached method we need to return a tuple, not a list
        # as the result of a cached method must be immutable
        return tuple(all_types)

    """
    (Semi-)stability
    """

    def canonical_stability_parameter(self, d):
        r"""
        Returns the canonical stability parameter for ``d``

        INPUT:

        - ``d``: dimension vector

        OUTPUT: canonical stability parameter

        The canonical stability parameter is given by
        :math:`\langle d,-\rangle - \langle -,d\rangle`.

        EXAMPLES:

        Our usual example of the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.canonical_stability_parameter((2, 3))
            (9, -6)

        For the 5-subspace quiver::

            sage: Q = SubspaceQuiver(5)
            sage: Q.canonical_stability_parameter((1, 1, 1, 1, 1, 2))
            (2, 2, 2, 2, 2, -5)

        It takes vertex labels (if present) into account::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q.canonical_stability_parameter((2, 3))
            {'bar': -6, 'foo': 9}

        EXAMPLES:

        Canonical stability parameter for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q, d = GeneralizedKroneckerQuiver(3), (2, 3)
            sage: Q.canonical_stability_parameter(d)
            (9, -6)

        This method also works with vertex labels::

            sage: from quivercombinatorics import *
            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: d = {"foo": 2, "bar": 3}
            sage: Q.canonical_stability_parameter(d)
            {'bar': -6, 'foo': 9}
        """
        d = self._coerce_dimension_vector(d)
        theta = vector(d) * (-self.euler_matrix().transpose() + self.euler_matrix())

        if self.__has_vertex_labels():
            return dict(zip(self.vertices(), theta))
        return theta

    def has_semistable_representation(self, d, theta=None, denom=sum):
        r"""Checks if there is a ``theta``-semistable of dimension vector ``d``

        INPUT:

        - ``d``: dimension vector

        - ``theta`` (default: canonical stability parameter): stability parameter

        OUTPUT: whether there is a ``theta``-semistable of dimension vector ``d``

        By MR1162487_ a dimension vector `d` admits a :math:`\theta`-semi-stable
        representation if and only if :math:`\mu_{\theta}(e) \leq \mu_{\theta}(d)` for
        all general subdimension vectors `e` of `d`.

        .. _MR1162487: https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487

        EXAMPLES:

        Semistables for the :math:`\mathrm{A}_2` quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: Q.has_semistable_representation((1, 1), (1, -1))
            True
            sage: Q.has_semistable_representation((2, 2), (1, -1))
            True
            sage: Q.has_semistable_representation((1, 2), (1, -1))
            False
            sage: Q.has_semistable_representation((0, 0), (1, -1))
            True

        Semistables for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.has_semistable_representation((2, 3))
            True
            sage: Q.has_semistable_representation((1, 4), (-3, 2))
            False

        """
        if theta is None:
            theta = self.canonical_stability_parameter(d)

        d = self._coerce_dimension_vector(d)
        theta = self._coerce_vector(theta)

        return all(
            self.slope(e, theta, denom=denom) <= self.slope(d, theta, denom=denom)
            for e in self.all_general_subdimension_vectors(d, nonzero=True)
        )

    def has_stable_representation(self, d, theta=None, denom=sum):
        r"""
        Checks if there is a ``theta``-stable representation of ``d``

        INPUT:

        - ``d``: dimension vector

        - ``theta`` (default: canonical stability parameter): stability parameter

        OUTPUT: whether there is a ``theta``-stable of dimension vector ``d``

        By MR1162487_ `d` admits a theta-stable representation if and only if
        :math:`\mu_{\theta}(e) < \mu_{\theta}(d)` for all proper general subdimension
        vectors :math:`e` of :math:`d`.

        .. _MR1162487: https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487

        EXAMPLES:

        Stables for the :math:`\mathrm{A}_2` quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: theta = (1, -1)
            sage: Q.has_stable_representation((1, 1), theta)
            True
            sage: Q.has_stable_representation((2, 2), theta)
            False
            sage: Q.has_stable_representation((0, 0), theta)
            False

        Stables for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = (2, 3)
            sage: theta = Q.canonical_stability_parameter(d)
            sage: Q.has_stable_representation(d, theta)
            True
            sage: Q.has_stable_representation(d)
            True

        """
        if theta is None:
            theta = self.canonical_stability_parameter(d)

        d = self._coerce_dimension_vector(d)
        theta = self._coerce_vector(theta)

        if d == self._coerce_dimension_vector(self.zero_vector()):
            return False

        return all(
            self.slope(e, theta, denom=denom) < self.slope(d, theta, denom=denom)
            for e in self.all_general_subdimension_vectors(d, proper=True, nonzero=True)
        )

    """
    Canonical decomposition
    """

    def canonical_decomposition(self, d):
        r"""
        Computes the canonical decomposition of a dimension vector.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: canonical decomposition as list of dimension vectors

        The canonical decomposition of a dimension vector `d` is the unique
        decomposition :math:`d = e_1 + e_2 + ... + e_k` such that
        :math:`e_1, e_2, ..., e_k` are such that for all
        :math:`i \neq j, \mathrm{ext}(e_i, e_j) = \mathrm{ext}(e_j, e_i) = 0`.

        The general representation of dimension vector `d` is isomorphic to the direct
        sum of representations of dimension vectors :math:`e_1, e_2, ..., e_k`.

        EXAMPLES:

        Canonical decomposition of the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.canonical_decomposition((2, 3))
            ((2, 3),)
            sage: for d in Q.all_subdimension_vectors((5, 5)):
            ....:     print(Q.canonical_decomposition(d))
            ((0, 0),)
            ((0, 1),)
            ((0, 1), (0, 1))
            ((0, 1), (0, 1), (0, 1))
            ((0, 1), (0, 1), (0, 1), (0, 1))
            ((0, 1), (0, 1), (0, 1), (0, 1), (0, 1))
            ((1, 0),)
            ((1, 1),)
            ((1, 2),)
            ((1, 3),)
            ((0, 1), (1, 3))
            ((0, 1), (0, 1), (1, 3))
            ((1, 0), (1, 0))
            ((2, 1),)
            ((2, 2),)
            ((2, 3),)
            ((2, 4),)
            ((2, 5),)
            ((1, 0), (1, 0), (1, 0))
            ((3, 1),)
            ((3, 2),)
            ((3, 3),)
            ((3, 4),)
            ((3, 5),)
            ((1, 0), (1, 0), (1, 0), (1, 0))
            ((1, 0), (3, 1))
            ((4, 2),)
            ((4, 3),)
            ((4, 4),)
            ((4, 5),)
            ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0))
            ((1, 0), (1, 0), (3, 1))
            ((5, 2),)
            ((5, 3),)
            ((5, 4),)
            ((5, 5),)
        """
        return self.__canonical_decomposition(d)

    @cached_method(key=lambda self, d: self._coerce_vector(d))
    def __canonical_decomposition(self, d):
        r"""The cacheable implementation of :meth:`Quiver.canonical_decomposition`

        EXAMPLES:

        Canonical decomposition of `(5, 3)` for the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = GeneralizedKroneckerQuiver(2)
            sage: Q._Quiver__canonical_decomposition((5, 5))
            ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1))
        """
        d = self._coerce_dimension_vector(d)

        ds = self.all_general_subdimension_vectors(d, proper=True, nonzero=True)
        for e in ds:
            if d - e in ds:
                return self.canonical_decomposition(e) + self.canonical_decomposition(
                    d - e
                )

        # because it is a cached method we need to return a tuple, not a list
        # as the result of a cached method must be immutable
        return (d,)

    def dimension_nullcone(self, d):
        r"""
        Returns the dimension of the nullcone

        The nullcone is the set of all nilpotent representations.

        INPUT:

        - ``d`` -- dimension vector

        OUTPUT: dimension of the nullcone

        EXAMPLES:

        The usual example of the 3-Kronecker quiver::

            sage: from quivercombinatorics import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.dimension_nullcone((2, 3))
            18

        """
        d = self._coerce_dimension_vector(d)

        if self.is_acyclic():
            return d * self.adjacency_matrix() * d
        else:
            raise NotImplementedError()

    def first_hochschild_cohomology(self):
        r"""
        Compute the dimension of the first Hochschild cohomology

        This uses the formula of Happel from Proposition 1.6 in MR1035222_.
        One needs the quiver to be acyclic for this, otherwise it is not necessarily
        finite-dimensional.

        EXAMPLES:

        The first Hochschild cohomology of the `m`-th generalized Kronecker quiver
        is the dimension of :math:`\mathrm{PGL}_{m+1}`::

            sage: from quivercombinatorics import *
            sage: GeneralizedKroneckerQuiver(3).first_hochschild_cohomology()
            8

        The first Hochschild cohomology vanishes if and only if the quiver is a tree::

            sage: from quivercombinatorics import *
            sage: SubspaceQuiver(7).first_hochschild_cohomology()
            0

        .. _MR1035222: https://mathscinet.ams.org/mathscinet/relay-station?mr=1035222
        """
        assert self.is_acyclic(), "the quiver needs to be acyclic"

        return (
            1
            - self.number_of_vertices()
            + sum(
                len(self.graph().all_paths(a[0], a[1], use_multiedges=True))
                for a in self.graph().edges()
            )
        )
