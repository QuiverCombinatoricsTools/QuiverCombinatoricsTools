==================
QuiverCombinatoricsTools
==================

`QuiverCombinatoricsTools` is a SageMath package that adds combinatorial functions to `QuiverTools` to calculate symplectic leaves of quiver varieties, available here `https://github.com/QuiverCombinatoricsTools/QuiverCombinatoricsTools <https://github.com/QuiverCombinatoricsTools/QuiverCombinatoricsTools>`_. It adds on to the `QuiverTools` package written by Pieter Belmans, Hans Franzen and Gianni Petrella, as seen here `https://sage.quiver.tools/ <https://sage.quiver.tools/>`_ and here `https://github.com/QuiverTools/QuiverTools <https://github.com/QuiverTools/QuiverTools>`_, so consult their documentation when needed.

To install it, make sure you have both `QuiverTools` and `QuiverCombinatoricsTools`

.. code-block::

   sage --pip install git+https://github.com/QuiverTools/QuiverTools.git
   sage --pip install git+https://github.com/QuiverCombinatoricsTools/QuiverCombinatoricsTools.git

and then you can simply run

.. code-block:: python

   from quiver import *
   from quivercombinatorics import *

to get started.

**Authors**

* Tudor-Ioan Caba (University of Edinburgh)
* Mia Lam (University of Edinburgh)
* Emanuel Roth (University of Edinburgh)

We were supervised by Gwyn Bellamy (University of Glasgow), as part of an `AGQ <https://www.agq-cdt.org/>`_ computing project.

Generating quivers
==================

Here are some functions that help generate quivers to test examples.

.. autofunction:: quivercombinatorics.quiver_from_cartan_matrix
.. autofunction:: quivercombinatorics.random_quiver

We notate quivers by :math:`Q`, with vertices in :math:`Q_0` and edges in :math:`Q_1`.

Constructing :math:`\Sigma_{\lambda}`
=====================================

We follow William Crawley-Boevey in `this paper <https://link.springer.com/article/10.1023/A:1017558904030>`__. Given a quiver :math:`Q` with :math:`n` vertices, and :math:`{\lambda}\in\mathbb{Z}^n`, :math:`\Sigma_{\lambda}` is the set of :math:`\alpha\in\mathbb{N}^n` such that :math:`\alpha` is a positive root of :math:`Q`, :math:`\alpha\cdot\lambda = 0`, and

.. math::

   p(\alpha)>\sum_{t=1}^r p(\beta^{(t)})

for any decomposition :math:`\alpha=\beta^{(1)}+\dots+\beta^{(r)}` with :math:`r\geq 2` and :math:`\beta^{(t)}` a positive root of :math:`Q` with :math:`\lambda\cdot\beta^{(t)}=0` for all `t`. 

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.p_function

We write :math:`R_\lambda^+` to denote the set of positive roots :math:`\alpha` with :math:`\alpha\cdot\lambda=0` and :math:`\mathbb{N}R_\lambda^+` for the set of sums of elements of :math:`R_\lambda^+`. Using `Theorem 5.6 of Crawley-Boevey's paper <https://link.springer.com/article/10.1023/A:1017558904030>`_.

.. admonition:: Theorem (Crawley-Boevey, 2001)

   For :math:`\alpha\in\mathbb{N}^n`, then :math:`\alpha\in\Sigma_\lambda` if and only if :math:`0\neq\alpha\in\mathbb{N}R_\lambda^+` and :math:`(\beta,\alpha-\beta)\leq -2` whenever :math:`\beta\in\mathbb{N}^n` and :math:`0<\alpha<\beta`.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.R_lambda_plus

We also need the helper function ``N_set``.

.. autofunction:: quivercombinatorics.N_set

With these functions, we can define :math:`\Sigma_{\lambda}`.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.sigma_lambda

CB-decompositions
==================

A *representation type* of :math:`x` of a quiver :math:`Q` is a sequence :math:`\tau=(\beta^{(1)},n_1;\dots;\beta^{(k)},n_k)` with :math:`\beta^{(i)}\in\Sigma_\lambda` and :math:`n_i\in\mathbb{N}` such that

.. math::

   x=\sum_{i=1}^k n_i\beta^{(i)}

We allow :math:`\beta^{(i)}` to occur multiple times if it is an imaginary root, i.e. :math:`p(\beta^{(i)})>0`. The *CB-decomposition* of :math:`x` is defined to be the unique representation type of :math:`x` for which

.. math::

   p(\tau):=\sum_{i=1}^k p(\beta^{(i)})

is maximal. In order to construct CB-decompositions (called canonical decompositions in `this paper <https://link.springer.com/article/10.1023/A:1017558904030>`__, but are *not* the same as ``canonical_decomposition`` from `QuiverTools`). We first need to find all representation types of :math:`x`, up to a bound :math:`v`, for which we need the following helper functions.

.. autofunction:: quivercombinatorics.vector_decomposition

.. autofunction:: quivercombinatorics.small_decomposition

With these functions, we can determine all representation types.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.all_representation_types

Representation types of :math:`x` correspond to a symplectic leaf of the quiver variety :math:`M_0(Q,x)`, and the symplectic leaf dimension is :math:`2p` applied to the representation type, where :math:`p` is the ``p_function``.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.symplectic_leaf_dimension

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.CB_decomposition

Knowing the CB-decomposition, it is then easy to calculate the dimension of the quiver variety :math:`M_0(Q,x)`, as it gives the largest (open) symplectic leaf inside :math:`M_0(Q,x)`.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.quiver_variety_dimension

It's often important to know the codimension 2 leaves of the quiver variety :math:`M_0(Q,x)`, so we code this here.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.codimension_two_leaves

:math:`\mathrm{ext}`-quivers
============================

Let :math:`\tau=(\beta^{(1)},n_1;\dots;\beta^{(k)},n_k)` be a representation type of :math:`Q`. The :math:`\mathrm{ext}`\ *-quiver* :math:`\tilde{Q}` associated to the symplectic leaf is the quiver with :math:`k` vertices where the number of edges between vertex :math:`i` and vertex :math:`j` is :math:`-(\beta^{(i)},\beta^{(j)})` and the number of loops at vertex :math:`i` is :math:`p(\beta^{(i)})`. The corresponding dimension vector of :math:`\tilde{Q}` is :math:`\mathbf{n}=(n_1,\dots,n_k)`, and we take :math:`\tilde{\lambda}=0`.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.ext_quiver

.. autofunction:: quivercombinatorics.ext_dimension_vector

Classification of minimal degenerations of symplectic leaves
============================================================

For a quiver :math:`Q`, a positive root :math:`\alpha` is *minimal imaginary* if it is imaginary and given any positive root :math:`\beta`, :math:`\alpha>\beta` implies :math:`\beta` is real.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.is_minimal_imaginary_root

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.all_minimal_imaginary_positive_roots

Given a quiver :math:`Q`, a *subquiver* :math:`T = (T_0, T_1)` is specified by a subset :math:`T_0\subset Q_0` of the vertices. Then, for :math:`i, j \in T_0 \subset Q_0` the number of edges between vertex :math:`i` and :math:`j` in :math:`T` equals the number of edges between vertex :math:`i` and :math:`j` in :math:`Q`. In particular, if :math:`i\in T_0` then the number of loops at :math:`i` in :math:`T` is the same as the number of loops at :math:`i` in :math:`Q`. The support of a dimension vector :math:`v` is the subset :math:`T_0 \subset Q_0` of all vertices :math:`i` such that :math:`v_i\neq 0`. We think of the support of :math:`v` as a subquiver of :math:`Q`.

We say that the subquiver :math:`T` is of *affine ADE type* if it is one of the graphs appearing in the classification of simply-laced affine Dynkin diagrams :math:`\tilde{A}_l`, :math:`\tilde{D}_l` or :math:`\tilde{E}_6`, :math:`\tilde{E}_7`, :math:`\tilde{E}_8`. Then there is a unique imaginary root :math:`\delta`, with :math:`p(\delta)=1`, whose support is :math:`T`. The key result is the classification of minimal imaginary roots is the following (Bellamy-Schedler 2026).

.. admonition:: Theorem (Bellamy-Schedler, 2026)

      If :math:`M` is a closed leaf in :math:`M_0(Q,x)` and :math:`M\subset\bar{L}` is a minimal degeneration, then :math:`\bar{L}\cong M\times S`, where :math:`S` is isomorphic to one of the following isolated singularities:

      1. A Kleinian singularity :math:`(\mathbb{C}^2/\Gamma, 0)`.

      2. The type A minimal nilpotent orbit closure :math:`(\mathcal{O}_{\min},0)` in :math:`\mathfrak{sl}_r(\mathbb{C})`.

      3. :math:`(\mathbb{C}^{2g}/\mathbb{Z}_2,0)`.

      4. :math:`\operatorname{Spec}(\mathbb{C}[x\mid\deg x\geq 2])`.

We define:

   .. math::
      
      \tilde{\tau}(\beta,n):=(\beta,n;e_1,\alpha_1-n\beta_1;\dots;e_r,\alpha_r-n\beta_r)

In each of the above cases, the leaf :math:`L` from this theorem corresponds to representation types:

   1. :math:`\tilde{\tau}(\delta,n)` for :math:`\delta` a minimal imaginary root on an affine ADE subquiver.

   2. :math:`\tilde{\tau}((1,1),n)` for a two-vertex subquiver with :math:`t\geq 3` edges between the vertices and no loops at the vertices.

   3.
      .. math::
      
         \tau=(e_i,a;e_i,a;e_1,\alpha_1;\dots;e_{i-1},\alpha_{i-1};e_{i+1},\alpha_{i+1};\dots)
   
   where :math:`i` is a vertex with :math:`g\geq 1` loops and :math:`2a=\alpha_i`.

   4.
      .. math::
      
         \tau=(e_i,a;e_i,b;e_1,\alpha_1;\dots;e_{i-1},\alpha_{i-1};e_{i+1},\alpha_{i+1};\dots)
      
   where :math:`i` is a vertex with :math:`g\geq 1` loops and :math:`0<a\neq b<\alpha_i` with :math:`a+b=\alpha_i`.

We assign the following labels to each subminimal representation type (i.e., representation types corresponding to a minimal degeneration):

   1. :math:`A_l,D_l` or :math:`E_6,E_7,E_8` for affine Dynkin diagram of type :math:`\tilde{A}_l,\tilde{D}_l` or :math:`\tilde{E}_6,\tilde{E}_7,\tilde{E}_8`.

   2. :math:`a_{r-1}`.

   3. :math:`c_g`.

   4. :math:`m_g`.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.all_subminimal_representation_types

Plotting minimal degenerations of symplectic leaves in Hasse diagrams
============================================

To compute these Hasse diagrams, make sure you have ``dot2tex`` and ``graphviz`` installed.

.. code-block::

   sage --pip install dot2tex
   sage --pip install graphviz

We visualize the poset of symplectic leaves (equivalently of representation types) as a Hasse diagram: The vertices in this diagram are the representation types and the minimal degenerations, represented by an edge connecting two representation types, correspond to the situation where :math:`L_1\leq L_2` (equivalently, :math:`\tau_1\leq\tau_2`) but there does not exist :math:`L_3` such that :math:`L_1 < L_3 < L_2`. We wish to compute the Hasse diagram, for which there are two methods.

.. admonition:: Method 1

   Traverse the Hasse diagram from the bottom, starting with the lowest leaves. When we reach a leaf :math:`L_i`, we compute the :math:`\mathrm{ext}`-quiver associated to :math:`L_i`. On this :math:`\mathrm{ext}`-quiver, we compute subminimal representation types of :math:`\mathbf{n}`. For each such representation type :math:`\tilde{\tau}=(\gamma^{(1)},m_1;\dots;\gamma^{(r)},m_r)`, define:

   .. math::

      D(\tilde{\tau}):=\left(D(\gamma^{(1)}),m_1;\dots;D(\gamma^{(r)}),m_r\right)

   By (Bellamy-Schedler 2026), it is known that :math:`D(\tilde{\tau})` is a representation type for :math:`v`.

   In this way, we have a rule that assigns to each subminimal representation type of :math:`\tilde{Q}` a representation type :math:`D(\tilde{\tau}(\alpha))` of :math:`v`. It is known that :math:`\tau < D(\tilde{\tau}(\alpha))` is a minimal degeneration and they all occur in this way. Proceeding in this way allows one to build the Hasse diagram from the bottom up.

For Method 1, we need the following helper functions.

.. autofunction:: quivercombinatorics.D_map

.. autofunction:: quivercombinatorics.D_lifting

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.minimal_degenerations

The second method to compute this Hasse diagram is as follows.

.. admonition:: Method 2

   Consider the larger set :math:`\mathcal{D}` of all decompositions of a given dimension vector :math:`v`. Here, a decomposition is simply a way of writing :math:`v` as a sum of smaller dimension vectors with multiplicities:

   .. math::
   
      \left(n_1,\beta^{(1)};\dots;n_k,\beta^{(k)}\right),
   
   where the :math:`\beta^{(i)}` are dimension vectors satisfying :math:`\beta^{(i)}\leq v`. Following `this paper <http://matrix.uantwerpen.be/lieven.lebruyn/b2hd-LeBruyn1988b.html>`__, we say that:

   .. math::
   
      \varepsilon=(p_1,\gamma^{(1)};\dots;p_l,\gamma^{(l)})>\rho=(m_1,\alpha^{(1)};\dots;m_r,\alpha^{(r)})
   
   is a *direct successor* (:math:`\varepsilon` is the successor of :math:`\rho`) if either:

      1. :math:`l=r+1` and :math:`\varepsilon` can be ordered such that for all :math:`1\leq i\leq r-1` we have :math:`(m_i,\alpha^{(i)})=(p_i,\gamma^{(i)})` and :math:`m_r=p_r=p_{r+1},\gamma^{(r)}+\gamma^{(r+1)}=\alpha^{(r)}`.

      2. :math:`l=r-1` and :math:`\varepsilon` can be ordered such that for all :math:`1\leq i\leq r-2` we have :math:`(m_i,\alpha^{(i)})=(p_i,\gamma^{(i)})` and :math:`\alpha^{(r)}=\alpha^{(r-1)}=\gamma^{(r-1)},m_{r-1}+m_r=p_{r-1}`.

The symplectic leaves form a subset of this set of decompositions, and by restriction, they inherit a sub-poset structure. It is shown in (Bellamy-Schedler, 2026) that this is the partial ordering given by leaf closure. We need the following helper functions for Method 2.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.all_decompositions

.. autofunction:: quivercombinatorics.is_direct_successor

We are now able to obtain and plot the Hasse diagram for Methods 1 and 2. Method 1 seems to be much faster, so it is enabled by default.

.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.get_Hasse_diagram
   
.. automethod:: quivercombinatorics.quivercombinatorics.Quiver.plot_Hasse_diagram

The example above should export the following diagram:

.. image:: _static/output.svg
   :alt: Poset diagram
   :align: center
   :class: only-light
   :width: 85%

.. image:: _static/output-dark.svg
   :alt: Poset diagram
   :align: center
   :class: only-dark
   :width: 85%