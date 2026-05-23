***********
QuiverCombinatoricsTools
***********

`QuiverCombinatoricsTools` is a SageMath package that adds combinatorial functions to `QuiverTools` to calculate symplectic leaves of quiver varieties, available here `https://github.com/emanuel-roth/QuiverCombinatoricsTools <https://github.com/emanuel-roth/QuiverCombinatoricsTools>`_. It adds on to the `QuiverTools` package written by Pieter Belmans, Hans Franzen and Gianni Petrella, as seen here `https://sage.quiver.tools/ <https://sage.quiver.tools/>`_ and here `https://github.com/QuiverTools/QuiverTools <https://github.com/QuiverTools/QuiverTools>`_, so consult their documentation when needed.

To install it, make sure you have both QuiverTools and QuiverCombinatoricsTools

.. code-block::

   sage --pip install git+https://github.com/QuiverTools/QuiverTools.git
   sage --pip install git+https://github.com/emanuel-roth/QuiverCombinatoricsTools.git

and then you can simply run

.. code-block:: python

   from quiver import *
   from quivercombinatorics import *

to get started.

**Authors**

* Tudor-Ioan Caba (University of Edinburgh, AGQ)
* Mia Lam (University of Edinburgh, AGQ)
* Emanuel Roth (University of Edinburgh, AGQ)

We were supervised by Gwyn Bellamy (University of Glasgow), as part of an AGQ computing project.

Generating quivers
=======
.. autofunction:: quivercombinatorics.quiver_from_cartan_matrix
.. autofunction:: quivercombinatorics.random_quiver


The construction of :math:`\Sigma_{\lambda}`
=======
.. autofunction:: quivercombinatorics.N_set
.. autoclass:: quivercombinatorics.Quiver
    :members:
    :special-members:
    :member-order: bysource

Constructing quivers