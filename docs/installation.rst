.. _installation_setup:

Installation and setup
======================

.. _vanilla_installation:

Vanilla installation
--------------------

This package requires Python >= 3.5. Assuming you have the correct version of
Python installed, you can install ``abagen`` by opening a terminal and running
the following:

.. code-block:: bash

   git clone https://github.com/rmarkello/abagen.git
   cd abagen
   python setup.py install


Alternatively, you can install from PyPi with:

.. code-block:: bash

    pip install abagen

.. _full_installation:

Full installation
-----------------

The data supplied by the Allen Human Brain Atlas is quite large---on the order
of ~4GB for all six donors. Because loading these datasets into memory can be
quite time-consuming, ``abagen`` has integrated support for `parquet <https://
parquet.apache.org/>`_ and can do some on-the-fly conversions to speed things
up. However, using ``parquet`` is completely optional, and therefore support
for such is not installed when using the "vanilla" installation procedures
detailed above.

If you would like to enable parquet support, you will need to install some
additional dependencies, which can be done using ``pip``:

.. code-block:: bash

   git clone https://github.com/rmarkello/abagen.git
   cd abagen
   pip install .[io]

.. note::

    If you are using a Linux OS, you will need to install the ``libsnappy-dev``
    (Debian) or ``snappy-devel`` (Fedora) package before running the above
    installation command!

You can also perform the full installation directly from PyPi with:

.. code-block:: bash

    pip install abagen[io]
