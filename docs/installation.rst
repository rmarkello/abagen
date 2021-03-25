.. _installation_setup:

----------------------
Installation and setup
----------------------

.. note::

    Using the instructions for the :ref:`io_installation` is highly
    recommended! It may take a little more set-up depending on your operating
    system, but the benefits during processing are noticeable!

.. _basic_installation:

Basic installation
==================

This package requires Python 3.6+. Assuming you have the correct version of
Python installed, you can install ``abagen`` by opening a terminal and running
the following:

.. code-block:: bash

    pip install abagen

Alternatively, you can install the most up-to-date version of from GitHub:

.. code-block:: bash

   git clone https://github.com/rmarkello/abagen.git
   cd abagen
   pip install .

.. important::

   If you are going to install the version directly from GitHub make sure that
   you are using the most `up-to-date documentation
   <https://abagen.readthedocs.io/en/latest/>`_!

.. _io_installation:

IO installation
===============

The data supplied by the Allen Human Brain Atlas is quite large---on the order
of ~4GB for all six donors. Because loading these datasets into memory can be
quite time-consuming, ``abagen`` has integrated support for `parquet <https://
parquet.apache.org/>`_ and can do some on-the-fly conversions to speed things
up. However, using parquet is completely optional, and therefore support for it
is not installed when using the "vanilla" installation procedures.

If you would like to enable parquet support, you will need to install some
additional dependencies. This can be done using ``pip``:

.. code-block:: bash

    pip install abagen[io]

.. note::

    If you are using a Linux OS, you will need to install the ``libsnappy-dev``
    (Debian) or ``snappy-devel`` (Fedora) package before running the above
    installation command!

You can also install these extra packages from the GitHub source:

.. code-block:: bash

   git clone https://github.com/rmarkello/abagen.git
   cd abagen
   pip install .[io]

.. important::

    If you are going to install the version directly from GitHub make sure that
    you are using the most `up-to-date documentation
    <https://abagen.readthedocs.io/en/latest/>`_!
