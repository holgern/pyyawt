PyYAWT - Yet Another Wavelet Toolbox in Python
=================================================

PyYAWT is a free Open Source wavelet toolbox for Python_
programming language. 

This toolbox is aimed to mimic matlab wavelet toolbox. Most of the functions are
similiar to their counterparts in Matlab equivalents.

  .. sourcecode:: python

    >>> import pyyawt
    >>> cA, cD = pyyawt.dwt([1, 2, 3, 4], 'db1')


Main features
-------------

The main features of PyYAWT are:

  * 1D, 2D and 3D Forward and Inverse Discrete Wavelet Transform (DWT and IDWT)
  * 1D and 2D Stationary Wavelet Transform (Undecimated Wavelet Transform)
  * Continuous wavelet transform
  * Dualtree real and complex wavelet transform
  * Double precision calculations
  * Results compatible with Matlab Wavelet Toolbox (TM)

Requirements
------------

PyYAWT is a package for the Python programming language. It requires:

 - Python_ 2.7 or >=3.3
 - Numpy_ >= 1.6.2

Download
--------

The most recent *development* version can be found on GitHub at
https://github.com/holgern/pyyawt.

Latest release, including source and binary package for Windows, is available
for download from the `Python Package Index`_ or on the `Releases Page`_.

Install
-------

In order to build PyYAWT from source, a working C compiler (GCC or MSVC)
and a recent version of Cython_ is required.

 - Install PyYAWT with ``pip install pyyawt``.

 - To build and install from source, navigate to downloaded PyYAWT source
   code directory and type ``python setup.py install``.

Prebuilt Windows binaries and source code packages are also
available from `Python Package Index`_.


.. seealso::  :ref:`Development notes <dev-index>` section contains more
              information on building and installing from source code.

Documentation
-------------

Documentation with detailed examples and links to more resources is available
online at http://pyyawt.readthedocs.org.

For more usage examples see the `demo`_ directory in the source package.

State of development & Contributing
-----------------------------------

PyYAWT started in 2009 as a scilab toolbox and was maintained until 2009 by its `original developer`_.  
Authors were
* Professor Mei, Supervisor
* Roger Liu 
* Isaac Zhi 
* Jason Huang
* Du HuiQian

In 2010, maintenance was taken over in a `new repo <http://forge.scilab.org/index.php/p/swt/>`_)
by Holger Nahrstaedt. 
Daubechies wavelets coefficents DB2 - DB50 were calculated  by Bob Strunz - University of Limerick, Ireland


Finally, all the c-source files from the SWT-Toolbox are forked into this python toolbox  `pyawt <https://github.com/holgern/pyyawt/>`_)

Contributions recarding bug reports, bug fixes and new features are welcome.  


Python 3
--------

Python 3.x is fully supported from release v0.0.1 on.

Contact
-------

Use `GitHub Issues`_ to post your comments or questions.

License
-------

PyYAWT is a free Open Source software released under the GPL license.

Contents
--------
.. toctree::
   :maxdepth: 2
  
   dev/index
   resources
   
   modules


  


.. _Cython: http://cython.org/
.. _demo: https://github.com/holgern/pyyawt/tree/master/demo
.. _GitHub: https://github.com/holgern/pyyawt
.. _GitHub Issues: https://github.com/holgern/pyyawt/issues
.. _Numpy: http://www.numpy.org
.. _original developer: http://scwt.sourceforge.net/
.. _Python: http://python.org/
.. _Python Package Index: http://pypi.python.org/pypi/pyyawt/
.. _Releases Page: https://github.com/holgern/pyyawt/releases
