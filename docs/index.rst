.. Asynch documentation master file, created by
   sphinx-quickstart on Mon Nov 28 14:41:46 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. meta::
   :github_url: https://github.com/Iowa-Flood-Center/asynch

Welcome to Asynch's documentation!
==================================

A numerical library for solving differential equations with a tree structure. Emphasis is given to hillslope-link river basin models.

Introduction
------------

ASYNCH is a software package for solving large systems of ordinary differential equations with a tree structure. At its heart, ASYNCH is a library that
can be mixed with other software to allow users greater fexibility. It is written entirely in the C programming language. ASYNCH currently supports
an interface for programs written in C/C++ or Python.

Although ASYNCH is a library, it does come prepackaged with simple programs for performing basic simulations. These programs can be easily
modified to increase their fexibility. ASYNCH is designed for a distributed memory computer. The Message Passing Interface (MPI) is used for inter-process communication.
ASYNCH also contains routines for downloading and uploading data from and to PostgreSQL databases.

ASYNCH was designed with hydrological applications in mind. Much of the naming conventions and sample models reflect this. However, there is
no reason why ASYNCH cannot be applied to other problems with ordinary differential equations having tree structures.

The novelty of this solver is that it uses an asynchronous integrator to solve the differential equations.
Further, the communication between processes occurs asynchronously. Details about how this works can be found in
*Small, et. al. An Asynchronous Solver for Systems of ODEs Linked by a Directed Tree Structure, Advances in Water Resources, 53, March 2013, 23-32*.


Folder structure
================

::

    +---docs
    +---examples
    +---forecasters
    +---forecasters_js
    +---ide
    +---m4
    +---src
    \---tests

The main documentation for the site is organized into a couple sections:

* :ref:`user-docs`
* :ref:`feature-docs`
* :ref:`about-docs`

Information for developers is also available:

* :ref:`dev-docs`
* :ref:`design-docs`
* :ref:`ops-docs`

.. _user-docs:

.. toctree::
  :maxdepth: 2
  :caption: User Documentation

  terminology
  getting_started
  input_output
  builtin_options
  builtin_models

.. _dev-docs:

.. toctree::
  :maxdepth: 2
  :caption: Developer Documentation

  install
  changelog
  tests
  custom_models
  c_api
  python_api

.. _about-docs:

.. toctree::
   :maxdepth: 2
   :caption: About Asynch

   contribute
   team
