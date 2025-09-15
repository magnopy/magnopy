.. _user-guide_input_model-file:

*******************************
Specification of the model file
*******************************

The model file is a file that contains all the information about the spin Hamiltonian
It is intendant to describe the structure of the crystal (cell and atoms) and the
parameters (bilinear exchange, on-site anisotropy, etc.) of the spin Hamiltonian.

The model file of magnopy can be a .txt file or an .hdf5 file.
The .txt file is a plain text file with a specific format. The .hdf5 file is a binary file
with a specific structure. These two files are designed to be completely interchangeable in
the context of magnopy. Magnopy by itself offers a conversion from one format to another.

We suggest to use |myHDF5|_ web-site to browse the content of the .hdf5 file.

* One full .txt file corresponds to the group in .hdf5 file with the attribute "type"
  set to the value "SpinHamiltonian". In this guide we reference this root group as
  "spinham/"
* All dataset's, group's, and attributes's **names** are lowercase
  (i.e "units", "spinham", "cell", ...).

* Sections
  :ref:`user-guide_input_model-file_cell`,
  :ref:`user-guide_input_model-file_atoms`,
  are required. Other sections are optional.
* Section :ref:`user-guide_input_model-file_notation` is required if section
  :ref:`user-guide_input_model-file_exchange` or :ref:`user-guide_input_model-file_on-site`
  is present.



.. include:: cell.inc

.. include:: atoms.inc

.. include:: notation.inc

.. include:: exchange.inc

.. include:: on-site.inc


Full specification reference
============================

.. literalinclude:: full-hdf5.inc
  :language: text
