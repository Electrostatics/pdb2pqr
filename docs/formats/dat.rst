PDB2PQR DAT files
=================

A PDB2PQR ``.DAT`` file has the following whitespace-delimited columns:

1. Name of the residue (:py:class:`str`, e.g., RU, ARG, WAT, etc.)
2. Name of the atom (:py:class:`str`, e.g., HB3, CA, HA, etc.)
3. Charge of the atom (:py:class:`float`, in electrons)
4. Radius of the atom (:py:class:`float`, in Å)
5. (optional) Atom type (:py:class:`str`, e.g., HT, etc.)

Lines that start with ``#`` are treated as comments.
