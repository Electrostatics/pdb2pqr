"""PDB2PQR exceptions

This module represents errors specific to PDB2PQR. Exists mainly to allow us
to more easily distinguish between code errors and input errors.

Author:  Kyle Monson
"""
import inspect


class PDB2PQRError(Exception):
    """General PDB2PQR error."""
    def __init__(self, message):
        self.message = message
        self.line = inspect.currentframe().f_back.f_lineno
        self.filename = inspect.currentframe().f_back.f_code.co_filename

    def __str__(self):
        estr = "DEBUG INFO: {errorType} {filename}: {line} Error encountered: {message}"
        return estr.format(message=self.message, errorType=self.__class__.__name__,
                           line=self.line, filename=self.filename)


class PDBInternalError(PDB2PQRError):
    """Something strange happened inside the PDB file."""
    pass


class PDBInputError(PDB2PQRError):
    """Something strange happend with the PDB input."""
    pass


class PDB2PKAError(PDB2PQRError):
    """Error within PDB2PKA."""
    pass


class PDBFileParseError(PDB2PQRError):
    """Parsing error with PDB file."""
    def __init__(self, lineno, error_str):
        self.lineno = lineno
        self.error_str = error_str

    def __str__(self):
        estr = 'PDB file parsing error line {lineno}: {error_str}'
        return estr.format(lineno=self.lineno, error_str=self.error_str)
