"""PDB2PQR exceptions

This module represents errors specific to PDB2PQR. Exists mainly to allow us
to more easily distinguish between code errors and input errors. 

Author:  Kyle Monson
"""
import inspect


class PDB2PQRError(Exception):    
    def __init__(self, message):
        self.message = message
        self.line = inspect.currentframe().f_back.f_lineno 
        self.filename = inspect.currentframe().f_back.f_code.co_filename
        
    def __str__(self):
        return "DEBUG INFO: {errorType} {filename}: {line} Error encountered: {message}".format(message=self.message,
            errorType=self.__class__.__name__, line=self.line,
            filename=self.filename)

class PDBInternalError(PDB2PQRError):
    pass

class PDBInputError(PDB2PQRError):
    pass

class PDB2PKAError(PDB2PQRError):
    pass

class PDBFileParseError(PDB2PQRError):
    def __init__(self, lineno, errorStr):
        self.lineno = lineno
        self.errorStr = errorStr 
        
    def __str__(self):
        return 'PDB file parsing error line {lineno}: {errorStr}'.format(lineno=self.lineno, 
                                                                         errorStr=self.errorStr) 