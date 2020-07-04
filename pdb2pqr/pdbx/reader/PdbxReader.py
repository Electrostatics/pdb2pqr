##
# File: PdbxReader.py
# Date: 2012-01-09 Jdw Adapted from PdbxParser
#
# Updates:
#
# 2012-01-09 - (jdw) Separate reader and writer classes.
#
# 2012-09-02 - (jdw) Revise tokenizer to better handle embedded quoting.
#
##
"""
PDBx/mmCIF dictionary and data file parser.

Acknowledgements:

 The tokenizer used in this module is modeled after the clever parser design
 used in the PyMMLIB package.

 PyMMLib Development Group
 Authors: Ethan Merritt: merritt@u.washington.ed & Jay Painter: jay.painter@gmail.com
 See: http://pymmlib.sourceforge.net/

"""
import re
from .PdbxContainers import DataCategory, DefinitionContainer, DataContainer


class PdbxError(Exception):
    """ Class for catch general errors
    """
    pass


class SyntaxError(Exception):
    """ Class for catching syntax errors
    """
    def __init__(self, line_number, text):
        Exception.__init__(self)
        self.line_number = line_number
        self.text = text

    def __str__(self):
        return "%%ERROR - [at line: %d] %s" % (self.line_number, self.text)


class PdbxReader:
    """ PDBx reader for data files and dictionaries.

    """
    def __init__(self, ifh):
        """ ifh - input file handle returned by open()
        """
        #
        self.__cur_line_number = 0
        self.__ifh = ifh
        self.__state_dict = {
            "data": "ST_DATA_CONTAINER",
            "loop": "ST_TABLE",
            "global": "ST_GLOBAL_CONTAINER",
            "save": "ST_DEFINITION",
            "stop": "ST_STOP"
        }

    def read(self, container_list):
        """
        Appends to the input list of definition and data containers.

        """
        self.__cur_line_number = 0
        try:
            self.__parser(self.__tokenizer(self.__ifh), container_list)
        except StopIteration:
            pass
        else:
            raise PdbxError()

    def __syntaxError(self, err_text):
        raise SyntaxError(self.__cur_line_number, err_text)

    def __getContainerName(self, in_word):
        """ Returns the name of the data_ or save_ container
        """
        return str(in_word[5:]).strip()

    def __getState(self, in_word):
        """Identifies reserved syntax elements and assigns an associated state.

           Returns: (reserved word, state)
           where -
              reserved word - is one of CIF syntax elements:
                               data_, loop_, global_, save_, stop_
              state - the parser state required to process this next section.
        """
        i = in_word.find("_")
        if i == -1:
            return None, "ST_UNKNOWN"

        try:
            r_word = in_word[:i].lower()
            return r_word, self.__state_dict[r_word]
        except:
            return None, "ST_UNKNOWN"

    def __parser(self, tokenizer, container_list):
        """ Parser for PDBx data files and dictionaries.

            Input - tokenizer() reentrant method recognizing
                    data item names (_category.attribute)
                    quoted strings (single, double and multi-line semi-colon delimited),
                    and unquoted strings.

                    container_list - list-type container for data and definition objects parsed
                                     from from the input file.

            Return:
                    container_list - is appended with data and definition objects -
        """
        # Working container - data or definition
        cur_container = None
        #
        # Working category container
        category_index = {}
        cur_category = None
        #
        cur_row = None
        state = None

        # Find the first reserved word and begin capturing data.
        #
        while True:
            cur_cat_name, cur_att_name, cur_quoted_string, cur_word = next(tokenizer)
            if cur_word is None:
                continue
            reserved_word, state = self.__getState(cur_word)
            if reserved_word is not None:
                break

        while True:
            #
            # Set the current state -
            #
            # At this point in the processing cycle we are expecting a token containing
            # either a '_category.attribute' or a reserved word.
            #
            if cur_cat_name is not None:
                state = "ST_KEY_VALUE_PAIR"
            elif cur_word is not None:
                reserved_word, state = self.__getState(cur_word)
            else:
                self.__syntaxError("Miscellaneous syntax error")
                return

            #
            # Process _category.attribute value assignments
            #
            if state == "ST_KEY_VALUE_PAIR":
                try:
                    cur_category = category_index[cur_cat_name]
                except KeyError:
                    # A new category is encountered - create a container and add a row
                    cur_category = category_index[cur_cat_name] = DataCategory(cur_cat_name)

                    try:
                        cur_container.append(cur_category)
                    except AttributeError:
                        self.__syntaxError("Category cannot be added to data_ block")
                        return

                    cur_row = []
                    cur_category.append(cur_row)
                else:
                    # Recover the existing row from the category
                    try:
                        cur_row = cur_category[0]
                    except IndexError:
                        self.__syntaxError("Internal index error accessing category data")
                        return

                # Check for duplicate attributes and add attribute to table.
                if cur_att_name in cur_category.getAttributeList():
                    self.__syntaxError("Duplicate attribute encountered in category")
                    return
                else:
                    cur_category.appendAttribute(cur_att_name)

                # Get the data for this attribute from the next token
                t_cat, t_att, cur_quoted_string, cur_word = next(tokenizer)

                if t_cat is not None or (cur_quoted_string is None and cur_word is None):
                    self.__syntaxError("Missing data for item _%s.%s" %
                                       (cur_cat_name, cur_att_name))

                if cur_word is not None:
                    #
                    # Validation check token for misplaced reserved words -
                    #
                    reserved_word, state = self.__getState(cur_word)
                    if reserved_word is not None:
                        self.__syntaxError("Unexpected reserved word: %s" % (reserved_word))

                    cur_row.append(cur_word)

                elif cur_quoted_string is not None:
                    cur_row.append(cur_quoted_string)

                else:
                    self.__syntaxError("Missing value in item-value pair")

                cur_cat_name, cur_att_name, cur_quoted_string, cur_word = next(tokenizer)
                continue

            #
            # Process a loop_ declaration and associated data -
            #
            elif state == "ST_TABLE":

                # The category name in the next cur_cat_name, cur_att_name pair
                # defines the name of the category container.
                cur_cat_name, cur_att_name, cur_quoted_string, cur_word = next(tokenizer)

                if cur_cat_name is None or cur_att_name is None:
                    self.__syntaxError("Unexpected token in loop_ declaration")
                    return

                # Check for a previous category declaration.
                if cur_cat_name in category_index:
                    self.__syntaxError("Duplicate category declaration in loop_")
                    return

                cur_category = DataCategory(cur_cat_name)

                try:
                    cur_container.append(cur_category)
                except AttributeError:
                    self.__syntaxError("loop_ declaration outside of data_ block or save_ frame")
                    return

                cur_category.appendAttribute(cur_att_name)

                # Read the rest of the loop_ declaration
                while True:
                    cur_cat_name, cur_att_name, cur_quoted_string, cur_word = next(tokenizer)

                    if cur_cat_name is None:
                        break

                    if cur_cat_name != cur_category.getName():
                        self.__syntaxError("Changed category name in loop_ declaration")
                        return

                    cur_category.appendAttribute(cur_att_name)

                # If the next token is a 'word', check it for any reserved words -
                if cur_word is not None:
                    reserved_word, state = self.__getState(cur_word)
                    if reserved_word is not None:
                        if reserved_word == "stop":
                            return
                        else:
                            self.__syntaxError(
                                "Unexpected reserved word after loop declaration: %s" %
                                (reserved_word))

                # Read the table of data for this loop_ -
                while True:
                    cur_row = []
                    cur_category.append(cur_row)

                    for t_att in cur_category.getAttributeList():
                        if cur_word is not None:
                            cur_row.append(cur_word)
                        elif cur_quoted_string is not None:
                            cur_row.append(cur_quoted_string)

                        cur_cat_name, cur_att_name, cur_quoted_string, cur_word = next(tokenizer)

                    # loop_ data processing ends if -

                    # A new _category.attribute is encountered
                    if cur_cat_name is not None:
                        break

                    # A reserved word is encountered
                    if cur_word is not None:
                        reserved_word, state = self.__getState(cur_word)
                        if reserved_word is not None:
                            break

                continue

            elif state == "ST_DEFINITION":
                # Ignore trailing unnamed saveframe delimiters e.g. 'save_'
                s_name = self.__getContainerName(cur_word)
                if len(s_name) > 0:
                    cur_container = DefinitionContainer(s_name)
                    container_list.append(cur_container)
                    category_index = {}
                    cur_category = None

                cur_cat_name, cur_att_name, cur_quoted_string, cur_word = next(tokenizer)

            elif state == "ST_DATA_CONTAINER":
                #
                d_name = self.__getContainerName(cur_word)
                if len(d_name) == 0:
                    d_name = "unidentified"
                cur_container = DataContainer(d_name)
                container_list.append(cur_container)
                category_index = {}
                cur_category = None
                cur_cat_name, cur_att_name, cur_quoted_string, cur_word = next(tokenizer)

            elif state == "ST_STOP":
                return
            elif state == "ST_GLOBAL":
                cur_container = DataContainer("blank-global")
                cur_container.setGlobal()
                container_list.append(cur_container)
                category_index = {}
                cur_category = None
                cur_cat_name, cur_att_name, cur_quoted_string, cur_word = next(tokenizer)

            elif state == "ST_UNKNOWN":
                self.__syntaxError("Unrecogized syntax element: " + str(cur_word))
                return

    def __tokenizer(self, ifh):
        """ Tokenizer method for the mmCIF syntax file -

            Each return/yield from this method returns information about
            the next token in the form of a tuple with the following structure.

            (category name, attribute name, quoted strings, words w/o quotes or white space)

            Differentiated the reqular expression to the better handle embedded quotes.

        """
        #
        # Regex definition for mmCIF syntax - semi-colon delimited strings are handled
        # outside of this regex.
        mmcif_re = re.compile(
            r"(?:"
            r"(?:_(.+?)[.](\S+))" "|"  # _category.attribute
            r"(?:['](.*?)(?:[']\s|[']$))" "|"  # single quoted strings
            r"(?:[\"](.*?)(?:[\"]\s|[\"]$))" "|"  # double quoted strings
            r"(?:\s*#.*$)" "|"  # comments (dumped)
            r"(\S+)"  # unquoted words
            r")")

        file_iter = iter(ifh)

        # Tokenizer loop begins here ---
        while True:
            line = next(file_iter)
            self.__cur_line_number += 1

            # Dump comments
            if line.startswith("#"):
                continue

            # Gobble up the entire semi-colon/multi-line delimited string and
            # and stuff this into the string slot in the return tuple
            #
            if line.startswith(";"):
                ml_string = [line[1:]]
                while True:
                    line = next(file_iter)
                    self.__cur_line_number += 1
                    if line.startswith(";"):
                        break
                    ml_string.append(line)

                # remove trailing new-line that is part of the \n; delimiter
                ml_string[-1] = ml_string[-1].rstrip()
                #
                yield (None, None, "".join(ml_string), None)
                #
                # Need to process the remainder of the current line -
                line = line[1:]
                # continue

            # Apply regex to the current line consolidate the single/double
            # quoted within the quoted string category
            for itor in mmcif_re.finditer(line):
                tgroups = itor.groups()
                if tgroups != (None, None, None, None, None):
                    if tgroups[2] is not None:
                        qstr = tgroups[2]
                    elif tgroups[3] is not None:
                        qstr = tgroups[3]
                    else:
                        qstr = None
                    groups = (tgroups[0], tgroups[1], qstr, tgroups[4])
                    yield groups

    def __tokenizerOrg(self, ifh):
        """ Tokenizer method for the mmCIF syntax file -

            Each return/yield from this method returns information about
            the next token in the form of a tuple with the following structure.

            (category name, attribute name, quoted strings, words w/o quotes or white space)

        """
        #
        # Regex definition for mmCIF syntax - semi-colon delimited strings are handled
        # outside of this regex.
        mmcif_re = re.compile(
            r"(?:"
            r"(?:_(.+?)[.](\S+))" "|"  # _category.attribute
            r"(?:['\"](.*?)(?:['\"]\s|['\"]$))" "|"  # quoted strings
            r"(?:\s*#.*$)" "|"  # comments (dumped)
            r"(\S+)"  # unquoted words
            r")")

        file_iter = iter(ifh)

        # Tokenizer loop begins here ---
        while True:
            line = next(file_iter)
            self.__cur_line_number += 1

            # Dump comments
            if line.startswith("#"):
                continue

            # Gobble up the entire semi-colon/multi-line delimited string and
            # and stuff this into the string slot in the return tuple
            #
            if line.startswith(";"):
                ml_string = [line[1:]]
                while True:
                    line = next(file_iter)
                    self.__cur_line_number += 1
                    if line.startswith(";"):
                        break
                    ml_string.append(line)

                # remove trailing new-line that is part of the \n; delimiter
                ml_string[-1] = ml_string[-1].rstrip()
                #
                yield (None, None, "".join(ml_string), None)
                #
                # Need to process the remainder of the current line -
                line = line[1:]
                # continue

            # Apply regex to the current line
            for itor in mmcif_re.finditer(line):
                groups = itor.groups()
                if groups != (None, None, None, None):
                    yield groups
