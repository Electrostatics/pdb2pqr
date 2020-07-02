##
#
# File: PdbxContainers.py
# Original: 02-Feb-2009 jdw
#
# Update:
# 23-Mar-2011 jdw Added method to rename attributes in category containers.
# 05-Apr-2011 jdw Change cif writer to select double quoting as preferred
# quoting style where possible.
# 16-Jan-2012 jdw Create base class for DataCategory class
# 22-Mar-2012 jdw when append attributes to existing categories update
# existing rows with placeholder null values.
# 2-Sep-2012 jdw add option to avoid embedded quoting that might
# confuse simple parsers.
# 28-Jun-2013 jdw export remove method
# 29-Jun-2013 jdw export remove row method
##
"""

A collection of container classes supporting the PDBx/mmCIF storage model.

A base container class is defined which supports common features of
data and definition containers. PDBx data files are organized in
sections called data blocks which are mapped to data containers.
PDBx dictionaries contain definition sections and data sections
which are mapped to definition and data containes respectively.

Data in both PDBx data files and dictionaries are organized in
data categories. In the PDBx syntax individual items or data
identified by labels of the form '_categoryName.attribute_name'.
The terms category and attribute in PDBx jargon are analogous
table and column in relational data model, or class and attribute
in an object oriented data model.

The DataCategory class provides base storage container for instance
data and definition meta data.

"""
import re
import sys
import traceback

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


class CifName:
    ''' Class of utilities for CIF-style data names -
    '''
    def __init__(self):
        pass

    @staticmethod
    def categoryPart(name):
        tname = ""
        if name.startswith("_"):
            tname = name[1:]
        else:
            tname = name

        idx = tname.find(".")
        if idx == -1:
            return tname
        else:
            return tname[:idx]

    @staticmethod
    def attributePart(name):
        idx = name.find(".")
        if idx == -1:
            return None
        else:
            return name[idx + 1:]


class ContainerBase:
    ''' Container base class for data and definition objects.
    '''
    def __init__(self, name):
        # The enclosing scope of the data container (e.g. data_/save_)
        self.__name = name
        # List of category names within this container -
        self.__obj_name_list = []
        # dictionary of DataCategory objects keyed by category name.
        self.__obj_catalog = {}
        self.__type = None

    def getType(self):
        return self.__type

    def setType(self, type_name):
        self.__type = type_name

    def getName(self):
        return self.__name

    def setName(self, name):
        self.__name = name

    def exists(self, name):
        return name in self.__obj_catalog

    def getObj(self, name):
        if name in self.__obj_catalog:
            return self.__obj_catalog[name]
        else:
            return None

    def getObjNameList(self):
        return self.__obj_name_list

    def append(self, obj):
        """ Add the input object to the current object catalog. An existing object
            of the same name will be overwritten.
        """
        if obj.getName() is not None:
            if obj.getName() not in self.__obj_catalog:
                # self.__obj_name_list is keeping track of object order here --
                self.__obj_name_list.append(obj.getName())
            self.__obj_catalog[obj.getName()] = obj

    def replace(self, obj):
        """ Replace an existing object with the input object
        """
        if ((obj.getName() is not None) and (obj.getName() in self.__obj_catalog)):
            self.__obj_catalog[obj.getName()] = obj

    def printIt(self, wfh=sys.stdout, type_name="brief"):
        wfh.write("+ %s container: %30s contains %4d categories\n" %
                  (self.getType(), self.getName(), len(self.__obj_name_list)))
        for name in self.__obj_name_list:
            wfh.write("--------------------------------------------\n")
            wfh.write("Data category: %s\n" % name)
            if type_name == 'brief':
                self.__obj_catalog[name].printIt(wfh)
            else:
                self.__obj_catalog[name].dumpIt(wfh)

    def rename(self, cur_name, new_name):
        """ Change the name of an object in place -
        """
        try:
            idx = self.__obj_name_list.index(cur_name)
            self.__obj_name_list[idx] = new_name
            self.__obj_catalog[new_name] = self.__obj_catalog[cur_name]
            self.__obj_catalog[new_name].setName(new_name)
            return True
        except:
            return False

    def remove(self, cur_name):
        """ Revmove object by name. Return True on success or False otherwise.
        """
        try:
            if cur_name in self.__obj_catalog:
                del self.__obj_catalog[cur_name]
                idx = self.__obj_name_list.index(cur_name)
                del self.__obj_name_list[idx]
                return True
            else:
                return False
        except:
            pass

        return False


class DefinitionContainer(ContainerBase):
    def __init__(self, name):
        super(DefinitionContainer, self).__init__(name)
        self.setType('definition')

    def isCategory(self):
        if self.exists('category'):
            return True
        return False

    def isAttribute(self):
        if self.exists('item'):
            return True
        return False

    def printIt(self, wfh=sys.stdout, type_name="brief"):
        wfh.write("Definition container: %30s contains %4d categories\n" %
                  (self.getName(), len(self.getObjNameList())))
        if self.isCategory():
            wfh.write("Definition type: category\n")
        elif self.isAttribute():
            wfh.write("Definition type: item\n")
        else:
            wfh.write("Definition type: undefined\n")

        for name in self.getObjNameList():
            wfh.write("--------------------------------------------\n")
            wfh.write("Definition category: %s\n" % name)
            if type == 'brief':
                self.getObj(name).printIt(wfh)
            else:
                self.getObj(name).dumpId(wfh)


class DataContainer(ContainerBase):
    ''' Container class for DataCategory objects.
    '''
    def __init__(self, name):
        super(DataContainer, self).__init__(name)
        self.setType('data')
        self.__global_flag = False
        self.__row = 0

    def invokeDataBlockMethod(self, type_name, method, db_name):
        self.__row = 1
        exec(method.getInline())

    def setGlobal(self):
        self.__global_flag = True

    def getGlobal(self):
        return self.__global_flag


class DataCategoryBase:
    """ Base object definition for a data category -
    """
    def __init__(self, name, attribute_name_list=None, row_list=None):
        self.__name = name
        #
        if row_list is not None:
            self.__row_list = row_list
        else:
            self.__row_list = []

        if attribute_name_list is not None:
            self.__attribute_name_list = attribute_name_list
        else:
            self.__attribute_name_list = []
        #
        # Derived class data -
        #
        self.__catalog = {}
        self.__num_attributes = 0
        #
        self.__setup()

    def __setup(self):
        self.__num_attributes = len(self.__attribute_name_list)
        self.__catalog = {}
        for attribute_name in self.__attribute_name_list:
            attribute_name_lc = attribute_name.lower()
            self.__catalog[attribute_name_lc] = attribute_name

    def setRowList(self, row_list):
        self.__row_list = row_list

    def setAttributeNameList(self, attribute_name_list):
        self.__attribute_name_list = attribute_name_list
        self.__setup()

    def setName(self, name):
        self.__name = name

    def get(self):
        return (self.__name, self.__attribute_name_list, self.__row_list)


class DataCategory(DataCategoryBase):
    """ Methods for creating, accessing, and formatting PDBx cif data categories.
    """
    def __init__(self, name, attribute_name_list=None, row_list=None):
        super(DataCategory, self).__init__(name, attribute_name_list, row_list)
        #
        self.__lfh = sys.stdout

        self.__current_row_index = 0
        self.__current_attribute = None
        #
        self.__avoid_embedded_quoting = False
        #
        # --------------------------------------------------------------------
        # any whitespace
        self.__ws_re = re.compile(r"\s")
        self.__ws_and_quotes_re = re.compile(r"[\s'\"]")
        # any newline or carriage control
        self.__nl_re = re.compile(r"[\n\r]")
        #
        # single quote
        self.__sq_re = re.compile(r"[']")
        #
        self.__sq_ws_re = re.compile(r"('\s)|(\s')")

        # double quote
        self.__dq_re = re.compile(r'["]')
        self.__dq_ws_re = re.compile(r'("\s)|(\s")')
        #
        self.__int_re = re.compile(r'^[0-9]+$')
        self.__float_re = re.compile(
            r'^-?(([0-9]+)[.]?|([0-9]*[.][0-9]+))([(][0-9]+[)])?([eE][+-]?[0-9]+)?$')
        #
        self.__data_type_list = [
            'DT_NULL_VALUE', 'DT_INTEGER', 'DT_FLOAT', 'DT_UNQUOTED_STRING', 'DT_ITEM_NAME',
            'DT_DOUBLE_QUOTED_STRING', 'DT_SINGLE_QUOTED_STRING', 'DT_MULTI_LINE_STRING',
        ]
        self.__format_type_list = [
            'FT_NULL_VALUE', 'FT_NUMBER', 'FT_NUMBER', 'FT_UNQUOTED_STRING',
            'FT_QUOTED_STRING', 'FT_QUOTED_STRING', 'FT_QUOTED_STRING', 'FT_MULTI_LINE_STRING'
        ]

    def __getitem__(self, name):
        """ Implements list-type functionality -
             Implements op[name] for some special cases -
                name = integer - returns the row in category (normal list behavior)
                name = string - returns the value of attribute 'name' in first row.
        """
        if isinstance(name, int):
            return self.__row_list[name]

        elif isinstance(name, str):
            try:
                return self.__row_list[0][self.getAttributeIndex(name)]
            except (IndexError, KeyError):
                raise KeyError
        raise TypeError(name)

    def getCurrentAttribute(self):
        return self.__current_attribute

    def getRowIndex(self):
        return self.__current_row_index

    def getRowList(self):
        return self.__row_list

    def getRowCount(self):
        return len(self.__row_list)

    def getRow(self, index):
        try:
            return self.__row_list[index]
        except:
            return []

    def removeRow(self, index):
        try:
            if ((index >= 0) and (index < len(self.__row_list))):
                del self.__row_list[index]
                if self.__current_row_index >= len(self.__row_list):
                    self.__current_row_index = len(self.__row_list) - 1
                return True
            else:
                pass
        except:
            pass

        return False

    def getFullRow(self, index):
        """ Return a full row based on the length of the the attribute list.
        """
        try:
            if len(self.__row_list[index] < self.__num_attributes):
                # TODO: 2020/07/02 intendo - range should be over 2 values or defaults
                # from (0 to (x - y) - 1)
                for idx in range(self.__num_attributes - len(self.__row_list[index])):
                    self.__row_list[index].append('?')
            return self.__row_list[index]
        except:
            return ['?' for idx in range(self.__num_attributes)]

    def getName(self):
        return self.__name

    def getAttributeList(self):
        return self.__attribute_name_list

    def getAttributeCount(self):
        return len(self.__attribute_name_list)

    def getAttributeListWithOrder(self):
        ordered_list = []
        for idx, att in enumerate(self.__attribute_name_list):
            ordered_list.append((att, idx))
        return ordered_list

    def getAttributeIndex(self, attribute_name):
        try:
            return self.__attribute_name_list.index(attribute_name)
        except:
            return -1

    def hasAttribute(self, attribute_name):
        return attribute_name in self.__attribute_name_list

    def getIndex(self, attribute_name):
        try:
            return self.__attribute_name_list.index(attribute_name)
        except:
            return -1

    def getItemNameList(self):
        item_name_list = []
        for att in self.__attribute_name_list:
            item_name_list.append("_" + self.__name + "." + att)
        return item_name_list

    def append(self, row):
        # self.__lfh.write("PdbxContainer(append) category %s row %r\n" % (self.__name, row))
        self.__row_list.append(row)

    def appendAttribute(self, attribute_name):
        attribute_name_lc = attribute_name.lower()
        if attribute_name_lc in self.__catalog:
            idx = self.__attribute_name_list.index(self.__catalog[attribute_name_lc])
            self.__attribute_name_list[idx] = attribute_name
            self.__catalog[attribute_name_lc] = attribute_name
            # self.__lfh.write("Appending existing attribute %s\n" % attribute_name)
        else:
            # self.__lfh.write("Appending existing attribute %s\n" % attribute_name)
            self.__attribute_name_list.append(attribute_name)
            self.__catalog[attribute_name_lc] = attribute_name
            #
        self.__num_attributes = len(self.__attribute_name_list)

    def appendAttributeExtendRows(self, attribute_name):
        attribute_name_lc = attribute_name.lower()
        if attribute_name_lc in self.__catalog:
            idx = self.__attribute_name_list.index(self.__catalog[attribute_name_lc])
            self.__attribute_name_list[idx] = attribute_name
            self.__catalog[attribute_name_lc] = attribute_name
            self.__lfh.write("Appending existing attribute %s\n" % attribute_name)
        else:
            self.__attribute_name_list.append(attribute_name)
            self.__catalog[attribute_name_lc] = attribute_name
            # add a placeholder to any existing rows for the new attribute.
            if len(self.__row_list) > 0:
                for row in self.__row_list:
                    row.append("?")
            #
        self.__num_attributes = len(self.__attribute_name_list)

    def getValue(self, attribute_name=None, row_index=None):
        if attribute_name is None:
            attribute = self.__current_attribute
        else:
            attribute = attribute_name
        if row_index is None:
            row_i = self.__current_row_index
        else:
            row_i = row_index

        if isinstance(attribute, str) and isinstance(row_i, int):
            try:
                return self.__row_list[row_i][self.__attribute_name_list.index(attribute)]
            except IndexError:
                raise IndexError
        raise IndexError(attribute)

    def setValue(self, value, attribute_name=None, row_index=None):
        if attribute_name is None:
            attribute = self.__current_attribute
        else:
            attribute = attribute_name

        if row_index is None:
            row_i = self.__current_row_index
        else:
            row_i = row_index

        if isinstance(attribute, str) and isinstance(row_i, int):
            try:
                # if row index is out of range - add the rows -
                for idx in range(row_i + 1 - len(self.__row_list)):
                    self.__row_list.append(self.__emptyRow())
                # self.__row_list[row_i][attribute] = value
                my_ll = len(self.__row_list[row_i])
                ind = self.__attribute_name_list.index(attribute)

                # extend the list if needed -
                if ind >= my_ll:
                    self.__row_list[row_i].extend([None for idx in range(2 * ind - my_ll)])
                self.__row_list[row_i][ind] = value
            except IndexError:
                self.__lfh.write(
                    "DataCategory(setvalue) index error category"
                    " %s attribute %s index %d value %r\n" %
                    (self.__name, attribute, row_i, value))
                traceback.print_exc(file=self.__lfh)
                # raise IndexError
            except ValueError:
                self.__lfh.write(
                    "DataCategory(setvalue) value error category"
                    " %s attribute %s index %d value %r\n" %
                    (self.__name, attribute, row_i, value))
                traceback.print_exc(file=self.__lfh)
                # raise ValueError

    def __emptyRow(self):
        return [None for idx in range(len(self.__attribute_name_list))]

    def replaceValue(self, old_value, new_value, attribute_name):
        num_replace = 0
        if attribute_name not in self.__attribute_name_list:
            return num_replace
        ind = self.__attribute_name_list.index(attribute_name)
        for row in self.__row_list:
            if row[ind] == old_value:
                row[ind] = new_value
                num_replace += 1
        return num_replace

    def replaceSubstring(self, old_value, new_value, attribute_name):
        check = False
        if attribute_name not in self.__attribute_name_list:
            return False
        ind = self.__attribute_name_list.index(attribute_name)
        for row in self.__row_list:
            val = row[ind]
            row[ind] = val.replace(old_value, new_value)
            if val != row[ind]:
                check = True
        return check

    def invokeAttributeMethod(self, attribute_name, type_name, method, db_name):
        self.__current_row_index = 0
        self.__current_attribute = attribute_name
        self.appendAttribute(attribute_name)
        # currentRowIndex = self.__current_row_index # assigned by never used
        #
        ind = self.__attribute_name_list.index(attribute_name)
        if len(self.__row_list) == 0:
            row = [None for idx in range(len(self.__attribute_name_list) * 2)]
            row[ind] = None
            self.__row_list.append(row)

        for row in self.__row_list:
            num = len(row)
            if ind >= num:
                row.extend([None for idx in range(2 * ind - num)])
                row[ind] = None
            exec(method.getInline())
            self.__current_row_index += 1
            # currentRowIndex = self.__current_row_index # assigned by never used

    def invokeCategoryMethod(self, type_name, method, db_name):
        self.__current_row_index = 0
        exec(method.getInline())

    def getAttributeLengthMaximumList(self):
        m_list = [0 for idx in len(self.__attribute_name_list)]
        for row in self.__row_list:
            for indx, val in enumerate(row):
                m_list[indx] = max(m_list[indx], len(val))
        return m_list

    def renameAttribute(self, cur_attribute_name, new_attribute_name):
        """ Change the name of an attribute in place -
        """
        try:
            idx = self.__attribute_name_list.index(cur_attribute_name)
            self.__attribute_name_list[idx] = new_attribute_name
            del self.__catalog[cur_attribute_name.lower()]
            self.__catalog[new_attribute_name.lower()] = new_attribute_name
            return True
        except:
            return False

    def printIt(self, wfh=sys.stdout):
        wfh.write("--------------------------------------------\n")
        wfh.write(" Category: %s attribute list length: %d\n" %
                  (self.__name, len(self.__attribute_name_list)))
        for attr in self.__attribute_name_list:
            wfh.write(" Category: %s attribute: %s\n" % (self.__name, attr))

        wfh.write(" Row value list length: %d\n" % len(self.__row_list))
        #
        for row in self.__row_list[:2]:
            #
            if len(row) == len(self.__attribute_name_list):
                for idx, val in enumerate(row):
                    wfh.write(" %30s: %s ...\n" % (self.__attribute_name_list[idx], str(val)[:30]))
            else:
                wfh.write("+WARNING - %s data length %d attribute name length %s mismatched\n" %
                          (self.__name, len(row), len(self.__attribute_name_list)))

    def dumpIt(self, wfh=sys.stdout):
        wfh.write("--------------------------------------------\n")
        wfh.write(" Category: %s attribute list length: %d\n" %
                  (self.__name, len(self.__attribute_name_list)))
        for attr in self.__attribute_name_list:
            wfh.write(" Category: %s attribute: %s\n" % (self.__name, attr))

        wfh.write(" Value list length: %d\n" % len(self.__row_list))
        for row in self.__row_list:
            for idx, val in enumerate(row):
                wfh.write(" %30s: %s\n" % (self.__attribute_name_list[idx], val))

    def __formatPdbx(self, inp):
        """ Format input data following PDBx quoting rules -
        """
        try:
            if inp is None:
                return ("?", 'DT_NULL_VALUE')

            # pure numerical values are returned as unquoted strings
            if (isinstance(inp, int) or self.__int_re.search(str(inp))):
                return ([str(inp)], 'DT_INTEGER')

            if (isinstance(inp, float) or self.__float_re.search(str(inp))):
                return ([str(inp)], 'DT_FLOAT')

            # null value handling -

            if (inp == "." or inp == "?"):
                return ([inp], 'DT_NULL_VALUE')

            if inp == "":
                return (["."], 'DT_NULL_VALUE')

            # Contains white space or quotes ?
            if not self.__ws_and_quotes_re.search(inp):
                if inp.startswith("_"):
                    return (self.__doubleQuotedList(inp), 'DT_ITEM_NAME')
                else:
                    return ([str(inp)], 'DT_UNQUOTED_STRING')
            else:
                if self.__nl_re.search(inp):
                    return (self.__semiColonQuotedList(inp), 'DT_MULTI_LINE_STRING')
                else:
                    if self.__avoid_embedded_quoting:
                        # change priority to choose double quoting where possible.
                        if not self.__dq_re.search(inp) and not self.__sq_ws_re.search(inp):
                            return (self.__doubleQuotedList(inp), 'DT_DOUBLE_QUOTED_STRING')
                        elif not self.__sq_re.search(inp) and not self.__dq_ws_re.search(inp):
                            return (self.__singleQuotedList(inp), 'DT_SINGLE_QUOTED_STRING')
                        else:
                            return (self.__semiColonQuotedList(inp), 'DT_MULTI_LINE_STRING')
                    else:
                        # change priority to choose double quoting where possible.
                        if not self.__dq_re.search(inp):
                            return (self.__doubleQuotedList(inp), 'DT_DOUBLE_QUOTED_STRING')
                        elif not self.__sq_re.search(inp):
                            return (self.__singleQuotedList(inp), 'DT_SINGLE_QUOTED_STRING')
                        else:
                            return (self.__semiColonQuotedList(inp), 'DT_MULTI_LINE_STRING')
        except:
            traceback.print_exc(file=self.__lfh)

    def __dataTypePdbx(self, inp):
        """ Detect the PDBx data type -
        """
        if inp is None:
            return 'DT_NULL_VALUE'

        # pure numerical values are returned as unquoted strings
        if isinstance(inp, int) or self.__int_re.search(str(inp)):
            return 'DT_INTEGER'

        if isinstance(inp, float) or self.__float_re.search(str(inp)):
            return 'DT_FLOAT'

        # null value handling -

        if (inp == "." or inp == "?"):
            return 'DT_NULL_VALUE'

        if inp == "":
            return 'DT_NULL_VALUE'

        # Contains white space or quotes ?
        if not self.__ws_and_quotes_re.search(inp):
            if inp.startswith("_"):
                return 'DT_ITEM_NAME'
            else:
                return 'DT_UNQUOTED_STRING'
        else:
            if self.__nl_re.search(inp):
                return 'DT_MULTI_LINE_STRING'
            else:
                if self.__avoid_embedded_quoting:
                    if not self.__sq_re.search(inp) and not self.__dq_ws_re.search(inp):
                        return 'DT_DOUBLE_QUOTED_STRING'
                    elif not self.__dq_re.search(inp) and not self.__sq_ws_re.search(inp):
                        return 'DT_SINGLE_QUOTED_STRING'
                    else:
                        return 'DT_MULTI_LINE_STRING'
                else:
                    if not self.__sq_re.search(inp):
                        return 'DT_DOUBLE_QUOTED_STRING'
                    elif not self.__dq_re.search(inp):
                        return 'DT_SINGLE_QUOTED_STRING'
                    else:
                        return 'DT_MULTI_LINE_STRING'

    def __singleQuotedList(self, inp):
        myl = []
        myl.append("'")
        myl.append(inp)
        myl.append("'")
        return myl

    def __doubleQuotedList(self, inp):
        myl = []
        myl.append('"')
        myl.append(inp)
        myl.append('"')
        return myl

    def __semiColonQuotedList(self, inp):
        myl = []
        myl.append("\n")
        if inp[-1] == '\n':
            myl.append(";")
            myl.append(inp)
            myl.append(";")
            myl.append("\n")
        else:
            myl.append(";")
            myl.append(inp)
            myl.append("\n")
            myl.append(";")
            myl.append("\n")
        return myl

    def getValueFormatted(self, attribute_name=None, row_index=None):
        if attribute_name is None:
            attribute = self.__current_attribute
        else:
            attribute = attribute_name

        if row_index is None:
            row_i = self.__current_row_index
        else:
            row_i = row_index

        if isinstance(attribute, str) and isinstance(row_i, int):
            try:
                list_name, type_name = self.__formatPdbx(
                    self.__row_list[row_i][self.__attribute_name_list.index(attribute)])
                return "".join(list_name)
            except IndexError:
                self.__lfh.write(
                    "attribute_name %s row_i %r rowdata %r\n" %
                    (attribute_name, row_i, self.__row_list[row_i]))
                raise IndexError
        raise TypeError(attribute)

    def getValueFormattedByIndex(self, attributeIndex, row_index):
        try:
            list_name, type_name = self.__formatPdbx(self.__row_list[row_index][attributeIndex])
            return "".join(list_name)
        except IndexError:
            raise IndexError

    def getAttributeValueMaxLengthList(self, steps=1):
        mList = [0 for i in range(len(self.__attribute_name_list))]
        for row in self.__row_list[::steps]:
            for indx in range(len(self.__attribute_name_list)):
                val = row[indx]
                mList[indx] = max(mList[indx], len(str(val)))
        return mList

    def getFormatTypeList(self, steps=1):
        try:
            cur_data_type_list = ['DT_NULL_VALUE' for i in range(len(self.__attribute_name_list))]
            for row in self.__row_list[::steps]:
                for indx in range(len(self.__attribute_name_list)):
                    val = row[indx]
                    # print "index ", indx, " val ", val
                    d_type = self.__dataTypePdbx(val)
                    d_indx = self.__data_type_list.index(d_type)
                    # print "d type", dType, " d type index ", dIndx
                    c_type = cur_data_type_list[indx]
                    c_indx = self.__data_type_list.index(c_type)
                    c_indx = max(c_indx, d_indx)
                    cur_data_type_list[indx] = self.__data_type_list[c_indx]

            # Map the format types to the data types
            cur_format_type_list = []
            for value in cur_data_type_list:
                idx = self.__data_type_list.index(value)
                cur_format_type_list.append(self.__format_type_list[idx])
        except:
            self.__lfh.write(
                "PdbxDataCategory(getFormatTypeList) ++Index error at index %d in row %r\n" %
                (indx, row))

        return cur_format_type_list, cur_data_type_list

    def getFormatTypeListX(self):
        cur_data_type_list = ['DT_NULL_VALUE' for i in range(len(self.__attribute_name_list))]
        for row in self.__row_list:
            for indx in range(len(self.__attribute_name_list)):
                val = row[indx]
                # print "index ", indx, " val ", val
                d_type = self.__dataTypePdbx(val)
                d_indx = self.__data_type_list.index(d_type)
                # print "d type", dType, " d type index ", dIndx

                c_type = cur_data_type_list[indx]
                c_indx = self.__data_type_list.index(c_type)
                c_indx = max(c_indx, d_indx)
                cur_data_type_list[indx] = self.__data_type_list[c_indx]

        # Map the format types to the data types
        cur_format_type_list = []
        for dt in cur_data_type_list:
            idx = self.__data_type_list.index(dt)
            cur_format_type_list.append(self.__format_type_list[idx])
        return cur_format_type_list, cur_data_type_list
