##
# File: PdbxWriter.py
# Date: 2011-10-09 Jdw Adapted from PdbxParser.py
#
# Updates:
# 5-Apr-2011 jdw Using the double quote format preference
# 23-Oct-2012 jdw update path details and reorganize.
#
###
"""
Classes for writing data and dictionary containers in PDBx/mmCIF format.

"""
import sys
from pdbx.reader.PdbxContainers import DefinitionContainer, DataContainer

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


class PdbxError(Exception):
    """ Class for catch general errors
    """
    pass


class PdbxWriter:
    """Write PDBx data files or dictionaries using the input container
       or container list.
    """
    def __init__(self, ofh=sys.stdout):
        self.__ofh = ofh
        self.__container_list = []
        #self.__MAXIMUM_LINE_LENGTH = 2048
        self.__spacing = 2
        self.__indent_definition = 3
        self.__indent_space = " " * self.__indent_definition
        self.__do_definition_indent = False
        # Maximum number of rows checked for value length and format
        self.__row_partition = None

    def setRowPartition(self, num_rows):
        ''' Maximum number of rows checked for value length and format
        '''
        self.__row_partition = num_rows

    def write(self, container_list):
        self.__container_list = container_list
        for container in self.__container_list:
            self.writeContainer(container)

    def writeContainer(self, container):
        ind_s = " " * self.__indent_definition
        if isinstance(container, DefinitionContainer):
            self.__write("save_%s\n" % container.getName())
            self.__do_definition_indent = True
            self.__write(ind_s + "#\n")
        elif isinstance(container, DataContainer):
            if container.getGlobal():
                self.__write("global_\n")
                self.__do_definition_indent = False
                self.__write("\n")
            else:
                self.__write("data_%s\n" % container.getName())
                self.__do_definition_indent = False
                self.__write("#\n")

        for name in container.getObjNameList():
            obj = container.getObj(name)
            obj_l = obj.getRowList()

            # Skip empty objects
            if len(obj_l) == 0:
                continue

            # Item - value formattting
            elif len(obj_l) == 1:
                self.__writeItemValueFormat(obj)

            # Table formatting -
            elif len(obj_l) > 1 and len(obj.getAttributeList()) > 0:
                self.__writeTableFormat(obj)
            else:
                raise PdbxError()

            if self.__do_definition_indent:
                self.__write(ind_s + "#")
            else:
                self.__write("#")

        # Add a trailing saveframe reserved word
        if isinstance(container, DefinitionContainer):
            self.__write("\nsave_\n")
        self.__write("#\n")

    def __write(self, value):
        self.__ofh.write(value)

    def __writeItemValueFormat(self, my_category):

        # Compute the maximum item name length within this category -
        attribute_name_length_max = 0
        for attribute_name in my_category.getAttributeList():
            attribute_name_length_max = max(attribute_name_length_max, len(attribute_name))
        item_name_length_max = self.__spacing + len(my_category.getName())
        item_name_length_max += attribute_name_length_max + 2
        #
        line_list = []
        line_list.append("#\n")
        for attribute_name, i_pos in my_category.getAttributeListWithOrder():
            if self.__do_definition_indent:
                # - add indent --
                line_list.append(self.__indent_space)

            item_name = "_%s.%s" % (my_category.getName(), attribute_name)
            line_list.append(item_name.ljust(item_name_length_max))

            line_list.append(my_category.getValueFormatted(attribute_name, 0))
            line_list.append("\n")

        self.__write("".join(line_list))

    def __writeTableFormat(self, my_category):

        # Write the declaration of the loop_
        #
        line_list = []
        line_list.append('#\n')
        if self.__do_definition_indent:
            line_list.append(self.__indent_space)
        line_list.append("loop_")
        for attribute_name in my_category.getAttributeList():
            line_list.append('\n')
            if self.__do_definition_indent:
                line_list.append(self.__indent_space)
            item_name = "_%s.%s" % (my_category.getName(), attribute_name)
            line_list.append(item_name)
        self.__write("".join(line_list))

        #
        # Write the data in tabular format -
        #
        # print my_category.getName()
        # print my_category.getAttributeList()

        # For speed make the following evaluation on a portion of the table
        num_steps = 1
        if self.__row_partition is not None:
            num_steps = max(1, my_category.getRowCount() / self.__row_partition)

        format_type_list, data_type_list = my_category.getFormatTypeList(steps=num_steps)
        max_length_list = my_category.getAttributeValueMaxLengthList(steps=num_steps)
        spacing = " " * self.__spacing
        #

        # print format_type_list
        # print data_type_list
        # print max_length_list
        #
        for i_row in range(my_category.getRowCount()):
            line_list = []
            line_list.append('\n')
            if self.__do_definition_indent:
                line_list.append(self.__indent_space + " ")

            for i_at in range(my_category.getAttributeCount()):
                format_type = format_type_list[i_at]
                max_length = max_length_list[i_at]

                if (format_type == 'FT_UNQUOTED_STRING' or format_type == 'FT_NULL_VALUE'):
                    val = my_category.getValueFormattedByIndex(i_at, i_row)
                    line_list.append(val.ljust(max_length))

                elif format_type == 'FT_NUMBER':
                    val = my_category.getValueFormattedByIndex(i_at, i_row)
                    line_list.append(val.rjust(max_length))

                elif format_type == 'FT_QUOTED_STRING':
                    val = my_category.getValueFormattedByIndex(i_at, i_row)

                    line_list.append(val.ljust(max_length + 2))

                elif format_type == "FT_MULTI_LINE_STRING":
                    val = my_category.getValueFormattedByIndex(i_at, i_row)
                    line_list.append(val)

                line_list.append(spacing)

            self.__write("".join(line_list))
        self.__write("\n")
