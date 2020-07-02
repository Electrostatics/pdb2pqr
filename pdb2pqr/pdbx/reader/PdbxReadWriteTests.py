##
# File: PdbxReadWriteTests.py
# Author: jdw
# Date: 9-Oct-2011
# Version: 0.001
#
# Updated:
# 24-Oct-2012 jdw update path details and reorganize.
#
##
""" Various tests caess for PDBx/mmCIF data file and dictionary reader and writer.
"""
import sys
import unittest
import traceback

from pdbx.reader.PdbxReader import PdbxReader
from pdbx.writer.PdbxWriter import PdbxWriter
from pdbx.reader.PdbxContainers import DataContainer, DataCategory

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


class PdbxReadWriteTests(unittest.TestCase):
    def setUp(self):
        self.lfh = sys.stdout
        self.verbose = False
        self.path_pdbx_data_file = "../tests/1kip.cif"
        self.path_output_file = "testOutputDataFile.cif"

    def tearDown(self):
        pass

    def test_simple_initialization(self):
        """Test case - Simple initialization of a data category and data block
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            filename = "test-simple.cif"
            attribute_name_list = [
                'aOne', 'aTwo', 'aThree', 'aFour', 'aFive', 'aSix', 'aSeven', 'aEight', 'aNine',
                'aTen'
            ]
            row_list = [
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            ]
            name_cat = 'myCategory'
            #
            #
            cur_container = DataContainer("myblock")
            a_cat = DataCategory(name_cat, attribute_name_list, row_list)
            a_cat.printIt()
            cur_container.append(a_cat)
            cur_container.printIt()
            #
            my_container_list = []
            my_container_list.append(cur_container)
            ofh = open(filename, "w")
            pdbx_w = PdbxWriter(ofh)
            pdbx_w.write(my_container_list)
            ofh.close()

            my_container_list = []
            ifh = open(filename, "r")
            prd = PdbxReader(ifh)
            prd.read(my_container_list)
            ifh.close()
            for container in my_container_list:
                for obj_name in container.getObjNameList():
                    name, a_list, r_list = container.getObj(obj_name).get()
                    self.lfh.write("Recovered data category %s\n" % name)
                    self.lfh.write("Attribute list %r\n" % repr(a_list))
                    self.lfh.write("Row list %r\n" % repr(r_list))
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()

    def test_write_data_file(self):
        """Test case - write data file
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            my_data_list = []
            ofh = open("test-output.cif", "w")
            cur_container = DataContainer("myblock")
            a_cat = DataCategory("pdbx_seqtool_mapping_ref")
            a_cat.appendAttribute("ordinal")
            a_cat.appendAttribute("entity_id")
            a_cat.appendAttribute("auth_mon_id")
            a_cat.appendAttribute("auth_mon_num")
            a_cat.appendAttribute("pdb_chain_id")
            a_cat.appendAttribute("ref_mon_id")
            a_cat.appendAttribute("ref_mon_num")
            a_cat.append([1, 2, 3, 4, 5, 6, 7])
            a_cat.append([1, 2, 3, 4, 5, 6, 7])
            a_cat.append([1, 2, 3, 4, 5, 6, 7])
            a_cat.append([1, 2, 3, 4, 5, 6, 7])
            a_cat.append([7, 6, 5, 4, 3, 2, 1])
            a_cat.printIt()
            cur_container.append(a_cat)
            cur_container.printIt()
            #
            my_data_list.append(cur_container)
            pdbx_w = PdbxWriter(ofh)
            pdbx_w.write(my_data_list)
            ofh.close()
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()

    def test_update_data_file(self):
        """Test case - update data file
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            # Create a initial data file --
            #
            my_data_list = []

            cur_container = DataContainer("myblock")
            a_cat = DataCategory("pdbx_seqtool_mapping_ref")
            a_cat.appendAttribute("ordinal")
            a_cat.appendAttribute("entity_id")
            a_cat.appendAttribute("auth_mon_id")
            a_cat.appendAttribute("auth_mon_num")
            a_cat.appendAttribute("pdb_chain_id")
            a_cat.appendAttribute("ref_mon_id")
            a_cat.appendAttribute("ref_mon_num")
            a_cat.append([9, 2, 3, 4, 5, 6, 7])
            a_cat.append([10, 2, 3, 4, 5, 6, 7])
            a_cat.append([11, 2, 3, 4, 5, 6, 7])
            a_cat.append([12, 2, 3, 4, 5, 6, 7])

            # self.lfh.write("Assigned data category state-----------------\n")
            # aCat.dumpIt(fh = self.lfh)

            cur_container.append(a_cat)
            my_data_list.append(curContainer)
            ofh = open("test-output-1.cif", "w")
            pdbx_w = PdbxWriter(ofh)
            pdbx_w.write(my_data_list)
            ofh.close()
            #
            #
            # Read and update the data -
            #
            my_data_list = []
            ifh = open("test-output-1.cif", "r")
            pRd = PdbxReader(ifh)
            pRd.read(my_data_list)
            ifh.close()
            #
            my_block = my_data_list[0]
            my_block.printIt()
            my_cat = my_block.getObj('pdbx_seqtool_mapping_ref')
            my_cat.printIt()
            for i_row in range(0, my_cat.getRowCount()):
                my_cat.setValue('some value', 'ref_mon_id', i_row)
                my_cat.setValue(100, 'ref_mon_num', i_row)
            ofh = open("test-output-2.cif", "w")
            pdbx_w = PdbxWriter(ofh)
            pdbx_w.write(my_data_list)
            ofh.close()
            #

        except:
            traceback.print_exc(file=self.lfh)
            self.fail()

    def test_read_data_file(self):
        """Test case - read data file
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            my_data_list = []
            ifh = open(self.path_pdbx_data_file, "r")
            prd = PdbxReader(ifh)
            prd.read(my_data_list)
            ifh.close()
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()

    def test_read_write_data_file(self):
        """Test case - data file read write test
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            my_data_list = []
            ifh = open(self.path_pdbx_data_file, "r")
            prd = PdbxReader(ifh)
            prd.read(my_data_list)
            ifh.close()

            ofh = open(self.path_output_file, "w")
            pwr = PdbxWriter(ofh)
            pwr.write(my_data_list)
            ofh.close()
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()


def simple_suite():
    suite_select = unittest.TestSuite()
    suite_select.addTest(PdbxReadWriteTests("test_simple_initialization"))
    suite_select.addTest(PdbxReadWriteTests("test_update_data_file"))
    suite_select.addTest(PdbxReadWriteTests("test_read_write_data_file"))
    return suite_select


if __name__ == '__main__':
    #
    MYSUITE = simple_suite()
    unittest.TextTestRunner(verbosity=2).run(MYSUITE)
