# File: PdbxWriterTests.py
# Author: jdw
# Date: 3-November-2009
# Version: 0.001
#
# Update:
# 5-Apr-2011 jdw Using the double quote format preference
# 24-Oct-2012 jdw Update path and examples.
##
"""
Test implementing PDBx/mmCIF write and formatting operations.

"""
import sys
import traceback
import unittest

from pdbx.reader.PdbxReader import PdbxReader
from pdbx.writer.PdbxWriter import PdbxWriter
from pdbx.reader.PdbxContainers import DataContainer, DataCategory


__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


class PdbxWriterTests(unittest.TestCase):
    def setUp(self):
        self.lfh = sys.stderr
        self.verbose = False
        self.path_pdbx_data_file = "../tests/1kip.cif"
        self.path_output_file = "testOutputDataFile.cif"

    def tearDown(self):
        pass

    def testWriteDataFile(self):
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
            a_cat.append((1, 2, 3, 4, 5, 6, 7))
            a_cat.append((1, 2, 3, 4, 5, 6, 7))
            a_cat.append((1, 2, 3, 4, 5, 6, 7))
            a_cat.append((1, 2, 3, 4, 5, 6, 7))
            cur_container.append(a_cat)
            my_data_list.append(cur_container)
            pdbx_w = PdbxWriter(ofh)
            pdbx_w.write(my_data_list)
            ofh.close()
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def testUpdateDataFile(self):
        """Test case - write data file
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            # Create a initial data file --
            #
            my_data_list = []
            ofh = open("test-output-1.cif", "w")
            cur_container = DataContainer("myblock")
            a_cat = DataCategory("pdbx_seqtool_mapping_ref")
            a_cat.appendAttribute("ordinal")
            a_cat.appendAttribute("entity_id")
            a_cat.appendAttribute("auth_mon_id")
            a_cat.appendAttribute("auth_mon_num")
            a_cat.appendAttribute("pdb_chain_id")
            a_cat.appendAttribute("ref_mon_id")
            a_cat.appendAttribute("ref_mon_num")
            a_cat.append((1, 2, 3, 4, 5, 6, 7))
            a_cat.append((1, 2, 3, 4, 5, 6, 7))
            a_cat.append((1, 2, 3, 4, 5, 6, 7))
            a_cat.append((1, 2, 3, 4, 5, 6, 7))
            cur_container.append(a_cat)
            my_data_list.append(cur_container)
            pdbx_w = PdbxWriter(ofh)
            pdbx_w.write(my_data_list)
            ofh.close()
            #
            # Read and update the data -
            #
            my_data_list = []
            ifh = open("test-output-1.cif", "r")
            prd = PdbxReader(ifh)
            prd.read(my_data_list)
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

        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def testReadDataFile(self):
        """Test case - read data file
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            my_data_list = []
            ifh = open(self.path_pdbx_data_file, "r")
            prd = PdbxReader(ifh)
            prd.read(my_data_list)
            ifh.close()
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def testReadWriteDataFile(self):
        """Test case - data file read write test
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
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
            traceback.print_exc(file=sys.stderr)
            self.fail()


def suite():
    return unittest.makeSuite(PdbxWriterTests, 'test')


if __name__ == '__main__':
    unittest.main()
