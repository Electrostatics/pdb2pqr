##
# File: PdbxReaderTests.py
# Author: jdw
# Date: 9-Jan-2012
# Version: 0.001
#
# Update:
# 27-Sep-2012 jdw add test case for reading PDBx structure factor file
#
#  TODO: 2020/07/02 intendo - Update to Python3 (e.g. use "with open" for files)
##
"""
Test cases for reading PDBx/mmCIF data files PdbxReader class -

"""
import sys
import unittest
import traceback

from pdbx.reader.PdbxReader import PdbxReader


class PdbxReaderTests(unittest.TestCase):
    def setUp(self):
        self.lfh = sys.stderr
        self.verbose = False
        self.path_pdbx_data_file = "../tests/1kip.cif"
        self.path_big_pdbx_data_file = "../tests/1ffk.cif"
        self.path_sf_data_file = "../tests/1kip-sf.cif"

    def tearDown(self):
        pass

    def test_read_small_data_file(self):
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

    def test_read_big_data_file(self):
        """Test case - read data file
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            my_data_list = []
            ifh = open(self.path_big_pdbx_data_file, "r")
            prd = PdbxReader(ifh)
            prd.read(my_data_list)
            ifh.close()
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def test_read_sf_data_file(self):
        """Test case - read PDB structure factor data file and compute statistics on f/sig(f).
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            my_container_list = []
            ifh = open(self.path_sf_data_file, "r")
            prd = PdbxReader(ifh)
            prd.read(my_container_list)
            first_item = my_container_list[0]
            #
            cat_obj = first_item.getObj("refln")
            if cat_obj is None:
                # TODO: 2020/07/02 intendo - Shouldn't this just be self.fail()?
                return False

            # nRows = catObj.getRowCount()  # assigned but never used
            #
            # Get column name index.
            #
            it_dict = {}
            it_name_list = cat_obj.getItemNameList()
            for idx_it, it_name in enumerate(it_name_list):
                it_dict[str(it_name).lower()] = idx_it
                #
            idf = it_dict['_refln.f_meas_au']
            idsigf = it_dict['_refln.f_meas_sigma_au']
            min_r = 100
            max_r = -1
            sum_r = 0
            icount = 0
            for row in cat_obj.getRowList():
                try:
                    f_val = float(row[idf])
                    sigf = float(row[idsigf])
                    ratio = sigf / f_val
                    # self.lfh.write(" %f %f %f\n" % (f, sigf, ratio))
                    max_r = max(max_r, ratio)
                    min_r = min(min_r, ratio)
                    sum_r += ratio
                    icount += 1
                except:
                    # TODO: 2020/07/02 intendo - Try/Except should not be used this way
                    continue

            ifh.close()
            self.lfh.write(
                "f/sig(f) min %f max %f avg %f count %d\n" %
                (min_r, max_r, sum_r / icount, icount))
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()


def simple_suite():
    suite_select = unittest.TestSuite()
    suite_select.addTest(PdbxReaderTests("test_read_big_data_file"))
    suite_select.addTest(PdbxReaderTests("test_read_small_data_file"))
    suite_select.addTest(PdbxReaderTests("test_read_sf_data_file"))
    return suite_select


if __name__ == '__main__':
    MY_SUITE = simple_suite()
    unittest.TextTestRunner(verbosity=2).run(MY_SUITE)
