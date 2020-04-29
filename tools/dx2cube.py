"""Utility to translate a DX file to CUBE format"""

import argparse

PARSER = argparse.ArgumentParser(description='Translate a dx file to cube format')
PARSER.add_argument('dx_input',
                    help='Name of the dx_input file (required arguement)')
PARSER.add_argument('pqr_input',
                    help='Name of the pqr_input file (required arguement)')
PARSER.add_argument('output',
                    help='Name of the output file (required arguement)')

ARGS = PARSER.parse_args()

#DX STUFF
if ARGS.dx_input.endswith('.dx'):
    print("Success")
else:
    print("Error converting file")

try:
    with open(ARGS.dx_input, 'r') as in_f,\
         open(ARGS.output, 'w') as out_f,\
         open(ARGS.pqr_input, 'r') as in_pqr:
        out_f.write("CPMD CUBE FILE.\n"
                    "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")

        # Discard comments at top of file.
        LINE = in_f.readline()
        NEWLINE = in_pqr.readline()
        while LINE.startswith('#'):
            LINE = in_f.readline()
        while NEWLINE.startswith('REMARK'):
            NEWLINE = in_pqr.readline()

        SPLIT_LINE = LINE.split()
        GRID_SIZES = [int(x)*-1 for x in SPLIT_LINE[-3:]]

        SPLIT_LINE = in_f.readline().split()

        ORIGIN = [float(x) for x in SPLIT_LINE[-3:]]

        PARAMETER_FMT = "{:>4} {:>11.6f} {:>11.6f} {:>11.6f}\n"
        try:
            while NEWLINE.startswith('ATOM'):
                PQR_LINE = in_pqr.readline()
                NEW_SPLIT_LINE = PQR_LINE.split()
                ATOM_NUM = NEW_SPLIT_LINE[1]
        except IndexError:
            pass
        in_pqr.seek(0)
        NEWLINE = in_pqr.readline()
        while NEWLINE.startswith('REMARK'):
            NEWLINE = in_pqr.readline()

        ORIGIN_LINE = PARAMETER_FMT.format(ATOM_NUM, *ORIGIN)
        out_f.write(ORIGIN_LINE)


        for x in range(3):
            SPLIT_LINE = in_f.readline().split()
            grid_dims = [float(item) for item in SPLIT_LINE[-3:]]

            dim_lin = PARAMETER_FMT.format(GRID_SIZES[x], *grid_dims)
            out_f.write(dim_lin)

        ATOMS_PARAMETER_FMT = "{:>4} {:>11.6f} {:>11.6f} {:>11.6f} {:>11.6f}\n"
        A = True
        XREAL_CENTER = []
        YREAL_CENTER = []
        ZREAL_CENTER = []
        try:
            while A is True:
                if not NEWLINE.startswith('TER'):
                    NEW_SPLIT_LINE = NEWLINE.split()
                    RADIUS = NEW_SPLIT_LINE[-1]
                    XYZ = NEW_SPLIT_LINE[-5:-2]
                    LINE_ATOM_NUM = NEW_SPLIT_LINE[1]
                    ATOM_RADIUS = NEW_SPLIT_LINE[-1]
                    PQR_LIN = ATOMS_PARAMETER_FMT.format(
                        int(LINE_ATOM_NUM),
                        float(NEW_SPLIT_LINE[-2]),
                        float(XYZ[0]), float(XYZ[1]),
                        float(XYZ[2]))
                    out_f.write(PQR_LIN)
                    NEWLINE = in_pqr.readline()
                    XREAL_CENTER.append(float(XYZ[0]))
                    YREAL_CENTER.append(float(XYZ[1]))
                    ZREAL_CENTER.append(float(XYZ[2]))
                else:
                    A = False
        except IndexError:
            pass

        X_AVG = sum(XREAL_CENTER)/float(ATOM_NUM)
        Y_AVG = sum(YREAL_CENTER)/float(ATOM_NUM)
        Z_AVG = sum(ZREAL_CENTER)/float(ATOM_NUM)
        print(X_AVG, Y_AVG, Z_AVG)

        # print origin
        # new_origin = []
        # for item in origin:
        #    newitem = item/0.529177
        #    new_new = item/2 + newitem/2
        #    new_origin.append(newitem)
        # print new_origin

        # Consume unneeded object lines.
        in_f.readline()
        in_f.readline()

        ## TODO_OLD: put atoms here - This NOTE is over 5 years old

        VALUE_FORMAT = ["{:< 13.5E}"]
        VALUE_FORMAT = ' '.join(VALUE_FORMAT * 6) + '\n'
        print(VALUE_FORMAT)
        GROUP = []
        LINE = in_f.readline()
        while not LINE.startswith('attribute'):
            VALUES = [float(ITEM) for ITEM in LINE.split()]
            GROUP.extend(VALUES)

            if len(GROUP) >= 6:
                out_f.write(VALUE_FORMAT.format(*GROUP))
                GROUP = []

            LINE = in_f.readline()

        if GROUP:
            GROUP_STRS = ["{:< 13.5E}".format(ITEM) for ITEM in GROUP]
            out_f.write(' '.join(GROUP_STRS))



except IOError:
    print("file doesn't exist")
