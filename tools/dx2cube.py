import argparse

parser = argparse.ArgumentParser(description='Translate a dx file to cube format')
parser.add_argument('dx_input',
                   help='Name of the dx_input file (required arguement)')
parser.add_argument('pqr_input',
                   help='Name of the pqr_input file (required arguement)')
parser.add_argument('output',
                   help='Name of the output file (required arguement)')

args = parser.parse_args()  
counter = 0                 

#DX STUFF
if args.dx_input.endswith('.dx'):
    print "Success"
else:
    print "Error converting file"

try:
    with open(args.dx_input, 'r') as in_f, open(args.output, 'w') as out_f, open(args.pqr_input, 'r') as in_pqr:
        out_f.write("CPMD CUBE FILE.\n"
                    "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
        
        #Discard comments at top of file.
        line = in_f.readline()
        newline = in_pqr.readline()
        while line.startswith('#'):
            line = in_f.readline()
        while newline.startswith('REMARK'):
            newline = in_pqr.readline()
        
        split_line = line.split()        
        grid_sizes = [int(x)*-1 for x in split_line[-3:]]
                    
        split_line = in_f.readline().split()
        
        origin = [float(x) for x in split_line[-3:]]
        
        parameter_fmt = "{:>4} {:>11.6f} {:>11.6f} {:>11.6f}\n"
        try:    
            while newline.startswith('ATOM'):
                pqr_line =  in_pqr.readline()
                new_split_line = pqr_line.split()
                atom_num = new_split_line[1]
        except IndexError:
            pass
        in_pqr.seek(0)
        newline = in_pqr.readline()
        while newline.startswith('REMARK'):
            newline = in_pqr.readline()
        
        origin_line = parameter_fmt.format(atom_num, *origin)
        out_f.write(origin_line)
        
        
        for x in xrange(3):
            split_line = in_f.readline().split()
            grid_dims = [float(item) for item in split_line[-3:]]
            
            dim_lin = parameter_fmt.format(grid_sizes[x], *grid_dims)
            out_f.write(dim_lin)
            
        atoms_parameter_fmt = "{:>4} {:>11.6f} {:>11.6f} {:>11.6f} {:>11.6f}\n"
        a = True
        xreal_center = []
        yreal_center = []
        zreal_center = []
        try:
            while a == True:
                if not newline.startswith('TER'):
                    new_split_line = newline.split()
                    radius = new_split_line[-1]
                    xyz = new_split_line[-5:-2]
                    line_atom_num = new_split_line[1]
                    atom_radius = new_split_line[-1]
                    pqr_lin = atoms_parameter_fmt.format(int(line_atom_num), float(new_split_line[-2]), float(xyz[0]), float(xyz[1]), float(xyz[2]))
                    out_f.write(pqr_lin)
                    newline = in_pqr.readline()
                    xreal_center.append(float(xyz[0]))
                    yreal_center.append(float(xyz[1]))
                    zreal_center.append(float(xyz[2]))                    
                else:
                    a = False
        except IndexError:
            pass
            
        x_avg = sum(xreal_center)/float(atom_num)
        y_avg = sum(yreal_center)/float(atom_num)
        z_avg = sum(zreal_center)/float(atom_num)       
        print x_avg, y_avg, z_avg
            
        #print origin
        #new_origin = []
        #for item in origin:
        #    newitem = item/0.529177
            #new_new = item/2 + newitem/2
        #    new_origin.append(newitem)
        #print new_origin

        #Consume unneeded object lines.
        in_f.readline()
        in_f.readline()
        
        ##TODO: put atoms here
        
        value_format = ["{:< 13.5E}"]
        value_format = ' '.join(value_format * 6) + '\n'
        print value_format 
        group = []
        line = in_f.readline()
        while not line.startswith('attribute'):
            values = [float(item) for item in line.split()]
            group.extend(values)
            
            if len(group) >= 6:
                out_f.write(value_format.format(*group))
                group = []
                
            line = in_f.readline()
            
        if group:
            group_strs = ["{:< 13.5E}".format(item) for item in group]
            out_f.write(' '.join(group_strs))
        
            

except IOError:
    print "file doesn't exist"
'''
for line in dx_data:
    dx_data_2 = dx_data
    data_point = line.strip("\n")
    dx_data_2.pop(0)
    dx_data._2append(data_point)
    

for line in dx_data_2:
    if line.startswith('#'):
        dx_data_2.remove(line)
        continue
    elif "gridpositions" in line:
        new_line = line.split()
        grid_values = new_line[-3:]
        #grid values for cube file, will need grid_valies[0]
        break
    elif line.startswith('delta'):
        line.insert('0')
        line.remove(delta)



#PQR STUFF
if args.pqr_input.endswith('.pqr'):
    print "yay"
else:
    print "nay"

try:
    with open(args.pqr_input, 'r') as f:
        pqr_data = f.readlines()
        
except IOError:
    print "file doesn't exist"
    
for line in pqr_data:
    data_point = line.strip("\n")
    pqr_data.pop(0)
    pqr_data.append(data_point)

for line in pqr_data:
    if line.startswith('REMARK'):
        pqr_data.remove(line)
        continue
    counter +=1
    #the number next to the origin in cube format
'''