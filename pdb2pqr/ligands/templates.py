"""There aren't many templates here..."""
templates = {}
## templates['CarboxyGroup']={'O10':{'sybyl_type':'O.co2','neighbours':['C3'],'alreadyvisited':False},
##                            'O20':{'sybyl_type':'O.co2','neighbours':['C3'],'alreadyvisited':False},
##                            'C3': {'sybyl_type':'C.2',  'neighbours':['O10','O20'],'alreadyvisited':False}}
templates['AceticAcid']={'O1' :{'sybyl_type':'O.co2','neighbours':['C'],'alreadyvisited':False},
                         'O2' :{'sybyl_type':'O.co2','neighbours':['C'],'alreadyvisited':False},
                         'C'  :{'sybyl_type':'C.2'  ,'neighbours':['O1','O2','C40'],'alreadyvisited':False},
                         'C40':{'sybyl_type':'C.3'  ,'neighbours':['C'],'alreadyvisited':False}}
templates['IDD594']={'O1':  {'sybyl_type':'O.co2','neighbours':['C'],'alreadyvisited':False},
                     'O2':  {'sybyl_type':'O.co2','neighbours':['C'],'alreadyvisited':False},
                     'C' :  {'sybyl_type':'C.2',  'neighbours':['O1','O2','C40'],'alreadyvisited':False},
                     'C40': {'sybyl_type':'C.3',  'neighbours':['C','O50'],'alreadyvisited':False},
                     'O50': {'sybyl_type':'O.3',  'neighbours':['C40','C60'],'alreadyvisited':False},
                     'C60': {'sybyl_type':'C.ar', 'neighbours':['O50'],'alreadyvisited':False}}
templates['acetylsalicylicacid']={'O1' :{'sybyl_type':'O.co2','neighbours':['C'],'alreadyvisited':False},
                                   'O2' :{'sybyl_type':'O.co2','neighbours':['C'],'alreadyvisited':False},
                                   'C'  :{'sybyl_type':'C.2'  ,'neighbours':['O1','O2','C40'],'alreadyvisited':False},
                                   'C20'  :{'sybyl_type':'C.3'  ,'neighbours':['C','C40','C30'],'alreadyvisited':False},
                                   'C30'  :{'sybyl_type':'C.3'  ,'neighbours':['C20','C50'],'alreadyvisited':False},
                                   'C50'  :{'sybyl_type':'C.3'  ,'neighbours':['C30','C80'],'alreadyvisited':False},
                                   'C80'  :{'sybyl_type':'C.3'  ,'neighbours':['C50','C60'],'alreadyvisited':False},
                                   'C60'  :{'sybyl_type':'C.3'  ,'neighbours':['C80','C40'],'alreadyvisited':False},
                                   'C40'  :{'sybyl_type':'C.3'  ,'neighbours':['C','C60'],'alreadyvisited':False}}
#templates['CRAP']={'O10':{'sybyl_type':'O.co2','neighbours':['C30'],'alreadyvisited':False},
#                   'O20':{'sybyl_type':'O.co2','neighbours':['C30'],'alreadyvisited':False},
#                   'C30':{'sybyl_type':'C.2'  ,'neighbours':['O10','O20','C40'],'alreadyvisited':False},
#                   'C40':{'sybyl_type':'C.3'  ,'neighbours':['C30','O50'],'alreadyvisited':False},
#                   'O50':{'sybyl_type':'O.3',  'neighbours':['C40'],'alreadyvisited':False}}
## templates['PropanoicAcid']={'O10':{'sybyl_type':'O.co2','neighbours':['C30'],'alreadyvisited':False},
##                             'O20':{'sybyl_type':'O.co2','neighbours':['C30'],'alreadyvisited':False},
##                             'C30':{'sybyl_type':'C.2'  ,'neighbours':['O10','O20','C40'],'alreadyvisited':False},
##                             'C40':{'sybyl_type':'C.3'  ,'neighbours':['C30','C50'],'alreadyvisited':False},
##                             'C50':{'sybyl_type':'C.3'  ,'neighbours':['C40'],'alreadyvisited':False}}
## templates['piperidine']={'N1':{ 'sybyl_type':'N.4','neighbours':['C1','C5'],'alreadyvisited':False},
##                          'C1':{ 'sybyl_type':'C.3','neighbours':['N1','C2'],'alreadyvisited':False},
##                          'C2':{ 'sybyl_type':'C.3','neighbours':['C1','C3'],'alreadyvisited':False},
##                          'C3':{ 'sybyl_type':'C.3','neighbours':['C2','C4'],'alreadyvisited':False},
##                          'C4':{ 'sybyl_type':'C.3','neighbours':['C3','C5'],'alreadyvisited':False},
##                          'C5':{ 'sybyl_type':'C.3','neighbours':['N1','C4'],'alreadyvisited':False}}
## templates['imidazole']={'N10':{ 'sybyl_type':'N.pl3','neighbours':['C30','C10'],'alreadyvisited':False},
##                         'C10':{ 'sybyl_type':'C.2',  'neighbours':['N10','N20'],'alreadyvisited':False},
##                         'N20':{ 'sybyl_type':'N.2',  'neighbours':['C10','C20'],'alreadyvisited':False},
##                         'C20':{ 'sybyl_type':'C.2',  'neighbours':['N20','C30'],'alreadyvisited':False},
##                         'C30':{ 'sybyl_type':'C.2',  'neighbours':['C20','N10'],'alreadyvisited':False}}
## templates['Ethine']={'C10':{ 'sybyl_type':'C.1','neighbours':['C20'],'alreadyvisited':False},
##                      'C20':{ 'sybyl_type':'C.1','neighbours':['C10'],'alreadyvisited':False}}
