import numpy as np
# DOMAIN SETTINGS SPHERE
nSubdiv = 1
i_loadGrid = 0 #<=4

# TIME STEPS
dt = 1600/2**nSubdiv
ndays = 1.0
save_nth_hour = 12
nts = np.ceil((ndays*24*3600)/dt).astype(np.int)+1
save_nth_ts = np.ceil(save_nth_hour*3600/dt).astype(np.int)



omega = 7.2921E-5
#omega = 0
