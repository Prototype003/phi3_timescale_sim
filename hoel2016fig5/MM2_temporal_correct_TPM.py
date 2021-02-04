import pyphi
import numpy as np
from pprint import pprint
import scipy.io as sio


mm_state = (1,1,1,1)
m_state = (1,1)

tpm_mm = sio.loadmat('./MM2_TemporalExample_Noisy_Ors_Micro.mat')
tpm_m = sio.loadmat('./Noisy_Ors_1timestep_micro.mat')

# ============================================================================
if __name__ == "__main__":
    # micro results
    net_m = pyphi.Network(tpm_m['tpm'])
    mc_m = pyphi.compute.major_complex(net_m, m_state)
    print(mc_m)

    # Macro results
    net_mm = pyphi.Network(tpm_mm['tpm'])
    #M = pyphi.macro.emergence(net_mm, mm_state) #This does not work here, as there are many spatial coarse grainings that are not allowed temporally

    grains = list(pyphi.macro.all_coarse_grains(net_mm.node_indices))

    print(grains[41])    

    coarse_grain = grains[41]
    print(coarse_grain.macro_tpm(net_mm.tpm).reshape([4]+[2], order = 'F'))

    macro_subsystem = pyphi.macro.MacroSubsystem(net_mm, mm_state, coarse_grain=coarse_grain)

    macro_sia = pyphi.compute.sia(macro_subsystem)
    
    print(macro_sia)