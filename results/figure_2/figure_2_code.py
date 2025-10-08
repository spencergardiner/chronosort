import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
 
legend_mapping = {
    'seq51': 'WT',
    'seq42': 'A.02',
    'seq43': 'A.03',
    'seq44': 'A.04',
    'seq45': 'A.05',
    'seq46': 'A.06',
    'seq47': 'A.07',
    'seq48': 'A.08',
    'seq49': 'A.09',
    'seq410': 'A.10',
    'seq41': 'A.01',
    'seq52': 'B.02',
    'seq53': 'B.07',
    'seq54': 'B.03',
    'seq55': 'B.04',
    'seq56': 'B.05',
    'seq57': 'B.06',
    'seq58': 'B.08',
    'seq59': 'B.09',
    'seq510': 'B.10'
}
# Collect all seq*/run*/rmsf.xvg data
seq_dirs = sorted([d for d in os.listdir('.') if d.startswith('seq') and os.path.isdir(d)])
 
avg_rmsf_by_seq = {}
residue_indices = None
 
for seq in seq_dirs:
    run_rmsfs = []
    run_dirs = sorted([os.path.join(seq, d) for d in os.listdir(seq) if d.startswith('run') and os.path.isdir(os.path.join(seq, d))])
 
    for run in run_dirs:
        rmsf_path = os.path.join(run, 'rmsf.xvg')
        if os.path.isfile(rmsf_path):
            try:
                data = np.loadtxt(rmsf_path, comments=('@', '#'))
                if residue_indices is None:
                    residue_indices = data[:, 0]
                run_rmsfs.append(data[:, 1])
            except Exception as e:
                print(f"Failed to read {rmsf_path}: {e}")
 
 
    rmsf_array = np.array(run_rmsfs)
    avg_rmsf_by_seq[legend_mapping[seq]] = rmsf_array.mean(axis=0)
 
 
avg_rmsf=list(avg_rmsf_by_seq.values())
plt.rcParams.update({'font.size': 12})
 
 
plt.figure(figsize=(10, 6))
stats_df = pd.read_csv('stats.csv')
 
normalize=50
# Divide the values by 100 and add 2
processed_stats = (stats_df['Mean'] / normalize)
plt.plot(np.abs(processed_stats-2), color='black',linestyle='--', label='AF3 pLDDT')
 
 
stats_df = pd.read_csv('avg_rmsf.csv')
processed_stats = (stats_df['avg'])
 
plt.plot(processed_stats, color='black', label ='AF3 RMSF')
 
 
start=avg_rmsf[0]
for x in range(1,20):
    if x==7:
        start+=avg_rmsf[x][0:171]
    elif x==8:  
        start+=avg_rmsf[x][5:176]
    else:
        start+=avg_rmsf[x]
start=start/2
plt.plot(start,label='Mean RMSF')
ax = plt.gca()
# Add a secondary axis on the right
max_rmsf = 3.0  # adjust to match your RMSF plot scale
 
def rmsf_to_plddt(y):
    return 100 * (1 - y / max_rmsf)
 
def plddt_to_rmsf(y):
    return max_rmsf * (1 - y / 100)
 
secax = ax.secondary_yaxis('right', functions=(rmsf_to_plddt, plddt_to_rmsf))
secax.set_ylabel('pLDDT [%]')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.legend(frameon=False)
plt.xlabel('Residue Number')
plt.ylabel('Mean RMSF [Å]')
plt.ylim(0,2.8)
 
plt.savefig('si_rmsf.pdf', dpi=300)
 
r = np.corrcoef(start, processed_stats)[0, 1]
print("r =", r)
