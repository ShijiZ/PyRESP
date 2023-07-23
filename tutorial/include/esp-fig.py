import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import seaborn as sns
sns.set()

resp_filename = "resp/2nd.esp"
pgm_ind_filename = "pgm-ind/2nd.esp"
pgm_perm_filename = "pgm-perm/2nd.esp"
fig_name = "esp-fig"

nesp = 7775
qm = np.ndarray((nesp))
resp = np.ndarray((nesp))
pgm_ind = np.ndarray((nesp))
pgm_perm = np.ndarray((nesp))

resp_file = open(resp_filename, 'r')
for i in range(3):
    line = resp_file.readline().split()
for i in range(nesp):
    line = resp_file.readline().split()
    qm[i] = float(line[3])
    resp[i] = float(line[4])
resp_file.close()

pgm_ind_file = open(pgm_ind_filename, 'r')
for i in range(3):
    line = pgm_ind_file.readline().split()
for i in range(nesp):
	line = pgm_ind_file.readline().split()
	pgm_ind[i] = float(line[4])
pgm_ind_file.close()

pgm_perm_file = open(pgm_perm_filename, 'r')
for i in range(3):
    line = pgm_perm_file.readline().split()
for i in range(nesp):
	line = pgm_perm_file.readline().split()
	pgm_perm[i] = float(line[4])
pgm_perm_file.close()

resp_R = stats.pearsonr(qm, resp)[0]
pgm_ind_R = stats.pearsonr(qm, pgm_ind)[0]
pgm_perm_R = stats.pearsonr(qm, pgm_perm)[0]
print(resp_R)
print(pgm_ind_R)
print(pgm_perm_R)

esp_list = [resp, pgm_ind, pgm_perm]
ylabel_list = ['RESP Model', 'pGM-ind Model', 'pGM-perm Model']
R_list = [resp_R, pgm_ind_R, pgm_perm_R]
color_list = ['#1f77b4', '#2ca02c', '#d62728']
plt.figure(figsize=(18, 5), dpi=350)
for i in range(1,4):
    plt.subplot(1,3,i)
    
    plt.scatter(qm, esp_list[i-1], alpha=0.5, color=color_list[i-1])
    plt.plot([], label="R = %f"%(R_list[i-1]))
    
    min_val = min(min(qm), min(esp_list[i-1]))
    max_val = max(max(qm), max(esp_list[i-1]))
    plt.xlim(min_val-0.006, max_val+0.006)
    plt.ylim(min_val-0.006, max_val+0.006)

    xpoints = ypoints = plt.xlim()
    plt.plot(xpoints, ypoints, linestyle='--', color='k', lw=2, scalex=False, scaley=False)
    
    plt.legend(handlelength=0, prop={'size':'large','weight':'bold'})
    plt.xlabel('QM ESP (a.u.)')
    plt.ylabel(ylabel_list[i-1]+' ESP (a.u.)')
    
plt.savefig(fig_name,bbox_inches='tight')