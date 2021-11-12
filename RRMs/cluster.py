import seaborn as sns 
#sns.set(color_codes=True)

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy

#qw=np.loadtxt('qw.data')
qw=np.loadtxt('qw_inter.data')

q=[]
for i in range(0,len(qw)):
    for j in range(i+1,len(qw)):
        q.append(1-qw[i][j])



#row_linkage = hierarchy.linkage(q, method='single')
row_linkage = hierarchy.linkage(q, method='average')
#row_linkage = hierarchy.linkage(q, method='average', optimal_ordering=True)
cluster = hierarchy.fcluster(row_linkage, t=0.4, criterion='distance')
np.savetxt('cluster.log', cluster, fmt="%d")

plt.title('Hierarchical Clustering Dendrogram (truncated)')
#plt.xlabel('cluster size')
plt.ylabel('distance')
#C = hierarchy.dendrogram(row_linkage, truncate_mode='lastp', p=200, leaf_rotation=90., leaf_font_size=8., show_contracted=True, color_threshold=0.4)
C = hierarchy.dendrogram(row_linkage, truncate_mode='lastp', p=100, no_labels=True, color_threshold=0.4)
plt.savefig('cluster_size.png')
plt.close()

T=sns.clustermap(qw,col_linkage=row_linkage,row_linkage=row_linkage)
plt.savefig('cluster_qw_inter.png')

index = T.dendrogram_col.dendrogram['leaves']
np.savetxt('col.dat',row_linkage, fmt='%f')
np.savetxt('cluster.index',index, fmt='%d')
plt.close()
