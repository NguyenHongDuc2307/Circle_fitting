from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt



def set_axes_equal_3d(ax):
    limits = array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    spans = abs(limits[:,0] - limits[:,1])
    centers = mean(limits, axis=1)
    radius = 0.5 * max(spans)
    ax.set_xlim3d([centers[0]-radius, centers[0]+radius])
    ax.set_ylim3d([centers[1]-radius, centers[1]+radius])
    ax.set_zlim3d([centers[2]-radius, centers[2]+radius])


Fitcircle_list = []
Cluster_list = []
f1 = open("fitting_circle.txt", "r")
f2 = open("data.txt","r")

while f1:
	string = f1.readline()
	if string == "":
		break
	list1 = list(string.split())
	list2=list(map(float,list1))
	Fitcircle_list.append(list2)


while f2:
	string  = f2.readline()
	if string == "":
		break
	list1 = list(string.split())
	#print(list1)
	list2 = list(map(float,list1))
	Cluster_list.append(list2)


P_fitcircle = array(Fitcircle_list)
P = array(Cluster_list)

fig = figure(figsize=(15,15))
ax = fig.add_subplot(1,1,1,projection='3d')

#--- Cluster points
ax.plot(*P.T, ls='', marker='o', alpha=0.5, label='Cluster points P')


#--- Plot fitting circle
ax.plot(*P_fitcircle, color='k', ls='--', lw=2, label='Fitting circle')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()

ax.set_aspect('auto', 'datalim')
set_axes_equal_3d(ax)
plt.show()

f1.close() 
f2.close() 