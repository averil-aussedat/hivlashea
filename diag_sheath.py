#python diag_peSW.py
import h5py
import sys
import numpy as np
from pylab import * 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import colors
#0 0.27 0.34 0.43 0.63 1
# coarse graining
# 100

np.random.seed(19680801)
Nr = 2
Nc = 2
cmap = "cool"


#x1_mesh = h5py.File('cartesian_mesh-f-x1-128x128.h5','r')
#x2_mesh = h5py.File('cartesian_mesh-f-x2-128x128.h5','r')
#x1 = x1_mesh['x1'] 
#Nv = x1.shape[0]-1
#Nx = x1.shape[1]-1
#print("Nx=",Nx);
#print("Nv=",Nv);
#f = h5py.File('f0001-values.h5','r')
#f = np.transpose(f["values"])
#f = fi[7649:8738,:]
#print("f:",np.min(f),np.max(f))
#im = plt.imshow(f,cmap='jet',origin='lower',aspect=0.1875)
#im = plt.imshow(f,cmap='jet',origin='lower',aspect=1,interpolation='spline36')#2.823529412)
#plt.clim(0,0.4)
#plt.axis('off')
#plt.colorbar()


#fig, axs = plt.subplots(Nr, Nc,constrained_layout=False,figsize=(1.1*18,1.1*7.))
#fig, axs = plt.subplots(Nr, Nc,constrained_layout=False,figsize=(1.1*18,1.1*7.))

fig, axs = plt.subplots(Nr, Nc,constrained_layout=False,figsize=(15,7.5))

left = 0.06  # the left side of the subplots of the figure
right = 1.   # the right side of the subplots of the figure
bottom = 0. #1  # the bottom of the subplots of the figure
top = 0.97     # the top of the subplots of the figure
wspace = 0.01  # the amount of width reserved for space between subplots,
              # expressed as a fraction of the average axis width
hspace = 0.01  # the amount of height reserved for space between subplots,
              # expressed as a fraction of the average axis height



#fig.suptitle('f at time T=1000')
fig.subplots_adjust(hspace=0.01, wspace=0.01)

plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

#outer_grid = fig.add_gridspec(4, 4, wspace=0.0, hspace=0.0)
#fig.set_figheight(15)
#fig.set_figwidth(15)

images = []
aspect0=2. #0.5 #25
simu = [["f0e0001-values.h5" for i in range(Nr)] for j in range(Nc)]
#print(simu)

for i in range(Nr):
	for j in range(Nc):
		simu[i][j]="f0e0001-values.h5"

simu[0][0]="f0e0001-values.h5"
simu[1][0]="f0i0001-values.h5"
simu[0][1]="fe0001-values.h5"
simu[1][1]="fi0001-values.h5"



# simu[0][0]="fpeSW_N32_splT400.h5"
# simu[1][0]="fpeSW_N32_wenoT400.h5"
# simu[2][0]="fpeSW_N32_d2limT400.h5"
# simu[3][0]="fpeSW_N32_d2DaTeT400.h5"
# simu[4][0]="fpeSW_N32_d2nolimT400.h5"
# simu[5][0]="fpeSW_N32_d4limT400.h5"
# simu[6][0]="fpeSW_N32_d4DaTeT400.h5"
# simu[7][0]="fpeSW_N32_d4nolimT400.h5"

# simu[0][0]="fpeSW_N32_splT400.h5"
# simu[1][0]="fpeSW_N32_wenoT400.h5"
# simu[2][0]="fpeSW_N32_d2limT400.h5"
# simu[3][0]="fpeSW_N32_d2DaTeT400.h5"
# simu[4][0]="fpeSW_N32_d2nolimT400.h5"
# simu[5][0]="fpeSW_N32_d4limT400.h5"
# simu[6][0]="fpeSW_N32_d4DaTeT400.h5"
# simu[7][0]="fpeSW_N32_d4nolimT400.h5"
# 
# simu[0][1]="fpeSW_N64_splT400.h5"
# simu[1][1]="fpeSW_N64_wenoT400.h5"
# simu[2][1]="fpeSW_N64_d2limT400.h5"
# simu[3][1]="fpeSW_N64_d2DaTeT400.h5"
# simu[4][1]="fpeSW_N64_d2nolimT400.h5"
# simu[5][1]="fpeSW_N64_d4limT400.h5"
# simu[6][1]="fpeSW_N64_d4DaTeT400.h5"
# simu[7][1]="fpeSW_N64_d4nolimT400.h5"
# 
# simu[0][2]="fpeSW_N128_splT400.h5"
# simu[1][2]="fpeSW_N128_wenoT400.h5"
# simu[2][2]="fpeSW_N128_d2limT400.h5"
# simu[3][2]="fpeSW_N128_d2DaTeT400.h5"
# simu[4][2]="fpeSW_N128_d2nolimT400.h5"
# simu[5][2]="fpeSW_N128_d4limT400.h5"
# simu[6][2]="fpeSW_N128_d4DaTeT400.h5"
# simu[7][2]="fpeSW_N128_d4nolimT400.h5"
# 
# simu[0][3]="fpeSW_N256_splT400.h5"
# simu[1][3]="fpeSW_N256_wenoT400.h5"
# simu[2][3]="fpeSW_N256_d2limT400.h5"
# simu[3][3]="fpeSW_N256_d2DaTeT400.h5"
# simu[4][3]="fpeSW_N256_d2nolimT400.h5"
# simu[5][3]="fpeSW_N256_d4limT400.h5"
# simu[6][3]="fpeSW_N256_d4DaTeT400.h5"
# simu[7][3]="fpeSW_N256_d4nolimT400.h5"


for i in range(Nr):
	for j in range(Nc):	
		f = h5py.File(simu[i][j],'r')
		f = np.transpose(f["values"])
		#if(j==0):
		#	f = np.transpose(f["values"][60//4:(256-60)//4,75//4:182//4])
		#if(j==1):
		#	f = np.transpose(f["values"][60//2:(256-60)//2,75//2:182//2])
		#if(j==2):
		#	f = np.transpose(f["values"][60:(256-60),75:182])
		#if(j==3):
		#	f = np.transpose(f["values"][60*2:(256-60)*2,75*2:182*2])

		#print("fd1_lim_128:",np.min(f),np.max(f))
		images.append(axs[i, j].imshow(f,cmap='jet',origin='lower',aspect=aspect0,interpolation='spline36'))#2.823529412)
		#axs[i, j].label_outer()
		axs[i, j].axis('off')
		if(i==0):
			if(j==0):
				axs[i,j].set_title('init time')
			if(j==1):
				axs[i,j].set_title('final time')
		if(j==0):
			if(i==0):
				axs[i,j].text(-30,30,'electrons')
			if(i==1):
				axs[i,j].text(-30,30,'ions')		
# 			if(i==2):
# 				axs[i,j].text(-8,10,'d=2 \nlim')		
# 			if(i==3):
# 				axs[i,j].text(-8,10,'d=2 \nDaTe')		
# 			if(i==4):
# 				axs[i,j].text(-8,10,'d=2 \nno lim')		
# 			if(i==5):
# 				axs[i,j].text(-8,10,'d=4 \nlim')		
# 			if(i==6):
# 				axs[i,j].text(-8,10,'d=4 \nDaTe')		
# 			if(i==7):
# 				axs[i,j].text(-8,10,'d=4 \nno lim')		






#for i in range(Nr):
#    for j in range(Nc):
#        images.append(axs[i, j].imshow(f,cmap='jet',origin='lower',aspect=1,interpolation='spline36'))#2.823529412)
#        axs[i, j].label_outer()
#        axs[i, j].axis('off')


# Find the min and max of all colors for use in setting the color scale.
vmin = min(image.get_array().min() for image in images)
vmax = max(image.get_array().max() for image in images)
#norm = colors.Normalize(vmin=vmin, vmax=vmax)
norm = colors.Normalize(vmin=0., vmax=0.04) #0.15915494309189533576)
#(3/2)*exp(-1)*sqrt(2/Pi)
for im in images:
    im.set_norm(norm)

#fig.colorbar(images[0], ax=axs, orientation='horizontal', fraction=.1)


# Make images respond to changes in the norm of other images (e.g. via the
# "edit axis, curves and images parameters" GUI on Qt), but be careful not to
# recurse infinitely!
def update(changed_image):
    for im in images:
        if (changed_image.get_cmap() != im.get_cmap()
                or changed_image.get_clim() != im.get_clim()):
            im.set_cmap(changed_image.get_cmap())
            im.set_clim(changed_image.get_clim())


for im in images:
    im.callbacksSM.connect('changed', update)

plt.savefig('f.png')

plt.show()
sys.exit()
x1_mesh = h5py.File('cartesian_mesh-f-x1.h5','r')
x2_mesh = h5py.File('cartesian_mesh-f-x2.h5','r')
x1 = x1_mesh['x1'] 
Nv = x1.shape[0]-1
Nx = x1.shape[1]-1
print("Nx=",Nx);
print("Nv=",Nv);
f = h5py.File('f0001-values.h5','r')
f = np.transpose(f["values"])
#f = fi[7649:8738,:]
print("f:",np.min(f),np.max(f))
#im = plt.imshow(f,cmap='jet',origin='lower',aspect=0.1875)
im = plt.imshow(f,cmap='jet',origin='lower',aspect=1,interpolation='spline36')#2.823529412)
#plt.clim(0,0.4)
#plt.axis('off')
plt.colorbar()
#plt.pcolormesh(f)
pl.show()
sys.exit()
x2_mesh = h5py.File('cartesian_mesh-x2.h5','r')

#for f in x1_mesh:
#  print(f)

x1 = x1_mesh['x1'] 
x2 = x2_mesh['x2']

Nv = x1.shape[0]
Nx = x1.shape[1]-1

#print(x2[:,0])

#search x2_min
x2 = x2[:,0]
x2_diff = np.diff(x2[:])
min_val = np.min(x2_diff)
max_val = np.max(x2_diff)

print (min_val,max_val)

i=0
while(x2_diff[i]> (1+1e-4)*min_val):
  i=i+1     
#print(x2[i],x2[i+1]-x2[i],x2[i]-x2[i-1])
i_fine_min = i
while(x2_diff[i]< (1+1e-4)*min_val):
  i=i+1     
#print(x2[i],x2[i+1]-x2[i],x2[i]-x2[i-1])
i_fine_max = i

x2_fine = x2[i_fine_min:i_fine_max+1]
x2_fine_diff = np.diff(x2_fine[:])

if(np.max(x2_fine_diff)-np.min(x2_fine_diff)>1e-14):
  print(np.max(x2_fine_diff)-np.min(x2_fine_diff))
  print("bad detection of fine region")
  sys.exit()

#look for f_min and f_max
compute_f_bounds = False
#compute_f_bounds = True
f_min = 0.
f_max = 0.
if(compute_f_bounds):
  for i0 in range(0,2):
    for i1 in range(0,10):
      for i2 in range(0,10):
        for i3 in range(0,10):
          if((i0==0)and (i1==0) and(i2==0)and(i3<2)):
            continue
          f1 = h5py.File('f'+str(i0)+str(i1)+str(i2)+str(i3)+'-values.h5','r')
          f1 = f1["values"][i_fine_min:i_fine_max+1,:]
          f_min = min([np.min(f1),f_min])
          f_max = max([np.max(f1),f_max])
  print(f_min,f_max)
  sys.exit()

f_min = 0.
f_max = 0.3813095405431341
#deltaf_min = -0.21061378232176298
#deltaf_max = 0.21061378232176298
num_df = 100

         
f1 = h5py.File('f1600-values.h5','r')
f2 = h5py.File('f1601-values.h5','r')

f1 = f1["values"][i_fine_min:i_fine_max+1,:]
f2 = f2["values"][i_fine_min:i_fine_max+1,:]

Nv = f1.shape[0]
Nx = f1.shape[1]-1

#normalization to be between 0 and num_df
f1 = num_df*(f1-f_min)/(f_max-f_min)
f2 = num_df*(f2-f_min)/(f_max-f_min)

#first we do coarse graining
f1 = np.round(f1)
f2_store = f2
f2 = np.round(f2)


colors1 = np.zeros(num_df+1,dtype=integer)
for i in range(i_fine_max+1-i_fine_min):
  for j in range(Nx+1):
    colors1[int(f1[i,j])] += 1 
colors2 = np.zeros(num_df+1,dtype=integer)
for i in range(i_fine_max+1-i_fine_min):
  for j in range(Nx+1):
    colors2[int(f2[i,j])] += 1 

print(colors1)
print(colors2)

val = 0.5

val = int(val*num_df)
 
f1_val = np.zeros((i_fine_max+1-i_fine_min,Nx+1))
for i in range(i_fine_max+1-i_fine_min):
  for j in range(Nx+1):
    if int(f1[i,j])==val: 
      f1_val[i,j] = 1

f2_val = np.zeros((i_fine_max+1-i_fine_min,Nx+1))
for i in range(i_fine_max+1-i_fine_min):
  for j in range(Nx+1):
    if int(f2[i,j])==val: 
      f2_val[i,j] = 1

print(val,colors1[val],colors2[val])


#store all the elements for one specific val
index1_val = np.zeros((colors1[val],2),dtype=integer)
index2_val = np.zeros((colors2[val],2),dtype=integer)

s=0
for i in range(i_fine_max+1-i_fine_min):
  for j in range(Nx+1):
    if int(f1[i,j])==val: 
      index1_val[s,0] = i
      index1_val[s,1] = j
      s+=1
N1_val = s
#print(s)      
      
s=0
for i in range(i_fine_max+1-i_fine_min):
  for j in range(Nx+1):
    if int(f2[i,j])==val: 
      index2_val[s,0] = i
      index2_val[s,1] = j
      s+=1
N2_val = s
#print(s)
#print index1_val[:,0]
#print index1_val[:,1]
#print index2_val[:,0]
#print index2_val[:,1]

direction_val = np.zeros((colors1[val],2),dtype=integer)
exist_val = np.zeros(colors1[val],dtype=integer)

stencil1 = 3
stencil2 = 3

loc_val = np.zeros(((2*stencil1+1)*(2*stencil2+1),2),dtype=integer)


num=0

dx = x1_mesh['x1'][i_fine_min,1]-x1_mesh['x1'][i_fine_min,0]
dv = x2_mesh['x2'][i_fine_min+1,0]-x2_mesh['x2'][i_fine_min,0]

print('dx,dv=',dx,dv)

##we look for the nearest value of f for each point




for s in range(N1_val):
  i = index1_val[s,0]
  j = index1_val[s,1]
  s_loc = 0
  for ii in range(-stencil1,stencil1+1):
    for jj in range(-stencil2,stencil2+1):
      ii1 = (i+ii)
      if(ii1>i_fine_max-i_fine_min):
        ii1 = i_fine_max-i_fine_min
      ii2 = (j+jj)%Nx        
      if(int(f2[ii1,ii2])==val):
        exist_val[s] += 1
        loc_val[s_loc,0] = ii 
        loc_val[s_loc,1] = jj 
        s_loc +=1
  if(exist_val[s]>0):
    num+=1
    dist_min = (loc_val[0,1]*dx)**2+(loc_val[0,0]*dv)**2
    s_loc_min = 0
    for ss in range(s_loc):
      dist = (loc_val[ss,1]*dx)**2+(loc_val[ss,0]*dv)**2
      if(dist<dist_min):
        dist_min = dist
        s_loc_min = ss  
    direction_val[s,0] = loc_val[s_loc_min,0]
    direction_val[s,1] = loc_val[s_loc_min,1]
    

print("num=",num,max(exist_val))


# check = np.zeros(N1_val,dtype=integer)
# for s in range(N1_val):
#   che

ax = axes()
for s in range(N1_val):
  if(exist_val[s]>0):
    i = index1_val[s,0]
    j = index1_val[s,1]
    x0 = x1_mesh['x1'][i_fine_min+i,j] 
    v0 = x2_mesh['x2'][i_fine_min+i,j]
    x1 = x0+direction_val[s,1]*dx 
    v1 = v0+direction_val[s,0]*dv 
    #arrow(x0, v0, x1, v1) #, head_width=0.05, head_length=0.1, fc='k', ec='k')
    plot([x0,x1],[v0,v1])


show()
        
        
          
sys.exit()

X = x1_mesh['x1'][i_fine_min:i_fine_max+1,:]
V = x2_mesh['x2'][i_fine_min:i_fine_max+1,:]
subplot(211)
contourf(X,V,f1_val,interpolation='none')
subplot(212)
contourf(X,V,f2_val,interpolation='none')


show()


#fig = plt.figure()
#ax = Axes3D(fig)
#ax.set_xlabel('x')
#ax.set_ylabel('v')
#ax.plot_surface(T, OMEGA, W)
#ax.plot_surface(X, V, f1)
#contourf(X,V,f2)



#print(np.min(f1),np.max(f1))
#print(np.min(f2),np.max(f2))


