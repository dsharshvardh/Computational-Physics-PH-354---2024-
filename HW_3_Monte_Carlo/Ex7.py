
# Question 7(a) and 7(b)   

import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import os

def random_walk(N): # N = width of grid 
    random_walk_points=[] # Shows the complete path 
    random_walk_final_points=[[100,100]] # # Shows the final picture and here I am just adding a random point outside the grid because i need it for my code, in the end i will remove it
    f=0
    while random_walk_final_points[len(random_walk_final_points)-1]!=[0,0]:
        f=f+1
       
        anchored_points=[[0,0]]
        pos=[0,0] # first element is x and second element is y
        while True:
            k=random.randint(0,3) # 0 = left, 1 = right, 2 = up, 3 = down
            if k==0:
                newpos=[pos[0]-1,pos[1]]
                
                l=0
                for i in random_walk_final_points:
                    if i==newpos:
                        l=l+1
                    else:
                        continue 
                if newpos[0] > (N-1)/2 or newpos[1] > (N-1)/2 or newpos[0] < -(N-1)/2 or newpos[1] < -(N-1)/2 or l > 0:
                    break
                else:
                    pos = newpos
                


            elif k==1:
                newpos=[pos[0]+1,pos[1]]
                l=0
                for i in random_walk_final_points:
                    if i==newpos:
                        l=l+1
                    else:
                        continue 
                if newpos[0] > (N-1)/2 or newpos[1] > (N-1)/2 or newpos[0] < -(N-1)/2 or newpos[1] < -(N-1)/2 or l > 0:
                    break
                else:
                    pos = newpos

            
            elif k==2:
                newpos=[pos[0],pos[1]+1]
                l=0
                for i in random_walk_final_points:
                    if i==newpos:
                        l=l+1
                    else:
                        continue 
                if newpos[0] > (N-1)/2 or newpos[1] > (N-1)/2 or newpos[0] < -(N-1)/2 or newpos[1] < -(N-1)/2 or l > 0:
                    break
                else:
                    pos = newpos

            
            elif k==3:
                newpos=[pos[0],pos[1]-1]
                l=0
                for i in random_walk_final_points:
                    if i==newpos:
                        l=l+1
                    else:
                        continue 
                if newpos[0] > (N-1)/2 or newpos[1] > (N-1)/2 or newpos[0] < -(N-1)/2 or newpos[1] < -(N-1)/2 or l > 0:
                    break
                else:
                    pos = newpos
            anchored_points.append(pos)
        random_walk_final_points.append(pos)
        random_walk_points.append(anchored_points)
    random_walk_final_points.pop(0) # Removing the first point od this list as we took that initially 
    return random_walk_final_points,random_walk_points



random_walk_final_points_1,random_walk_points_1=random_walk(101)


'''
For overall picture we get in the end of the whole process. virdis colormap is used to show age of the particle
'''

x_coordinates = [i[0] for i in random_walk_final_points_1]
y_coordinates = [j[1] for j in random_walk_final_points_1]

indices = np.arange(len(random_walk_final_points_1))
plt.scatter(x_coordinates, y_coordinates, c=indices, marker='s', cmap='viridis', label='Points',s=15)

plt.xlim(-50, 50)
plt.ylim(-50, 50)

plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Diffusion limited aggregation')

cbar = plt.colorbar()
cbar.set_label('Particle age (0 means oldest particle)')

plt.show()




'''
For animation of anchored particles
'''

fig, ax = plt.subplots()
particles_scatter = ax.scatter([], [], color='red', marker='s', s=10, alpha=0.8)
ax.set_xlim(-50, 50)
ax.set_ylim(-50, 50)

def update(frame):
    particles = random_walk_final_points_1[:frame+1][::-1] 
    x_coordinates = [i[0] for i in particles]
    y_coordinates = [j[1] for j in particles]    
    particles_scatter.set_offsets(np.column_stack((x_coordinates, y_coordinates)))
    particles_scatter.set_array(np.array([1.0] * len(x_coordinates)))  
    
    return particles_scatter,

total_frames = len(random_walk_final_points_1)
ani = FuncAnimation(fig, update, frames=total_frames, blit=True, interval=5)

plt.show()




'''
For snapshots of path of the particle. I am showing this in smaller lattice such that particles will anchored earlier and you will not have to watch too many steps.
Note that sometimes you may see that the particle is touching the boundary but it is not anchored, it is because in the code if we particle is on the boundary and the next step
should take it outside the boundary then only it will anchor. Below commented off code is for this purpose.
'''


# random_walk_final_points_1,random_walk_points_1=random_walk(11)

# # Create a folder to save plots
# output_folder = 'HW_3_Fig_7_snapshots_of_trajectory_of_particles'
# os.makedirs(output_folder, exist_ok=True)

# particles = []
# steps = 0
# for i in range(len(random_walk_points_1)):
#     for j in range(len(random_walk_points_1[i])):
#         particles.append([random_walk_points_1[i][j][0], random_walk_points_1[i][j][1]])

#         x_coordinates = [p[0] for p in particles]
#         y_coordinates = [p[1] for p in particles]

#         plt.scatter(x_coordinates, y_coordinates, color='red', marker='s', label='Points', s=100)

#         plt.xlim(-5, 5)
#         plt.ylim(-5, 5)
#         filename = os.path.join(output_folder, f'step_{steps}.png')
#         plt.savefig(filename)

#         plt.close()  
#         if j != len(random_walk_points_1[i]) - 1:
#             particles = particles[:-1]  

#         steps += 1
#         if steps == 100:
#             break

#     if steps == 100:
#         break







# Question 7(c)


import math

def random_circle(r): # By this function I am taking a random point on circumference of the circle where r is the radius of the circle
    theta=random.uniform(0, 2*math.pi)

    # Now Converting polar coordinates to Cartesian coordinates
    x = r*math.cos(theta)
    y = r*math.sin(theta)
    return x, y

def point_out_circle(n): # This function will be used to take a lattice point out of the circle when it is not giving a lattice point on the circumference
    if n==int(n):
        r=n
    elif n>0:
        r=int(n)+1
    elif n<0:
        r=int(n)-1
    return r

def distance(pos): # for distance of elements of pos from the origin
    distance = math.sqrt(pos[0]**2+pos[1]**2)
    return distance


# Now let start making the code for our problem

def DLA(N): # N is the length of x axis and y axis of square lattice
    random_walk_points_1=[] # Shows the complete path 
    random_walk_final_points_1=[[0,0]]
    r=1
    while r<=(N-1)/4:
        x,y=random_circle(r)
        pos=[point_out_circle(x),point_out_circle(y)]
        point_path=[pos]
        
        while True:
            k=random.randint(0,3) # 0 = left, 1 = right, 2 = up, 3 = down
            if k==0:
                newpos=[pos[0]-1,pos[1]]                              

            elif k==1:
                newpos=[pos[0]+1,pos[1]]
               
            elif k==2:
                newpos=[pos[0],pos[1]+1]

            elif k==3:
                newpos=[pos[0],pos[1]-1]


            l = 0
            for i in random_walk_final_points_1:
                if i == newpos:
                    l = l+1

            if l > 0:
                break
            elif distance(newpos) > 2 * r:
                a, b = random_circle(r)
                pos = [point_out_circle(a), point_out_circle(b)]
                point_path = [pos]
            else:
                pos = newpos
                point_path.append(pos)
            

        random_walk_points_1.append(point_path)
        random_walk_final_points_1.append(pos)
        
        
        m = distance(pos)
        if m>r:
            r=m
        else:
            r=r    

    return random_walk_final_points_1




random_walk_final_points_1 = DLA(101)




import matplotlib.pyplot as plt
import numpy as np

x_coordinates = [i[0] for i in random_walk_final_points_1]
y_coordinates = [j[1] for j in random_walk_final_points_1]

indices = np.arange(len(random_walk_final_points_1))
plt.scatter(x_coordinates, y_coordinates, c=indices, marker='s', cmap='viridis', label='Points',s=10)
plt.xlim(-50, 50)
plt.ylim(-50, 50)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Diffusion limited aggregation')

cbar = plt.colorbar()
cbar.set_label('Particle age (0 means oldest particle)')

plt.show()