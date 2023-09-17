import matplotlib.pyplot as plt
import numpy as np
# from mpl_toolkits.mplot3d import Axes3D

message_count = sum(1 for line in open("./pos_out.txt"))

pos = np.zeros([message_count, 2])

index = 0

with open("./pos_out.txt") as f:
    for line in f:
        pos[index][0], pos[index][1] = line.split()
        index += 1

plt.plot(pos[:, 0], pos[:, 1], 'black')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Estimated path using ekf')
plt.legend("ekf")
plt.legend("odometry")
plt.grid()
plt.show()
