# -- Import the required libraries --
import rosbag
from sys import stdout
import numpy as np
import transformations as tf


# -- Constants and Definitions --
POSITION = 'position'
ORIENTATION = 'orientation'
LINEAR_A = 'linear_a'
ANGULAR_V = 'angular_v'
LINEAR_T = 'linear_t'
ANGULAR_T = 'angular_t'
POSE_COV = 'pose_cov'
TWIST_COV = 'twist_cov'
ORI_COV = 'ori_cov'
TIME = 'time'
ODOM = '/odom'
IMU = '/imu/data'

ori_cov_size = [3, 3]
pose_cov_size = [6, 6]
twist_cov_size = [6, 6]


# -- Ros Database class --
class RosDB:
    # -- Constructor --
    def __init__(self, rosbag_path, enable_debug):
        self.rosbag_path = rosbag_path  # -- set rosbag path


    # -- Get Odometry values in a numpy matrix format --
    def get_odom_values(self):
        message_count = sum(1 for line in open(self.rosbag_path))
        # -- Initialize the matrix
        odometry = {
            POSITION: np.zeros([message_count, 3]),
            ORIENTATION: np.zeros([message_count, 3]),  # It's quaternion
            LINEAR_T: np.zeros([message_count, 3]),
            ANGULAR_T: np.zeros([message_count, 3]),  # It's quaternion
            TIME: np.zeros([message_count]),
            TWIST_COV: [],
            POSE_COV: [],
            ORI_COV: np.zeros([3, 3])}

        index = 0
        start_time = None
        # portion = (message_count / 10)

        # -- Get the data --
        with open(self.rosbag_path) as f:
            for line in f:
                ori = np.zeros(4)
                lt = np.zeros(3)
                at = np.zeros(3)

                pose_cov = np.zeros(36)
                ori_cov = np.zeros(9)
                twist_cov = np.zeros(36)

                ori[0], ori[1], ori[2], ori[3], lt[0], lt[1], lt[2], at[0], at[1], at[2], pose_cov[0], pose_cov[1], pose_cov[2], pose_cov[3], pose_cov[4], pose_cov[5], pose_cov[6], pose_cov[7], pose_cov[8], pose_cov[9], pose_cov[10], pose_cov[11], pose_cov[12], pose_cov[13], pose_cov[14], pose_cov[15], pose_cov[16], pose_cov[17], pose_cov[18], pose_cov[19], pose_cov[20], pose_cov[21], pose_cov[22], pose_cov[23],pose_cov[24], pose_cov[25], pose_cov[26], pose_cov[27], pose_cov[28],pose_cov[29], pose_cov[30], pose_cov[31], pose_cov[32], pose_cov[33], pose_cov[34], pose_cov[35], ori_cov[0], ori_cov[1], ori_cov[2], ori_cov[3], ori_cov[4], ori_cov[5], ori_cov[6], ori_cov[7], ori_cov[8], pose_cov[0], pose_cov[1], pose_cov[2], pose_cov[3], pose_cov[4], pose_cov[5], pose_cov[6], pose_cov[7], pose_cov[8], pose_cov[9], pose_cov[10], pose_cov[11],pose_cov[12], pose_cov[13], pose_cov[14], pose_cov[15],pose_cov[16], pose_cov[17], pose_cov[18], pose_cov[19], pose_cov[20], pose_cov[21], pose_cov[22], pose_cov[23], pose_cov[24], pose_cov[25], pose_cov[26], pose_cov[27], pose_cov[28], pose_cov[29], pose_cov[30], pose_cov[31], pose_cov[32], pose_cov[33], pose_cov[34], pose_cov[35], time= line.split(' ')

                # -- Convert quaternion to euler angles --
                quaternion = (ori[1], ori[2], ori[3], ori[0])


                euler = tf.euler_from_quaternion(quaternion)

                # if index < 10:
                #     print(euler)

                # -- Set initial T and get covariance
                if start_time is None:
                    start_time = int(time)
                    odometry[POSE_COV] = np.matrix(pose_cov).reshape(pose_cov_size)
                    odometry[TWIST_COV] = np.matrix(twist_cov).reshape(twist_cov_size)

                # -- Add values --
                # odometry[POSITION][index] = [pos[0], pos[1], pos[2]]
                odometry[ORIENTATION][index] = [euler[0], euler[1], euler[2]]
                odometry[LINEAR_T][index] = [lt[0], lt[1], lt[2]]
                odometry[ANGULAR_T][index] = [at[0], at[1], at[2]]
                odometry[TIME][index] = int(time) - start_time

                ####################################################################
                if start_time is None:
                    start_time = t
                    odometry[ORI_COV] = np.matrix(ori_cov).reshape(ori_cov_size)


                index += 1

        return odometry
