#!/usr/bin/python

import readchar
import rospy
from geometry_msgs.msg import PoseStamped, Vector3Stamped, QuaternionStamped, Pose
from std_msgs.msg import Bool
# from relaxed_ik.msg import EEPoseGoals
# import RelaxedIK.Utils.transformations as T

rospy.init_node('keyboard_ikgoal_driver')

ik_goal_r_pub = rospy.Publisher('/ik_goal',PoseStamped,queue_size=5)
ik_goal_l_pub = rospy.Publisher('/ik_goal_l',PoseStamped,queue_size=5)
goal_pos_pub = rospy.Publisher('vive_position', Vector3Stamped)
goal_quat_pub = rospy.Publisher('vive_quaternion', QuaternionStamped)
# ee_pose_goals_pub = rospy.Publisher('/relaxed_ik/ee_pose_goals', EEPoseGoals, queue_size=5)
#quit_pub = rospy.Publisher('/relaxed_ik/quit',Bool,queue_size=5)

pos_stride = 0.015
rot_stride = 0.055

position_r = [0,0,0]
rotation_r = [1,0,0,0]

position_l = [0,0,0]
rotation_l = [1,0,0,0]

seq = 1
rate = rospy.Rate(1000)
while not rospy.is_shutdown():
    pose = PoseStamped()
    pose.pose.position.x = position_r[0]
    pose.pose.position.y = position_r[1]
    pose.pose.position.z = position_r[2]

    pose.pose.orientation.w = rotation_r[0]
    pose.pose.orientation.x = rotation_r[1]
    pose.pose.orientation.y = rotation_r[2]
    pose.pose.orientation.z = rotation_r[3]
    ik_goal_r_pub.publish(pose)

    pose = PoseStamped()
    pose.pose.position.x = position_l[0]
    pose.pose.position.y = position_l[1]
    pose.pose.position.z = position_l[2]

    pose.pose.orientation.w = rotation_l[0]
    pose.pose.orientation.x = rotation_l[1]
    pose.pose.orientation.y = rotation_l[2]
    pose.pose.orientation.z = rotation_l[3]
    ik_goal_l_pub.publish(pose)

    key = readchar.readkey()
    if key == 'w':
        position_r[0] += pos_stride
    elif key == 'x':
        position_r[0] -= pos_stride
    elif key == 'a':
        position_r[1] += pos_stride
    elif key == 'd':
        position_r[1] -= pos_stride
    elif key == 'q':
        position_r[2] += pos_stride
    elif key == 'z':
        position_r[2] -= pos_stride
    elif key == 'c':
        rospy.signal_shutdown()


    pose = PoseStamped()
    pose.pose.position.x = position_r[0]
    pose.pose.position.y = position_r[1]
    pose.pose.position.z = position_r[2]

    pose.pose.orientation.w = rotation_r[0]
    pose.pose.orientation.x = rotation_r[1]
    pose.pose.orientation.y = rotation_r[2]
    pose.pose.orientation.z = rotation_r[3]
    ik_goal_r_pub.publish(pose)


    pose = PoseStamped()
    pose.pose.position.x = position_l[0]
    pose.pose.position.y = position_l[1]
    pose.pose.position.z = position_l[2]

    pose.pose.orientation.w = rotation_l[0]
    pose.pose.orientation.x = rotation_l[1]
    pose.pose.orientation.y = rotation_l[2]
    pose.pose.orientation.z = rotation_l[3]
    ik_goal_l_pub.publish(pose)

    pos_goal = Vector3Stamped()
    pos_goal.vector.x = position_r[0]
    pos_goal.vector.y = position_r[1]
    pos_goal.vector.z = position_r[2]
    goal_pos_pub.publish(pos_goal)

    quat_goal = QuaternionStamped()
    quat_goal.quaternion.w = rotation_r[0]
    quat_goal.quaternion.x = rotation_r[1]
    quat_goal.quaternion.y = rotation_r[2]
    quat_goal.quaternion.z = rotation_r[3]
    goal_quat_pub.publish(quat_goal)

    
    pose_r = Pose()
    pose_r.position.x = position_r[0]
    pose_r.position.y = position_r[1]
    pose_r.position.z = position_r[2]

    pose_r.orientation.w = rotation_r[0]
    pose_r.orientation.x = rotation_r[1]
    pose_r.orientation.y = rotation_r[2]
    pose_r.orientation.z = rotation_r[3]

    pose_l = Pose()
    pose_l.position.x = position_l[0]
    pose_l.position.y = position_l[1]
    pose_l.position.z = position_l[2]

    pose_l.orientation.w = rotation_l[0]
    pose_l.orientation.x = rotation_l[1]
    pose_l.orientation.y = rotation_l[2]
    pose_l.orientation.z = rotation_l[3]
  
    seq += 1

