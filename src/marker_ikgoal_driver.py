#!/usr/bin/python
# 'Borrowed' from the relaxed IK github

import rospy
import os
from interactive_markers.interactive_marker_server import *
from visualization_msgs.msg import *
# from RelaxedIK.relaxedIK import get_relaxedIK_from_info_file
from interactive_marker_utils import InteractiveMarkerFeedbackUtil, InteractiveMarkerUtil, InteractiveMarkerServerUtil
# from relaxed_ik.msg import EEPoseGoals
from geometry_msgs.msg import PoseStamped, Vector3Stamped, QuaternionStamped, Pose
# import RelaxedIK.Utils.transformations as T

rospy.init_node('marker_ikgoal_driver')

# path_to_src = os.path.dirname(__file__)

# relaxedIK = get_relaxedIK_from_info_file(path_to_src)
# num_chains = relaxedIK.vars.robot.numChains

# init_ee_positions =  relaxedIK.vars.init_ee_positions
# init_ee_quats =  relaxedIK.vars.init_ee_quats
# print(init_ee_positions)
# print(init_ee_quats)
server = InteractiveMarkerServer("simple_marker")
ik_goal_pub = rospy.Publisher('/ik_goal',PoseStamped,queue_size=5)
rospy.sleep(0.2)

int_markers = []

int_marker = InteractiveMarkerUtil(init_pos=[0,0,0], init_quat=[0,0,0,1])
int_marker.add_6dof_controls()
int_markers.append(int_marker)

server.insert(int_marker.interactive_marker, int_marker.feedback_util.feedback_handler)

server.applyChanges()


rate = rospy.Rate(40)
while not rospy.is_shutdown():
    # ikgoal = PoseStamped()

    # for i in xrange(num_chains):
    i=0
    if not int_markers[i].feedback_util.active:
        pose = PoseStamped()
        pose.pose.position.x = 0.0
        pose.pose.position.y = 0.0
        pose.pose.position.z = 0.0

        pose.pose.orientation.w = 1.0
        pose.pose.orientation.x = 0.0
        pose.pose.orientation.y = 0.0
        pose.pose.orientation.z = 0.0
    else:
        pose = PoseStamped()

        pose.pose.position.x = int_markers[i].feedback_util.feedback.pose.position.x
        pose.pose.position.y = int_markers[i].feedback_util.feedback.pose.position.y 
        pose.pose.position.z = int_markers[i].feedback_util.feedback.pose.position.z

        pose.pose.orientation.w = int_markers[i].feedback_util.feedback.pose.orientation.w
        pose.pose.orientation.x = int_markers[i].feedback_util.feedback.pose.orientation.x
        pose.pose.orientation.y = int_markers[i].feedback_util.feedback.pose.orientation.y
        pose.pose.orientation.z = int_markers[i].feedback_util.feedback.pose.orientation.z

    # ikgoal.append(pose)

    ik_goal_pub.publish(pose)

    rate.sleep()
