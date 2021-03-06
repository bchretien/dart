## DART 4.0

### Version 4.0.0 (2014-06-02)

1. Added implicit joint spring force and damping force
1. Added planar joint
1. Added soft body dynamics
1. Added computation of velocity and acceleration of COM
1. Added bullet collision detector
  * [Pull request #156](https://github.com/dartsim/dart/pull/156)
1. Improved performance of forward dynamics algorithm
  * [Pull request #188](https://github.com/dartsim/dart/pull/188)
1. Improved dynamics API for Skeleton and Joint
  * [Pull request #161](https://github.com/dartsim/dart/pull/161)
  * [Pull request #192](https://github.com/dartsim/dart/pull/192)
  * [Pull request #193](https://github.com/dartsim/dart/pull/193)
1. Improved constraint dynamics solver
  * [Pull request #184](https://github.com/dartsim/dart/pull/184)
1. Improved calculation of equations of motion using Featherstone algorithm
  * [Issue #85](https://github.com/dartsim/dart/issues/87)
1. Improved optimizer interface and added nlopt solver
  * [Pull request #152](https://github.com/dartsim/dart/pull/152)
1. Fixed self collision bug
  * [Issue #125](https://github.com/dartsim/dart/issues/125)
1. Fixed incorrect integration of BallJoint and FreeJoint
  * [Issue #122](https://github.com/dartsim/dart/issues/122)
  * [Pull request #168](https://github.com/dartsim/dart/pull/168)

### Version 3.0 (2013-11-04)

1. Removed Transformation classes. Their functionality is now included in joint classes.
1. Added Featherstone algorithm. Can currently only be used without collision handling. The old algortihm is still present and used for that case.
1. Removed kinematics namespace. Functionality is moved to dynamics classes.
1. Added dart root namespace
1. A lot of function and variable renames
1. Added constraint namespace
1. Added "common" namespace

### Version 2.6 (2013-09-07)

1. Clean-up of build system:
  * Renamed DART_INCLUDEDIR to the standard-compliant DART_INCLUDE_DIRS in CMake files. Users need to adapt their CMake files for this change.
  * Users no longer need to call find_package(DARTExt) in the CMake files. A call to find_package(DART) also finds its dependencies now.
  * Allow user to overwrite installation prefix
  * Add possibility to include DART header files as '#include \<dart/dynamics/Skeleton.h\>' in addition to '#include \<dynamics/Skeleton.h\>'
  * Allow out-of-source builds
1. URDF loader:
  * Major clean-up
  * Consider mesh scaling factor

### Version 2.5 (2013-07-16)

1. Replaced robotics::World with simulation::World
1. Removed robotics::Robot
1. Added simulation::SimWindow
1. Some speed-up of Eigen calculations
1. Added abstract trajectory interface
1. ConstraintDynamics handles contact, joint limit and other constraint forces simultaneously
1. Improved Lemke algorithm for solving LCP
1. Renamed skeletonDynamics::getQDotVector() to getPoseVelocity()
1. Added abstract CollisionDetector interface allowing for multiple different collision detector implementations.
1. Created math namespace
1. Added System class as base class to Skeleton and Joint
1. URDF loader: Removed ability to load nonstandard URDF files with an object tag
1. Added support for multiple shapes per BodyNode
1. Made urdfdom a dependency instead of including it in the DART source
1. Added function to CollisionDetector to let user check a specific pair of BodyNodes for collision

### Version 2.4 (2013-03-05)

1. Mass and inertia are no longer stored in Shape but in BodyNode.
1. Different shapes for collision and visualization (not just different meshes)
1. Shapes are no longer centered at the COM but can be transformed independently relative to the link frame.
1. Improved URDF support
  * Support for non-mesh shapes
  * Does not create dummy root nodes anymore
  * Support for continuous joints
  * Support for arbitrary joint axes for revolute joints (but not for prismatic joints) instead of only axis-aligned joint axes
  * Support for relative mesh paths even if the robot and world URDF files are in different directories
  * All supported joint types can be root joints
1. Clean-up of the Robot class
1. Removed Object class
1. More robust build and installation process on Linux

