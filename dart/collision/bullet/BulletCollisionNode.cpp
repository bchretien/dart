/*
 * Copyright (c) 2014, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jeongseok Lee <jslee02@gmail.com>
 *
 * Georgia Tech Graphics Lab and Humanoid Robotics Lab
 *
 * Directed by Prof. C. Karen Liu and Prof. Mike Stilman
 * <karenliu@cc.gatech.edu> <mstilman@cc.gatech.edu>
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#include "dart/collision/bullet/BulletCollisionNode.h"

#include <iostream>

#include <assimp/scene.h>

#include "dart/collision/bullet/BulletTypes.h"
#include "dart/dynamics/BodyNode.h"
#include "dart/dynamics/BoxShape.h"
#include "dart/dynamics/EllipsoidShape.h"
#include "dart/dynamics/CylinderShape.h"
#include "dart/dynamics/PlaneShape.h"
#include "dart/dynamics/MeshShape.h"

namespace dart {
namespace collision {

//==============================================================================
BulletCollisionNode::BulletCollisionNode(dynamics::BodyNode* _bodyNode)
  : CollisionNode(_bodyNode)
{
  for (int i = 0; i < _bodyNode->getNumCollisionShapes(); i++)
  {
    dynamics::Shape* shape = _bodyNode->getCollisionShape(i);
    switch (shape->getShapeType())
    {
      case dynamics::Shape::BOX:
      {
        dynamics::BoxShape* box = static_cast<dynamics::BoxShape*>(shape);

        btBoxShape* btBox = new btBoxShape(btVector3(box->getSize()[0]*0.5,
                                                     box->getSize()[1]*0.5,
                                                     box->getSize()[2]*0.5));
        btBox->setMargin(0.0);
        btCollisionObject* btCollObj = new btCollisionObject();
        btCollObj->setCollisionShape(btBox);
        BulletUserData* userData = new BulletUserData;
        userData->bodyNode = _bodyNode;
        userData->shape = shape;
        userData->btCollNode = this;
        btCollObj->setUserPointer(userData);
        mbtCollsionObjects.push_back(btCollObj);

        break;
      }
      case dynamics::Shape::ELLIPSOID:
      {
        dynamics::EllipsoidShape* ellipsoid =
            static_cast<dynamics::EllipsoidShape*>(shape);

        if (ellipsoid->isSphere())
        {
          btSphereShape* btSphere = new btSphereShape(ellipsoid->getSize()[0] *
                                                      0.5);
          btSphere->setMargin(0.0);
          btCollisionObject* btCollObj = new btCollisionObject();
          btCollObj->setCollisionShape(btSphere);
          BulletUserData* userData = new BulletUserData;
          userData->bodyNode = _bodyNode;
          userData->shape = shape;
          userData->btCollNode = this;
          btCollObj->setUserPointer(userData);
          mbtCollsionObjects.push_back(btCollObj);
        }
        else
        {
          // TODO(JS): Add mesh for ellipsoid
        }

        break;
      }
      case dynamics::Shape::CYLINDER:
      {
        dynamics::CylinderShape* cylinder =
            static_cast<dynamics::CylinderShape*>(shape);

        btCylinderShapeZ* btCylinder =
            new btCylinderShapeZ(btVector3(cylinder->getRadius(),
                                           cylinder->getRadius(),
                                           cylinder->getHeight() * 0.5));
        btCylinder->setMargin(0.0);
        btCollisionObject* btCollObj = new btCollisionObject();
        btCollObj->setCollisionShape(btCylinder);
        BulletUserData* userData = new BulletUserData;
        userData->bodyNode = _bodyNode;
        userData->shape = shape;
        userData->btCollNode = this;
        btCollObj->setUserPointer(userData);
        mbtCollsionObjects.push_back(btCollObj);

        break;
      }
      case dynamics::Shape::PLANE:
      {
        dynamics::PlaneShape* plane =
            static_cast<dynamics::PlaneShape*>(shape);

        btScalar d = plane->getNormal().dot(plane->getPoint())
                   / plane->getNormal().squaredNorm();

        btStaticPlaneShape* btStaticPlane =
            new btStaticPlaneShape(convertVector3(plane->getNormal()), d);
        btStaticPlane->setMargin(0.0);
        btCollisionObject* btCollObj = new btCollisionObject();
        btCollObj->setCollisionShape(btStaticPlane);
        BulletUserData* userData = new BulletUserData;
        userData->bodyNode = _bodyNode;
        userData->shape = shape;
        userData->btCollNode = this;
        btCollObj->setUserPointer(userData);
        mbtCollsionObjects.push_back(btCollObj);

        break;
      }
      case dynamics::Shape::MESH:
      {
        dynamics::MeshShape* shapeMesh
            = static_cast<dynamics::MeshShape*>(shape);
        btConvexTriangleMeshShape* btMesh = _createMesh(shapeMesh->getScale(),
                                                        shapeMesh->getMesh());
        btMesh->setMargin(0.0);
        btCollisionObject* btCollObj = new btCollisionObject();

        // Add user data
        btCollObj->setCollisionShape(btMesh);
        BulletUserData* userData = new BulletUserData;
        userData->bodyNode   = _bodyNode;
        userData->shape      = shape;
        userData->btCollNode = this;
        btCollObj->setUserPointer(userData);

        //
        mbtCollsionObjects.push_back(btCollObj);

        break;
      }
      default:
      {
        std::cout << "ERROR: Collision checking does not support "
                  << _bodyNode->getName()
                  << "'s Shape type\n";
        break;
      }
    }
  }
}

//==============================================================================
BulletCollisionNode::~BulletCollisionNode()
{
}

//==============================================================================
void BulletCollisionNode::updateBulletCollisionObjects()
{
  for (int i = 0; i < mbtCollsionObjects.size(); ++i)
  {
    BulletUserData* userData =
        static_cast<BulletUserData*>(mbtCollsionObjects[i]->getUserPointer());
    dynamics::Shape* shape = userData->shape;
    btTransform T = convertTransform(mBodyNode->getTransform() *
                                     shape->getLocalTransform());
    mbtCollsionObjects[i]->setWorldTransform(T);
  }
}

//==============================================================================
int BulletCollisionNode::getNumBulletCollisionObjects() const
{
  return mbtCollsionObjects.size();
}

//==============================================================================
btCollisionObject*BulletCollisionNode::getBulletCollisionObject(int _i)
{
  return mbtCollsionObjects[_i];
}

//==============================================================================
btConvexTriangleMeshShape*_createMesh(const Eigen::Vector3d& _scale,
                                      const aiScene* _mesh)
{
  btTriangleMesh* btMesh = new btTriangleMesh();

  for (unsigned int i = 0; i < _mesh->mNumMeshes; i++)
  {
    for (unsigned int j = 0; j < _mesh->mMeshes[i]->mNumFaces; j++)
    {
      btVector3 vertices[3];
      for (unsigned int k = 0; k < 3; k++)
      {
        const aiVector3D& vertex = _mesh->mMeshes[i]->mVertices[
                                   _mesh->mMeshes[i]->mFaces[j].mIndices[k]];
        vertices[k] = btVector3(vertex.x * _scale[0],
                                vertex.y * _scale[1],
                                vertex.z * _scale[2]);
      }
      btMesh->addTriangle(vertices[0], vertices[1], vertices[2]);
    }
  }

  btConvexTriangleMeshShape* btMeshShape =
      new btConvexTriangleMeshShape(btMesh);

  return btMeshShape;
}

}  // namespace collision
}  // namespace dart
