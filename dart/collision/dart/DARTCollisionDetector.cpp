/*
 * Copyright (c) 2013-2014, Georgia Tech Research Corporation
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

#include "dart/collision/dart/DARTCollisionDetector.h"

#include <vector>

#include "dart/dynamics/Shape.h"
#include "dart/dynamics/BodyNode.h"
#include "dart/dynamics/BoxShape.h"
#include "dart/dynamics/EllipsoidShape.h"
#include "dart/dynamics/CylinderShape.h"
#include "dart/collision/dart/DARTCollide.h"
#include "dart/collision/dart/NarrowPhaseAlgorithm.h"

namespace dart {
namespace collision {

//==============================================================================
DARTCollisionDetector::DARTCollisionDetector()
  : CollisionDetector()
{
}

//==============================================================================
DARTCollisionDetector::~DARTCollisionDetector()
{
}

//==============================================================================
CollisionNode* DARTCollisionDetector::createCollisionNode(
    dynamics::BodyNode* _bodyNode)
{
  return new CollisionNode(_bodyNode);
}

//==============================================================================
bool DARTCollisionDetector::detectCollision(bool /*_checkAllCollisions*/,
                                            bool /*_calculateContactPoints*/)
{
  clearAllContacts();

  // Set all the body nodes are not in colliding
  for (int i = 0; i < mCollisionNodes.size(); i++)
    mCollisionNodes[i]->getBodyNode()->setColliding(false);

  std::vector<Contact> contacts;

  for (int i = 0; i < mCollisionNodes.size(); i++)
  {
    for (int j = i + 1; j < mCollisionNodes.size(); j++)
    {
      CollisionNode* collNode1 = mCollisionNodes[i];
      CollisionNode* collNode2 = mCollisionNodes[j];
      dynamics::BodyNode* BodyNode1 = collNode1->getBodyNode();
      dynamics::BodyNode* BodyNode2 = collNode2->getBodyNode();

      if (!isCollidable(collNode1, collNode2))
        continue;

      for (int k = 0; k < BodyNode1->getNumCollisionShapes(); k++)
      {
        for (int l = 0; l < BodyNode2->getNumCollisionShapes(); l++)
        {
          int currContactNum = mContacts.size();

          // Select the narrow phase algorithm according to two collision shapes
          NarrowPhaseAlgorithm* narrowPhaseAlgorithm
              = selectNarrowPhaseAlgorithm(BodyNode1->getCollisionShape(k),
                                           BodyNode2->getCollisionShape(l));

          contacts.clear();

          // Debug code
          if (narrowPhaseAlgorithm == NULL)
            continue;

          //
          narrowPhaseAlgorithm->collide(
                BodyNode1->getCollisionShape(k),
                BodyNode1->getTransform()
                * BodyNode1->getCollisionShape(k)->getLocalTransform(),
                BodyNode2->getCollisionShape(l),
                BodyNode2->getTransform()
                * BodyNode2->getCollisionShape(l)->getLocalTransform(),
                &contacts);

          size_t numContacts = contacts.size();

          for (size_t m = 0; m < numContacts; ++m)
          {
            Contact contactPair;
            contactPair = contacts[m];
            contactPair.bodyNode1 = BodyNode1;
            contactPair.bodyNode2 = BodyNode2;
            assert(contactPair.bodyNode1 != NULL);
            assert(contactPair.bodyNode2 != NULL);

            mContacts.push_back(contactPair);
          }

          std::vector<bool> markForDeletion(numContacts, false);

          for (int m = 0; m < numContacts; m++)
          {
            for (int n = m + 1; n < numContacts; n++)
            {
              Eigen::Vector3d diff = mContacts[currContactNum + m].point
                                     - mContacts[currContactNum + n].point;

              if (diff.dot(diff) < 1e-6)
              {
                markForDeletion[m] = true;
                break;
              }
            }
          }

          for (int m = numContacts - 1; m >= 0; m--)
          {
            if (markForDeletion[m])
              mContacts.erase(mContacts.begin() + currContactNum + m);
          }
        }
      }
    }
  }

  for (size_t i = 0; i < mContacts.size(); ++i)
  {
    // Set these two bodies are in colliding
    mContacts[i].bodyNode1->setColliding(true);
    mContacts[i].bodyNode2->setColliding(true);
  }

  return !mContacts.empty();
}

//==============================================================================
bool DARTCollisionDetector::detectCollision(CollisionNode* _collNode1,
                                            CollisionNode* _collNode2,
                                            bool /*_calculateContactPoints*/)
{
  std::vector<Contact> contacts;
  dynamics::BodyNode* BodyNode1 = _collNode1->getBodyNode();
  dynamics::BodyNode* BodyNode2 = _collNode2->getBodyNode();

  for (int i = 0; i < BodyNode1->getNumCollisionShapes(); i++)
  {
    for (int j = 0; j < BodyNode2->getNumCollisionShapes(); j++)
    {
      collide(BodyNode1->getCollisionShape(i),
              BodyNode1->getTransform()
              * BodyNode1->getCollisionShape(i)->getLocalTransform(),
              BodyNode2->getCollisionShape(j),
              BodyNode2->getTransform()
              * BodyNode2->getCollisionShape(j)->getLocalTransform(),
              &contacts);
    }
  }

  return contacts.size() > 0 ? true : false;
}

//==============================================================================
NarrowPhaseAlgorithm* DARTCollisionDetector::selectNarrowPhaseAlgorithm(
    const dynamics::Shape* _shape1, const dynamics::Shape* _shape2)
{
  dynamics::Shape::ShapeType type1 = _shape1->getShapeType();
  dynamics::Shape::ShapeType type2 = _shape2->getShapeType();

  switch(type1)
  {
    case dynamics::Shape::BOX:
    {
      const dynamics::BoxShape* box0 = static_cast<const dynamics::BoxShape*>(_shape1);

      switch(type2)
      {
        case dynamics::Shape::BOX:
        {
          return &mBoxBoxAlgorithm;
        }
        case dynamics::Shape::ELLIPSOID:
        {
          return NULL;
        }
        case dynamics::Shape::CYLINDER:
        {
          //----------------------------------------------------------
          // NOT SUPPORT CYLINDER
          //----------------------------------------------------------
          return NULL;
        }
        default:
        {
          return NULL;
        }
      }

      break;
    }
    case dynamics::Shape::ELLIPSOID:
    {
      const dynamics::EllipsoidShape* ellipsoid0 = static_cast<const dynamics::EllipsoidShape*>(_shape1);

      switch(type2)
      {
        case dynamics::Shape::BOX:
        {
          return NULL;
        }
        case dynamics::Shape::ELLIPSOID:
        {
          return NULL;
        }
        case dynamics::Shape::CYLINDER:
        {
          //----------------------------------------------------------
          // NOT SUPPORT CYLINDER
          //----------------------------------------------------------
          return NULL;
        }
        default:
          return false;

          break;
      }

      break;
    }
    case dynamics::Shape::CYLINDER:
    {
      //----------------------------------------------------------
      // NOT SUPPORT CYLINDER
      //----------------------------------------------------------
      const dynamics::CylinderShape* cylinder0 = static_cast<const dynamics::CylinderShape*>(_shape1);

      Eigen::Vector3d dimTemp0(cylinder0->getRadius() * sqrt(2.0),
                               cylinder0->getRadius() * sqrt(2.0),
                               cylinder0->getHeight());
      switch(type2)
      {
        case dynamics::Shape::BOX:
        {
          return NULL;
        }
        case dynamics::Shape::ELLIPSOID:
        {
          return NULL;
        }
        case dynamics::Shape::CYLINDER:
        {
          //----------------------------------------------------------
          // NOT SUPPORT CYLINDER
          //----------------------------------------------------------
          return NULL;
        }
        default:
        {
          return false;
        }
      }

      break;
    }
    default:
      return false;

      break;
  }
}

}  // namespace collision
}  // namespace dart
