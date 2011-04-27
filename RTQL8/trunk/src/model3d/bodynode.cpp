#include "primitive.h"
#include "joint.h"
#include "marker.h"
#include "dof.h"
#include "skeleton.h"

#include "utils/misc.h"
#include "utils/utils.h"
#include "utils/eigen_helper.h"

#include "bodynode.h"

extern VectorXd gravity;

namespace model3d {
  
  BodyNode::BodyNode(char *_name) {
    mJointOut.clear();
    mHandles.clear();

    mJointIn = NULL;
    mNodeIn = NULL;
    mPrimitive = NULL;

    mModelIndex = -1;

    if(_name == NULL) {
      strcpy(mName, "Node3D");
    } else {
      strcpy(mName, _name);
    }

    dependsOnDof = NULL;

    mOffset = VectorXd::Zero(3);
    mMass = 0;
  }

  BodyNode::~BodyNode() {
    for(size_t i = 0; i < mHandles.size(); ++i){
      delete mHandles[i];
    }
    mHandles.clear();

    if(mPrimitive != NULL) {
      delete mPrimitive;
    }
    mJointOut.clear();
    delete[] dependsOnDof;
    dependsOnDof = NULL;
  }

  void BodyNode::init() {
    mMass = mPrimitive->getMass();
    
    int nDofs = mSkel->getNumDofs();
    int localDofs = getNumDofs();

    mJWqQd = MatrixXd::Zero(3, nDofs);
    mJCqQd = MatrixXd::Zero(3, nDofs);

    T = MatrixXd::Identity(4, 4);
    Tq.resize(localDofs, MatrixXd::Identity(4, 4));
    Tqq.resize(localDofs, Tq);

    W = MatrixXd::Identity(4, 4);
    Wq.resize(localDofs, MatrixXd::Identity(4, 4));
    Wqq.resize(localDofs, Wq);

    mJC = MatrixXd::Zero(3, nDofs);
    mJW = MatrixXd::Zero(3, nDofs);
    mJCq.resize(nDofs, mJC);
    mJWq.resize(nDofs, mJW);
  }

  void BodyNode::update(VectorXd& state) {
    int nDofs = mSkel->getNumDofs();
    int localDofs = getNumDofs();

    T = getLocalTrans();
    for(int i=0; i<localDofs; i++){
      Tq.at(i) = getLocalDeriv(getDof(i));
    }

    if (mNodeIn) {
      W = mNodeIn->W*T;
      for (int i = 0; i < nDofs; ++i) {
        if (dependsOnDof[i]) {
          if (mNodeIn->dependsOnDof[i]) {
            Wq.at(i) = mNodeIn->Wq.at(i) * T;
          } else {
            Wq.at(i) = mNodeIn->W * Tq.at( mJointIn->getIndex(i) );
          }
        }
      }
    } else {

      W = T;
      for(int i = 0; i < localDofs; ++i) {
        int index1 = mJointIn->getDof(i)->getModelIndex();
        Wq.at(index1) = Tq.at(i);
      }
    }

    evalJC();
    evalJW();
	
    mCOM = utils::transform(W, mOffset);
  }

  void BodyNode::evalSecondOrder(const VectorXd& state) {
    int nDofs = mSkel->getNumDofs();
    int localDofs = getNumDofs();

    for (int i = 0; i < localDofs; ++i) {
      for (int j = 0; j < localDofs; ++j) {
        Tqq.at(i).at(j) = getLocalDeriv2(getDof(i),getDof(j));
      }
    }

    if (mNodeIn) {
      for (int i = 0; i < nDofs; ++i) {
        if (dependsOnDof[i]) {

          if (mNodeIn->dependsOnDof[i]) {
            for (int j = 0; j < nDofs; ++j) {
              if(dependsOnDof[j]) {
                if(mNodeIn->dependsOnDof[j]) {
                  Wqq.at(i).at(j) = mNodeIn->Wqq.at(i).at(j) * T;
                } else {
                  Wqq.at(i).at(j) = mNodeIn->Wq.at(i)*Tq.at(mJointIn->getIndex(j));
                }
              } // if(dependsOnDof[j]) {
            }
          } else { // if (mNodeIn->dependsOnDof[i]) {
            for(int j = 0; j < nDofs; ++j) {
              if (dependsOnDof[j]) {
                if (mNodeIn->dependsOnDof[j]) {
                  Wqq.at(i).at(j) = Wqq.at(j).at(i);
                } else {
                  Wqq.at(i).at(j) = mNodeIn->W*Tqq.at(mJointIn->getIndex(i))
                    .at(mJointIn->getIndex(j));
                }
              } // if (dependsOnDof[j]) {
            }
          } // if (mNodeIn->dependsOnDof[i]) {

        } // if (dependsOnDof[i]) {
      }

    } else { // if (mNodeIn) {
      for(int i = 0; i < localDofs; ++i) {
        int index1 = mJointIn->getDof(i)->getModelIndex();
        for(int j = 0; j < localDofs; ++j) {
          Wqq.at(index1).at(mJointIn->getDof(j)->getModelIndex()) = Tqq.at(i).at(j);
        }
      }
    }

    for(int i = 0; i < nDofs; ++i) {
      if(dependsOnDof[i]){
        evalJCq(i);
        evalJWq(i);
      }
    }	
  }

  MatrixXd BodyNode::evalM(VectorXd& state) {
    using eigenhelper::trans;

    Matrix3d inert = mPrimitive->getInertia();
    Matrix3d JM_w = trans(mJW) * inert;
    Matrix3d JM_c = mMass * trans(mJC);
    mK1 = JM_w * mJW + JM_c * mJC;

    return mK1;
  }

  VectorXd BodyNode::evalC(VectorXd& state) {
    using eigenhelper::trans;
    using eigenhelper::last;
    
    int nDofs = mSkel->getNumDofs();

    mJWqQd.setZero();
    mJCqQd.setZero();

    for(int i = 0; i < nDofs; ++i) {
      if(dependsOnDof[i]){
        mJWqQd += mJWq.at(i) * state[mSkel->getNumDofs() + i];
        mJCqQd += mJCq.at(i) * state[mSkel->getNumDofs() + i];
      }
    }
	
    Matrix3d inert = mPrimitive->getInertia();
    Matrix3d JM_w = trans(mJW) * inert;
    Matrix3d JM_c = mMass * trans(mJC);
    mK2Qd = JM_w*mJWqQd + JM_c*mJCqQd;

    VectorXd ret = VectorXd::Zero(mSkel->getNumDofs());
    ret = mK2Qd * last(state,mSkel->getNumDofs()) - mMass * gravity * mJC;
    return ret;

  }

  VectorXd BodyNode::evalWorldPos(VectorXd& lp) {
    VectorXd result = utils::transform(W,lp);
    return result;
  }

  MatrixXd BodyNode::evalDpDq(VectorXd& lp) {
    int nDofs = mSkel->getNumDofs();
    MatrixXd ret = MatrixXd::Zero(3, nDofs);
    
    for(int i = 0; i < nDofs; ++i) {
      if(dependsOnDof[i]) {
        ret.col(i) = utils::transform(Wq.at(i),lp);
      }
    }
    return ret;
  }

  void BodyNode::evalSecDpDq(VectorXd& lp, vector<MatrixXd>& H) {
    int nDofs = mSkel->getNumDofs();

    for (int i = 0; i < nDofs; ++i) {

      for (int j = 0; j <= i; ++j) {
        if (dependsOnDof[i] && dependsOnDof[j]) {
          H.at(i).col(j) = utils::transform(Wqq.at(i).at(j),lp);
          H.at(j).col(i) = H.at(i).col(j);
        }
      }

    }
    
  }

  VectorXd BodyNode::evalMomenta(VectorXd& state, VectorXd& COM) {
    using eigenhelper::first;
    using eigenhelper::last;
    using eigenhelper::sub;
    
    VectorXd L = VectorXd::Zero(6);
    VectorXd qDot = last(state, mSkel->getNumDofs());
	
    first(L, 3) = mMass * mJC * qDot;
	
    MatrixXd inertia = mPrimitive->getInertia();
    MatrixXd R = sub(W, 3, 3);
    last(L, 3) = (R * inertia) * (mJW * qDot);
    VectorXd r = evalCOM() - COM;
    last(L, 3) += utils::cross(r,first(L,3));
    return L;
  }

  MatrixXd BodyNode::evalP(VectorXd& COM) {
    using eigenhelper::sub;
    
    int nDofs = mSkel->getNumDofs();
    MatrixXd P = MatrixXd::Zero(6, nDofs);
	
    sub(P, 3, nDofs) = mMass*mJC;
	
    MatrixXd inertia = mPrimitive->getInertia();
    MatrixXd R = sub(W, 3, 3);
    VectorXd r = evalCOM() - COM;
    sub(P,3,0,3,nDofs) = (R*inertia)*mJW + mMass*utils::makeSkewSymmetric(r)*mJC;
	
    return P;
  }

  MatrixXd BodyNode::evaldLdq(VectorXd& state, VectorXd& COM) {
    using eigenhelper::last;
    using eigenhelper::sub;
    using eigenhelper::col;
    
    int nDofs = mSkel->getNumDofs();
    VectorXd qDot = last(state, nDofs);
    MatrixXd JCDot = MatrixXd::Zero(3, nDofs);
    MatrixXd JWDot = MatrixXd::Zero(3, nDofs);
	
    for (int i = 0; i < nDofs; ++i) {
      if(dependsOnDof[i]){
        JCDot += mJCq.at(i)*qDot(i);
        JWDot += mJWq.at(i)*qDot(i);
      }
    }

    MatrixXd dLdq = mMass * JCDot;
	
    MatrixXd R = sub(W, 3, 3);
    VectorXd omega = mJW * qDot;
    VectorXd v = mJC * qDot;
    MatrixXd inertia = mPrimitive->getInertia();
    VectorXd r = evalCOM() - COM;
    MatrixXd dPdq = MatrixXd::Zero(3, nDofs);
    
    for (int i = 0; i < nDofs; ++i) {
      if(dependsOnDof[i]){
        MatrixXd Rq = sub(Wq.at(i),3,3);
        VectorXd rq = col(mJC, i);
        col(dPdq,i) = (Rq*inertia)*omega + mMass*utils::makeSkewSymmetric(rq)*v;
      }
    }
    dPdq  += (R*inertia)*JWDot + mMass*utils::makeSkewSymmetric(r)*JCDot;

    MatrixXd ret = MatrixXd::Zero(6, nDofs);
    sub(ret,3,nDofs) = dLdq;
    sub(ret,3,0,3,nDofs) = dPdq;
    return ret;
  }

  double BodyNode::evalG()
  {
    VectorXd com = evalCOM();
    return mMass * gravity.dot(com); // dot(gravity, com);
  }

  VectorXd BodyNode::evaldGdq()
  {
    // return mMass*(gravity*mJC);
  }


  void BodyNode::evalJC() {
    mJC.setZero();
    int nDofs = mSkel->getNumDofs();
    for (int i = 0; i < nDofs; ++i) {
      if (dependsOnDof[i]) {
        // cout << "Wq = " << Wq.at(i).Rows() << " " << Wq.at(i).Cols() << endl;
        // cout << "mOffset = " << mOffset.Elts() << endl;
        VectorXd J = utils::transform(Wq.at(i), mOffset);
        mJC(0, i) = J[0];
        mJC(1, i) = J[1];
        mJC(2, i) = J[2];
      }
    }
  }

  void BodyNode::evalJCq(int dofIndex) {
    int nDofs = mSkel->getNumDofs();
    mJCq.at(dofIndex).setZero();
    for (int i = 0; i < nDofs; ++i) {
      if (dependsOnDof[i]) {
        VectorXd Jq = utils::transform(Wqq.at(dofIndex).at(i),mOffset);
        mJCq.at(dofIndex)(0,i) = Jq[0];
        mJCq.at(dofIndex)(1,i) = Jq[1];
        mJCq.at(dofIndex)(2,i) = Jq[2];
      }
    }
  }

  void BodyNode::evalJW() {
    int nDofs = mSkel->getNumDofs();
    mJW.setZero();
    
    for (int i = 0; i < nDofs; ++i) {
      if (dependsOnDof[i]) {
        MatrixXd transR = W.topLeftCorner(3, 3).transpose();
        MatrixXd dRdq = Wq.at(i).topLeftCorner(3, 3);
        MatrixXd omegaSkewSymmetric = dRdq * transR;
        VectorXd omega = utils::fromSkewSymmetric(omegaSkewSymmetric);
        omega = transR * omega;
			
        mJW(0, i) = omega[0];
        mJW(1, i) = omega[1];
        mJW(2, i) = omega[2];
      }
    }
  }

  void BodyNode::evalJWq(int dofIndex) {
    using eigenhelper::sub;
    using eigenhelper::trans;
    
    int nDofs = mSkel->getNumDofs();
    mJWq.at(dofIndex).setZero();

    MatrixXd transR = trans(sub(W, 3, 3));
    MatrixXd dRdq = sub(Wq.at(dofIndex), 3, 3);

    for(int i = 0; i < nDofs; ++i) {
      if(dependsOnDof[i]) {
        MatrixXd dRdqi = sub(Wq.at(i), 3, 3);
        MatrixXd dRdqq = sub(Wqq.at(dofIndex).at(i), 3, 3);
        MatrixXd omegaSkewSymmetric1 = dRdqi * transR;
        VectorXd omega1 = utils::fromSkewSymmetric(omegaSkewSymmetric1);
        MatrixXd omegaSkewSymmetric2 = dRdqq * transR + dRdqi*trans(dRdq);
        VectorXd omega2 = utils::fromSkewSymmetric(omegaSkewSymmetric2);
        VectorXd omega = (dRdq.transpose()) *omega1 + transR * omega2;
        mJWq.at(dofIndex)(0,i) = omega[0];
        mJWq.at(dofIndex)(1,i) = omega[1];
        mJWq.at(dofIndex)(2,i) = omega[2];
      }
    }
  }


  vector<Marker*> BodyNode::clearHandles() {
    vector<Marker*> ret = mHandles; 
    mHandles.clear(); 
    return ret;
  }

  void BodyNode::removeHandle(Marker *h) {
    for(int i = 0; i < mHandles.size(); ++i) {
      if(mHandles[i] == h) {
        mHandles.erase(mHandles.begin()+i);
        break;
      }
    }

  }

} // namespace model3d
