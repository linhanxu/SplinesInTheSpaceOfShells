#pragma once
#include "functions.h"

typedef OpenMesh::TriMesh_ArrayKernelT<> PGMesh;


int nE_ = 0;
int nF_ = 0;
int nV_ = 0;
vector<vector<double>> EdgeLen1;//已知mesh的各帧各条边的边长
vector<vector<double>> DihAngel1;//已知mesh的各帧各条边的二面角
vector<vector<double>> FaceAreas1;//已知mesh的各帧各个面的面积

Eigen::MatrixXd DihEdg1;//已知mesh的各帧各条边的二面角和边长


vector<vector<double>> EdgeLen2;//参考帧
vector<vector<double>> DihAngel2;//参考帧
vector<vector<double>> FaceAreas2;//参考帧
Eigen::MatrixXd EdgeLen;//所有mesh的各帧各条边的边长
Eigen::MatrixXd DihAngel;//所有mesh的各帧各条边的二面角
Eigen::MatrixXd FaceAreas;//所有mesh的各帧各个面的面积


Eigen::MatrixXd DiEdgeDataMatrix2;//参考mesh的各帧各条边的二面角和边长
Eigen::MatrixXd DiEdgeDataMatrix_all;//所有mesh的各帧各条边的二面角和边长
int meshNum = 5;//已知的关键帧个数
const int n_ = 4;//中间插n帧
int K = (meshNum - 1)*(n_ + 1);//P116页矩阵数组的大小,即共有0-K个mesh，即一共K+1个
Eigen::SparseMatrix<double> linSysMatrice(K + 1, K + 1);//P116页的矩阵
