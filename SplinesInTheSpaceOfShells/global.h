#pragma once
#include "functions.h"

typedef OpenMesh::TriMesh_ArrayKernelT<> PGMesh;


int nE_ = 0;
int nF_ = 0;
int nV_ = 0;
vector<vector<double>> EdgeLen1;//��֪mesh�ĸ�֡�����ߵı߳�
vector<vector<double>> DihAngel1;//��֪mesh�ĸ�֡�����ߵĶ����
vector<vector<double>> FaceAreas1;//��֪mesh�ĸ�֡����������

Eigen::MatrixXd DihEdg1;//��֪mesh�ĸ�֡�����ߵĶ���Ǻͱ߳�


vector<vector<double>> EdgeLen2;//�ο�֡
vector<vector<double>> DihAngel2;//�ο�֡
vector<vector<double>> FaceAreas2;//�ο�֡
Eigen::MatrixXd EdgeLen;//����mesh�ĸ�֡�����ߵı߳�
Eigen::MatrixXd DihAngel;//����mesh�ĸ�֡�����ߵĶ����
Eigen::MatrixXd FaceAreas;//����mesh�ĸ�֡����������


Eigen::MatrixXd DiEdgeDataMatrix2;//�ο�mesh�ĸ�֡�����ߵĶ���Ǻͱ߳�
Eigen::MatrixXd DiEdgeDataMatrix_all;//����mesh�ĸ�֡�����ߵĶ���Ǻͱ߳�
int meshNum = 5;//��֪�Ĺؼ�֡����
const int n_ = 4;//�м��n֡
int K = (meshNum - 1)*(n_ + 1);//P116ҳ��������Ĵ�С,������0-K��mesh����һ��K+1��
Eigen::SparseMatrix<double> linSysMatrice(K + 1, K + 1);//P116ҳ�ľ���