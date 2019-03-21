#pragma once
#include "functions.h"

typedef OpenMesh::TriMesh_ArrayKernelT<> PGMesh;

static const char MODEL[] = "man";
const int meshNum = 5;//��֪�Ĺؼ�֡����
const int n_ = 4;//�м��n֡

int nE_ = 0;
int nF_ = 0;
int nV_ = 0;


Eigen::MatrixXd DihEdg1;//��֪mesh�ĸ����߸�֡�Ķ���Ǻͱ߳�
Eigen::MatrixXd Dih1;//��֪mesh�ĸ����߸�֡�Ķ����
Eigen::MatrixXd Edg1;//��֪mesh�ĸ����߸�֡�ı߳�
Eigen::MatrixXd FaceAreas1;//��֪mesh�ĸ�֡����������
Eigen::MatrixXd Points1;//��֪mesh�ĸ������֡����ά����


vector<vector<double>> EdgeLen2;//�ο�֡
vector<vector<double>> DihAngel2;//�ο�֡
vector<vector<double>> FaceAreas2;//�ο�֡
Eigen::MatrixXd EdgeLen;//����mesh�ĸ�֡�����ߵı߳�
Eigen::MatrixXd DihAngel;//����mesh�ĸ�֡�����ߵĶ����
Eigen::MatrixXd FaceAreas;//����mesh�ĸ�֡����������


Eigen::MatrixXd DiEdgeDataMatrix2;//�ο�mesh�ĸ�֡�����ߵĶ���Ǻͱ߳�
Eigen::MatrixXd DiEdgeDataMatrix_all;//����mesh�ĸ�֡�����ߵĶ���Ǻͱ߳�
int K = (meshNum - 1)*(n_ + 1);//P116ҳ��������Ĵ�С,������0-K��mesh����һ��K+1��
Eigen::SparseMatrix<double> linSysMatrice(K + 1, K + 1);//P116ҳ�ľ���

Eigen::MatrixXd DihEdg_blend;//���Ի����֪mesh�ĸ����߸�֡�Ķ���Ǻͱ߳�
double w0=0,w1=0,w2=0,w3=0,w4=0;//���Ȩ��,��СΪmeshNum����ֵΪ0
