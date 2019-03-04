#pragma once
#include "functions.h"
//#pragma comment(lib,"tapkee.lib")

typedef OpenMesh::TriMesh_ArrayKernelT<> PGMesh;


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
int meshNum = 5;//��֪�Ĺؼ�֡����
const int n_ = 4;//�м��n֡
int K = (meshNum - 1)*(n_ + 1);//P116ҳ��������Ĵ�С,������0-K��mesh����һ��K+1��
Eigen::SparseMatrix<double> linSysMatrice(K + 1, K + 1);//P116ҳ�ľ���
