#pragma once
#include "global.h"
//#include "tapkee/tapkee.hpp"

using namespace std;
//using namespace tapkee;

void getNum_EFV()
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "inMesh\\body_0.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);

	functions *calcu = new functions(pg, meshNum);

	nE_ = calcu->nE;
	nF_ = calcu->nF;
	nV_ = calcu->nV;
}

void getDihEdgeArea()
{
	DihAngel1.resize(meshNum);
	EdgeLen1.resize(meshNum);
	FaceAreas1.resize(meshNum);

	/*1...DihEdg1 = Eigen::MatrixXd::Zero(nE_, 2);*/
	/*2...DihEdg1 = Eigen::MatrixXd::Zero(2*nE_, 1);*/
	DihEdg1 = Eigen::MatrixXd::Zero(2 * nE_, meshNum);

	cout << "load mesh..." << endl;
	cout << "calculate dihedral angle、edge length and face area..." << endl;

	char *DihEdgDataPath = new char[100];
	sprintf_s(DihEdgDataPath, 100, "data\\DihEdg\\body_DihEdg.txt");
	fstream opDE(DihEdgDataPath, ios::out);

	for (int i = 0; i < meshNum; i++)
	{
		char *meshpath = new char[100];
		sprintf_s(meshpath, 100, "inMesh\\body_%d.obj", i);
		PGMesh *pg = new PGMesh();
		OpenMesh::IO::read_mesh(*pg, meshpath);

		functions *calcu = new functions(pg, meshNum);
		calcu->calculateFacesNormals();
		calcu->calculateFacesFrame();
		calcu->computeDiEdge();

		//存二面角和边长到data/DihEdg文件夹
		char *DihEdgDataPath = new char[100];
		sprintf_s(DihEdgDataPath, 100, "data\\DihEdg\\body_DihEdg.txt");
		fstream opDE(DihEdgDataPath, ios::out);
		DihEdg1.col(i) = calcu->DiEdgeDataMatrix.row(i).transpose();
		opDE << DihEdg1 << endl;
		opDE.close();

		//DihAngel1[i].resize(nE_);
		//EdgeLen1[i].resize(nF_);
		//FaceAreas1[i].resize(nV_);

		////存面积到data/FacesAreas文件夹
		//calcu->calculateFacesAreas();
		//char *FaceAreasPath = new char[100];
		//sprintf_s(FaceAreasPath, 100, "data\\FaceAreas\\body_%d.txt", i);
		//fstream opFA(FaceAreasPath, ios::out);
		//for (int j = 0; j < calcu->nF; j++)
		//{
		//	opFA << calcu->FaceAreas[j][i] << endl;
		//	FaceAreas1[i][j] = calcu->FaceAreas[j][i];
		//}
		//opFA.close();

		////存二面角到data/DihedralAngle文件夹			
		//char *DihAngleDataPath = new char[100];
		//sprintf_s(DihAngleDataPath, 100, "data\\DihedralAngle\\body_%d.txt", i);
		//fstream opDA(DihAngleDataPath, ios::out);
		//
		////存边长到data/EdgeLength文件夹	
		//char *EdgeLengthDataPath = new char[100];
		//sprintf_s(EdgeLengthDataPath, 100, "data\\EdgeLength\\body_%d.txt", i);
		//fstream opEL(EdgeLengthDataPath, ios::out);

		//char *DihEdgDataPath = new char[100];
		//sprintf_s(DihEdgDataPath, 100, "data\\DihEdg\\body_%d.txt",i);
		//fstream opDE(DihEdgDataPath, ios::out);

		//int temp = 2 * (nE_);
		//for (int j = 0/*,k=0*/; j<temp/*,k<nE_*/; j++/*,k++*/)
		//{
		//	opDA << calcu->DiEdgeDataMatrix(0,j) << endl;
		//	DihAngel1[i][0.5*j] = calcu->DiEdgeDataMatrix(0,j);

		//	/*1...DihEdg1(k,0)= calcu->DiEdgeDataMatrix(0, k*2);
		//	DihEdg1(k,1) = calcu->DiEdgeDataMatrix(0, k*2+1);*/
		//	/*2...DihEdg1 = calcu->DiEdgeDataMatrix.row(0).transpose();*/				

		//	j++;
		//	opEL << calcu->DiEdgeDataMatrix(0,j) << endl;
		//	EdgeLen1[i][0.5*j - 0.5] = calcu->DiEdgeDataMatrix(0,j);

		//}	
		//opDA.close();
		//opEL.close();
		

		delete calcu;
		delete pg;
	}

}

//把2*nE_维降到2维
void dimReduction()
{
	const int ROW = 2*nE_;
	const int COL = meshNum;
	/*tapkee::DenseMatrix testMatrix(ROW, COL);

	TapkeeOutput output = initialize()
	.withParameters((method = KernelLocallyLinearEmbedding,
	target_dimension = 2))
	.embedUsing(matrix);

	cout << output.embedding.transpose() << endl;*/

}

void getReferData()//简化的样条无需此函数，非简化则需重写
{
	EdgeLen2.resize(K + 1);
	DihAngel2.resize(K + 1);
	double _n = 1 / (n_ + 1.0);

	for (int i = 0; i < K + 1; i++)
	{
		EdgeLen2[i].resize(nE_);
		DihAngel2[i].resize(nE_);
	}

	for(int p = 0; p < meshNum-1; p++)
	{
		for (int q = 0; q<n_ + 1; q++)
		{
			for (int j = 0; j < nE_; j++)
			{
				EdgeLen2[p*(n_ + 1) + q][j] = EdgeLen1[p][j] + q*_n*(EdgeLen1[p + 1][j] - EdgeLen1[p][j]);
				DihAngel2[p*(n_ + 1) + q][j] = DihAngel1[p][j] + q*_n*(DihAngel1[p + 1][j] - DihAngel1[p][j]);
			}
		}	
	}
	for (int j = 0; j < nE_; j++)
	{
		EdgeLen2[(meshNum - 1)*(n_ + 1)][j] = EdgeLen1[meshNum - 1][j];
		DihAngel2[(meshNum - 1)*(n_ + 1)][j] = DihAngel1[meshNum - 1][j];
	}

	DiEdgeDataMatrix2 = Eigen::MatrixXd(K + 1, 2 * nE_);
	for (int i = 0; i < K + 1; i++)
	{
		for (int j = 0; j < nE_; j++)
		{
			DiEdgeDataMatrix2(i, j * 2) = DihAngel2[i][j];
			DiEdgeDataMatrix2(i, j * 2 + 1) = EdgeLen2[i][j];
		}
	}

	cout << "get reference data ..." << endl;
	for (int i = 0; i < K + 1; i++)
	{
		//存参考二面角和边长分别到data/efer_DihedralAngle、data/refer_EdgeLength文件夹
		char *DihAngleDataPath = new char[100];
		sprintf_s(DihAngleDataPath, 100, "data\\refer_DihedralAngle\\refer_body_%d.txt", i);
		fstream opDA(DihAngleDataPath, ios::out);

		char *EdgeLengthDataPath = new char[100];
		sprintf_s(EdgeLengthDataPath, 100, "data\\refer_EdgeLength\\refer_body_%d.txt", i);
		fstream opEL(EdgeLengthDataPath, ios::out);
		
		for (int j = 0; j<nE_; j++)
		{
			opDA << DihAngel2[i][j]<< endl;		
			opEL << EdgeLen2[i][j] << endl;
		}
		opDA.close();
		opEL.close();
	}

}

void getlinSysMatrice(int n1)
{
	linSysMatrice.setZero();
	std::vector<Eigen::Triplet<double> > triple;

	triple.push_back(Eigen::Triplet<double>(0, 0, 1));
	triple.push_back(Eigen::Triplet<double>(0, 1, -2));
	triple.push_back(Eigen::Triplet<double>(0, 2, 1));

	triple.push_back(Eigen::Triplet<double>(1, 0, -2));
	triple.push_back(Eigen::Triplet<double>(1, 1, 5));
	triple.push_back(Eigen::Triplet<double>(1, 2, -4));
	triple.push_back(Eigen::Triplet<double>(1, 3, 1));

	for (int i = 2; i < n1 - 2; i++)
	{
		triple.push_back(Eigen::Triplet<double>(i, i - 2, 1));
		triple.push_back(Eigen::Triplet<double>(i, i - 1, -4));
		triple.push_back(Eigen::Triplet<double>(i, i, 6));
		triple.push_back(Eigen::Triplet<double>(i, i + 1, -4));
		triple.push_back(Eigen::Triplet<double>(i, i + 2, 1));
	}
	triple.push_back(Eigen::Triplet<double>(n1 - 2, n1 - 4, 1));
	triple.push_back(Eigen::Triplet<double>(n1 - 2, n1 - 3, -4));
	triple.push_back(Eigen::Triplet<double>(n1 - 2, n1 - 2, 5));
	triple.push_back(Eigen::Triplet<double>(n1 - 2, n1 - 1, -2));

	triple.push_back(Eigen::Triplet<double>(n1 - 1, n1 - 3, 1));
	triple.push_back(Eigen::Triplet<double>(n1 - 1, n1 - 2, -2));
	triple.push_back(Eigen::Triplet<double>(n1 - 1, n1 - 1, 1));

	linSysMatrice.setFromTriplets(triple.begin(), triple.end());
	triple.clear();
}
void solveLinearEquations()//需要重写
{
	vector<vector<double>> constZ;//存基础解系
	constZ.resize(2);
	constZ[0].resize(n_ + 2);
	constZ[1].resize(n_ + 2);
	for (int j = 0, m =n_+1; j <n_ + 2, m >= 0; j++, m--)
	{
		constZ[0][j] = m;
		constZ[1][j] = -m + 1;
	}

	/*cout<< "constZ[0][]"<<endl;
	for (int j = 0; j <n_ + 2; j++)
	{
		cout << constZ[0][j]<<endl;
		
	}
	cout << "constZ[1][]" << endl;
	for (int j = 0; j <n_ + 2; j++)
	{
		
		cout<<constZ[1][j] << endl;
	}*/

	//遍历每一边，求解线性方程组，得到所有mesh的边长EdgeLen
	EdgeLen = Eigen::MatrixXd::Zero(K + 1, nE_);
	DihAngel = Eigen::MatrixXd::Zero(K + 1, nE_);
	FaceAreas = Eigen::MatrixXd::Zero(K + 1, nF_);

	MatrixXd coefficientA = MatrixXd::Zero(2, constZ.size());
	for (int i = 0; i <2; i++)
	{
		coefficientA(i, 0) = constZ[0][i * (n_ + 1)];
		coefficientA(i, 1) = constZ[1][i * (n_ + 1)];
	}	
	
	for (int i=0;i<meshNum-1;i++) //若有3个keyPose,则需对每一边求解2组线性方程组
	{
		for (int j = 0; j < nE_; j++)
		{
			VectorXd constK_edge(constZ.size());
			VectorXd constK_dihAngel(constZ.size());
			VectorXd b_edge(2);
			VectorXd b_dihAngel(2);

			for (int m = 0,q=i; m <2; m++,q++)
			{
				b_edge(m) = EdgeLen1[q][j];
				b_dihAngel(m) = DihAngel1[q][j];				
			}
			constK_edge = coefficientA.fullPivLu().solve(b_edge);
			constK_dihAngel = coefficientA.fullPivLu().solve(b_dihAngel);
			//cout << "第" << j << "边的b_edge向量为:" << endl << b_edge << endl << endl << "常数K为:" << endl << constK_edge << endl << endl;
			//cout << "第" << j << "边的b_dihAngel向量为:" << endl << b_dihAngel << endl << endl << "常数K为:" << endl << constK_dihAngel << endl << endl;
			for (int m = i*(n_+1),q=0; m <= (i+1)*(n_ + 1),q<=(n_+1); m++,q++)
			{
				EdgeLen(m, j) = constK_edge(0)*constZ[0][q] + constK_edge(1)*constZ[1][q];
				DihAngel(m, j) = constK_dihAngel(0)*constZ[0][q] + constK_dihAngel(1)*constZ[1][q];
			}
		}
		
		for (int j = 0; j < nF_; j++)
		{
			VectorXd constK_faceAreas(constZ.size());
			VectorXd b_faceAreas(2);
			for (int m = 0, q = i; m <2; m++, q++)
			{
				b_faceAreas(m) = FaceAreas1[q][j];
			}
			constK_faceAreas = coefficientA.fullPivLu().solve(b_faceAreas);
			//cout << "第" << j << "边的b_faceAreas向量为:" << endl << b_faceAreas << endl << endl << "常数K为:" << endl << constK_faceAreas << endl << endl;
			for (int m = i*(n_ + 1), q = 0; m <= (i + 1)*(n_ + 1), q <= (n_ + 1); m++, q++)
			{
				FaceAreas(m, j) = constK_faceAreas(0)*constZ[0][q] + constK_faceAreas(1)*constZ[1][q];
			}
		}
	}

	DiEdgeDataMatrix_all = Eigen::MatrixXd(K + 1, 2 * nE_);
	for (int i = 0; i < K + 1; i++)
	{
		for (int j = 0; j < nE_; j++)
		{
			DiEdgeDataMatrix_all(i, j * 2) = DihAngel(i, j);
			DiEdgeDataMatrix_all(i, j * 2 + 1) = EdgeLen(i, j);
		}
	}
	
}

void getRecst()
{
	cout << endl << "getRecst..." << endl;		
	char *meshpath = new char[100];
	char *outMesh_refer = new char[100];
	char *outMesh_recst = new char[100];

	for (int m=0;m<meshNum-1;m++) 
	{
		//读取keyPose并new对象		
		sprintf_s(meshpath, 100, "inMesh\\body_%d.obj", m);
		PGMesh *pg = new PGMesh();
		OpenMesh::IO::read_mesh(*pg, meshpath);
		functions *calcu_key = new functions(pg, meshNum);

		//重构keyPose
		Eigen::VectorXd LK1(3 * nV_);
		sprintf_s(outMesh_recst, 100, "outMesh\\recst_body_%d.obj", m * (n_+1));
		Eigen::VectorXd LK0 = DiEdgeDataMatrix_all.row(m * (n_ + 1));
		calcu_key->presolve();
		calcu_key->reconstructionFromDiEdges(LK0, LK1, 0, 0.5, 0);
		calcu_key->outMesh(LK1, outMesh_recst);		
		
		for (int i = 1; i <= n_; i++)
		{
			//重构referPose
			Eigen::VectorXd LR1(3 * nV_);
			sprintf_s(outMesh_refer, 100, "outMeshRefer\\refer_body_%d.obj", m * (n_ + 1) + i);
			Eigen::VectorXd LR0 = DiEdgeDataMatrix2.row(m * (n_ + 1) + i);			
			calcu_key->reconstructionFromDiEdges(LR0, LR1, 0, 0.5, 0);
			calcu_key->outMesh(LR1, outMesh_refer);

			//读取referPose并new对象		
			PGMesh *pg = new PGMesh();
			OpenMesh::IO::read_mesh(*pg, outMesh_refer);
			functions *calcu_refer = new functions(pg, meshNum);

			//重构otherPose
			Eigen::VectorXd L1(3 * nV_);
			sprintf_s(outMesh_recst, 100, "outMesh\\recst_body_%d.obj", m * (n_ + 1) + i);
			Eigen::VectorXd L0 = DiEdgeDataMatrix_all.row(m * (n_ + 1) + i);
			calcu_refer->presolve();
			calcu_refer->reconstructionFromDiEdges(L0, L1, 0, 0.5, 0);
			calcu_refer->outMesh(L1, outMesh_recst);

			delete calcu_refer;
			delete pg;
		}

		delete calcu_key;
		delete pg;
	}
	//读取最后一个keyPose并new对象		
	sprintf_s(meshpath, 100, "inMesh\\body_%d.obj", (meshNum - 1));
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	functions *calcu = new functions(pg, meshNum);

	//重构最后一个keyPose
	Eigen::VectorXd L1(3 * calcu->nV);
	sprintf_s(outMesh_recst, 100, "outMesh\\recst_body_%d.obj", (meshNum-1) * (n_ + 1));
	Eigen::VectorXd L0 = DiEdgeDataMatrix_all.row((meshNum - 1) * (n_ + 1));
	calcu->presolve();
	calcu->reconstructionFromDiEdges(L0, L1, 0, 0.5, 0);
	calcu->outMesh(L1, outMesh_recst);

	delete calcu;
	delete pg;

}

void getAllData()
{	
	cout << "get all data from reconstruction..." << endl;
	for (int i = 0; i < K+1; i++)
	{
		//存面积到data/recst_FaceAreas文件夹
		char *FaceAreasPath = new char[100];
		sprintf_s(FaceAreasPath, 100, "data\\recst_FaceAreas\\recst_body_%d.txt", i);
		fstream opFA(FaceAreasPath, ios::out);
		for (int j = 0; j <nF_; j++)
		{
			opFA << FaceAreas(i,j)<< endl;			
		}
		opFA.close();

		//存二面角和边长分别到data/recst_DihedralAngle、data/recst_EdgeLength文件夹
		char *DihAngleDataPath = new char[100];
		sprintf_s(DihAngleDataPath, 100, "data\\recst_DihedralAngle\\recst_body_%d.txt", i);
		fstream opDA(DihAngleDataPath, ios::out);

		char *EdgeLengthDataPath = new char[100];
		sprintf_s(EdgeLengthDataPath, 100, "data\\recst_EdgeLength\\recst_body_%d.txt", i);
		fstream opEL(EdgeLengthDataPath, ios::out);

		/*int temp = 2 * (nE_);
		for (int j = 0; j<temp; j++)
		{
			opDA << DiEdgeDataMatrix_all(i,j) << endl;	
			j++;
			opEL << DiEdgeDataMatrix_all(i,j) << endl;			

		}*/
		for (int j = 0; j < nE_; j++)
		{
			opDA << DihAngel(i, j) << endl;
			opEL << EdgeLen(i, j) << endl;
		}
		opDA.close();
		opEL.close();
	}	
	
}


int main()
{
	getNum_EFV();
	getDihEdgeArea();	
	/*getlinSysMatrice(K + 1);*/

	/*dimReduction();*/
	/*solveLinearEquations();
	getReferData();
	getRecst();
	getAllData();*/


	//system("pause");
	cout << "请输入任意字符..." << endl;
	getchar();
	return 0;
}