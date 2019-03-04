#pragma once
#include "global.h"
#include <tapkee/tapkee.hpp>

using namespace std;
using namespace tapkee;

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
	DihEdg1 = Eigen::MatrixXd::Zero(2 * nE_, meshNum);
	Points1 = Eigen::MatrixXd::Zero(3 * nV_, meshNum);
	Dih1 = Eigen::MatrixXd::Zero(nE_, meshNum);
	Edg1 = Eigen::MatrixXd::Zero(nE_, meshNum);
	FaceAreas1 = Eigen::MatrixXd::Zero(nF_, meshNum);

	cout << "load mesh..." << endl;
	cout << "calculate dihedral angle、edge length and face area..." << endl;

	//定义文件流，存二面角和边长到data/DihEdg文件夹	
	char *DihEdgDataPath = new char[100];
	sprintf_s(DihEdgDataPath, 100, "data\\DihEdg\\body_DihEdg.txt");
	fstream opDE(DihEdgDataPath, ios::out);
	//定义文件流，存顶点三维坐标到data/Points文件夹	
	char *PointsPath = new char[100];
	sprintf_s(PointsPath, 100, "data\\Points\\body_Points.txt");
	fstream opPoints(PointsPath, ios::out);
	//定义文件流，分别存二面角和边长到data/DihEdg文件夹	
	char *DihDataPath = new char[100];
	sprintf_s(DihDataPath, 100, "data\\DihEdg\\body_Dih.txt");
	fstream opD(DihDataPath, ios::out);
	char *EdgDataPath = new char[100];
	sprintf_s(EdgDataPath, 100, "data\\DihEdg\\body_Edg.txt");
	fstream opE(EdgDataPath, ios::out);
	/*//定义文件流，存面积到data/FacesAreas文件夹			
	char *FaceAreasPath = new char[100];
	sprintf_s(FaceAreasPath, 100, "data\\FaceAreas\\body_FaceAreas.txt");
	fstream opFA(FaceAreasPath, ios::out);*/

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
		DihEdg1.col(i) = calcu->DiEdgeDataMatrix.row(i).transpose();
		//存顶点三维坐标到data/Points文件夹			
		for (int j = 0,k=0; j <nV_ ,k<3 * nV_; j++)
		{
			Points1(k, i) = calcu->SpatialTemData[j][0][0];
			Points1(k+1, i) = calcu->SpatialTemData[j][0][1];
			Points1(k+2, i) = calcu->SpatialTemData[j][0][2];

			k = k + 3;
		}	
		//分别存二面角和边长到data/DihEdg文件夹			
		for (int j = 0; j<nE_; j++)
		{
			Dih1(j, i) = calcu->DiEdgeDataMatrix(0, 2 * j);
			Edg1(j, i) = calcu->DiEdgeDataMatrix(0, 2 * j + 1);
		}		
		/*//存面积到data/FacesAreas文件夹	
		calcu->calculateFacesAreas();
		for (int j = 0; j < calcu->nF; j++)
		{			
			FaceAreas1(j,i) = calcu->FaceAreas[j][i];
		}*/		

		delete calcu;
		delete pg;
	}
	opDE << DihEdg1 << endl;
	opDE.close();
	opPoints << Points1 << endl;
	opPoints.close();
	opD << Dih1 << endl;
	opD.close();
	opE << Edg1 << endl;
	opE.close();
	/*opFA << FaceAreas1;
	opFA.close();*/
}

//把2*nE_维降到2维
tapkee::DenseMatrix dimReduction()
{
	const int ROW = 2*nE_;
	const int COL = meshNum;
	tapkee::DenseMatrix testMatrix(ROW, COL);

	ifstream infile;//定义读取文件流，相对于程序来说是in
	infile.open("data\\DihEdg\\body_DihEdg.txt");//打开文件								   
	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			infile >> testMatrix(i, j);
		}
	}
	infile.close();

	TapkeeOutput output = tapkee::initialize()
		/*.withParameters((method = KernelLocallyLinearEmbedding,*/
		/*.withParameters((method = KernelLocalTangentSpaceAlignment,*/
		.withParameters((method = Isomap,
			num_neighbors = 4,
			target_dimension = 2))
		.embedUsing(testMatrix);
	return output.embedding.transpose();
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
				EdgeLen2[p*(n_ + 1) + q][j] = Edg1(j,p) + q*_n*(Edg1(j,p + 1) - Edg1(p,j));
				DihAngel2[p*(n_ + 1) + q][j] = Dih1(j, p) + q*_n*(Dih1(j, p + 1) - Dih1(p, j));
			}
		}	
	}
	for (int j = 0; j < nE_; j++)
	{
		EdgeLen2[(meshNum - 1)*(n_ + 1)][j] = Edg1(j, meshNum - 1);
		DihAngel2[(meshNum - 1)*(n_ + 1)][j] = Dih1(j, meshNum - 1);
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
				b_edge(m) = Edg1(j, q);
				b_dihAngel(m) = Dih1(j,q);
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
				b_faceAreas(m) = FaceAreas1(q,j);
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
	tapkee::DenseMatrix matrix_after_dimReduction = dimReduction();
	cout << matrix_after_dimReduction << endl;
	
		
	/*getlinSysMatrice(K + 1);*/

	/*dimReduction();*/
	/*solveLinearEquations();
	getReferData();
	getRecst();
	getAllData();*/

	cout << "请输入任意字符..." << endl;
	getchar();
	return 0;
}