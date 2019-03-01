#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/Cholesky>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include<cmath>

using namespace std;
using namespace Eigen;

class functions
{
	typedef OpenMesh::TriMesh_ArrayKernelT<> PGMesh;
	typedef Eigen::VectorXd VectorXd;
	typedef PGMesh::EdgeHandle EdgeHandle;
	typedef PGMesh::VertexHandle VertexHandle;

protected:
	struct faceanchor
	{
		int frameNo;
		int idx;
		Eigen::Matrix3d localframe;
	};
	struct vertexanchor
	{
		int frameNo;
		int idx;
		Eigen::VectorXd cord;
	};


public:
	/**初始化*/
	functions(PGMesh *_mesh, int meshNum);
	~functions();

	PGMesh *mesh_;//参考网格
	int nSeq;//序列的帧数
	int nV;//点的个数
	int nF;//面的个数
	int nE;//边的个数	
	
	Eigen::MatrixXd DiEdgeDataMatrix;//各个面各帧的二面角和边长
	vector<vector<double>> FaceAreas;//各个面各帧的面积
	//vector<vector<double>> de;//de[边][帧]，de=(1/3)*(a1+a2)，以边e的为公共边的两个面面积之和的三分之一
	vector<vector<Eigen::VectorXd>> SpatialTemData;//各个顶点各帧的坐标
	vector<vector<Eigen::Vector3d>> FaceNormals;//各个面各帧的法向量	
	vector<vector<Eigen::Matrix3d>> FaceLocalFrames;//各个面各帧的局部标架

	std::vector<unsigned int> edgeanchors_;//边锚点向量
	std::vector<unsigned int> Vertexanchors_;//顶点锚点向量
	std::vector<unsigned int> Faceanchors_;//面锚点向量
	int anchorPointidx;//锚点的index
	int anchorEidx;//锚点的邻边的index
	std::vector<std::vector<faceanchor>> faceAnchorFrames;//保存每帧面锚点的标架
	std::vector<std::vector<vertexanchor>> vertexAnchorCords;//保存每帧点锚点的坐标
	SimplicialCholesky<SparseMatrix<double>> presolver;
	Eigen::SparseMatrix<double> vertexCoff;

	VertexHandle POINT0_HANDLE(EdgeHandle edgeHandle) { return mesh_->from_vertex_handle(mesh_->halfedge_handle(edgeHandle, 0)); }
	VertexHandle POINT1_HANDLE(EdgeHandle edgeHandle) { return mesh_->to_vertex_handle(mesh_->halfedge_handle(edgeHandle, 0)); }
	VectorXd POINT(VertexHandle vertexhandle, int frame) { return SpatialTemData[vertexhandle.idx()][frame]; }
	VectorXd EDGE_POINT0(EdgeHandle edgeHandle, int frame) { return POINT(POINT0_HANDLE(edgeHandle), frame); }
	VectorXd EDGE_POINT1(EdgeHandle edgeHandle, int frame) { return POINT(POINT1_HANDLE(edgeHandle), frame); }
	VectorXd VECTOR(EdgeHandle edgeHandle, int frame) { return EDGE_POINT1(edgeHandle, frame) - EDGE_POINT0(edgeHandle, frame); }
	Eigen::Vector3d computeCross(Eigen::Vector3d x1, Eigen::Vector3d x2)//叉乘
	{
		Eigen::Vector3d x;
		x(0) = x1(1) * x2(2) - x1(2) * x2(1);
		x(1) = -x1(0) * x2(2) + x2(0) * x1(2);
		x(2) = x1(0) * x2(1) - x1(1) * x2(0);
		return x;
	}
	/** 设置锚点*/
	void functions::setAnchor(unsigned int _idx);
	/** 计算各个面各帧的面积以及每个面的法向量*/
	void calculateFacesAreas();
	/** 计算各个面各帧的法向量*/
	void calculateFacesNormals();
	/** 计算单帧单个面的法向量*/
	Eigen::Vector3d calculateFaceNormal(int faceNo, int frameNo);
	/** 计算单帧单个面的局部标架*/
	faceanchor calculateFaceFrame(int faceNo, int frameNo);
	/** 计算各个面各帧的局部标架*/
	void calculateFacesFrame();
	/** 计算各条边各帧的二面角和边长*/
	void computeDiEdge();
	/** 根据二面角和边长重构网格*/
	void functions::reconstructionFromDiEdges(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
	/** 根据点坐标导出obj网格*/
	void functions::outMesh(Eigen::VectorXd cordComponent, char* path);
	/** 根据点坐标导出off网格*/
	void functions::outMesh1(Eigen::VectorXd cordComponent, string path);
	/** 先解presolver*/
	void presolve();
};