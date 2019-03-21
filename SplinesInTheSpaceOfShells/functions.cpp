#pragma once
#include "functions.h"

/**初始化构造函数*/
functions::functions(PGMesh *_mesh, int meshNum)
{
	mesh_ = _mesh;
	nV = mesh_->n_vertices();
	nF = mesh_->n_faces();
	nE = mesh_->n_edges();
	nSeq = meshNum;

	SpatialTemData.resize(nV);
	for (int i = 0; i < nV; i++)
	{
		SpatialTemData[i].resize(nSeq);
		PGMesh::VertexHandle vh = _mesh->vertex_handle(i);
		PGMesh::Point &p = _mesh->point(vh);
		Eigen::VectorXd x(3);
		x << p[0], p[1], p[2];
		for (int j = 0; j < nSeq; j++)
		{
			SpatialTemData[i][j] = x;
		}
	}

	for (int i = 4817; i < 4820; i++)//模型man右手位置基本不动，所以可以把锚点设置在那里，每个模型都不一样
	{
		setAnchor(i);
	}
	faceAnchorFrames.resize(nSeq);
	vertexAnchorCords.resize(nSeq);

	for (int i = 0; i < nSeq; i++)
	{
		faceAnchorFrames[i].resize(Faceanchors_.size());
		vertexAnchorCords[i].resize(Vertexanchors_.size());
		for (int j = 0; j < Faceanchors_.size(); j++)
		{
			faceAnchorFrames[i][j] = calculateFaceFrame(Faceanchors_[j], i);
		}
		for (int j = 0; j < Vertexanchors_.size(); j++)
		{
			vertexAnchorCords[i][j].idx = Vertexanchors_[j];
			vertexAnchorCords[i][j].frameNo = i;
			vertexAnchorCords[i][j].cord = SpatialTemData[Vertexanchors_[j]][i];
		}
	}	
}
functions::~functions() {}

/** 设置锚点*/
void functions::setAnchor(unsigned int _idx)
{
	OpenMesh::VertexHandle vh = mesh_->vertex_handle(_idx);
	anchorPointidx = vh.idx();
	anchorEidx = mesh_->ve_begin(vh)->idx();
	int faceNo;
	faceNo = mesh_->vf_begin(vh)->idx();
	edgeanchors_.push_back(anchorEidx);
	Vertexanchors_.push_back(anchorPointidx);
	Faceanchors_.push_back(faceNo);
}
/** 计算单帧单个面的法向量*/
Eigen::Vector3d functions::calculateFaceNormal(int faceNo, int frameNo)
{
	PGMesh::FaceHandle fit = mesh_->face_handle(faceNo);
	PGMesh::VertexHandle v0;
	PGMesh::VertexHandle v1;
	PGMesh::VertexHandle v2;
	PGMesh::FVIter fvit = mesh_->fv_begin(fit);
	v0 = *fvit;
	++fvit;
	v1 = *fvit;
	++fvit;
	v2 = *fvit;
	Eigen::Vector3d e0 = POINT(v1, frameNo) - POINT(v0, frameNo);
	Eigen::Vector3d e1 = POINT(v2, frameNo) - POINT(v1, frameNo);
	Eigen::Vector3d Normal = computeCross(e0, e1);
	Normal = Normal / Normal.norm();
	return Normal;
}

/** 计算单帧单个面的标架*/
functions::faceanchor functions::calculateFaceFrame(int faceNo, int frameNo)
{
	faceanchor fanchor;
	PGMesh::FaceHandle fit = mesh_->face_handle(faceNo);
	PGMesh::VertexHandle v0;
	PGMesh::VertexHandle v1;
	PGMesh::VertexHandle v2;
	PGMesh::FVIter fvit = mesh_->fv_begin(fit);
	PGMesh::HalfedgeHandle he0;
	PGMesh::FHIter fheit = mesh_->fh_begin(fit);
	he0 = *fheit;
	v0 = mesh_->from_vertex_handle(he0);
	v1 = mesh_->to_vertex_handle(he0);
	//计算第一条边
	Eigen::Vector3d x_axis, y_axis, z_axis;
	x_axis = POINT(v1, frameNo) - POINT(v0, frameNo);
	x_axis = x_axis / x_axis.norm();
	z_axis = calculateFaceNormal(faceNo, frameNo);
	y_axis = computeCross(z_axis, x_axis);
	y_axis = y_axis / y_axis.norm();
	fanchor.localframe.col(0) = x_axis;
	fanchor.localframe.col(1) = y_axis;
	fanchor.localframe.col(2) = z_axis;
	fanchor.frameNo = frameNo;
	fanchor.idx = faceNo;
	return fanchor;
}

/** 计算各个面各帧的面积以及每个面的法向量*/
void functions::calculateFacesAreas()
{
	FaceAreas.resize(nF);
	PGMesh::ConstFaceIter cfIt = mesh_->faces_begin();
	PGMesh::ConstFaceIter cfItEnd = mesh_->faces_end();
	for (; cfIt != cfItEnd; ++cfIt)
	{
		FaceAreas[cfIt->idx()].resize(nSeq);
		for (int j = 0; j < nSeq; j++)
		{
			PGMesh::ConstFaceEdgeIter cfeIt = mesh_->cfe_iter(*cfIt);

			Eigen::VectorXd edgeVec1 = EDGE_POINT1(*cfeIt, j) - EDGE_POINT0(*cfeIt, j);
			++cfeIt;
			Eigen::VectorXd edgeVec2 = EDGE_POINT1(*cfeIt, j) - EDGE_POINT0(*cfeIt, j);
			double area = (computeCross(edgeVec1, edgeVec2)).norm() * 0.5L;
			Eigen::Vector3d facenormal = computeCross(edgeVec1, edgeVec2) / (computeCross(edgeVec1, edgeVec2)).norm();
			FaceAreas[cfIt->idx()][j] = area;
			FaceNormals[cfIt->idx()][j] = facenormal;
		}
	}
}

/** 计算各个面各帧的法向量*/
void functions::calculateFacesNormals()
{
	/*初始化FaceNormals*/
	FaceNormals.resize(nF);
	for (int i = 0; i < nF; i++)
	{
		FaceNormals[i].resize(nSeq);
	}
	/*遍历每个面的每一帧，计算出其面的法向量*/
	for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		int faceNo = fit->idx();
		PGMesh::VertexHandle v0;
		PGMesh::VertexHandle v1;
		PGMesh::VertexHandle v2;
		PGMesh::FVIter fvit = mesh_->fv_begin(*fit);
		v0 = *fvit;
		++fvit;
		v1 = *fvit;
		++fvit;
		v2 = *fvit;
		for (int frameNo = 0; frameNo < nSeq; frameNo++)
		{
			Eigen::Vector3d e0 = POINT(v1, frameNo) - POINT(v0, frameNo);
			Eigen::Vector3d e1 = POINT(v2, frameNo) - POINT(v1, frameNo);
			Eigen::Vector3d Normal = computeCross(e0, e1);
			Normal = Normal / Normal.norm();
			FaceNormals[faceNo][frameNo] = Normal;
		}
	}
}

/** 计算各个面各帧的局部标架*/
/*以列向量存储*/
void functions::calculateFacesFrame()
{
	/*初始化FaceLocalFrames*/
	FaceLocalFrames.resize(nF);
	for (int i = 0; i < nF; i++)
	{
		FaceLocalFrames[i].resize(nSeq);
	}

	/*遍历每个面的每一帧，得到面的局部标架*/
	for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		int faceNo = fit->idx();
		PGMesh::VertexHandle v0;
		PGMesh::VertexHandle v1;
		PGMesh::VertexHandle v2;
		PGMesh::FVIter fvit = mesh_->fv_begin(*fit);
		PGMesh::HalfedgeHandle he0;
		PGMesh::FHIter fheit = mesh_->fh_begin(*fit);
		he0 = *fheit;
		v0 = mesh_->from_vertex_handle(he0);
		v1 = mesh_->to_vertex_handle(he0);
		for (int frameNo = 0; frameNo < nSeq; frameNo++)
		{
			Eigen::Matrix3d faceFrame;
			Eigen::Vector3d x_axis, y_axis, z_axis;
			x_axis = POINT(v1, frameNo) - POINT(v0, frameNo);
			x_axis = x_axis / x_axis.norm();
			z_axis = FaceNormals[faceNo][frameNo];
			y_axis = computeCross(z_axis, x_axis);
			y_axis = y_axis / y_axis.norm();
			faceFrame.col(0) = x_axis;
			faceFrame.col(1) = y_axis;
			faceFrame.col(2) = z_axis;
			FaceLocalFrames[faceNo][frameNo] = faceFrame;
		}
	}
	//nSeq = SpatialTemData[0].size();
}

/** 计算各条边各帧的二面角和边长*/
void functions::computeDiEdge()
{	
	/*初始化DiEdgeData矩阵*/
	DiEdgeDataMatrix = Eigen::MatrixXd(nSeq, 2 * nE);	

	/*遍历每条边，计算二面角和边长*/
	for (PGMesh::EdgeIter eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++eit)
	{
		int edgeNo = eit->idx();
		PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle(*eit, 0);
		PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(*eit, 1);
		PGMesh::VertexHandle vj = mesh_->from_vertex_handle(heh0);
		PGMesh::VertexHandle vi = mesh_->to_vertex_handle(heh0);
		PGMesh::FaceHandle fh0 = mesh_->face_handle(heh0);
		PGMesh::FaceHandle fh1 = mesh_->face_handle(heh1);
		//cout << "经过边" << edgeNo << "的：" << endl<<"第一个面是："<<fh0<< endl<<"第二个面是：" << fh1 << endl<<endl;

		/*找出fh0中的三条半边的handle*/
		PGMesh::HalfedgeHandle he0;
		PGMesh::HalfedgeHandle he1;
		PGMesh::HalfedgeHandle he2;
		PGMesh::EdgeHandle e0;
		PGMesh::EdgeHandle e1;
		PGMesh::EdgeHandle e2;
		PGMesh::FHIter fheit = mesh_->fh_begin(fh0);
		he0 = *fheit;
		e0 = mesh_->edge_handle(he0);
		++fheit;
		he1 = *fheit;
		e1 = mesh_->edge_handle(he1);
		++fheit;
		he2 = *fheit;
		e2 = mesh_->edge_handle(he2);

		for (int frameNo = 0; frameNo < nSeq; frameNo++)
		{
			//cout << "帧" << frameNo << "以" << edgeNo << "为公共边的："<<endl<<"面"<<fh0.idx()<<"的面积是:" << FaceAreas[fh0.idx()][frameNo] << endl<<"面"<<fh1.idx()<<"的面积是:"  << FaceAreas[fh1.idx()][frameNo] << endl<<endl;
			//计算de
			//de[edgeNo][frameNo]=(1/3)*(FaceAreas[fh0.idx()][frameNo]+ FaceAreas[fh1.idx()][frameNo]);

			double length = VECTOR(*eit, frameNo).norm();//第eit条边，第frameNo帧的边长
														 /*若为边缘边，则标架无变化，因此两角度分别为0*/
			if (!(fh0.is_valid() && fh1.is_valid()))
			{
				Eigen::Matrix3d trans;
				trans << 1, 0, 0,
					0, 1, 0,
					0, 0, 1;
				Eigen::Vector2d Iic(0, length);
				DiEdgeDataMatrix.block<1, 2>(frameNo, 2 * edgeNo) = Iic.transpose();//从原矩阵中第(frameNo, 2 * edgeNo)开始，获取一个1行2列的子矩阵				
			}

			else
			{
				Eigen::Matrix3d frame1 = FaceLocalFrames[fh0.idx()][frameNo];
				Eigen::Matrix3d frame2 = FaceLocalFrames[fh1.idx()][frameNo];				
				/*先求出二面角，即phi*/
				double phi;
				Eigen::Vector3d normal1 = frame1.col(2);//提取矩阵frame1的第3列
				Eigen::Vector3d normal2 = frame2.col(2);
				Eigen::Matrix3d trans = frame1.transpose() * frame2;
				Eigen::Vector3d e_ij = POINT(vj, frameNo) - POINT(vi, frameNo);
				e_ij = e_ij / (e_ij.norm());
				Eigen::Vector3d crossnormal = computeCross(normal1, normal2);
				double dotcross = e_ij.dot(crossnormal);
				if (dotcross >= 0)
				{
					double dotNorm1Norm2 = normal1.dot(normal2);
					if (dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
					{
						phi = acos(dotNorm1Norm2);
					}
					else if (dotNorm1Norm2 >= 1 - 1e-6)
					{
						phi = 0;
					}
					else
					{
						phi = 3.1415926;
					}
				}
				else
				{
					double dotNorm1Norm2 = normal1.dot(normal2);
					if (dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
					{
						phi = -acos(dotNorm1Norm2);
					}
					else if (dotNorm1Norm2 >= 1 - 1e-6)
					{
						phi = 0;
					}
					else
					{
						phi = -3.1415926;
					}
				}
				Eigen::Vector2d Iic(phi, length);//第eit条边，第frameNo帧中的共享该边的二面角及该边边长						
				DiEdgeDataMatrix.block<1, 2>(frameNo, 2 * edgeNo) = Iic.transpose();
			}
		}
	}
}

/** 根据二面角和边长重构网格*/
void functions::reconstructionFromDiEdges(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp)
{

	if (edgeComponent.size() != 2 * nE)
		return;
	Eigen::SparseMatrix<double> Coefficients(9 * nE + Faceanchors_.size() * 9, 9 * nF);
	Coefficients.setZero();
	Eigen::VectorXd frameCor(9 * nE + Faceanchors_.size() * 9);
	frameCor.setZero();
	std::vector<Eigen::Triplet<double> > triple;//记录一个非零元素的行、列、值，填充一个稀疏矩阵只需要将所有表示非零元素的Triplet放在一个vector中传入即可。
	triple.clear();

	for (PGMesh::EdgeIter eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++eit)
	{
		if (mesh_->is_boundary(*eit))
			continue;//终止本次循环的执行，即跳过当前这次循环中continue语句后尚未执行的语句，接着进行下一次循环条件的判断。
		int edgeNo = eit->idx();
		PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle(*eit, 0);
		PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(*eit, 1);
		PGMesh::VertexHandle vj = mesh_->from_vertex_handle(heh0);
		PGMesh::VertexHandle vi = mesh_->to_vertex_handle(heh0);
		PGMesh::FaceHandle fh0 = mesh_->face_handle(heh0);
		PGMesh::FaceHandle fh1 = mesh_->face_handle(heh1);
		int fromFaceNo = fh0.idx();
		int toFaceNo = fh1.idx();		
		double theta = edgeComponent(2 * edgeNo + 0);

		/*找出fh0中的三条半边的handle*/
		PGMesh::HalfedgeHandle he0;
		PGMesh::HalfedgeHandle he1;
		PGMesh::HalfedgeHandle he2;
		PGMesh::EdgeHandle e0;
		PGMesh::EdgeHandle e1;
		PGMesh::EdgeHandle e2;
		PGMesh::FHIter fheit = mesh_->fh_begin(fh0);
		he0 = *fheit;
		e0 = mesh_->edge_handle(he0);
		++fheit;
		he1 = *fheit;
		e1 = mesh_->edge_handle(he1);
		++fheit;
		he2 = *fheit;
		e2 = mesh_->edge_handle(he2);
		
		//开始求axis
		Eigen::Vector3d axis;
		int judge;
		if (heh0 == he0)
		{
			axis << -1, 0, 0;
			judge = 1;
		}
		else if (heh0 == he1)
		{
			double a = edgeComponent(2 * e0.idx() + 1);
			double b = edgeComponent(2 * e1.idx() + 1);
			double c = edgeComponent(2 * e2.idx() + 1);
			/*double a = edgelengths [e0.idx()];
			double b = edgelengths [e1.idx()];
			double c = edgelengths [e2.idx()];*/

			double sqcos = (a * a + b * b - c * c) / (2 * a * b);//cosC＝(a²＋b²－c²)／2ab
			double squaresin = 1 - sqcos * sqcos;//sin²C
			axis << sqcos, -sqrt(squaresin), 0;//cosC,-sinC,0
			if (squaresin < 0)
				axis(1) = 0;
			if (sqcos > 1)
				axis(0) = 1;
			if (sqcos < -1)
				axis(0) = -1;
			judge = 2;
		}
		else
		{
			double a = edgeComponent(2 * e0.idx() + 1);
			double b = edgeComponent(2 * e2.idx() + 1);
			double c = edgeComponent(2 * e1.idx() + 1);
			/*double a = edgelengths [e0.idx()];
			double b = edgelengths [e2.idx()];
			double c = edgelengths [e1.idx()];*/

			double sqcos = (a * a + b * b - c * c) / (2 * a * b);
			double squaresin = 1 - sqcos * sqcos;
			axis << sqcos, sqrt(squaresin), 0;
			if (squaresin < 0)
				axis(1) = 0;
			if (sqcos > 1)
				axis(0) = 1;
			if (sqcos < -1)
				axis(0) = -1;
			judge = 3;
		}
		//结束求axis
		//开始求expedge
		Eigen::Matrix3d matrixedge;
		matrixedge << 0, -axis(2), axis(1),
			axis(2), 0, -axis(0),
			-axis(1), axis(0), 0;
		Eigen::Matrix3d expedge = (theta * matrixedge).exp();
		////结束求expedge

		//开始求matrixNormal（axis对应的第二个坐标轴）
		Eigen::Vector3d axisY = computeCross(Eigen::Vector3d(0, 0, 1), axis);
		Eigen::Matrix3d matrixNormal;
		matrixNormal << axis(0), axisY(0), 0,
			axis(1), axisY(1), 0,
			axis(2), axisY(2), 1;
		//结束求matrixNormal

		//开始求matrixNormal1（以axis作为坐标轴，求出第二个三角形的第一条边在该坐标轴下的坐标）
		/*找出fh1中的三条半边的handle*/
		PGMesh::HalfedgeHandle he3;
		PGMesh::HalfedgeHandle he4;
		PGMesh::HalfedgeHandle he5;
		PGMesh::EdgeHandle e3;
		PGMesh::EdgeHandle e4;
		PGMesh::EdgeHandle e5;
		fheit = mesh_->fh_begin(fh1);
		he3 = *fheit;
		e3 = mesh_->edge_handle(he3);
		++fheit;
		he4 = *fheit;
		e4 = mesh_->edge_handle(he4);
		++fheit;
		he5 = *fheit;
		e5 = mesh_->edge_handle(he5);
		//开始求axis1
		Eigen::Vector3d axis1;
		if (heh1 == he3)
		{
			axis1 << 1, 0, 0;
		}
		else if (heh1 == he4)
		{
			double a = edgeComponent(2 * e3.idx() + 1);
			double b = edgeComponent(2 * e4.idx() + 1);
			double c = edgeComponent(2 * e5.idx() + 1);
			double sqcos = (a * a + b * b - c * c) / (2 * a * b);
			double squaresin = 1 - sqcos * sqcos;
			axis1 << -sqcos, -sqrt(squaresin), 0;
			if (squaresin < 0)
				axis1(1) = 0;
			if (sqcos > 1)
				axis1(0) = -1;
			if (sqcos < -1)
				axis1(0) = 1;
		}
		else
		{
			double a = edgeComponent(2 * e3.idx() + 1);
			double b = edgeComponent(2 * e5.idx() + 1);
			double c = edgeComponent(2 * e4.idx() + 1);
			double sqcos = (a * a + b * b - c * c) / (2 * a * b);
			double squaresin = 1 - sqcos * sqcos;
			axis1 << -sqcos, sqrt(squaresin), 0;
			if (squaresin < 0)
				axis1(1) = 0;
			if (sqcos > 1)
				axis1(0) = -1;
			if (sqcos < -1)
				axis1(0) = 1;
		}
		//结束求axis1

		//开始求matrixNormal1
		Eigen::Vector3d axisY1 = computeCross(Eigen::Vector3d(0, 0, 1), axis1);
		Eigen::Matrix3d matrixNormal1;
		matrixNormal1 << axis1(0), axisY1(0), 0,
			axis1(1), axisY1(1), 0,
			axis1(2), axisY1(2), 1;
		//结束求matrixNormal1

		Eigen::Matrix3d rot = expedge * matrixNormal* matrixNormal1;//提取旋转分量

		if ((fh0.is_valid() && fh1.is_valid()))
		{			
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 0, rot(0, 0)));//triple.push_back(Triplet<double>(i, j, A[i][j]));//ij处的值是A[i][j],相同下标重复插入，则表示相加
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 1, rot(1, 0)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 2, rot(2, 0)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * toFaceNo + 0, -1));         
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 3, rot(0, 0)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 4, rot(1, 0)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 5, rot(2, 0)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * toFaceNo + 3, -1));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 6, rot(0, 0)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 7, rot(1, 0)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 8, rot(2, 0)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * toFaceNo + 6, -1));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 0, rot(0, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 1, rot(1, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 2, rot(2, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * toFaceNo + 1, -1));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 3, rot(0, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 4, rot(1, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 5, rot(2, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * toFaceNo + 4, -1));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 6, rot(0, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 7, rot(1, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 8, rot(2, 1)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * toFaceNo + 7, -1));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 0, rot(0, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 1, rot(1, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 2, rot(2, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * toFaceNo + 2, -1));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 3, rot(0, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 4, rot(1, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 5, rot(2, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * toFaceNo + 5, -1));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 6, rot(0, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 7, rot(1, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 8, rot(2, 2)));
			triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * toFaceNo + 8, -1));
		}
	}
	//设置锚点
	for (int i = 0; i < Faceanchors_.size(); i++)
	{
		Eigen::Matrix3d anchorFrame = faceAnchorFrames[anchorIthSeq][i].localframe;
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 0, 9 * Faceanchors_[i] + 0, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 1, 9 * Faceanchors_[i] + 1, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 2, 9 * Faceanchors_[i] + 2, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 3, 9 * Faceanchors_[i] + 3, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 4, 9 * Faceanchors_[i] + 4, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 5, 9 * Faceanchors_[i] + 5, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 6, 9 * Faceanchors_[i] + 6, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 7, 9 * Faceanchors_[i] + 7, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nE + i * 9 + 8, 9 * Faceanchors_[i] + 8, 1));
		frameCor(9 * nE + i * 9 + 0) = anchorFrame(0, 0);
		frameCor(9 * nE + i * 9 + 1) = anchorFrame(0, 1);
		frameCor(9 * nE + i * 9 + 2) = anchorFrame(0, 2);
		frameCor(9 * nE + i * 9 + 3) = anchorFrame(1, 0);
		frameCor(9 * nE + i * 9 + 4) = anchorFrame(1, 1);
		frameCor(9 * nE + i * 9 + 5) = anchorFrame(1, 2);
		frameCor(9 * nE + i * 9 + 6) = anchorFrame(2, 0);
		frameCor(9 * nE + i * 9 + 7) = anchorFrame(2, 1);
		frameCor(9 * nE + i * 9 + 8) = anchorFrame(2, 2);
	}
	Coefficients.setFromTriplets(triple.begin(), triple.end());//得到系数矩阵
	triple.clear();
	SimplicialCholesky<SparseMatrix<double>> solver;
	solver.compute(Coefficients.transpose() * Coefficients);
	if (solver.info() != Success)
	{		
		std::cout << "Decomposition failed" << std::endl;
		return;
	}
	Eigen::VectorXd x_up;
	x_up = solver.solve(Coefficients.transpose() * frameCor);
	
	/*给定所有面的标架，重构顶点坐标*/
	Eigen::VectorXd fvertex(9 * nF + 3 * vertexAnchorCords[0].size());
	Eigen::VectorXd y_up(3 * nV);

	for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		PGMesh::FaceHandle fh = *fit;
		int faceNo = fh.idx();
		/*找出fh0中的三条半边的handle*/
		PGMesh::HalfedgeHandle he0;
		PGMesh::HalfedgeHandle he1;
		PGMesh::HalfedgeHandle he2;
		PGMesh::EdgeHandle e0;
		PGMesh::EdgeHandle e1;
		PGMesh::EdgeHandle e2;
		PGMesh::FHIter fheit = mesh_->fh_begin(fh);
		he0 = *fheit;
		e0 = mesh_->edge_handle(he0);
		++fheit;
		he1 = *fheit;
		e1 = mesh_->edge_handle(he1);
		++fheit;
		he2 = *fheit;
		e2 = mesh_->edge_handle(he2);
		/*找出三条边的长度*/
		double a, b, c;
		a = edgeComponent(2 * e0.idx() + 1);
		b = edgeComponent(2 * e1.idx() + 1);
		c = edgeComponent(2 * e2.idx() + 1);		

		/*找出fh1中的三个点的handle*/
		PGMesh::VertexHandle v0;
		PGMesh::VertexHandle v1;
		PGMesh::VertexHandle v2;
		v0 = mesh_->from_vertex_handle(he0);
		v1 = mesh_->to_vertex_handle(he0);
		v2 = mesh_->to_vertex_handle(he1);

		/*提取出该面的标架中的a，b轴*/
		Eigen::Matrix3d frameFace;
		frameFace << x_up(9 * faceNo + 0), x_up(9 * faceNo + 1), x_up(9 * faceNo + 2),
			x_up(9 * faceNo + 3), x_up(9 * faceNo + 4), x_up(9 * faceNo + 5),
			x_up(9 * faceNo + 6), x_up(9 * faceNo + 7), x_up(9 * faceNo + 8);
		Eigen::JacobiSVD<MatrixXd> svd(frameFace, ComputeThinU | ComputeThinV);
		Eigen::MatrixXd u = svd.matrixU();
		Eigen::MatrixXd v = svd.matrixV();
		Eigen::Matrix3d U = u * v.transpose();
		Eigen::Vector3d a_axis, b_axis;
		a_axis = U.col(0);
		b_axis = U.col(1);
		a_axis = a_axis / a_axis.norm();
		b_axis = b_axis / b_axis.norm();

		/*提取出角度 alpha, beta*/
		double cosalpha, sinalpha, cosbeta, sinbeta;
		cosalpha = (a * a + c * c - b * b) / (2 * a * c);
		sinalpha = sqrt(1 - cosalpha * cosalpha);
		cosbeta = (a * a + b * b - c * c) / (2 * a * b);
		sinbeta = sqrt(1 - cosbeta * cosbeta);

		if (cosalpha > 1)
		{
			cosalpha = 1;
			sinalpha = 0;
		}
		if (cosalpha < -1)
		{
			cosalpha = -1;
			sinalpha = 0;
		}

		if (cosbeta > 1)
		{
			cosbeta = 1;
			sinbeta = 0;
		}
		if (cosbeta < -1)
		{
			cosbeta = -1;
			sinbeta = 0;
		}

		/*得到各边向量*/
		Eigen::Vector3d v0v1, v1v2, v2v0;
		v0v1 = a * a_axis;
		v1v2 = -b * (cosbeta * a_axis - sinbeta * b_axis);
		v2v0 = -c * (cosalpha * a_axis + sinalpha * b_axis);

		/*填入矩阵和向量元素*/
		fvertex(faceNo * 9 + 0) = v0v1(0);
		fvertex(faceNo * 9 + 1) = v0v1(1);
		fvertex(faceNo * 9 + 2) = v0v1(2);
		fvertex(faceNo * 9 + 3) = v1v2(0);
		fvertex(faceNo * 9 + 4) = v1v2(1);		
		fvertex(faceNo * 9 + 5) = v1v2(2);		
		fvertex(faceNo * 9 + 6) = v2v0(0);		
		fvertex(faceNo * 9 + 7) = v2v0(1);		
		fvertex(faceNo * 9 + 8) = v2v0(2);
	}

	/*添加锚点*/
	for (int j = 0; j < vertexAnchorCords[0].size(); j++)
	{		
		fvertex(9 * nF + 3 * j + 0) = vertexAnchorCords[anchorIthSeq][j].cord(0);
		fvertex(9 * nF + 3 * j + 1) = vertexAnchorCords[anchorIthSeq][j].cord(1);		
		fvertex(9 * nF + 3 * j + 2) = vertexAnchorCords[anchorIthSeq][j].cord(2);
	}	
	y_up = presolver.solve(vertexCoff.transpose() * fvertex);
	cordComponent = y_up;
}
/** 先解presolver*/
void functions::presolve()
{
	vertexCoff = Eigen::SparseMatrix<double>(9 * nF + 3 * vertexAnchorCords[0].size(), 3 * nV);
	std::vector<Eigen::Triplet<double> > triple;
	triple.clear();
	for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		PGMesh::FaceHandle fh = *fit;
		int faceNo = fh.idx();
		/*找出fh0中的三条半边的handle*/
		PGMesh::HalfedgeHandle he0;
		PGMesh::HalfedgeHandle he1;
		PGMesh::HalfedgeHandle he2;
		PGMesh::EdgeHandle e0;
		PGMesh::EdgeHandle e1;
		PGMesh::EdgeHandle e2;
		PGMesh::FHIter fheit = mesh_->fh_begin(fh);
		he0 = *fheit;
		e0 = mesh_->edge_handle(he0);
		++fheit;
		he1 = *fheit;
		e1 = mesh_->edge_handle(he1);
		++fheit;
		he2 = *fheit;
		e2 = mesh_->edge_handle(he2);
		PGMesh::VertexHandle v0;
		PGMesh::VertexHandle v1;
		PGMesh::VertexHandle v2;
		v0 = mesh_->from_vertex_handle(he0);
		v1 = mesh_->to_vertex_handle(he0);
		v2 = mesh_->to_vertex_handle(he1);

		/*填入矩阵和向量元素*/
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v1.idx() + 0, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v0.idx() + 0, -1));		
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v1.idx() + 1, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v0.idx() + 1, -1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v1.idx() + 2, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v0.idx() + 2, -1));	
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v2.idx() + 0, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v1.idx() + 0, -1));		
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v2.idx() + 1, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v1.idx() + 1, -1));	
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v2.idx() + 2, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v1.idx() + 2, -1));	
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v0.idx() + 0, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v2.idx() + 0, -1));	
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v0.idx() + 1, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v2.idx() + 1, -1));	
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v0.idx() + 2, 1));
		triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v2.idx() + 2, -1));	
	}
	/*添加锚点*/
	for (int j = 0; j < vertexAnchorCords[0].size(); j++)
	{
		triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 0, 3 * vertexAnchorCords[0][j].idx + 0, 1));	
		triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 1, 3 * vertexAnchorCords[0][j].idx + 1, 1));
		triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 2, 3 * vertexAnchorCords[0][j].idx + 2, 1));
		vertexCoff.setFromTriplets(triple.begin(), triple.end());
	}
	triple.clear();
	presolver.compute(vertexCoff.transpose() * vertexCoff);
	if (presolver.info() != Success)
	{
		std::cout << "Decomposition failed." << std::endl;
		return;
	}
}

/** 根据点坐标导出obj网格*/
void functions::outMesh(Eigen::VectorXd cordComponent, char* path)
{
	ofstream offfile(path);

	for (int i = 0; i< nV; i++)
	{
		offfile << "v" << " " << cordComponent(3 * i + 0) << " " << cordComponent(3 * i + 1) << " " << cordComponent(3 * i + 2) << endl;
	}
	for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		offfile << "f" << " ";
		for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++fvit)
		{
			offfile << fvit->idx() + 1 << " ";
		}
		offfile << endl;
	}
	offfile.close();
}

/** 根据点坐标导出off网格*/
void functions::outMesh1(Eigen::VectorXd cordComponent, string path)
{
	ofstream offfile(path);
	offfile << "OFF" << endl;
	offfile << nV << " " << nF << " " << 0 << endl;

	for (int i = 0; i< nV; i++)
	{
		offfile << cordComponent(3 * i + 0) << " " << cordComponent(3 * i + 1) << " " << cordComponent(3 * i + 2) << endl;
	}
	for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		offfile << "3" << " ";
		for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++fvit)
		{
			offfile << fvit->idx() << " ";
		}
		offfile << endl;
	}
	offfile.close();
}








