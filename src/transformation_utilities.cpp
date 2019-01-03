#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include "transformation_utilities.hpp"

///////////////////////////////////////////////////////////

Eigen::Matrix4d rtf::hom_T(Eigen::Vector3d t, Eigen::Matrix3d r)
{
    Eigen::Matrix4d T;
    T.block(0,3,3,1) << t;
    T.block(0,0,3,3) << r;
    T.block(3,0,1,4) << 0,0,0,1;
    return T;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::apply_transformation(Eigen::MatrixXd data, Eigen::Matrix4d T_mat)
{
    //NOTE: Homogeneous Tranformation Matrix (4x4)

    // putting data in [x, y, z, 1]' format
    Eigen::MatrixXd data_with_fourth_row(data.cols()+1,data.rows());
    Eigen::VectorXd ones_vec = Eigen::VectorXd::Constant(data.rows(),1);
    data_with_fourth_row.block(0,0,data.cols(),data.rows()) = data.transpose();
    data_with_fourth_row.block(data.cols(),0,1,data.rows()) = ones_vec.transpose();
    Eigen::MatrixXd transformed_data = T_mat*data_with_fourth_row;
    Eigen::MatrixXd transformed_data_mat(transformed_data.rows()-1,transformed_data.cols());
    transformed_data_mat = transformed_data.block(0,0,transformed_data.rows()-1,transformed_data.cols());
    return transformed_data_mat.transpose();
}

////////////////////////////////////////////////////////////

std::string rtf::validate_seq(std::string seq)
{
	if(seq =="")
		seq = "ZYX";	
	bool invalid_flag = false;
	if(seq.size()>0 && seq.size()<3)
	{
		invalid_flag = true;
	}
	for (int i =0;i<3;++i)
		if(seq[i]!='X' && seq[i]!='Y' && seq[i]!='Z')
		{
			invalid_flag = true; 
			break;
		}
	if(invalid_flag)
	{
		std::cerr << "ERROR: Invalid Rotations Sequence: " << seq << std::endl;
		std::terminate();		
	}
	return seq;
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::eul2rot(Eigen::MatrixXd eul_angles, std::string seq)
{
	seq = rtf::validate_seq(seq);
	Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity();
	for (int i=0; i<3; ++i)
	{
		if(seq[i]=='X')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitX());
		else if(seq[i]=='Y')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitY());			
		else if(seq[i]=='Z')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitZ());					
	}
	return rot_mat; 
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::rot2eul(Eigen::Matrix3d rot_mat, std::string seq)
{
	seq = rtf::validate_seq(seq);
	int rot_idx[3];
	for (int i=0; i<3; ++i)
	{
		if(seq[i]=='X')
			rot_idx[i] = 0;
		else if(seq[i]=='Y')
			rot_idx[i] = 1;
		else if(seq[i]=='Z')
			rot_idx[i] = 2;
	}	
	Eigen::MatrixXd eul_angles(1,3);
	Eigen::Vector3d eul_angles_vec;
	eul_angles_vec = rot_mat.eulerAngles(rot_idx[0], rot_idx[1], rot_idx[2]);
	eul_angles(0,0) = eul_angles_vec[0];
	eul_angles(0,1) = eul_angles_vec[1];
	eul_angles(0,2) = eul_angles_vec[2];
	return eul_angles;
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::qt2rot(Eigen::MatrixXd quat)
{
	Eigen::Quaterniond q;
	q.x() = quat(0,0);
	q.y() = quat(0,1);
	q.z() = quat(0,2);
	q.w() = quat(0,3);	
	return q.normalized().toRotationMatrix();
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::rot2qt(Eigen::Matrix3d rot_mat)
{
	Eigen::MatrixXd quat(1,4);
	Eigen::Quaterniond q(rot_mat);
	quat(0,0) = q.x();
	quat(0,1) = q.y();
	quat(0,2) = q.z();
	quat(0,3) = q.w();
	return quat;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::eul2qt(Eigen::MatrixXd eul_angles,std::string seq)
{
	seq = rtf::validate_seq(seq);
	Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity();
	for (int i=0; i<3; ++i)
	{
		if(seq[i]=='X')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitX());
		else if(seq[i]=='Y')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitY());			
		else if(seq[i]=='Z')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitZ());					
	}
	Eigen::MatrixXd quat(1,4);
	Eigen::Quaterniond q(rot_mat);
	quat(0,0) = q.x();
	quat(0,1) = q.y();
	quat(0,2) = q.z();
	quat(0,3) = q.w();
    return quat; 
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::qt2eul(Eigen::MatrixXd quat, std::string seq)
{
	seq = rtf::validate_seq(seq);
	Eigen::Matrix3d rot_mat = rtf::qt2rot(quat);
	int rot_idx[3];
	for (int i=0; i<3; ++i)
	{
		if(seq[i]=='X')
			rot_idx[i] = 0;
		else if(seq[i]=='Y')
			rot_idx[i] = 1;
		else if(seq[i]=='Z')
			rot_idx[i] = 2;
	}	
	Eigen::MatrixXd eul_angles(1,3);
	Eigen::Vector3d eul_angles_vec;
	eul_angles_vec = rot_mat.eulerAngles(rot_idx[0], rot_idx[1], rot_idx[2]);
	eul_angles(0,0) = eul_angles_vec[0];
	eul_angles(0,1) = eul_angles_vec[1];
	eul_angles(0,2) = eul_angles_vec[2];	
	return eul_angles;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::get_rob_T_part(Eigen::MatrixXd part_pts, Eigen::MatrixXd rob_pts)
{
// part_pts: co-ordinates with respect to CAD part frame
// rob_pts: points with respect to robot base frame (or world frame if it matches with robot base frame)
// input: part_pts = [x1, y1, z1;
//                             x2, y2, z2;
//                                :
//                                :
//                             xn, yn, zn]
// input: rob_pts = [x1, y1, z1;
//                            x2, y2, z2;
//                                :
//                                :
//                            xn, yn, zn]
    Eigen::MatrixXd centroid_part_pts(1,part_pts.cols());
    Eigen::MatrixXd centroid_rob_pts(1,rob_pts.cols());
    Eigen::MatrixXd shifted_part_pts(part_pts.rows(),part_pts.cols());
    Eigen::MatrixXd shifted_rob_pts(rob_pts.rows(),rob_pts.cols());
    Eigen::MatrixXd cros_cov_mat(part_pts.cols(),rob_pts.cols());
    Eigen::Matrix3d R;
    Eigen::Matrix3d U_T;
    Eigen::Matrix3d V;
    Eigen::Vector3d T;
    Eigen::Matrix3d M = Eigen::Matrix3d::Identity();
    Eigen::Matrix4d transform_mat = Eigen::Matrix4d::Constant(0);
    if (part_pts.rows()==rob_pts.rows())
    { 
        centroid_part_pts = rtf::mean(part_pts);
        centroid_rob_pts = rtf::mean(rob_pts);
        for (int i=0;i<part_pts.rows();++i)
        {
            shifted_part_pts.block(i,0,1,shifted_part_pts.cols()) = part_pts.block(i,0,1,part_pts.cols()) - centroid_part_pts;
            shifted_rob_pts.block(i,0,1,shifted_rob_pts.cols()) = rob_pts.block(i,0,1,rob_pts.cols()) - centroid_rob_pts;   
        }
        cros_cov_mat = shifted_part_pts.transpose()*shifted_rob_pts;
    // Singular Value Decomposition
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(cros_cov_mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // Take care of reflection case due to negative eigen vectors
        U_T = svd.matrixU().transpose();    V = svd.matrixV();
        M(2,2) = (V*U_T).determinant();
        R = V*M*U_T;
        if (R.determinant()>0)
        {
            T = -R*centroid_part_pts.transpose() + centroid_rob_pts.transpose();
            transform_mat.block(0,0,3,3) = R;
            transform_mat.block(0,3,3,1) = T;
            transform_mat(3,3) = 1; 
        }
        else
        {
            std::cerr << "ERROR: Determinant of rotation matrix is negative..." << std::endl;
        }   
    }
    else
    {
        std::cerr << "ERROR: FUNCTION ERROR: For correspondance, number of rows of both matrices should be same..." << std::endl;
    }
    return transform_mat;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::mean(Eigen::MatrixXd mat)
{
    Eigen::VectorXd vec(mat.cols());
    for (int i=0;i<mat.cols();++i)
    {
        vec(i) =  mat.block(0,i,mat.rows(),1).sum()/mat.rows();
    }
    return vec.transpose();
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::rot_x(double t)
{
	Eigen::Matrix3d rx;
	rx << 	1,		0,		0,
			0, cos(t),-sin(t),
			0, sin(t), cos(t);
	return rx;  
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::rot_y(double t)
{
	Eigen::Matrix3d ry;
	ry << cos(t),	0, sin(t),
			   0,	1,		0,
		 -sin(t),	0, cos(t);
	return ry;  
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::rot_z(double t)
{
	Eigen::Matrix3d rz;
	rz << cos(t),-sin(t),	0,
		  sin(t), cos(t),	0,
			   0, 	   0,	1;
	return rz;  
}

///////////////////////////////////////////////////////////
