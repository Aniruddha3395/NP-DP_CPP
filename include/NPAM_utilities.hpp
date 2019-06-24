#ifndef NPAM_UTILITIES_HPP
#define NPAM_UTILITIES_HPP

#include <Eigen/Dense>
#include <vector>
#include <string>

class NPAM
{
public:

	// apply rotation to points about their centroid
	static Eigen::MatrixXd rotate_pts(Eigen::MatrixXd, double, double, double);
	
	// find number of layers
	static int number_of_layers(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, double);
	
	// generate uniform grid points on a plane 
	static Eigen::MatrixXd generate_grid_points(double, double, double, double, double, double, double);
	
	// identify bottom layer from stl file data
	static Eigen::MatrixXd identify_bottom_layer(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd);

	// identify top layer from stl file data
	static Eigen::MatrixXd identify_top_layer(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd);

	// project data on the bottom plane
	static Eigen::MatrixXd project_grid_points(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, double, double, double);

	// generate the hatching path with points
	static Eigen::MatrixXd Infill_Path(Eigen::MatrixXd, bool, double, double, double, double);
	static Eigen::MatrixXd Infill_Path_with_Normals(Eigen::MatrixXd, bool, double, double, double, double);
	static Eigen::MatrixXd Infill_Path_with_bxbybz(Eigen::MatrixXd, bool, double, double, double, double);
	static Eigen::MatrixXd Infill_Path_with_euler(Eigen::MatrixXd, bool, double, double, double, double);
};

#endif