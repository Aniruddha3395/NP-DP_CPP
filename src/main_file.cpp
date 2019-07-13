#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <igl/readSTL.h>
#include "stdlib.h"
#include "utilities.hpp"
#include "file_rw.hpp"
#include "transformation_utilities.hpp"
#include "NPAM_utilities.hpp"
#include <cstdlib>
#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h"
#include "plot_utilities.hpp"

int main()
{
    ///////////////////////////////////////////////////////
    ////////////////////// INPUTS /////////////////////////
	///////////////////////////////////////////////////////

	// if true generate all layers, else generate only one layer
	bool generate_all_layers = false;
	bool generate_cont_toolpath = true;		// add start and end point for each layer

	// Gap between 2 hatching lines -
	double pathgap_x = 1.0;
	double pathgap_y = 1.0;

	//Gap between 2 layers -
	double pathgap_z = 0.5;
	
	// grid aadition (NOTE: higher number for higher aspect ratio of part in xy plane)
	double grid_addition = 70.0;

	// hatching angle properties
	double start_hatch_angle = 0; //3.14159/4;
	double hatch_angle_change = 3.14159/4;	// consecutive layers will have this much change in hatching angle (degrees)
	
	// Path to generate -
	// 1 - Boundary
	// 2 - Hatching
	int Path_Number = 2;    // put number from above
	int FlipTravel = 0;     // 1 for yes, 0 for No
	int space = 2;          // number of points to skip for smoother path
	int skip_lines = 0;     // number of points to skip along the tool path
    ///////////////////////////
    // Bottom layer generation
	///////////////////////////

	Eigen::MatrixXd v;
	Eigen::MatrixXd f;
	Eigen::MatrixXd n;
	std::string stl_file_name = "/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/test_part4.STL";
	igl::readSTL(stl_file_name, v, f, n);
	f = f.array()+1;	//NOTE : adding 1 make make start index as 1... this is compatible for rest of the code 
    Eigen::MatrixXd fillpts;

	//////////////////////////////////////
    // xmax and ymax to get the grid size
	//////////////////////////////////////

    double xmax = v.block(0,0,v.rows(),1).maxCoeff() + grid_addition;
	double ymax = v.block(0,1,v.rows(),1).maxCoeff() + grid_addition;
	double xmin = v.block(0,0,v.rows(),1).minCoeff() - grid_addition;
    double ymin = v.block(0,1,v.rows(),1).minCoeff() - grid_addition;

	//////////////////////////////////////////////////////////////////////////////
    // finding thickness of shell part and number of layers required for printing
	//////////////////////////////////////////////////////////////////////////////

	int num_of_layers = 0;
    if (generate_all_layers)
    {
		num_of_layers = NPAM::number_of_layers(v,f,n,pathgap_z);
	}
	else
	{
		num_of_layers = 1;
	}	

	//////////////////////////////////////////////////////
    // identifying the co-ordinates of bottom layer first
	//////////////////////////////////////////////////////
	
	Eigen::MatrixXd fnew = NPAM::identify_bottom_layer(v,f,n);

	////////////////////////////
    // generating the tool path
	////////////////////////////

	Eigen::MatrixXd	pts;
	double hatch_angle;
	double x_avg;
	double y_avg;
	long int tool_path_strt_idx = 0;
	Eigen::MatrixXd tool_path(1,1);

	for (int layer=1;layer<=num_of_layers;++layer)
	{
		hatch_angle = start_hatch_angle + (layer-1)*hatch_angle_change;
		pts = NPAM::generate_grid_points(pathgap_x,pathgap_y,xmin,ymin,xmax,ymax,hatch_angle);

		// apply rotation to points
		x_avg = pts.block(0,0,pts.rows(),1).sum()/pts.rows();
		y_avg = pts.block(0,1,pts.rows(),1).sum()/pts.rows();
		Eigen::MatrixXd rotated_pts = NPAM::rotate_pts(pts,hatch_angle,x_avg,y_avg);

		////////////////////////////
		// project grid points		
		////////////////////////////
		fillpts = NPAM::project_grid_points(fnew, v, rotated_pts, hatch_angle, x_avg, y_avg);
	
		////////////////////////////
	    // infill path
		////////////////////////////
		Eigen::MatrixXd	layer_tool_path = NPAM::Infill_Path_with_bxbybz(fillpts, FlipTravel, space, hatch_angle, x_avg, y_avg, skip_lines);
		layer_tool_path.col(2) = layer_tool_path.col(2).array() + ((layer-1)*pathgap_z);
		
		if (generate_cont_toolpath)
		{	
			////////////////////////////////
			// Making continunous tool_path:
			// adding start and end points
			////////////////////////////////
			tool_path.conservativeResize(layer_tool_path.rows()+tool_path_strt_idx+2,layer_tool_path.cols());
			tool_path.block(tool_path_strt_idx,0,(layer_tool_path.rows()+2),layer_tool_path.cols()) = Eigen::MatrixXd::Constant((layer_tool_path.rows()+2),layer_tool_path.cols(),0);
			Eigen::MatrixXd start_pt(1,layer_tool_path.cols());
			start_pt = layer_tool_path.row(0);
	    	start_pt(0,2) = start_pt(0,2)+50;
	    	Eigen::MatrixXd end_pt(1,layer_tool_path.cols());
	    	end_pt = layer_tool_path.row(layer_tool_path.rows()-1);
	    	end_pt(0,2) = end_pt(0,2)+50;
	    	tool_path.block(tool_path_strt_idx,0,1,layer_tool_path.cols()) = start_pt;
	    	++tool_path_strt_idx;
	    	tool_path.block(tool_path_strt_idx,0,layer_tool_path.rows(),layer_tool_path.cols()) = layer_tool_path;
	    	tool_path_strt_idx = tool_path_strt_idx + layer_tool_path.rows();
	    	tool_path.block(tool_path_strt_idx,0,1,layer_tool_path.cols()) = end_pt;
	    	++tool_path_strt_idx;
		}
		else
		{	
			////////////////////////////////
			// Making toolpath layer by layer 
			////////////////////////////////
			tool_path.conservativeResize(layer_tool_path.rows()+tool_path_strt_idx,layer_tool_path.cols());
			tool_path.block(tool_path_strt_idx,0,(layer_tool_path.rows()),layer_tool_path.cols()) = Eigen::MatrixXd::Constant((layer_tool_path.rows()),layer_tool_path.cols(),0);
	    	tool_path.block(tool_path_strt_idx,0,layer_tool_path.rows(),layer_tool_path.cols()) = layer_tool_path;
	    	tool_path_strt_idx = tool_path_strt_idx + layer_tool_path.rows();
		}
	}

	file_rw::file_write("tool_path.csv",tool_path);
	std::cout << "Tool Path Generated!" << std::endl;
	
	// plotting final toolpath
	plt * plotptr = new plt();
	plotptr->xlabel("X-axis");
	plotptr->ylabel("Y-axis");
	plotptr->zlabel("Z-axis");
	plotptr->grid();
	plotptr->autoscale("ALL");
	plotptr->title("toolpath plot");
	plotptr->ActivePlot3d(tool_path, "blue", 1, "ToolPath");

    return 0;
}
