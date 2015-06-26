/*
*	@author seymayucer	
*	@date 26.06.2015
*   used opencv 3.0.0 beta
*	---------------------------------------------------------	
* 
*	Homography Estimation using Direct Linear Transformation
*
*	---------------------------------------------------------
*/
#ifndef DLT
#define DLT

#include "Libs.h"

//!
//	Loads 3D points from txt files
//!
bool loadPoints(string file1, bool isp1);
//!
// Matrix estimation using DLT method
// http://kwon3d.com/theory/dlt/dlt.html
//!
Mat DLTHomographyEstimation();
//!
//Gets avarage milimeter point error 
//!
double getHomograpyMatrixError(Mat H);

vector<Vec3d> pointSet1;//real kinect2 3d point for current files
vector<Vec3d> pointSet2;//mirror kinect2 3d point for current files can initiliaze different places

#if 1
int main()
{
	if (loadPoints("Points/PointSet1.txt", true) && loadPoints("Points/PointSet2.txt", false))
	{
		cout << "\nPoints are loaded. +" << pointSet1.size()<<" +"<<pointSet2.size();
		Mat H=DLTHomographyEstimation();
		getHomograpyMatrixError(H);
	}
	return 0;
}
#endif
bool loadPoints(string filename,bool isp1){
	
	ifstream file;
	file.open(filename);
	vector<double> all_el;
	Vec3d vecElem;
	double element;

	cout << endl << filename;
	if (!file.is_open()){
		cout << "\nfile not open! " << filename;
		return false;
	}
	while (!file.eof()){
		file >> element;
		//cout <<"\n"<< element;
		if (!file.eof())
			all_el.push_back(element);
	}
	file.close();
	int k = 0;	
	while (k < all_el.size())
	{
		vecElem[0] = all_el[k];
		vecElem[1] = all_el[k + 1];
		vecElem[2] = all_el[k + 2];
		k = k + 3;
		if (isp1) pointSet1.push_back(vecElem);
		else pointSet2.push_back(vecElem);
		
	}

}
Mat DLTHomographyEstimation(){
	Mat H;
	//Direct Linear Transform calculations
	if (pointSet1.size() > 0){
		int rowSize = pointSet1.size() * 4;// 3d - homogen world translation - 4d world
		Mat A = Mat::zeros(rowSize, 16, CV_64F);
		Mat b = Mat::zeros(rowSize, 1, CV_64F);
		int k; int l = -1;
		Mat aRow;

		for (int i = 0; i < rowSize; ++i){
			k = i % 4;

			switch (k){
			case 0:
				++l;
				aRow = (Mat_<double>(1, 16) << pointSet1[l][0], pointSet1[l][1], pointSet1[l][2], 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
				b.at<double>(i, 0) = pointSet2[l][0];
				for (int j = 0; j < 16; ++j)
					A.at<double>(i, j) = aRow.at<double>(0, j);
				break;
			case 1:
				aRow = (Mat_<double>(1, 16) << 0, 0, 0, 0, pointSet1[l][0], pointSet1[l][1], pointSet1[l][2], 1, 0, 0, 0, 0, 0, 0, 0, 0);
				b.at<double>(i, 0) = pointSet2[l][1];
				for (int j = 0; j < 16; ++j)
					A.at<double>(i, j) = aRow.at<double>(0, j);
				break;
			case 2:
				aRow = (Mat_<double>(1, 16) << 0, 0, 0, 0, 0, 0, 0, 0, pointSet1[l][0], pointSet1[l][1], pointSet1[l][2], 1, 0, 0, 0, 0);
				b.at<double>(i, 0) = pointSet2[l][2];
				for (int j = 0; j < 16; ++j)
					A.at<double>(i, j) = aRow.at<double>(0, j);
				break;
			case 3:
				aRow = (Mat_<double>(1, 16) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pointSet1[l][0], pointSet1[l][1], pointSet1[l][2], 1);
				b.at<double>(i, 0) = 1;
				for (int j = 0; j < 16; ++j)
					A.at<double>(i, j) = aRow.at<double>(0, j);
				break;

			}
		}

		/*cout << "\nA size:" << A.size()<<"\nA transpose size"<<A.t().size() << "\nb size:" << b.size();*/
		invert(A, A, DECOMP_SVD);

		H = A*b;

		H = (Mat_<double>(4, 4) << H.at<double>(0, 0), H.at<double>(1, 0), H.at<double>(2, 0), H.at<double>(3, 0),
			H.at<double>(4, 0), H.at<double>(5, 0), H.at<double>(6, 0), H.at<double>(7, 0),
			H.at<double>(8, 0), H.at<double>(9, 0), H.at<double>(10, 0), H.at<double>(11, 0),
			0, 0, 0, 1);//if distortion is important for spesific point set should be used H.at<double>(12,13,14,15, 0)....

		cout << "\nThis homography matrix is calculated with " << pointSet1.size() << " points.";
		cout << "\nH: " << H;

		return H;
	}
	
}
double getHomograpyMatrixError(Mat H){ 
	double totalError = 0;
	
	int ps = 0;

	Mat real_point = (Mat_<double>(4, 1));
	Mat mirror_point = (Mat_<double>(4, 1));

	double D = 0;

	for (int i = 0; i < pointSet1.size(); ++i){
		real_point = (Mat_<double>(4, 1) << pointSet1[i][0], pointSet1[i][1], pointSet1[i][2], 1);
		//cout << endl << real_point;
		mirror_point = H * real_point;

		Vec3d a;
		a[0] = mirror_point.at<double>(0, 0) / mirror_point.at<double>(3, 0);
		a[1] = mirror_point.at<double>(1, 0) / mirror_point.at<double>(3, 0);
		a[2] = mirror_point.at<double>(2, 0) / mirror_point.at<double>(3, 0);

		//	results.push_back(a);


		D = sqrt(pow(pointSet2[i][0] - a[0], 2) + pow(pointSet2[i][1] - a[1], 2) + pow(pointSet2[i][2] - a[2], 2));
	
		//per_point_error.push_back(D);
		totalError = totalError + D;
		//cout << "\ninf test:" << i << " " << mirror_corners[i][0];
	}
	//cout << "\n point size  " << pointSet1.size();
	cout << "\n error from estimation: (mm) " << totalError / pointSet1.size();

	return totalError;
}

#endif