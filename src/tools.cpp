#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
	rmse << 0,0,0,0;

    // check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size	
	try
    {
        if ((estimations.size()==0) |(estimations.size()!=ground_truth.size()) ) throw 20;
    }
        catch (int e)
    {
        cout<<"estimation vector size should not be zero OR estimation vector size should be equal to ground truth vector size"<<"\n";
    }

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        //cout<<"Estimation:\n"<<estimations[i]<<"\n";
        //cout<<"Ground_truth:\n"<<ground_truth[i]<<"\n";
        VectorXd residual=estimations[i]-ground_truth[i];
        //cout<<"Residual: "<<residual<<"\n";
        residual = residual.array()*residual.array();
		//cout<<"Residual squared: "<<residual<<"\n";
		rmse += residual;
		//cout<<"RMSE: "<<rmse<<"\n";
        
	}

	//calculate the mean
	int n = estimations.size();
	rmse = rmse/n;

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);
    
	//check division by zero
	try
    {
        if(c1<0.0001)
        {

        	throw 20;
        }
        else 
        {
        	//compute the Jacobian matrix
        	Hj << (px/c2), (py/c2), 0, 0,
        		  -(py/c1), (px/c1), 0, 0,
				  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
            
        }
            
    }
    catch (int e)
    {
        cout << "CalculateJacobian() - Error - Division by Zero " << '\n';
    }
	
	
	return Hj;
}
