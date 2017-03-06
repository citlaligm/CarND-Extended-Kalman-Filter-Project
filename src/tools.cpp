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
        if (estimations.size()==0 |estimations.size()!=ground_truth.size() ) throw 20;
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
    
	//check division by zero
	try
    {
        if(px == 0 | py == 0) throw 20;
        else 
        {
            float px2 = px*px;
            float py2 = py*py;
            float p_00 = px/sqrt(px2 +py2);
            float p_01 = py/sqrt(px2 +py2);
            float p_10 = -py/(px2 +py2);
            float p_11 = px/(px2 +py2);
            float p_20 = (py*(vx*py-vy*px))/pow((px2 +py2),1.5);
            float p_21 = (px*(vy*px-vx*py))/pow((px2 +py2),1.5);
            float p_22 = px/sqrt(px2 +py2);
            float p_23 = px/sqrt(px2 +py2);
            
            //compute the Jacobian matrix
            Hj << p_00, p_01, 0, 0,
			  p_10, p_11, 0, 0,
			  p_20, p_21, p_22, p_23;
            
        }
            
    }
    catch (int e)
    {
        cout << "CalculateJacobian() - Error - Division by Zero " << '\n';
    }
	
	
	return Hj;
}
