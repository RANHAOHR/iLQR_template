/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2021 Case Western Reserve University
 *
 *    Ran Hao <rxh349@case.edu>
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Case Western Reserve University, nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef ILQR_H
#define ILQR_H

#include <vector>
#include <cstdio>
#include <iostream>

#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <string>
#include <cstring>
#include <random>

#include "matrixOperation.h"
#include "matplotlibcpp.h"


#define PI 3.14159265

class ILQR{
protected:
    int D_state;
    int D_control;
    int itr;
    int N_path;

    bool converged;
    double stopTol = 1e-7;
    // Matrix size is in the order of col vector, first size is row size, 2nd size is col size
    std::vector< std::vector<double> > Q_T; // terminal state cost
    std::vector< std::vector<double> > Q; // state cost
    std::vector< std::vector<double> > R; // control cost

    double lambda;
    double lambda_factor;
    double max_lambda;

    std::vector<double> x_0;
    std::vector<double> x_d;

public:
    //outputs
    std::vector<double> final_u_t;
    std::vector<double> final_x_t;
    double final_cost;

    ILQR();

    void initialize_traj( std::vector< std::vector<double> > &u_t, std::vector< std::vector<double> > &x_t);
    void planning_main();

    void value_iteration(std::vector< std::vector<double> > &x_t_stack, std::vector< std::vector<double> > &u_t_stack);

    std::vector<double> cost_function(std::vector<double>  &x_t, std::vector<double>  &u_t, std::vector< std::vector<double> > &Q_, std::vector< std::vector<double> > &R_);
    std::vector<double> quadratic_cost(std::vector<double>  &input_t, std::vector<double>  &input_d, std::vector< std::vector<double> > &Q_);
    std::vector<double> quadratic_cost(std::vector<double>  &input_t, std::vector< std::vector<double> > &Q_);

    std::vector<double> dynamic_function(std::vector<double>  &x_t, std::vector<double>  &u_t); //wrapper function of the dynamics
    std::vector<double> car_dynamic(std::vector<double>  &x_t, std::vector<double>  &u_t);
    std::vector<double> linear_dynamic_function(std::vector<double>  &x_t, std::vector<double>  &u_t);

    using obj_func_ = std::vector<double> (ILQR::*)(std::vector<double>&, std::vector<double>&);
    using obj_func_v = std::vector<double> (ILQR::*)(std::vector<double>&, std::vector< std::vector<double> >&);

    std::vector< std::vector<double> > jacobian_forward_difference(obj_func_ func, auto &param_1, auto &param_2, auto  *variable_);
    std::vector< std::vector<double> > jacobian_forward_difference(obj_func_v func, auto &param_1, auto &param_2, auto  *variable_);
    std::vector< std::vector<double> > jacobian_forward_difference(obj_func_v func, auto &param_1, auto &param_2, auto &param_3, auto &param_4,  auto *variable_);
    std::vector< std::vector<double> > heissian_cost_func(obj_func_v func, auto &param_1, auto &param_2);

    void plot_convergence(std::vector<std::vector<double> > &x_stack, std::string& title);
    void plot_convergence(std::vector<std::vector<double> > &x_stack);
};
#endif //ILQR_H
