#include "ILQR.h"
#include <iostream>

std::default_random_engine generator;
std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);

namespace plt = matplotlibcpp;

ILQR::ILQR() {
    D_state = 3;
    D_control = 2;
    itr = 100;
    N_path = 15;

    x_d.resize(D_state);
    x_d[0] = 15; x_d[1] = 1;

    Q_T.resize(D_state);
    Q.resize(D_state);
    for (int i = 0; i < D_state; ++i) {
        Q_T[i].resize(D_state);
        Q[i].resize(D_state);
    }

    Q_T[0][0] = 15000; Q_T[1][1] = 10000; Q_T[2][2] = 500000;
    Q[0][0] = 1; Q[1][1]  = 1; Q[2][2]  = 7;

    R.resize(D_control);
    for (int i = 0; i < D_control; ++i) {
        R[i].resize(D_control);
    }
    R[0][0] = 1; R[1][1] = 50;

    x_0.resize(D_state);
    x_0[0] = 0; x_0[1] = 0;

    lambda = 1;
    lambda_factor = 1.6;
    max_lambda = 1000000;

    converged = false;
}

void ILQR::initialize_traj(std::vector< std::vector<double> > &u_t, std::vector< std::vector<double> > &x_t){

    u_t.resize(N_path-1);
    for (int n_ = 0; n_ < N_path-1; ++n_) {
        u_t[n_].resize(D_control);
        for (int i = 0; i < D_control; ++i) {
            double rand_= uniform_distribution(generator) * 1;
            u_t[n_][i] = rand_;
        }
    }

    x_t.resize(N_path);
    x_t[0] = x_0;
    for (int n = 0; n < N_path-1; ++n) {
        x_t[n+1] = dynamic_function(x_t[n], u_t[n]);
    }

};
void ILQR::planning_main(){
    std::vector< std::vector<double> > x_t_stack;
    std::vector< std::vector<double> > u_t_stack;
    initialize_traj(u_t_stack, x_t_stack);
    //    for (int i = 0; i < u_t_stack.size(); ++i) {
//        std::cout << "u_t: " << u_t_stack[i][0] <<  " " << u_t_stack[i][1] << std::endl;
//    }
//    for (int i = 0; i < x_t_stack.size(); ++i) {
//        std::cout << "x_t: " << x_t_stack[i][0] <<  " " << x_t_stack[i][1] << " " << x_t_stack[i][2] << std::endl;
//    }
    value_iteration(x_t_stack, u_t_stack);

    final_cost = 0;
    for (int i = 0; i < N_path-1; ++i) {
        std::vector<double> cost_temp = cost_function(x_t_stack[i], u_t_stack[i],Q, R);
        final_cost += cost_temp[0];
    }
    std::vector<double> cost_T = quadratic_cost(x_t_stack[N_path-1], x_d, Q_T);
    final_cost += cost_T[0];
    double scale_ = norm(Q_T) + norm(Q) + norm(R);
    final_cost /= scale_;
    std::cout << "Final cost: " << final_cost << std::endl;
}

void ILQR::value_iteration(std::vector< std::vector<double> > &x_t_stack, std::vector< std::vector<double> > &u_t_stack){

    std::vector< std::vector<double> > K_1; //K_1: [Dim_control, N_path]
    std::vector< std::vector< std::vector<double> > > K_2; //[Dim_control, Dim_state, N_path]
    K_1.resize(D_control);
    for (int i = 0; i < D_control; ++i) {
        K_1[i].resize(N_path);
    }

    K_2.resize(D_control);
    for (int i = 0; i < D_control; ++i) {
        K_2[i].resize(D_state);
        for (int j = 0; j < D_state; ++j) {
            K_2[i][j].resize(N_path);
        }
    }

    std::string title_ = "Convergence plot";

    for (int cnt = 0; cnt < itr; ++cnt) {

        std::vector<double> d_x = diff_vector(x_t_stack[N_path-1], x_d );
        std::vector<double> Vx = left_Mult_AB(Q_T, d_x);
        std::vector< std::vector<double> > Vxx = Q_T;

        for (int ipath = N_path-2; ipath >=0 ; --ipath) {
            std::vector<double> x_t_temp = x_t_stack[ipath];
            std::vector<double> u_t_temp = u_t_stack[ipath];
            std::vector<double> *variable_x;   variable_x = &x_t_temp;
            std::vector<std::vector<double>> fx = jacobian_forward_difference(&ILQR::dynamic_function, x_t_temp, u_t_temp, variable_x);

//            std::cout << "fx: " << std::endl;
//            for (int j = 0; j < fx.size(); ++j) {
//                std::cout << fx[j][0] << " " << fx[j][1] << " " << fx[j][2] <<  std::endl;
//            }

            x_t_temp = x_t_stack[ipath];
            u_t_temp = u_t_stack[ipath];
            std::vector<double> *variable_u;   variable_u = &u_t_temp;
            std::vector< std::vector<double> > fu = jacobian_forward_difference(&ILQR::dynamic_function, x_t_temp, u_t_temp, variable_u);

//            std::cout << "fu: " << std::endl;
//            for (int j = 0; j < fu.size(); ++j) {
//                std::cout << fu[j][0] << " " << fu[j][1]  <<  std::endl;
//            }

            //Differentiate cost functions
            x_t_temp = x_t_stack[ipath];
            u_t_temp = u_t_stack[ipath];
            std::vector<double> *variable_cx;   variable_cx = &x_t_temp;
            std::vector<std::vector<double> > cx_ = jacobian_forward_difference(&ILQR::quadratic_cost, x_t_temp, Q, variable_cx);
            std::vector<double> cx = cx_[0];

            x_t_temp = x_t_stack[ipath];
            u_t_temp = u_t_stack[ipath];
            std::vector<double> *variable_cu;   variable_cu = &u_t_temp;
            std::vector<std::vector<double> > cu_ = jacobian_forward_difference(&ILQR::quadratic_cost, u_t_temp, R, variable_cu);
            std::vector<double> cu = cu_[0];

            x_t_temp = x_t_stack[ipath];

            std::vector< std::vector<double> > cxx = heissian_cost_func(&ILQR::quadratic_cost, x_t_temp, Q);

            u_t_temp = u_t_stack[ipath];

            std::vector< std::vector<double> > cuu = heissian_cost_func(&ILQR::quadratic_cost, u_t_temp, R);


            std::vector<double> Qx_temp = right_Mult_AB(Vx, fx);
            std::vector<double> Qx = add_vector(cx, Qx_temp);

            std::vector<double> Qu_temp = right_Mult_AB(Vx, fu);

            std::vector<double> Qu = add_vector(cu, Qu_temp);

            std::vector< std::vector<double> > fx_T = transpose_matrix(fx);
            std::vector< std::vector<double> > Qxx_temp = mult_ABC(fx_T, Vxx, fx );
            std::vector< std::vector<double> > Qxx = add_matrix(cxx, Qxx_temp);

            std::vector< std::vector<double> > fu_T = transpose_matrix(fu);
            std::vector< std::vector<double> > Quu_temp = mult_ABC(fu_T, Vxx, fu );

            std::vector< std::vector<double> > Quu = add_matrix(cuu, Quu_temp);


            std::vector< std::vector<double> > Qux = mult_ABC(fu_T, Vxx, fx ); //set cux to zero

            std::vector<std::vector<double> > inv_Quu;
            bool proceed_ = inverse(Quu, inv_Quu);

            if (proceed_){
                std::vector<double> k_1 = left_Mult_AB(inv_Quu, Qu);
                std::vector<double> k_1_vec; k_1_vec.resize(D_control);
                for (int ik = 0; ik < D_control; ++ik) {
                    K_1[ik][ipath] = -k_1[ik];
                    k_1_vec[ik] =  -k_1[ik];
                }

                std::vector< std::vector<double> > k_2 = mult_AB(inv_Quu, Qux);
                std::vector< std::vector<double> > k_2_mat; k_2_mat.resize(D_control);
                for (int i = 0; i < D_control; ++i) {
                    k_2_mat[i].resize(D_state);
                    for (int j = 0; j < D_state; ++j) {
                        K_2[i][j][ipath] = -k_2[i][j];
                        k_2_mat[i][j] = -k_2[i][j];
                    }
                }

                //iLQR value iteration
                std::vector< double> temp_ = left_Mult_AB(Quu, k_1_vec);
                std::vector< std::vector<double> > k_2_T = transpose_matrix(k_2_mat);
                std::vector< double > Vx_1 = left_Mult_AB(k_2_T, temp_);
                std::vector< double > Vx_2 = right_Mult_AB(Qu, k_2_mat);

                std::vector< double > Vx_3 = right_Mult_AB(k_1_vec, Qux);

                for (int i = 0; i < D_state; ++i) {
                    Vx[i] = Qx[i] + Vx_1[i] + Vx_2[i] + Vx_3[i];
                }

                std::vector<std::vector<double>> Vxx_1 = mult_ABC(k_2_T, Quu, k_2_mat);
                std::vector<std::vector<double>> Vxx_2 = mult_AB(k_2_T, Qux);

                std::vector<std::vector<double>> Q_ux_T = transpose_matrix(Qux);
                std::vector<std::vector<double>> Vxx_3 = mult_AB(Q_ux_T, k_2_mat);
                for (int i = 0; i < D_state; ++i) {
                    for (int j = 0; j < D_state; ++j) {
                        Vxx[i][j] = Qxx[i][j] + Vxx_1[i][j] + Vxx_2[i][j] + Vxx_3[i][j];
                    }
                }

                std::vector<std::vector<double>> Vxx_T = transpose_matrix(Vxx);
                std::vector<std::vector<double>> Vxx_temp = add_matrix(Vxx, Vxx_T);
                for (int i = 0; i < D_state; ++i) {
                    for (int j = 0; j < D_state; ++j) {
                        Vxx[i][j] = 0.5*Vxx_temp[i][j];
                    }
                }
//                std::cout << "Vxx: "<< std::endl;
//                for (int i = 0; i < Vxx.size(); ++i) {
//                    std::cout << Vxx[i][0] << " " << Vxx[i][1] << " " << Vxx[i][2]  << std::endl;
//                }
            }else{
                std::cout << "Singular Quu, Quiting..." << std::endl;
                return;
            }
        }

        //forward propagation
        double diff_cost = 0;
        std::vector<double> x_new = x_0;
        for (int ipath = 0; ipath < N_path-1; ++ipath) {
            std::vector< std::vector<double> > k_2_mat; k_2_mat.resize(D_control);
            for (int i = 0; i < D_control; ++i) {
                k_2_mat[i].resize(D_state);
                for (int j = 0; j < D_state; ++j) {
                    k_2_mat[i][j] = K_2[i][j][ipath];
                }
            }

            std::vector<double> dx = diff_vector(x_new, x_t_stack[ipath]);
            std::vector<double> temp = left_Mult_AB(k_2_mat, dx);
            for (int i = 0; i < D_control; ++i) {
                u_t_stack[ipath][i] = u_t_stack[ipath][i]+ K_1[i][ipath] + temp[i];
            }

            std::vector<double> x_nextstep = dynamic_function(x_new, u_t_stack[ipath]) ;

            diff_cost += norm(x_nextstep, x_t_stack[ipath+1]);
            x_t_stack[ipath+1] = x_nextstep;
            x_new = x_nextstep;
            if(diff_cost < stopTol)
                converged = true;
        }
        if (converged)
            break;

        plot_convergence(x_t_stack, title_);
        std::cout << "itr: " << cnt << std::endl;
    }

    plt::figure(2);
    std::string res_title = "Final plot";
    plot_convergence(x_t_stack);

//    std::vector<double>  x_t_temp = x_t_stack[0];
//    std::cout << "Heissian: " << std::endl;
//    std::vector< std::vector<double> > Heissian_ = heissian_cost_func(&ILQR::quadratic_cost, x_t_temp, Q_T);
//    for (int i = 0; i < Heissian_.size(); ++i) {
//        std::cout << Heissian_[i][0] << " " << Heissian_[i][1] << " " << Heissian_[i][2] << std::endl;
//    }
}

std::vector<double> ILQR::dynamic_function(std::vector<double>  &x_t, std::vector<double>  &u_t){
    std::vector<double> x_new = car_dynamic(x_t, u_t);
    return x_new;
}

std::vector<double> ILQR::linear_dynamic_function(std::vector<double>  &x_t, std::vector<double>  &u_t){
    std::vector< std::vector<double> > A;
    A.resize(D_state);
    for (int i = 0; i < D_state; ++i) {
        A[i].resize(D_state);
    }
    A[0][0] = 15; A[1][1] = 2; A[2][2] = 10;

    std::vector< std::vector<double> > B;
    B.resize(D_state);
    for (int i = 0; i < D_state; ++i) {
        B[i].resize(D_control);
    }
    B[0][0] = 24; B[1][1] = 10;B[2][0] = 15;B[2][1] = 7;

    std::vector<double> x_new1 = left_Mult_AB(A, x_t);
    std::vector<double> x_new2 = left_Mult_AB(B, u_t);
    std::vector<double> x_new; x_new.resize(D_state);
    for (int i = 0; i < D_state; ++i) {
        x_new[i] = x_new1[i] + x_new2[i];
    }

    return x_new;
}

void ILQR::plot_convergence(std::vector<std::vector<double> > &x_stack, std::string& title){
    int xlim = x_stack.size();

    std::vector<double> x_vec(xlim), y_vec(xlim), theta_vec(xlim), path_(xlim);

    for (int ipath = 0; ipath < xlim; ++ipath) {
        path_[ipath] = ipath;
        x_vec[ipath] = x_stack[ipath][0];
        y_vec[ipath] = x_stack[ipath][1];
        theta_vec[ipath] = x_stack[ipath][2];
    }
    //PLOT x
    plt::suptitle(title);
    plt::subplot(1, 3, 1);

    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(path_, x_vec, "r-");
    //plot goal data
    double end_ = path_.at(xlim-1);
    std::vector<double> x_(1, end_), y_(1,x_d[0]);
    plt::scatter(x_, y_, 50);
    plt::xlim(path_[0], end_);

    //PLOT y
    plt::subplot(1, 3, 2);
    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(path_, y_vec, "g-");
    //plot goal data
    std::vector<double> yy_(1,x_d[1]);
    plt::scatter(x_, yy_, 50);
    plt::xlim(path_[0], end_);

    //PLOT theta
    plt::subplot(1, 3, 3);
    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(path_, theta_vec, "k-");
    //plot goal data
    std::vector<double> ty_(1,x_d[2]);
    plt::scatter(x_, ty_, 50);
    plt::xlim(path_[0], end_);

    plt::draw();
    plt::pause(0.05);

}

void ILQR::plot_convergence(std::vector<std::vector<double> > &x_stack){

    int xlim = x_stack.size();

    std::vector<double> x_vec(xlim), y_vec(xlim), theta_vec(xlim), path_(xlim);

    for (int ipath = 0; ipath < xlim; ++ipath) {
        path_[ipath] = ipath;
        x_vec[ipath] = x_stack[ipath][0];
        y_vec[ipath] = x_stack[ipath][1];
        theta_vec[ipath] = x_stack[ipath][2];
    }
    //PLOT x
    plt::suptitle("Final plot");
    plt::subplot(1, 3, 1);

    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(path_, x_vec, "r-");
    //plot goal data
    double end_ = path_.at(xlim-1);
    std::vector<double> x_(1, end_), y_(1,x_d[0]);
    plt::scatter(x_, y_, 50);
    plt::xlim(path_[0], end_);

    //PLOT y
    plt::subplot(1, 3, 2);
    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(path_, y_vec, "g-");
    //plot goal data
    std::vector<double> yy_(1,x_d[1]);
    plt::scatter(x_, yy_, 50);
    plt::xlim(path_[0], end_);


    //PLOT theta
    plt::subplot(1, 3, 3);
    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(path_, theta_vec, "k-");
    //plot goal data
    std::vector<double> ty_(1,x_d[2]);
    plt::scatter(x_, ty_, 50);
    plt::xlim(path_[0], end_);

    plt::show();
}

std::vector<double> ILQR::car_dynamic(std::vector<double>  &x_t, std::vector<double>  &u_t){
//    if (D_state != 3){ printf("Wrong model...");}
    double delta_t = 0.01;
    double angle = x_t[2];
    double v_ = u_t[0];
    double w_ = u_t[1];
    double x_ = v_ * cos(angle);
    double y_ = v_ * sin(angle);
    std::vector<double> x_t_new; x_t_new.resize(D_state);
    x_t_new[0] = x_*delta_t + x_t[0];  x_t_new[1] = y_*delta_t + x_t[1]; x_t_new[2] = w_*delta_t + x_t[2];

    return x_t_new;
};

std::vector<double> ILQR::cost_function(std::vector<double>  &x_t, std::vector<double>  &u_t, std::vector< std::vector<double> > &Q_, std::vector< std::vector<double> > &R_){
    std::vector<double> x_cost = quadratic_cost(x_t, Q_);
    std::vector<double> u_cost = quadratic_cost(u_t, R_);
    std::vector<double> cost_(1);
    for (int i = 0; i < cost_.size(); ++i) {
        cost_[i] = x_cost[i] + u_cost[i];
    }

    return cost_;
}


std::vector<double> ILQR::quadratic_cost(std::vector<double>  &input_t, std::vector<double>  &input_d, std::vector< std::vector<double> > &Q_){
    std::vector<double> d_x; d_x.resize(input_t.size());
    for (int x = 0; x < d_x.size(); ++x) {
        d_x[x] = input_t[x] - input_d[x];
    }

    std::vector<double> left_mul_vec = left_Mult_AB(Q_, d_x);
    std::vector<double> cost_; cost_.resize(1);
    cost_[0] = 0.5*dot_product ( d_x, left_mul_vec);

    return cost_;
};

std::vector<double>  ILQR::quadratic_cost(std::vector<double>  &input_t, std::vector< std::vector<double> > &Q_){
    std::vector<double> left_mul_vec = left_Mult_AB(Q_, input_t);
    std::vector<double> cost_; cost_.resize(1);
    cost_[0] = 0.5*dot_product ( input_t, left_mul_vec);
    return cost_;
};

std::vector< std::vector<double> >  ILQR::jacobian_forward_difference(obj_func_ func, auto &param_1, auto &param_2, auto  *variable_ ){
    std::vector< std::vector<double> > Df_;
    std::vector<double> f_0 = (this->*func)(param_1, param_2);
    double h = 1e-08;
    int f_size = f_0.size();

    std::vector<double> x_0 = *variable_;

    int var_size = (*variable_).size();
    Df_.resize(f_size);
    for (int i = 0; i < f_size; ++i) {
        Df_[i].resize(var_size);
    }
    for (int j = 0; j < var_size; ++j) {
        std::vector<double> x_h = x_0;  x_h[j] = x_0[j] + h;
        *variable_ = x_h;
        std::vector<double> f_h = (this->*func)(param_1, param_2);
        for (int k = 0; k < f_size; ++k) {
            Df_[k][j] = (f_h[k] - f_0[k]) / h;
        }
    }
    return Df_;
}

std::vector< std::vector<double> > ILQR::jacobian_forward_difference(obj_func_v func, auto &param_1, auto &param_2, auto  *variable_){
    std::vector< std::vector<double> > Df_;
    std::vector<double> f_0 = (this->*func)(param_1, param_2);
    double h = 1e-8;
    int f_size = f_0.size();

    std::vector<double> x_0 = *variable_;
    int var_size = (*variable_).size();
    Df_.resize(f_size);
    for (int i = 0; i < f_size; ++i) {
        Df_[i].resize(var_size);
    }
    for (int j = 0; j < var_size; ++j) {
        std::vector<double> x_h = x_0;  x_h[j] = x_0[j] + h;
        *variable_ = x_h;
        std::vector<double> f_h = (this->*func)(param_1, param_2);
        for (int k = 0; k < f_size; ++k) {
            Df_[k][j] = (f_h[k] - f_0[k]) / h;
        }
    }
    return Df_;
}

std::vector< std::vector<double> > ILQR::jacobian_forward_difference(obj_func_v func, auto &param_1, auto &param_2, auto &param_3, auto &param_4,  auto *variable_){
    std::vector< std::vector<double> > Df_;
    std::vector<double> f_0 = (this->*func)(param_1, param_2, param_3, param_4);
    double h = 1e-8;
    int f_size = f_0.size();

    std::vector<double> x_0 = *variable_;
    int var_size = (*variable_).size();
    Df_.resize(f_size);
    for (int i = 0; i < f_size; ++i) {
        Df_[i].resize(var_size);
    }
    for (int j = 0; j < var_size; ++j) {
        std::vector<double> x_h = x_0;  x_h[j] = x_0[j] + h;
        *variable_ = x_h;
        std::vector<double> f_h = (this->*func)(param_1, param_2, param_3, param_4);
        for (int k = 0; k < f_size; ++k) {
            Df_[k][j] = (f_h[k] - f_0[k]) / h;
        }
    }
    return Df_;
}

/*
 * Computing Heissian of the cost function numerically
 */
std::vector< std::vector<double> > ILQR::heissian_cost_func(obj_func_v func, auto &param_1, auto &param_2){
    double h = 2^-8;
    std::vector<double> x_0 = param_1;

    int var_size = param_1.size();
    std::vector< std::vector<double> > DDF; DDF.resize(var_size);
    for (int i = 0; i < var_size; ++i) {
        DDF[i].resize(var_size);
    }

    for (int j = 0; j < var_size; ++j) {
        std::vector<double> x_1 = x_0 ; x_1[j] = x_0[j] - h;
        std::vector<double> x_2 = x_0 ; x_2[j] = x_0[j] + h;

        std::vector<double> *variable_1;   variable_1 = &x_1;
        std::vector<double> *variable_2;   variable_2 = &x_2;

        std::vector< std::vector<double> > f_h1 = jacobian_forward_difference(func, x_1, param_2, variable_1);
        std::vector< std::vector<double> > f_h2 = jacobian_forward_difference(func, x_2, param_2, variable_2);
        for (int k = 0; k < var_size; ++k) {
            DDF[k][j] = (f_h2[0][k] - f_h1[0][k]) / (2* h);
        }
    }

    return DDF;

};

